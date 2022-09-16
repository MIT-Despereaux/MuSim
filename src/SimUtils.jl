### This part contains functions used for simulations

using FLoops
using SparseArrays
using DataStructures
using Distributions
using BangBang, MicroCollections
using Random

"""
Returns the relative direction and error going from center of T1 to center of T2, in spherical coordinates.
"""
function _relativedir(T1::RectBox{T}, T2::RectBox{T}) where {T<:Real}
    center_1 = T1.position
    center_2 = T2.position
    dir = center_2 - center_1
    (θ, φ) = _cart2unitsph(dir...)
    θ_max = θ
    θ_min = θ
    φ_max = φ
    φ_min = φ
    for x_sign in (-1, 1)
        for y_sign in (-1, 1)
            for z_sign in (-1, 1)
                v = (x_sign * T2.delta_x / 2, y_sign * T2.delta_y / 2, z_sign * T2.delta_z / 2) |> SVector{3,Float64}
                v = _rotate(v, T2.orientation..., 0)
                dir2 = center_2 - center_1 + v
                (θ2, φ2) = _cart2unitsph(dir2...)
                θ_max = max(θ_max, θ2)
                θ_min = min(θ_min, θ2)
                φ_max = max(φ_max, φ2)
                φ_min = min(φ_min, φ2)
            end
        end
    end
    Δθ = θ_max - θ_min
    Δφ = φ_max - φ_min
    return (θ, φ, Δθ, Δφ)
end

"""
Runs the simulation with a hemispherical generating surface. 
Optional parameters θ_range and φ_range denotes the range of solid angles
for the Rays to spawn. The solid angle vector extends from the center of
the hemisphere to its surface. 
By default, the simulation will center on the first detector.
"""
function runhemisim(n_sim::Int, detectors::Vector{T},
    R::Real, ℓ::Real;
    exec=ThreadedEx(), center::Union{NTuple{3,Real},Nothing}=nothing) where {T<:LabObject{<:Real}}
    if center === nothing
        center = detectors[1].position |> Tuple
        optimize = true
    else
        optimize = false
    end
    if optimize
        # Find the relative direction between other detectors and the first detector
        rel_dirs = []
        for i in 2:length(detectors)
            push!(rel_dirs, _relativedir(detectors[1], detectors[i]))
        end
        mix_dist_θ = MixtureModel([Normal(θ > π / 2 ? π - θ : θ, Δθ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_θ = truncated(mix_dist_θ, 0, π / 2)
        mix_dist_φ = MixtureModel([Normal(φ, Δφ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_φ = truncated(mix_dist_φ, 0, 2π)
    else
        mix_dist_θ = nothing
        mix_dist_φ = nothing
    end
    let optimize = optimize, mix_dist_θ = mix_dist_θ, mix_dist_φ = mix_dist_φ, center = center
        @floop exec for i = 1:n_sim
            # Private mutable variables
            @init begin
                crosses = SortedDict{Float64,SVector{3,Float64}}()
                hit_vec = falses(length(detectors))
                i_vec = zeros(length(detectors))
                j_vec = zeros(length(detectors))
                crx = Vector{Union{Missing,Dict}}(missing, length(detectors))
                dir = Vector{Union{Missing,Dict}}(missing, length(detectors))
                angles = zeros(2)
            end
            # Set the hit_vec to false to prepare for a new ray
            hit_vec .*= false
            # Generate a random ray
            if optimize
                θ_sample = rand(mix_dist_θ)
                φ_sample = rand(mix_dist_φ)
                angles = (θ_sample, φ_sample)
                θs = π - θ_sample
                φs = π + φ_sample
            else
                θs = nothing
                φs = nothing
            end
            ray = Ray(R, center, ℓ; θ=θs, φ=φs)
            # Loop through each dectector to fill the pre-allocated vectors
            for (j, d) in enumerate(detectors)
                hit = isthrough!(ray, d, crosses)
                if hit
                    i_vec[j] = i
                    j_vec[j] = j
                    hit_vec[j] = hit
                    points = hcat(values(crosses)...)'
                    # Note the points are vertically concatenated
                    crx[j] = Dict(i => points)
                    dir[j] = Dict(i => Float64[ray.azimuth_agl, ray.polar_agl])
                    empty!(crosses)
                end
            end
            # Hit indices to be filled
            ii = i_vec[hit_vec]
            jj = j_vec[hit_vec]
            # Reduce the accumulators
            # Note that for threaded executions, the initial values will be
            # assigned to the variables (e.g. ii) at the end of the loop
            # This means e.g. ii will be Int[] at some time
            @reduce() do (i_list = Float64[]; ii),
            (j_list = Float64[]; jj),
            (crx_pts = [Dict{Int,Matrix{Float64}}() for i = 1:length(detectors)]; crx),
            (ray_dirs = [Dict{Int,Vector{Float64}}() for i = 1:length(detectors)]; dir),
            (ang_tuple = Float64[]; angles)
                append!(i_list, ii)
                append!(j_list, jj)
                append!!(ang_tuple, angles)
                for (j, c) in enumerate(crx)
                    if !ismissing(c) && !isempty(c)
                        merge!(crx_pts[j], crx[j])
                        merge!(ray_dirs[j], dir[j])
                    end
                end
            end
        end
        results = sparse(i_list, j_list, trues(length(i_list)), n_sim, length(detectors))
        if optimize
            return (results, crx_pts, ray_dirs, mix_dist_θ, mix_dist_φ, ang_tuple)
        else
            return (results, crx_pts, ray_dirs)
        end
    end
end


"""
Runs the simulation with a hemispherical generating surface, outputting only the sparse matrix.
    """
function runhemisimlite(n_sim::Int, detectors::Vector{T},
    R::Real, ℓ::Real;
    exec=ThreadedEx(), center::Union{NTuple{3,Real},Nothing}=nothing) where {T<:LabObject{<:Real}}
    if center === nothing
        center = detectors[1].position |> Tuple
        optimize = true
    else
        optimize = false
    end
    if optimize
        # Find the relative direction between other detectors and the first detector
        rel_dirs = []
        for i in 2:length(detectors)
            push!(rel_dirs, _relativedir(detectors[1], detectors[i]))
        end
        mix_dist_θ = MixtureModel([Normal(θ > π / 2 ? π - θ : θ, Δθ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_θ = truncated(mix_dist_θ, 0, π / 2)
        mix_dist_φ = MixtureModel([Normal(φ, Δφ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_φ = truncated(mix_dist_φ, -π, π)
    else
        mix_dist_θ = nothing
        mix_dist_φ = nothing
    end
    # See https://juliafolds.github.io/FLoops.jl/dev/howto/avoid-box/#avoid-box for the "let" phrase.
    let optimize = optimize, mix_dist_θ = mix_dist_θ, mix_dist_φ = mix_dist_φ, center = center
        @floop exec for i = 1:n_sim
            # Private mutable variables
            @init begin
                crosses = SortedDict{Float64,SVector{3,Float64}}()
                hit_vec = falses(length(detectors))
                i_vec = zeros(length(detectors))
                j_vec = zeros(length(detectors))
                angles = zeros(2)
            end
            # Set the hit_vec to false to prepare for a new ray
            hit_vec .*= false
            # Generate a random ray
            if optimize
                θ_sample = rand(mix_dist_θ)
                φ_sample = rand(mix_dist_φ)
                angles = (θ_sample, φ_sample)
                θs = π - θ_sample
                φs = π + φ_sample
            else
                θs = nothing
                φs = nothing
            end
            ray = Ray(R, center, ℓ; θ=θs, φ=φs)
            # Loop through each dectector to fill the pre-allocated vectors
            for (j, d) in enumerate(detectors)
                hit = isthrough!(ray, d, crosses)
                if hit
                    i_vec[j] = i
                    j_vec[j] = j
                    hit_vec[j] = hit
                    empty!(crosses)
                end
            end
            # Hit indices to be filled
            ii = i_vec[hit_vec]
            jj = j_vec[hit_vec]
            # Reduce the accumulators
            @reduce() do (i_list = Float64[]; ii),
            (j_list = Float64[]; jj),
            (ang_tuple = Float64[]; angles)
                append!!(i_list, ii)
                append!!(j_list, jj)
                append!!(ang_tuple, angles)
            end
        end
        results = sparse(i_list, j_list, trues(length(i_list)), n_sim, length(detectors))
        if optimize
            return (results, mix_dist_θ, mix_dist_φ, ang_tuple)
        else
            return (results)
        end
    end
end
