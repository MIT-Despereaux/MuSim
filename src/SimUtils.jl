### This part contains functions used for simulations

using FLoops
using SparseArrays
using DataStructures


"""
Runs simulation given a list of detectors, and spawns the rays with 2max_x × 2max_x area @ center.
Note that the detectors are a vector of types being the subtypes of the LabObject (general method).
The sparse boolean array is stored with dimension (n_sim × n_detector).
The information for each hit event is stored in dictionary for each detector.
    """
function runsim(n_sim::Int, detectors::Vector{T},
    max_x::Real, center::NTuple{3,Real}; exec = ThreadedEx()) where {T<:LabObject{<:Real}}
    @floop exec for i = 1:n_sim
        # Private mutable variables
        @init begin
            crosses = SortedDict{Float64,SVector{3,Float64}}()
            hit_vec = falses(length(detectors))
            i_vec = zeros(length(detectors))
            j_vec = zeros(length(detectors))
            crx = Vector{Union{Missing,Dict}}(missing, length(detectors))
            dir = Vector{Union{Missing,Dict}}(missing, length(detectors))
        end
        # Set the hit_vec to false to prepare for a new ray
        hit_vec .*= false
        ray = Ray((max_x, max_x), center)
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
        @reduce() do (I = Int[]; ii),
        (J = Int[]; jj),
        (crx_pts = [Dict{Int,Matrix{Float64}}() for i = 1:length(detectors)]; crx),
        (ray_dirs = [Dict{Int,Vector{Float64}}() for i = 1:length(detectors)]; dir)
            append!(I, ii)
            append!(J, jj)
            for (j, c) in enumerate(crx)
                if !ismissing(c) && !isempty(c)
                    merge!(crx_pts[j], crx[j])
                    merge!(ray_dirs[j], dir[j])
                end
            end
        end
    end
    results = sparse(I, J, trues(length(I)), n_sim, length(detectors))
    return results, crx_pts, ray_dirs
end


"""
Runs the simulation a hemispherical generating surface. 
    """
function runhemisim(n_sim::Int, detectors::Vector{T},
    R::Real, center::NTuple{3,Real}, ℓ::Real;
    exec = ThreadedEx(), θ_range = (π / 2, π), φ_range = (0, 2π)) where {T<:LabObject{<:Real}}
    @floop exec for i = 1:n_sim
        # Private mutable variables
        @init begin
            crosses = SortedDict{Float64,SVector{3,Float64}}()
            hit_vec = falses(length(detectors))
            i_vec = zeros(length(detectors))
            j_vec = zeros(length(detectors))
            crx = Vector{Union{Missing,Dict}}(missing, length(detectors))
            dir = Vector{Union{Missing,Dict}}(missing, length(detectors))
        end
        # Set the hit_vec to false to prepare for a new ray
        hit_vec .*= false
        ray = Ray(R, center, ℓ; θ = θ_range, φ = φ_range)
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
        @reduce() do (I = Int[]; ii),
        (J = Int[]; jj),
        (crx_pts = [Dict{Int,Matrix{Float64}}() for i = 1:length(detectors)]; crx),
        (ray_dirs = [Dict{Int,Vector{Float64}}() for i = 1:length(detectors)]; dir)
            append!(I, ii)
            append!(J, jj)
            for (j, c) in enumerate(crx)
                if !ismissing(c) && !isempty(c)
                    merge!(crx_pts[j], crx[j])
                    merge!(ray_dirs[j], dir[j])
                end
            end
        end
    end
    results = sparse(I, J, trues(length(I)), n_sim, length(detectors))
    return results, crx_pts, ray_dirs
end


"""
Runs the simulation a hemispherical generating surface, outputting only the sparse matrix.
    """
function runhemisimlite(n_sim::Int, detectors::Vector{T},
    R::Real, center::NTuple{3,Real}, ℓ::Real;
    exec = ThreadedEx(), θ_range = (π / 2, π), φ_range = (0, 2π)) where {T<:LabObject{<:Real}}
    @floop exec for i = 1:n_sim
        # Private mutable variables
        @init begin
            crosses = SortedDict{Float64,SVector{3,Float64}}()
            hit_vec = falses(length(detectors))
            i_vec = zeros(length(detectors))
            j_vec = zeros(length(detectors))
        end
        # Set the hit_vec to false to prepare for a new ray
        hit_vec .*= false
        ray = Ray(R, center, ℓ; θ = θ_range, φ = φ_range)
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
        @reduce() do (I = Int[]; ii),
        (J = Int[]; jj)
            append!(I, ii)
            append!(J, jj)
        end
    end
    results = sparse(I, J, trues(length(I)), n_sim, length(detectors))
    return results
end
