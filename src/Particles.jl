### This part stores properties of the generated particles and their interactions.

import Base.@kwdef
using LinearAlgebra
using Random
using StaticArrays
using QuadGK
using Roots
using JLD2
using StatsBase

CACHE_DIR = joinpath(dirname(@__FILE__), "cache")
mkpath(CACHE_DIR)

"""
A mutable ray object (just a line) that could go through multiple LabObjects.
(θ, φ) denotes the unit direction vector, so (0, π) goes up.
Note: This is only a demo for the use of @kwdef macro.
    """
@kwdef mutable struct Ray{T<:Real}
    azimuth::T # the angle θ in spherical coordinates in radian
    polar::T # angle φ in spherical coordinates in radian
    start_position::MVector{3,Float64} = @MVector [0, 0, 0]
end

"""
Convenience constructor.
"""
function Ray(θ::T, φ::U; start_pos::NTuple{3,<:Real}=(0, 0, 0)) where {T,U}
    args = promote(θ, φ)
    return Ray(azimuth=args[1], polar=args[2], start_position=MVector{3,Float64}(start_pos))
end


"""
Returns a random number with pdf: Normalized(cos²(x)) ∈ [x_min, x_max). We use inverse sampling for speed.
!!!Note: sampling on a sphere requires taking into account of the Jacobian sin(θ)dθdφ = dΩ,
The constraint is within any solid angle dΩ, the random point counts are constant: so ∫f(θ,φ)dΩ = 1 => pdf = f(θ,φ)sin(θ).
"""
function _randcos2(x_min::Real=π / 2, x_max::Real=π)
    @assert 0 <= x_min <= x_max <= π
    # Note (-rand)^(1/3) doesn't work; it gives a complex root (not sure why Julia dev does this...)
    N = 1 / 3 * cos(x_min)^3 - 1 / 3 * cos(x_max)^3
    u = rand() * N
    res = (1 / 3 * cos(x_min)^3 - u)
    return acos(cbrt(3 * res))
end


"""
Returns the distribution f(θ) = cos³(θ)sin(θ). The extra cos factor is required for sampling on a horizontal plane with rays coming
non-vertically (for the angle between the flux vector and norm vector of the surface).
    """
function _randcos3(x_min::Real=π / 2, x_max::Real=π)
    @assert π / 2 <= x_min <= x_max <= π
    # For x_min smaller than π/2 cos(x)^4 will not be monotonic which causes trouble
    N = 1 / 4 * cos(x_min)^4 - 1 / 4 * cos(x_max)^4
    u = rand() * N
    res = (1 / 4 * cos(x_min)^4 - u)
    return acos(-(4 * res)^(1 / 4))
end


"""
Generates the uniform distribution on the sphere for variable θ. f(θ) = sin(θ). Note the range is [0, π).
"""
_randsin() = acos(1 - 2rand())


"""
Randomly modifies a ray with start uniform x, y ∈ (-max_x, max_x) U (-max_y, max_y), with max_bounds = (max_x, max_y)
The center is given as a tuple.
If not given, the θ and φ are drawn from cos³sin and uniform distributions respectively.
"""
function modifyray!(r::Ray{<:Real}, max_bounds::NTuple{2,Real}, center::NTuple{3,Real}; θ::Union{Nothing,Real}=nothing, φ::Union{Nothing,Real}=nothing)::Nothing
    x = rand() * 2max_bounds[1] - max_bounds[1] + center[1]
    y = rand() * 2max_bounds[2] - max_bounds[2] + center[2]
    if θ === nothing
        # Generate by the distribution
        θ = _randcos3()
    end
    if φ === nothing
        φ = rand() * 2π
    end
    r.azimuth = θ
    r.polar = φ
    (r.start_position.x, r.start_position.y, r.start_position.z) = (x, y, center[3])
    return nothing
end


"""
Randomly modifies a ray on the upper hemisphere with radius R. The direction is along the normal vector.
If not given, the θ and φ are drawn from cos²sin and uniform distributions respectively.
See http://arxiv.org/abs/1912.05462.
"""
function modifyray!(r::Ray{<:Real}, R::Real, center::NTuple{3,Real}, ℓ::Real;
    θ::Union{Nothing,Real}=nothing, φ::Union{Nothing,Real}=nothing)::Nothing
    if θ === nothing
        # Generate by the distribution
        θ = _randcos2()
    end
    if φ === nothing
        φ = rand() * 2π
    end
    start_θ = π - θ
    start_φ = φ + π
    x_vec = SA_F64[sin(start_θ)*cos(start_φ), sin(start_θ)*sin(start_φ), cos(start_θ)] * R
    pos = (@SVector(rand(2)) * ℓ) .- ℓ / 2
    θ_hat = SA_F64[cos(θ)*cos(φ), cos(θ)*sin(φ), -sin(θ)]
    φ_hat = SA_F64[-sin(φ), cos(φ), 0]
    dx_vec = pos[1] * θ_hat + pos[2] * φ_hat
    x_vec += (dx_vec + (center |> SVector))
    r.azimuth = θ
    r.polar = φ
    (r.start_position.x, r.start_position.y, r.start_position.z) = (x_vec.x, x_vec.y, x_vec.z)
    return nothing
end


"""
Returns the direction vector in Cartesian coordinates.
"""
function raydir(r::Ray{<:Real})::SVector{3,Float64}
    return _unitsph2cart(r.azimuth, r.polar)
end


"""
Returns the intersection point a ray and an inifinite plane. 
It's a tuple of intersection time and intersection position vector.
"""
function rayplaneint(dir::SVector{3,Float64}, r_pos::SVector{3,Float64}, plane_norm::SVector{3,Float64}, plane_pos::SVector{3,Float64})::Union{Tuple,Nothing}
    β = dir ⋅ plane_norm
    if β == 0
        return nothing
    else
        t = (plane_pos - r_pos) ⋅ plane_norm / β
        tup = (t, r_pos + t * dir)
        return tup
    end
end

# --- Scratch ---

# μ_dist_cache_fname = "mu_spawn_dist"
# μ_dist_cache_file = joinpath(CACHE_DIR, μ_dist_cache_fname * ".jld2")
# μ_dist_cache = Float64[]
# if !isfile(μ_dist_cache_file)
#     println("Generating cache...")
#     # save(μ_dist_cache_file, μ_dist_cache_fname, μ_dist_cache)
# else
#     println("Loading cache...")
#     μ_dist_cache = load(μ_dist_cache_file, μ_dist_cache_fname)
# end
# """
# Generates the muon energy given its ray angle.
# Note currently the energy distribution is completely independent of the angle distribution.
# """
# function μenergy!(θ::Real)
#     sample = unisampler!(x -> μpdf(x, θ=π - θ), low_bound=1.0, upp_bound=100, cache=μ_dist_cache)
#     save(μ_dist_cache_file, μ_dist_cache_fname, μ_dist_cache)
#     return sample
# end


"""
Muon flux pdf.
See https://www.worldscientific.com/doi/abs/10.1142/S0217751X18501750.
"""
function μpdf(E::Real; θ=0, E₀=4.29, ϵ=854, n=3.0, Rd_ratio=174)
    # Parameter units: E₀ [GeV], ϵ [GeV]
    D = √(Rd_ratio^2 * cos(θ)^2 + 2Rd_ratio + 1) - Rd_ratio * cos(θ)
    return (E₀ + E)^(-n) * (1 + E / ϵ)^(-1) * D
end


"""
A sampler for an arbitrary univariate distribution.
It takes in a function with possible lower and upper bounds and performs inverse sampling.
If provided, it also saves the results into the cache vector.
If the total number of provided cache exceeds (by default) 1e5 the function will draw from the saved cache.
"""
function unisampler!(f::Function; low_bound=-Inf, upp_bound=Inf, cache=nothing, cache_size=Int(1e5))
    if cache !== nothing && length(cache) > cache_size
        return rand(cache)
    end
    (total, err) = quadgk(f, low_bound, upp_bound)
    x = rand()
    eqn = y -> unicdf(f, y, low_bound=low_bound, upp_bound=upp_bound, total=total) - x
    res = find_zero(eqn, (low_bound, upp_bound), A42())
    if cache !== nothing
        push!(cache, res)
    end
    return res
end


"""
A sampler for an arbitrary univariate distribution using its cdf.
NOTE: THIS METHOD IS NOT STABLE. 
"""
function unisamplercdf!(ecdf; low_bound=-Inf, upp_bound=Inf, cache=nothing, cache_size=Int(1e5))
    if cache !== nothing && length(cache) > cache_size
        return rand(cache)
    end
    x = rand()
    eqn = y -> ecdf(y) - x
    res = find_zero(eqn, (low_bound, upp_bound), Bisection())
    if cache !== nothing
        push!(cache, res)
    end
    return res
end


"""
Return the cdf of a univariate distribution. Note this is not normalized.
"""
function unicdf(f::Function, x::Real; low_bound=-Inf, upp_bound=Inf, total=1.0)
    @assert low_bound <= x <= upp_bound "$x exceeding function support [$low_bound, $upp_bound]."
    (res, err) = quadgk(f, low_bound, x)
    res /= total
    return res
end
