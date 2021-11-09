### This script stores properties of the generated particles and their interactions.

import Base.@kwdef
using LinearAlgebra
using Random
using StaticArrays
using QuadGK
using Roots

"""
A mutable ray object (just a line) that could go through multiple LabObjects.
(θ, φ) denotes the unit direction vector, so (0, π) goes up.
    """
@kwdef struct Ray
    azimuth_agl::Real # the angle θ in spherical coordinates in radian
    polar_agl::Real # angle φ in spherical coordinates in radian
    # coord z defaults to 100 m
    start_position::SVector{3,Float64} = SA_F64[0, 0, 100]
end


"""
Returns a random number with pdf: 0 for x ∈ [0, π/2) and Normalized(cos²(x)) ∈ [π/2, π). We use inverse sampling for speed.
!!!Note: sampling on a sphere requires taking into account of the Jacobian sin(θ)dθdφ = dΩ,
The constraint is within any solid angle dΩ, the random point counts are constant: so ∫f(θ,φ)dΩ = 1 => pdf = f(θ,φ)sin(θ).
"""
function _randcos2()
    # Note (-rand)^(1/3) doesn't work; it gives a complex root (not sure why Julia dev does this...)
    u = cbrt(-rand())
    return acos(u)
end


"""
Returns the distribution f(θ) = cos³(θ)sin(θ). The extra cos factor is required for sampling on a horizontal plane with rays coming
non-vertically (for the angle between the flux vector and norm vector of the surface).
    """
function _randcos3()
    u = -rand()^(1 / 4)
    return acos(u)
end


"""
Generates the uniform distribution on the sphere for variable θ. f(θ) = sin(θ). Note the range is [0, π).
"""
_randsin() = acos(1 - 2rand())

"""
Returns a vector with the specified distribution function
"""
function _rand_vec!(dist::Function, v::Vector)
    for i in 1:length(v)
        v[i] = dist()
    end
end


"""
Generates the muon energy given its ray angle.
Note currently the energy distribution is completely independent of the angle distribution.
"""
function μenergy!(θ::Real; cache=nothing)
    return unisampler!(x -> μpdf(x, θ=π-θ), low_bound=1.0, upp_bound=100, cache=cache)
end


"""
Muon flux pdf.
See https://www.worldscientific.com/doi/abs/10.1142/S0217751X18501750.
"""
function μpdf(E::Real; θ=0, E₀=4.29, ϵ=854, n=3.0, Rd_ratio=174)
    # Parameter units: E₀ [GeV], ϵ [GeV]
    D = √(Rd_ratio^2 * cos(θ)^2 + 2Rd_ratio + 1) - Rd_ratio * cos(θ)
    return (E₀ + E)^(-n) * (1 + E/ϵ)^(-1) * D
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
Return the cdf of a univariate distribution. Note this is not normalized.
"""
function unicdf(f::Function, x::Real; low_bound=-Inf, upp_bound=Inf, total=1.0)
    @assert low_bound <= x <= upp_bound "$x exceeding function support [$low_bound, $upp_bound]."
    (res, err) = quadgk(f, low_bound, x)
    res /= total
    return res
end


# For optional constructors, it's best to keep the required fields before keywords, and use keywords for optinal fields
"""
Randomly generate a ray with start uniform x, y ∈ (-max_x, max_x) U (-max_y, max_y), with max_bounds = (max_x, max_y)
The center is given as a tuple.
If not given, the θ and φ are drawn from cos³sin and uniform distributions respectively.
"""
function Ray(max_bounds::NTuple{2,Real}, center::NTuple{3,Real}; θ=nothing, φ=nothing)
    x = rand() * 2max_bounds[1] - max_bounds[1] + center[1]
    y = rand() * 2max_bounds[2] - max_bounds[2] + center[2]
    if θ === nothing
        # Generate by the distribution
        θ = _randcos3()
    end
        if φ === nothing
        φ = rand() * 2pi
    end
    return Ray(azimuth_agl=θ, polar_agl=φ, start_position=SA_F64[x, y, center[3]])
end


"""
Randomly generate a ray on the upper hemisphere with radius R. The direction is along the normal vector.
See http://arxiv.org/abs/1912.05462.
"""
function Ray(R::Real, center::NTuple{3,Real}, side_len::Real; θ=nothing, φ=nothing)
    if θ === nothing
        # Generate by the distribution
        θ = _randcos2()
    end
        if φ === nothing
        φ = rand() * 2pi
    end
    start_θ = pi - θ
    start_φ = φ + pi
    (x, y, z) = SA_F64[sin(start_θ) * cos(start_φ), sin(start_θ) * sin(start_φ), cos(start_θ)] * R
    (a, b) = rand(2) * side_len .- side_len / 2
    θ_hat = SA_F64[cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
    φ_hat = SA_F64[-sin(φ), cos(φ), 0]
    (dx, dy, dz) = a * θ_hat + b * φ_hat
    x += center[1] + dx
    y += center[2] + dy
    z += center[3] + dz
    return Ray(azimuth_agl=θ, polar_agl=φ, start_position=SA_F64[x, y, z])
end


"""
Returns the direction vector in Cartesian coordinates.
"""
function raydir(r::Ray)::SVector{3,Float64}
    θ = r.azimuth_agl
    φ = r.polar_agl
    return SA_F64[sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
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
