### This script stores functions for simulations
### DEPRECATED due to failed attempt to decrease GC time

using Random
using Printf
using Folds, FLoops, BangBang
using Combinatorics
import Base.@kwdef

# Abstract lab object that can be used in coincidence
# It is NOT possible to subtype a concrete type
abstract type LabObject{T <: Real} end

# A "rectangle" box object with 8 corners
struct RectBox{T} <: LabObject{T}
    # The box is assumed to be at a position, with ±(Δx | Δy |Δz) / 2 coverage
    delta_x::T # x-direction
    delta_y::T # y-direction
    delta_z::T # z-direction
    efficiency::T # detector efficiency
    position::NTuple{3,Real} # position
    material::String # material
end

# This is a constructor that promotes type T to U if efficiency is not provided
# It also sets the default values
function RectBox(delta_x::T, delta_y::T, delta_z::T; efficiency::U=1.0, position=(0, 0, 0), material="Unknown") where {T,U}
    args = promote(delta_x, delta_y, delta_z, efficiency)
    return RectBox(args..., position, material)
end

# A ray object (just a line) that could go through multiple LabObjects
# Since we need to simulate a large number of rays, this is going to be mutable
@kwdef mutable struct Ray
    # (θ, φ) denotes the unit velocity vector, so (0, π) goes up
    azimuth_agl::Real # the angle θ in spherical coordinates in radian
    polar_agl::Real # angle φ in spherical coordinates in radian
    # coord z defaults to 100 m
    start_position::NTuple{3,Real} = (0, 0, 100)
end

# Print information about the RectBox
function Base.show(io::IO, box::RectBox)
    @printf(io, "\n--- Detector material: %s, with efficiency %.2f", box.material, box.efficiency)
    @printf(io, "\n--- Detector dimensions: %.2fx%.2fx%.2f cm^3", box.delta_x * 100, box.delta_y * 100, box.delta_z * 100)
    @printf(io, "\n--- Detector position: (%.0fcm, %.0fcm, %.0fcm)\n", (100 .* box.position)...)
end

"""
Returns a random number with pdf: 0 for x ∈ [0, π/2) and Normalized(cos²(x)) ∈ [π/2, π), using rejection sampling (and recursive...).
!!!Note: sampling on a sphere requires taking into account of the Jacobian sin(θ)dθdφ = dΩ
The constraint is within any solid angle dΩ, the random point counts are constant: so ∫f(θ,φ)dΩ = 1 => pdf = f(θ,φ)sin(θ)
"""
function rand_cos2()
    pts = rand(2)
    pts[1] *= pi / 2
    pts[1] += pi / 2
    return pts[2] <= sin(pts[1]) * cos(pts[1])^2 ? pts[1] : rand_cos2()
end

function rand_cos2!(v::Vector)
    for i in 1:length(v)
        v[i] = rand_cos2()
    end
    return v
end

# For optional constructors, it's best to keep the required fields before keywords, and use keywords for optinal fields
function Ray(max_x::Real, max_y::Real, z_pos::Real; azimuth_agl=nothing, polar_agl=nothing)
    # Randomly generate a ray with start uniform x, y ∈ (-max_x, max_x) U (-max_y, max_y)
    # If not given, the θ and φ are drawn from cos^2 and uniform distributions respectively
    x = rand() * 2max_x - max_x
    y = rand() * 2max_y - max_y
    if azimuth_agl === nothing
        # Generate by the distribution
        azimuth_agl = rand_cos2()
        end
    if polar_agl === nothing
    polar_agl = rand() * 2pi
    end
    return Ray(azimuth_agl=azimuth_agl, polar_agl=polar_agl, start_position=(x, y, z_pos))
end

# Sets the start position of the ray in place
function set_ray_startpos!(ray::Ray, position::NTuple{3,Real})
    ray.start_position = position
end

function set_ray_dir!(ray::Ray, θ::Real, φ::Real)
    ray.azimuth_agl = θ
    ray.polar_agl = φ
end

# Returns the direction vector in Cartesian coordinates
function ray_dir(r::Ray)::Vector{Real}
    θ = r.azimuth_agl
    φ = r.polar_agl
    return [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end

@enum Ori xaxis = 1 yaxis = 2 zaxis = 3
"""
Returns the intersection point a ray and an inifinite plane with only x, y, and z norm vector orientation for Axis Aligned Bounding Box (AABB).
(There is no difference in speed from used one).
"""
function _ray_plane_int_aabb(r::Ray, plane_ori::Ori, plane_pos::SVector{3,Float64})::Union{SVector{3,Float64},Nothing}
    dir = _ray_dir(r)
    xr0 = r.start_position[Int(plane_ori)]
    xp0 = plane_pos[Int(plane_ori)]
    β = dir[Int(plane_ori)]
    t = (xp0 - xr0) / β
    return β == 0 ? nothing : r.start_position + t * dir
end

# Tests if the ray goes through the box
# The ray can be infinitely long
function isthrough(r::Ray, box::RectBox)::Bool
    dir = ray_dir(r)
    # ray_endcoord = dir * raylength + collect(r.start_position)
    # Calculate the six end points of the box
    # z0 = z_pos - h / 2
    z0 = box.position[3] - box.delta_z / 2
    zc_t0 = (z0 - r.start_position[3]) / dir[3]
    zc_x0 = zc_t0 * dir[1] + r.start_position[1]
    zc_y0 = zc_t0 * dir[2] + r.start_position[2]
    z1 = box.position[3] + box.delta_z / 2
    zc_t1 = (z1 - r.start_position[3]) / dir[3]
    zc_x1 = zc_t1 * dir[1] + r.start_position[1]
    zc_y1 = zc_t1 * dir[2] + r.start_position[2]
    # y0
    y0 = box.position[2] - box.delta_y / 2
    yc_t0 = (y0 - r.start_position[2]) / dir[2]
    yc_z0 = yc_t0 * dir[3] + r.start_position[3]
    yc_x0 = yc_t0 * dir[1] + r.start_position[1]
    y1 = box.position[2] + box.delta_y / 2
    yc_t1 = (y1 - r.start_position[2]) / dir[2]
    yc_z1 = yc_t1 * dir[3] + r.start_position[3]
    yc_x1 = yc_t1 * dir[1] + r.start_position[1]
    # x0
    x0 = box.position[1] - box.delta_x / 2
    xc_t0 = (x0 - r.start_position[1]) / dir[1]
    xc_y0 = xc_t0 * dir[2] + r.start_position[2]
    xc_z0 = xc_t0 * dir[3] + r.start_position[3]
    x1 = box.position[1] + box.delta_x / 2
    xc_t1 = (x1 - r.start_position[1]) / dir[1]
    xc_y1 = xc_t1 * dir[2] + r.start_position[2]
    xc_z1 = xc_t1 * dir[3] + r.start_position[3]

    # If any of the six points are within the bounded square, then return True
    if ((x0 < zc_x0 < x1) && (y0 < zc_y0 < y1)) || 
        ((x0 < zc_x1 < x1) && (y0 < zc_y1 < y1)) ||
        ((x0 < yc_x0 < x1) && (z0 < yc_z0 < z1)) ||
        ((x0 < yc_x1 < x1) && (z0 < yc_z1 < z1)) ||
        ((y0 < xc_y0 < y1) && (z0 < xc_z0 < z1)) ||
        ((y0 < xc_y1 < y1) && (z0 < xc_z1 < z1))
        return rand() < box.efficiency
    else
        return false
    end
end

# Tests if the ray goes through the box and stores the enter and exit locations
function isthrough!(r::Ray, box::RectBox, cross1::NTuple{3,Real}, cross2::NTuple{3,Real})::Bool
    dir = ray_dir(r)
    enter = false
    # End points of the box
    z0 = box.position[3] - box.delta_z / 2
    z1 = box.position[3] + box.delta_z / 2
    y0 = box.position[2] - box.delta_y / 2
    y1 = box.position[2] + box.delta_y / 2
    x0 = box.position[1] - box.delta_x / 2
    x1 = box.position[1] + box.delta_x / 2
    # z-planes
    zc_t0 = (z0 - r.start_position[3]) / dir[3]
    zc_x0 = zc_t0 * dir[1] + r.start_position[1]
    zc_y0 = zc_t0 * dir[2] + r.start_position[2]
    zc_t1 = (z1 - r.start_position[3]) / dir[3]
    zc_x1 = zc_t1 * dir[1] + r.start_position[1]
    zc_y1 = zc_t1 * dir[2] + r.start_position[2]
    if ((x0 < zc_x0 < x1) && (y0 < zc_y0 < y1))
        enter = true
        cross1 = (zc_x0, zc_y0, z0)
    end
    if ((x0 < zc_x1 < x1) && (y0 < zc_y1 < y1))
        if enter
            cross2 = (zc_x1, zc_y1, z1)
            return true
        else
            enter = true
            cross1 = (zc_x1, zc_y1, z1)
        end
    end
    # y-planes
    yc_t0 = (y0 - r.start_position[2]) / dir[2]
    yc_z0 = yc_t0 * dir[3] + r.start_position[3]
    yc_x0 = yc_t0 * dir[1] + r.start_position[1]
    yc_t1 = (y1 - r.start_position[2]) / dir[2]
    yc_z1 = yc_t1 * dir[3] + r.start_position[3]
    yc_x1 = yc_t1 * dir[1] + r.start_position[1]
    if ((x0 < yc_x0 < x1) && (z0 < yc_z0 < z1))
        if enter
            cross2 = (yc_x0, y0, yc_z0)
            return true
        else
            enter = true
            cross1 = (yc_x0, y0, yc_z0)
        end
    end
    if ((x0 < yc_x1 < x1) && (z0 < yc_z1 < z1))
        if enter
            cross2 = (yc_x1, y1, yc_z1)
            return true
        else
            enter = true
            cross1 = (yc_x1, y1, yc_z1)
        end
    end
    # x-planes
    xc_t0 = (x0 - r.start_position[1]) / dir[1]
    xc_y0 = xc_t0 * dir[2] + r.start_position[2]
    xc_z0 = xc_t0 * dir[3] + r.start_position[3]
    xc_t1 = (x1 - r.start_position[1]) / dir[1]
    xc_y1 = xc_t1 * dir[2] + r.start_position[2]
    xc_z1 = xc_t1 * dir[3] + r.start_position[3]
    if ((y0 < xc_y0 < y1) && (z0 < xc_z0 < z1))
        if enter
            cross2 = (x0, xc_y0, xc_z0)
            return true
        else
            enter = true
            cross1 = (x0, xc_y0, xc_z0)
        end
    end
    if ((y0 < xc_y1 < y1) && (z0 < xc_z1 < z1))
        cross2 = (x1, xc_y1, xc_z1)
return true
    end
    return false
end

# Functions that modifies arguments in place should end with "!"
# Runs simulation given a list of detectors, and spawns the rays with 2x_max × 2x_max area @ z_pos
# Returns the populated array of results
# Note that the detectors are a vector of types being the subtypes of the LabObject (general method)
function _runsim!(detectors::Vector{T}, results::Array{Bool,2}, max_x::Real, z_pos::Real) where {T <: LabObject{<:Real}}
    ray = Ray(max_x, max_x, z_pos)
    for i in 1:size(results, 1)
        # Reseed the ray in place
        θ = rand_cos2()
        φ = rand() * 2pi
        (x, y) = rand(2) .* 2max_x .- max_x
        set_ray_startpos!(ray, (x, y, z_pos))
        set_ray_dir!(ray, θ, φ)
        for (det, res) in zip(detectors, eachcol(results))
        res[i] = isthrough(ray, det)
        end
    end
end

# This function is intended to store the crossing points
function _runsimx!(i::Int, detectors::Vector{T}, results::Array{Bool,2}, x_max::Real, z_pos::Real) where {T <: LabObject{<:Real}}
    ray = Ray(x_max, x_max, z_pos)
    # For now, do not store the crossing points.
    # cross1::NTuple{3,Real} = (0.0, 0.0, 0.0)
    # cross2::NTuple{3,Real} = (0.0, 0.0, 0.0)
    for i in 1:size(results, 1)
        # Reseed the ray in place
        θ = rand_cos2()
        φ = rand() * 2pi
        (x, y) = rand(2) .* 2max_x .- max_x
        set_ray_startpos!(ray, (x, y, z_pos))
        set_ray_dir!(ray, θ, φ)
        for (det, res) in zip(detectors, eachcol(results))
        res[i] = isthrough(ray, det)
        end
    end
end

# A special function for rectangular box detectors
# Looks like this is not necessary (doesn't speed up execution)
# function _runsim!(i::Int, detectors::Vector{RectBox{Float64}}, results::Array{Bool,2}, x_max::Real, z_pos::Real)
#     ray = Ray(x_max, x_max, z_pos)
#     # For now, do not store the crossing points.
#     cross1::NTuple{3,Real} = (0.0, 0.0, 0.0)
#     cross2::NTuple{3,Real} = (0.0, 0.0, 0.0)
#     for (det, res) in zip(detectors, eachcol(results))
#         res[i] = isthrough!(ray, det, cross1, cross2)
#     end
# end


function coincidence_main()
    # number of rays to be simulated
    sim_num = 1e6 # Actual simulation number depends on the number of threads
    # "Batches" are what being parallelized
    batch_num = Threads.nthreads()
    batch_sim = round(Int, sim_num / batch_num)

    # size of the square plane from which rays are generated [m]
    x_bound = 0.25 / 2
    
    # Set up the detectors
    α = 45
    y_pos1 = 0.2
    z_pos1 = cotd(α) * y_pos1
    y_pos2 = -0.4
    z_pos2 = cotd(α) * y_pos2

    box1 = RectBox(0.05, 0.05, 0.0005, position=(0.0, 0.0, 0.9), efficiency=1.0, material="BC408")
    box2 = RectBox(0.05, 0.05, 0.0005, position=(0.0, 0.0, 0.0), efficiency=1.0, material="BC408")
    # chip = RectBox(0.005, 350e-6, 0.005, position=(0.0, 0.0, 0.0), efficiency=1.0, material="Si")
    
    total_rate = 120 * (2x_bound)^2 # 1 second normalisation
    duration = sim_num / total_rate

    detectors = [box1, box2]
    z_pos_max = maximum([d.position[3] + d.delta_z / 2 for d in detectors]) # generate from plane above the highest detectors
    @printf("Number of rays: %.1e, duration: %.1fhours, threads: %d\n", sim_num, duration / 3600, Threads.nthreads())
    
    # WIP: The following section needs to parallelize on "batch" number
    # if Threads.nthreads() == 1
    #     @time collect(_runsim!(i, detectors, results, x_bound, z_pos_max) for i in 1:sim_num)
    # else
    #     @time Folds.collect(_runsim!(i, detectors, results, x_bound, z_pos_max) for i in 1:sim_num)
    # end
    @time begin
        # This is the thread-safe version of batch-simulation
        # Unfortunately it doesn't decrease the GC time...
        @floop ThreadedEx() for i in 1:batch_num
            results = Matrix{Bool}(undef, batch_sim, length(detectors))
            _runsim!(detectors, results, x_bound, z_pos_max)
            @reduce() do (total_res = Matrix{Bool}[]; results)
                push!!(total_res, results)
            end
        end
        total_res = vcat(total_res...)
    end

    # Print single detector rates:
    for (det, res) in zip(detectors, eachcol(total_res))
        print(det)
        counts = sum(res)
        relerr = 1 / √counts
        if counts >= 1
            @printf("\tEvent rate: 1 per %.3f min (1/%.1fs, %.2gHz), uncertainty: %.2g%%.\n", 
                                                                 duration / counts / 60,
                                                                 duration / counts,
                                                                 counts / duration, relerr * 100)
        else
            println("\tNo hits detected.")
        end
    end
    # # Print double coincidences
    # for c in combinations(zip(detectors, eachcol(results)), 2)
        #     double_coin = trues(sim_num)
    #     for (res) in c
    #         double_coin = double_coin .& res
    #     double_coin = sum(double_coin)
    #     @printf("\ndoublecoincidences: 1 per %.3f min (1/%.1fs, %.2gHz), uncertainty: %.2g%%.\n", 
    #                                                                 duration / double_coin / 60,
    #                                                                 duration / double_coin,
    #                                                                 double_coin / duration, 100 / √double_coin)
    #     println("Nr doublecoincidences: $(double_coin).")
    #     end
    # end
end