### This part stores functions related to ray tracing

import Base: +, -, ==, hash, show
using Random
using LinearAlgebra
using Printf
using StaticArrays
using DataStructures
using HCubature
using MCIntegration

# Abstract lab object that can be used in coincidence
# It is NOT possible to subtype a concrete type
abstract type LabObject{T<:Real} end


"""
A "rectangle" box object with 8 corners.
The box is assumed to be at a position, with ±(Δx | Δy | Δz) / 2 coverage.
"""
mutable struct RectBox{T} <: LabObject{T}
    name::String # name
    delta_x::T # x-direction
    delta_y::T # y-direction
    delta_z::T # z-direction
    efficiency::T # detector efficiency
    position::SVector{3,Float64} # position
    orientation::SVector{2,Float64} # orientation in θ and φ
    material::String # material
end


"""
This is a constructor that promotes type T to U if efficiency is not provided.
It also sets the default values.
"""
function RectBox(name::String, delta_x::T, delta_y::T, delta_z::T;
    efficiency::U=1.0, position=(0, 0, 0),
    orientation=(0, 0), material="Unknown") where {T,U}
    args = promote(delta_x, delta_y, delta_z, efficiency)
    return RectBox(name, args..., SA_F64[position...], SA_F64[orientation...], material)
end


"""
A "Sphere" object.
"""
mutable struct Sphere{T} <: LabObject{T}
    name::String
    radius::T
    efficiency::T
    position::SVector{3,Float64}
    material::String
end


function Sphere(name::String, radius::T;
    efficiency::U=1.0, position=(0, 0, 0),
    material="Unknown") where {T,U}
    args = promote(radius, efficiency)
    return Sphere(name, args..., SA_F64[position...], material)
end


"""
A cylindrical object with two caps at position ± height/2.
"""
mutable struct Cylinder{T} <: LabObject{T}
    name::String
    radius::T
    height::T
    efficiency::T
    position::SVector{3,Float64}
    orientation::SVector{2,Float64}
    material::String
end


function Cylinder(name::String, radius::T, height::T;
    efficiency::U=1.0, position=(0, 0, 0),
    orientation=(0, 0), material="Unknown") where {T,U}
    args = promote(radius, height, efficiency)
    return Cylinder(name, args..., SA_F64[position...], SA_F64[orientation...], material)
end


"""
Combined objects that could take on addition or subtraction of different shapes.
Used for exotic detector shapes.
NOTE: WORK IN PROGRSS...
"""
mutable struct CombinedObj{T} <: LabObject{T}
    dets::Vector{LabObject{T}}
    ops::Vector{Symbol} # The list of operations to be performed between detectors. Only "+" and "-" are supported. 
    function CombinedObj{T}(dets::Vector{LabObject{T}}, ops::Vector{Symbol}) where {T}
        @assert length(ops) == length(dets) - 1 "Number of operators not strictly less than number of detectors by 1"
        new(dets, ops)
    end
end


function (+)(x::LabObject{T}, y::LabObject{T}) where {T}
    return CombinedObj{T}(LabObject{T}[x, y], Symbol[:+])
end


function (+)(x::CombinedObj{T}, y::LabObject{T}) where {T}
    push!(x.dets, y)
    push!(x.ops, :+)
    return CombinedObj{T}(x.dets, x.ops)
end


function (-)(x::LabObject{T}, y::LabObject{T}) where {T}
    return CombinedObj{T}(LabObject{T}[x, y], Symbol[:-])
end


function (-)(x::CombinedObj{T}, y::LabObject{T}) where {T}
    push!(x.dets, y)
    push!(x.ops, :-)
    return CombinedObj{T}(x.dets, x.ops)
end

"""
For comparison, we need to implement Base.hash and Base.:==
This is also used when hashing a dictionary that contains the object.
"""
function hash(o::T, h::UInt) where {T<:LabObject{<:Real}}
    all_prop_names = propertynames(o)
    out_hash = h
    for p in all_prop_names
        out_hash = hash(getproperty(o, p), out_hash)
    end
    return out_hash
end


function hash(o::T) where {T<:LabObject{<:Real}}
    return hash(o, zero(UInt))
end


"""
Generic comparison of different LabObjects.
"""
function ==(x::T, y::T)::Bool where {T<:LabObject{<:Real}}
    return hash(x) == hash(y)
end


# Print information about the LabObject.
function show(io::IO, box::RectBox)
    @printf(io, "name: %s\n", box.name)
    @printf(io, "material: %s, with efficiency %.2f\n", box.material, box.efficiency)
    @printf(io, "dimensions: %.2fx%.2fx%.2f cm^3\n", box.delta_x * 100, box.delta_y * 100, box.delta_z * 100)
    @printf(io, "position: (%.0fcm, %.0fcm, %.0fcm)\n", (100 .* box.position)...)
    @printf(io, "orientation: (θ: %.2f°, φ: %.2f°)\n", (rad2deg.(box.orientation))...)

end


function show(io::IO, s::Sphere)
    @printf(io, "name: %s\n", s.name)
    @printf(io, "material: %s, with efficiency %.2f\n", s.material, s.efficiency)
    @printf(io, "radius: %.2f cm\n", s.radius * 100)
    @printf(io, "position: (%.0fcm, %.0fcm, %.0fcm)\n", (100 .* s.position)...)
end


function show(io::IO, cyl::Cylinder)
    @printf(io, "name: %s\n", cyl.name)
    @printf(io, "material: %s, with efficiency %.2f\n", cyl.material, cyl.efficiency)
    @printf(io, "radius: %.2f cm\n", cyl.radius * 100)
    @printf(io, "height: %.2f cm\n", cyl.height * 100)
    @printf(io, "position: (%.0fcm, %.0fcm, %.0fcm)\n", (100 .* cyl.position)...)
    @printf(io, "orientation: (θ: %.2f°, φ: %.2f°)\n", (rad2deg.(cyl.orientation))...)
end


"""
Returns object orientation in Cartesian coordinates.
"""
function objorient(o::LabObject{T})::SVector{3,Float64} where {T}
    return _unitsph2cart(o.orientation...)
end


function _unitsph2cart(θ::Real, φ::Real)::SVector{3,Float64}
    return SA_F64[sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ)]
end


function _cart2unitsph(x::Real, y::Real, z::Real)::SVector{2,Float64}
    len = √(x^2 + y^2 + z^2)
    return SA_F64[acos(z / len), atan(y, x)]
end


function _translate(v::Union{SVector{3,Float64},MVector{3,Float64}}, Δx::SVector{3,Float64})::SVector{3,Float64}
    return v + Δx
end


"""
Actively rotate any vector around the y+ -> z+ -> y+ axes for three given angles.
"""
function _rotate(v::SVector{3,Float64}, θ₁::Real, φ::Real, θ₂::Real)
    # Note this is an active rotation around the y+ axis
    if θ₁ == 0
        mat_θ1 = 1I
    else
        mat_θ1 = SA_F64[cos(θ₁) 0 sin(θ₁)
            0 1 0
            -sin(θ₁) 0 cos(θ₁)]
    end
    # Note this is an active rotation around the z+ axis
    if φ == 0
        mat_φ = 1I
    else
        mat_φ = SA_F64[cos(φ) -sin(φ) 0
            sin(φ) cos(φ) 0
            0 0 1]
    end
    if θ₂ == 0
        mat_θ2 = 1I
    else
        mat_θ2 = SA_F64[cos(θ₂) 0 sin(θ₂)
            0 1 0
            -sin(θ₂) 0 cos(θ₂)]
    end
    new_v = mat_θ2 * mat_φ * mat_θ1 * v
    return new_v
end

"""
Solid angle factor for a single side of a box.
"""
function _fΩ(θ)
    dΩ = x -> abs(cos(x[1]) * cos(θ) + cos(x[2]) * sin(x[1]) * sin(θ)) * cos(x[1])^2 * sin(x[1])
    res, _ = hcubature(dΩ, (0, 0), (π / 2, 2π), rtol=1e-6)
    return res
end

"""
Calculates the analytic rate (assuming unit I₀) of a box, given the parameters.
The result needs to be divided by (2π / 3) to be consistent with normalization.
Note this rate is inclusive.
"""
function analytic_R(box::RectBox)::Real
    # Find the orientations of three sides
    θ₁ = box.orientation[1]
    φ = box.orientation[2]
    x̂ = _rotate(SA_F64[1, 0, 0], θ₁, φ, 0)
    ŷ = _rotate(SA_F64[0, 1, 0], θ₁, φ, 0)
    ẑ = _rotate(SA_F64[0, 0, 1], θ₁, φ, 0)
    R = 0
    # Calculate the flux for each side
    θx = _cart2unitsph(x̂...)[1]
    θy = _cart2unitsph(ŷ...)[1]
    θz = _cart2unitsph(ẑ...)[1]
    R += _fΩ(θx) * (box.delta_y * box.delta_z)
    R += _fΩ(θy) * (box.delta_x * box.delta_z)
    R += _fΩ(θz) * (box.delta_x * box.delta_y)
    return R * box.efficiency
end


"""
Tests if the ray goes through the box, assuming the ray can be infinitely long.
"""
function isthrough!(r::Ray, box::RectBox, crosses::SortedDict{Float64,SVector{3,Float64}})::Bool
    vhit = rand() < box.efficiency
    if vhit
        len = length(crosses)
        # Translation
        Δx = box.position
        # Orientation
        (θ, φ) = box.orientation
        # Ray direction
        dir = raydir(r)
        # Ray start position
        r_pos = r.start_position

        # Translate back to origin
        r_pos_trans = _translate(r_pos, -Δx)
        # Rotate back to normal
        dir_rot = _rotate(dir, 0, -φ, -θ)
        r_pos_trans_rot = _rotate(r_pos_trans, 0, -φ, -θ)

        # Calculate the cross points with a centered box
        (t1, c1) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[1, 0, 0], SA_F64[-box.delta_x/2, 0, 0])
        (t2, c2) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[1, 0, 0], SA_F64[+box.delta_x/2, 0, 0])
        (t3, c3) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[0, 1, 0], SA_F64[0, -box.delta_y/2, 0])
        (t4, c4) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[0, 1, 0], SA_F64[0, +box.delta_y/2, 0])
        (t5, c5) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[0, 0, 1], SA_F64[0, 0, -box.delta_z/2])
        (t6, c6) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[0, 0, 1], SA_F64[0, 0, +box.delta_z/2])

        (x0, y0, z0) = SA_F64[-box.delta_x/2, -box.delta_y/2, -box.delta_z/2]
        (x1, y1, z1) = SA_F64[+box.delta_x/2, +box.delta_y/2, +box.delta_z/2]

        # If any of the six points are within the bounded square, then note the cross points
        if c1 !== nothing && (y0 < c1[2] < y1) && (z0 < c1[3] < z1)
            raw_c1 = _translate(_rotate(c1, θ, φ, 0), Δx)
            push!(crosses, t1 => raw_c1)
        end
        if c2 !== nothing && (y0 < c2[2] < y1) && (z0 < c2[3] < z1)
            raw_c2 = _translate(_rotate(c2, θ, φ, 0), Δx)
            push!(crosses, t2 => raw_c2)
        end
        if c3 !== nothing && (x0 < c3[1] < x1) && (z0 < c3[3] < z1)
            raw_c3 = _translate(_rotate(c3, θ, φ, 0), Δx)
            push!(crosses, t3 => raw_c3)
        end
        if c4 !== nothing && (x0 < c4[1] < x1) && (z0 < c4[3] < z1)
            raw_c4 = _translate(_rotate(c4, θ, φ, 0), Δx)
            push!(crosses, t4 => raw_c4)
        end
        if c5 !== nothing && (x0 < c5[1] < x1) && (y0 < c5[2] < y1)
            raw_c5 = _translate(_rotate(c5, θ, φ, 0), Δx)
            push!(crosses, t5 => raw_c5)
        end
        if c6 !== nothing && (x0 < c6[1] < x1) && (y0 < c6[2] < y1)
            raw_c6 = _translate(_rotate(c6, θ, φ, 0), Δx)
            push!(crosses, t6 => raw_c6)
        end
        return length(crosses) > len
    else
        return false
    end
end


"""
Tests if the ray goes through the sphere, assuming the ray can be infinitely long.
"""
function isthrough!(r::Ray, s::Sphere, crosses::SortedDict{Float64,SVector{3,Float64}})::Bool
    vhit = rand() < s.efficiency
    if vhit
        dir = raydir(r)
        # First translate the system to position = (0, 0, 0)
        Δx = s.position

        x0 = r.start_position - Δx
        a = dir ⋅ dir
        b = 2 * (x0 ⋅ dir)
        c = x0 ⋅ x0 - (s.radius)^2
        Δ = b^2 - 4 * a * c

        if Δ <= 0
            return false
        else
            t1 = (-b - √Δ) / (2a)
            t2 = (-b + √Δ) / (2a)
            push!(crosses, t1 => x0 + Δx + t1 * dir)
            push!(crosses, t2 => x0 + Δx + t2 * dir)
        end
        return true
    else
        return false
    end
end


# Functions that modifies arguments in place should end with "!"
"""
Tests if the ray goes through the cylinder, assuming the ray can be infinitely long.
"""
function isthrough!(r::Ray, cyl::Cylinder, crosses::SortedDict{Float64,SVector{3,Float64}})::Bool
    vhit = rand() < cyl.efficiency
    if vhit
        len = length(crosses)
        # Translation
        Δx = cyl.position
        # Orientation
        (θ, φ) = cyl.orientation
        # Ray direction
        dir = raydir(r)
        # Ray start position
        r_pos = r.start_position

        # Translate back to origin
        r_pos_trans = _translate(r_pos, -Δx)
        # Rotate back to normal
        dir_rot = _rotate(dir, 0, -φ, -θ)
        r_pos_trans_rot = _rotate(r_pos_trans, 0, -φ, -θ)

        # Calculate the "time" for the intersection for cylinder x² + y² = r²
        a = dir_rot[1:2] ⋅ dir_rot[1:2]
        b = 2 * (r_pos_trans_rot[1:2] ⋅ dir_rot[1:2])
        c = r_pos_trans_rot[1:2] ⋅ r_pos_trans_rot[1:2] - (cyl.radius)^2
        Δ = b^2 - 4 * a * c
        # No solution unless the ray goes parallel to the z-axis (not really possible)
        if a == 0
            return c < 0
        elseif Δ < 0
            return false
        else
            t1 = (-b - √Δ) / (2a)
            t2 = (-b + √Δ) / (2a)
            c1 = r_pos_trans_rot + t1 * dir_rot
            c2 = r_pos_trans_rot + t2 * dir_rot
            h = cyl.height
            if (-h / 2 < c1[3] < h / 2)
                raw_c1 = _translate(_rotate(c1, θ, φ, 0), Δx)
                push!(crosses, t1 => raw_c1)
            end
            if (-h / 2 < c2[3] < h / 2)
                raw_c2 = _translate(_rotate(c2, θ, φ, 0), Δx)
                push!(crosses, t2 => raw_c2)
            end
            (t3, c3) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[0, 0, 1], SA_F64[0, 0, +h/2])
            (t4, c4) = rayplaneint(dir_rot, r_pos_trans_rot, SA_F64[0, 0, 1], SA_F64[0, 0, -h/2])
            if (c3[1:2] ⋅ c3[1:2] < cyl.radius^2)
                raw_c3 = _translate(_rotate(c3, θ, φ, 0), Δx)
                push!(crosses, t3 => raw_c3)
            end
            if (c4[1:2] ⋅ c4[1:2] < cyl.radius^2)
                raw_c4 = _translate(_rotate(c4, θ, φ, 0), Δx)
                push!(crosses, t4 => raw_c4)
            end
        end
        return length(crosses) > len
    else
        return false
    end
end


"""
Calculate the analytic rate (assuming unit I₀) of a set of detectors by using MC Integration (hemisphere).
This function will determine if the detector vector has the chip in it (detector name being "Chip").
The result needs to be divided by (2π / 3) to be consistent with normalization.
If yes, then it will optimize the integration over the chip.
"""
function analytic_R(detectors::Vector{T};
    excluded_detectors::Union{Vector{T},Nothing}=nothing,
    config::Union{Configuration,Nothing}=nothing,
    seed::Union{Int,Nothing}=nothing) where {T<:LabObject{<:Real}}
    # Find center of the set of detectors:
    center = SA_F64[0, 0, 0]
    ℓ = 0
    chip = false
    for d in detectors
        if d.name == "Chip"
            chip = true
            center = d.position
            # Determine ℓ
            ℓ = max(d.delta_x, d.delta_y, d.delta_z) * 3.0
            break
        else
            center += d.position
            ℓ = max(max(d.delta_x, d.delta_y, d.delta_z) * 3.0, ℓ)
        end
    end
    if !chip
        center /= length(detectors)
    end
    println("Center: $(center)")
    println("ℓ: $(ℓ) (m)")
    # Determine R
    R = 100.0
    # for d1 in detectors
    #     for d2 in detectors
    #         Δs = norm(d1.position .- d2.position)
    #         R = max(Δs * 5.0, R)
    #     end
    # end
    println("R: $(R) (m)")

    if config === nothing
        config = Configuration()
        # Note these are ray directions
        θ = MCIntegration.Continuous(π / 2, 1π)
        φ = MCIntegration.Continuous(0.0, 2π)
        XY = MCIntegration.Continuous(-ℓ / 2, ℓ / 2, adapt=false)
        int_config = Configuration(var=(θ, φ, XY), dof=[[1, 1, 2]], userdata=Dict("Ray" => Ray(0.0, 0.0)))
    else
        int_config = config
    end
    if seed !== nothing
        println("Using seed: $seed")
        int_config.seed = seed
        int_config.rng = MersenneTwister(seed)
    end
    prob_integrand = (X, c) -> (sin(X[1][1]) * cos(X[1][1])^2)
    if excluded_detectors !== nothing
        hit_func = (X, c) -> begin
            r = c.userdata["Ray"]
            (θ̃, φ̃, x̃, ỹ) = (X[1][1], X[2][1], X[3][1], X[3][2])
            modifyray!(r, R, center, ℓ; θ=θ̃, φ=φ̃, x=x̃, y=ỹ)
            return all(d -> isthrough!(r, d, SortedDict{Float64,SVector{3,Float64}}()), detectors) &&
                   !any(d -> isthrough!(r, d, SortedDict{Float64,SVector{3,Float64}}()), excluded_detectors)
        end
    else
        hit_func = (X, c) -> begin
            r = c.userdata["Ray"]
            (θ̃, φ̃, x̃, ỹ) = (X[1][1], X[2][1], X[3][1], X[3][2])
            modifyray!(r, R, center, ℓ; θ=θ̃, φ=φ̃, x=x̃, y=ỹ)
            return all(d -> isthrough!(r, d, SortedDict{Float64,SVector{3,Float64}}()), detectors)
        end
    end
    hit_integrand = (X, c) -> begin
        return hit_func(X, c) ? prob_integrand(X, c) : 0.0
    end
    if chip
        neval = 5e6
    else
        neval = 2e6
    end
    println("neval = $(neval)")
    hit_prob = integrate(hit_integrand, config=int_config, print=-1, neval=neval, niter=15)
    return hit_prob
end


# --- Scratch ---
