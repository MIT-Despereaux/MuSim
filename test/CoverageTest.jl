### Tests for sampling geometry validations.

module CoverageTest

# %%
include(joinpath(dirname(@__FILE__), "testutils.jl"))

using Test
# The following code is necessary to fix VSCode julia local module linting
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MuSim.jl")
    using .MuSim
else
    using MuSim
end
using Distributions, FLoops
# using Plots

# Submodules
# import MuSim:_randvec!, _randcos2, _randcos3

# %%
"""
Test for basic coverage of one box. Returns beta and beta_err multiplied by simulation ℓ^2.
"""
function testcoverage1_hemisimlite(n::Int=50000)
    I₀ = 58
    (x, y, z) = (1, 2, 3)
    expected_rate = (x * y * π / 2 + 2 * x * z * π / 8 + 2 * y * z * π / 8) * I₀

    # Construct a X x Y x Z box (unit m)
    box = [RectBox("A", x, y, z)]
    r = max(x, y, z) * 3
    ℓ = max(x * y, x * z, y * z) * 2
    # println("r = $r, ℓ = $ℓ")
    results = runhemisimlite(n, box, r, ℓ; center=(0, 0, 0))

    β = sum(results) / n
    β_err = √(β * (1 - β) / n)
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    # println("rate = $rate ± $rate_err")
    # println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

# %%
# Basic coverage test
testcoverage1_hemisimlite()

# %%
"""
Test for coincidence between multiple detectors. 
"""
function testcoverage2_hemisimlite(n::Int=100000)
    I₀ = 58

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ = π / 4
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    dets = [box1, box2]

    r = d * 2
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    results = runhemisimlite(n, dets, r, ℓ; center=(0, 0, 0))
    res = view(results, :, 1) .& view(results, :, 2)

    β = sum(res) / n
    β_err = √(β * (1 - β) / n)
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

testcoverage2_hemisimlite()

# %%
"""
Test for optimization.
"""
function testcoverage3_hemisimlite(n::Int=100000)
    I₀ = 58

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ = 0.0
    θ2 = θ + deg2rad(3)
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    box3 = RectBox("C", w, w, 0.001; position=(sin(θ2) * d * 1.5, 0, cos(θ2) * d * 1.5), orientation=(θ2, 0))
    dets = [box1, box2, box3]

    r = d * 2
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    (results, dist_θ, dist_φ, angles) = runhemisimlite(n, dets, r, ℓ)
    angles = reshape(angles, 2, :)
    res = view(results, :, 1) .& view(results, :, 2)
    # Importance sampling
    β = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            β += 3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ))
        end
    end
    β /= n
    # Calculate the variance
    σβ = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            σβ += (3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ)) - β)^2
        else
            σβ += (β)^2
        end
    end
    β_err = (√σβ) / n
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

testcoverage3_hemisimlite()

# %%
end
