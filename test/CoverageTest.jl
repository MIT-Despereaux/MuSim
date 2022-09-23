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

# %%
"""
Test for basic coverage of one box. Returns beta and beta_err multiplied by simulation ℓ^2.
"""
function test_coverage1_hemisimlite(n_sim::Int=50000)
    I₀ = 58
    (x, y, z) = (1, 2, 3)
    expected_rate = (x * y * π / 2 + 2 * x * z * π / 8 + 2 * y * z * π / 8) * I₀

    # Construct a X x Y x Z box (unit m)
    box = [RectBox("A", x, y, z)]
    r = max(x, y, z) * 3
    ℓ = max(x * y, x * z, y * z) * 2
    # println("r = $r, ℓ = $ℓ")
    results = runhemisimlite(n_sim, box, r, ℓ; center=(0, 0, 0))

    β = sum(results) / n_sim
    β_err = √(β * (1 - β) / n_sim)
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    # println("rate = $rate ± $rate_err")
    # println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

# %%
# Basic coverage test
test_coverage1_hemisimlite()

# %%
"""
Test for coincidence between multiple detectors. 
"""
function test_coverage2_hemisimlite(n_sim::Int=1000000)
    I₀ = 58

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ = 27π / 60
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    dets = [box1, box2]

    r = d * 2
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    @time results = runhemisimlite(n_sim, dets, r, ℓ; center=(0, 0, 0))
    res = view(results, :, 1) .& view(results, :, 2)

    β = sum(res) / n_sim
    β_err = √(β * (1 - β) / n_sim)
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

test_coverage2_hemisimlite()

# %%
"""
Test for optimization.
"""
function test_coverage3_hemisimlite(n_sim::Int=1000000)
    I₀ = 1

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ1 = 20π / 60
    box1 = RectBox("A", w, w, 0.001; orientation=(θ1, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ1) * d, 0, cos(θ1) * d), orientation=(θ1, 0))
    dets = [box1, box2]

    r = 100
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    @time (results, dist_θ, dist_φ, angles) = runhemisimlite(n_sim, dets, r, ℓ)
    # println("dist_θ = $dist_θ, dist_φ = $dist_φ")
    angles = reshape(angles, 2, :)
    res = view(results, :, 1) .& view(results, :, 2)
    # p1 = histogram(angles[1, :], bins=500, label="θ")
    # p2 = histogram(angles[2, :], bins=500, label="φ")
    # display(p1)
    # display(p2)
    # Importance sampling
    β = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            β += 3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ))
        end
    end
    β /= n_sim
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
    β_err = (√σβ) / n_sim
    β *= ℓ^2
    β_err *= ℓ^2
    # println("β = $β ± $β_err")
    rate = β * 2π / 3 * I₀
    rate_err = β_err * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ1)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

test_coverage3_hemisimlite()

# %%
"""
Test for optimization for the full simulation.
"""
function test_coverage4_hemisim(n_sim::Int=100000)
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
    (results, _, _, dist_θ, dist_φ, angles) = runhemisim(n_sim, dets, r, ℓ)
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
    β /= n_sim
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
    β_err = (√σβ) / n_sim
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

test_coverage4_hemisim()

# %%
end
