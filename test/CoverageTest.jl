### Tests for sampling geometry validations.

module CoverageTest

# %%
include("testutils.jl")

using Test

# The following code is necessary to fix VSCode julia local module linting
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MuSim.jl")
    using .MuSim
else
    using MuSim
end

# using Plots
using MCIntegration
using Infiltrator

# %%
"""
Test for basic coverage of one box. Returns beta and beta_err multiplied by simulation ℓ^2.
"""
function test_coverage1_hemisimlite(n_sim::Int=800000)
    I₀ = 58
    (x, y, z) = (0.12, 0.19, 0.22)
    ori = deg2rad.((33.0, 112.5))

    # Construct a X x Y x Z box (unit m)
    box = [RectBox("A", x, y, z, orientation=ori)]
    R = max(x, y, z) * 3.0
    ℓ = max(x, y, z) * 2.0
    # println("R = $R, ℓ = $ℓ")
    (results,) = runhemisimlite(n_sim, box, R, ℓ; center=(0, 0, 0))

    expected_rate = analytic_R(box[1], I₀=I₀)

    @time (int_R, int_result) = analytic_R(box, R, ℓ, seed=1234, I₀=I₀)
    println("int_R = $int_R")

    (hit_prob_mean, hit_prob_std, _) = MCIntegration.average(int_result.iterations, init=5)
    # Note the integration already includes the normalization factor 2π/3 * ℓ^2
    int_rate = hit_prob_mean * I₀
    int_rate_err = hit_prob_std * I₀

    β = sum(results) / n_sim
    β_err = √(β * (1 - β) / n_sim)
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    println("rate = $rate ± $rate_err")
    println("integrated rate = $int_rate ± $int_rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
    @test int_rate ≈ expected_rate atol = 2int_rate_err
end

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
    @time (results,) = runhemisimlite(n_sim, dets, r, ℓ)
    res = view(results, :, 1) .& view(results, :, 2)

    β = sum(res) / n_sim
    β_err = √(β * (1 - β) / n_sim)
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    # println("rate = $rate ± $rate_err")
    # println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end


# %%
"""
Test for MC Integration.
"""
function test_coverage3_hemisimlite(n_sim::Int=10000000)
    I₀ = 1

    # Construct two planar surfaces oriented at θ
    w = 0.01
    d = 1.0 # R distance between centers (on a sphere)
    for θ in [11π / 60, 20π / 60, 27π / 60]
        box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
        box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
        dets = [box1, box2]

        R = d * 3
        ℓ = w * 1.3
        # println("R = $R, ℓ = $ℓ")
        @time (results,) = runhemisimlite(n_sim, dets, R, ℓ)
        res = view(results, :, 1) .& view(results, :, 2)

        β = sum(res) / n_sim
        β_err = √(β * (1 - β) / n_sim)
        rate = β * ℓ^2 * 2π / 3 * I₀
        rate_err = β_err * ℓ^2 * 2π / 3 * I₀

        @time (_, int_result) = analytic_R(dets, R, w * 2.5, seed=1234, I₀=I₀)

        (hit_prob_mean, hit_prob_std, _) = MCIntegration.average(int_result.iterations, init=5)
        # # Note the integration already includes the normalization factor 2π/3 * ℓ^2
        int_rate = hit_prob_mean * I₀
        int_rate_err = hit_prob_std * I₀

        expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
        println("rate = $rate ± $rate_err")
        println("integrated rate = $int_rate ± $int_rate_err")
        println("expected rate = $expected_rate")
        @test rate ≈ expected_rate atol = 2rate_err
        @test int_rate ≈ expected_rate atol = 2int_rate_err
    end
end


function main()
    initrand()
    test_coverage1_hemisimlite()
    test_coverage2_hemisimlite()
    test_coverage3_hemisimlite()
end

!isinteractive() && main()

end
