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
    results = runhemisimlite(n, box, r, (0, 0, 0), ℓ)

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
function testcoverage2_hemisimlite(n::Int=50000)
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
    results = runhemisimlite(n, dets, r, (0, 0, 0), ℓ)
    res = view(results, :, 1) .& view(results, :, 2)

    β = sum(res) / n
    β_err = √(β * (1 - β) / n)
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    # println("rate = $rate ± $rate_err")
    # println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

testcoverage2_hemisimlite()

# %%
end
