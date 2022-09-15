### Tests for sampling geometry validations.

module CoverageTest

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
I₀ = 58
total_rate = 2pi / 3 * I₀ * ℓ^2 # 1 second normalisation    

"""
Test for basic coverage of one box.
"""
function testcoverage(n::Int=300000)
    # Construct a 1 m x 2 m x 3 m box



    detectors = [RectBox("A", 0.05, 0.05, 0.001, position=(0.0, 0.0, z), efficiency=1.0) for z in z_pos]

    println("--- Hemispherical simulation started ---")
    results, _, _ = runhemisim(n, detectors, r, (0, 0, 0), ℓ)

    single_hits = sum(results, dims=1)
    single_hits_err = .√single_hits
    single_hits *= total_rate / n
    single_hits_err *= total_rate / n

    expected_rate = @. pi^2 * r^2 * I₀ * (r^2 / z_pos^2)
    return (abs.(z_pos), vec(single_hits), vec(single_hits_err), expected_rate)
end

# %%
# Test 1/R² rule
(_, rates, errs, expected) = testhemisphericalcoverage()
@test all(rates - 3.5errs <= expected <= rates + 3.5errs)

function testefficientcoverage(n::Int=300000)
    # Testing partial coverage of the hemisphere to increase efficiency of the simulation
    println("--- Partial Coverage Test started ---")

    # Generate a reference setup
    r = 0.70
    ℓ = 0.10

    # Set up the detectors
    det1 = RectBox("A", 0.05, 0.05, 0.01, position=(0.0, 0.0, 0.0), efficiency=1.0)
    det2 = RectBox("B", 0.05, 0.05, 0.01, position=(0.25, 0.0, 0.25 * sqrt(3.0)), efficiency=1.0)

    total_rate = 2pi / 3 * I₀ * ℓ^2 # 1 second normalisation    

    detectors = [det1, det2]

    results, _, _ = runhemisim(n, detectors, r, (0, 0, 0), ℓ)

    coin_hits = sum(results[:, 1] .& results[:, 2])
    println("Total hits: $coin_hits")
    coin_hits_err = √coin_hits
    coin_hits *= total_rate / n
    coin_hits_err *= total_rate / n


    # Efficient simulation
    θs = deg2rad.((25, 35))
    φs = deg2rad.((-10, 10))
    results, _, _ = runhemisim(n, detectors, r, (0, 0, 0), ℓ; θ_range=θs, φ_range=φs)
    total_rate = (φs[2] - φs[1]) * abs(1 / 3 * (cos(θs[2])^3 - cos(θs[1])^3)) * I₀ * ℓ^2 # 1 second normalisation    

    coin_hits_eff = sum(results[:, 1] .& results[:, 2])
    println("Total hits: $coin_hits_eff")
    coin_hits_eff_err = √coin_hits_eff
    coin_hits_eff *= total_rate / n
    coin_hits_eff_err *= total_rate / n


    @test abs(coin_hits - coin_hits_eff) ≈ 0 atol = √(coin_hits_err^2 + coin_hits_eff_err^2)
    return nothing
end

testefficientcoverage()

end
