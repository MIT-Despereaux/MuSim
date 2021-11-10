### Tests for sampling geometry validations

module CoverageTest

include(joinpath(dirname(@__FILE__), "testutils.jl"))

using Test
using MuSim

# Submodules
# Note: Use explicit import when needed
# import MuSim:_randvec!, _randcos2, _randcos3

# %%
I₀ = 70

"""
Helper function that generates the planar sampling and report the counts, errors, and expected values.
"""
function testplanarcoverage(n::Int=300000)
    # size of the square plane from which rays are generated [m]
    x_bound = 0.05 / 2
    
    # Set up the detectors
    z_pos = range(-1.0, -3.0, step=-0.25)
    
    total_rate = pi / 2 * I₀ * (2x_bound)^2 # 1 second normalisation
    
    detectors = [RectBox(0.05, 0.05, 0.001, position=(0.0, 0.0, z), efficiency=1.0) for z in z_pos]
    
    println("--- Planar simulation started ---")
    results, _, _ = runsim(n, detectors, x_bound, (0, 0, 0))

    single_hits = sum(results, dims=1)
    single_hits_err = .√single_hits
    single_hits *= total_rate / n
    single_hits_err *= total_rate / n

    r = 0.0282
    expected_rate = @. pi^2 * r^2 * I₀ * (r^2 / z_pos^2)
    return (abs.(z_pos), vec(single_hits), vec(single_hits_err), expected_rate)
end


function testhemisphericalcoverage(n::Int=300000)
    # Radius is r (should have 25 cm^2 projected flat surface)
    r = 0.0282
    ℓ = 0.05
    
    # Set up the detectors
    z_pos = range(-1.0, -3.0, step=-0.25)
        
    total_rate = 2pi / 3 * I₀ * ℓ^2 # 1 second normalisation    
    
    detectors = [RectBox(0.05, 0.05, 0.001, position=(0.0, 0.0, z), efficiency=1.0) for z in z_pos]
    
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
(_, rates, errs, expected) = testplanarcoverage()
@test all(rates - 3.5errs <= expected <= rates + 3.5errs)
(_, rates, errs, expected) = testhemisphericalcoverage()
@test all(rates - 3.5errs <= expected <= rates + 3.5errs)

function testefficientcoverage(n::Int=300000)
    # Testing partial coverage of the hemisphere to increase efficiency of the simulation
    # Generate a reference setup
    r = 1.0
    ℓ = 0.45

    # Set up the detectors
    det1 = RectBox(0.05, 0.05, 0.01, position=(0.0, 0.0, 0.0), efficiency=1.0)
    det2 = RectBox(0.05, 0.05, 0.01, position=(0.25, 0.0, 0.25), efficiency=1.0)

    total_rate = 2pi / 3 * 70 * ℓ^2 # 1 second normalisation    

    detectors = [det1, det2]

    results, _, _ = runhemisim(n, detectors, r, (0, 0, 0), ℓ)

    coin_hits = sum(results[:, 1] .& results[:, 2])
    coin_hits_err = √coin_hits
    coin_hits *= total_rate / n
    coin_hits_err *= total_rate / n

    # println(coin_hits)

    # Efficient simulation
    θs = deg2rad.((130, 140))
    φs = deg2rad.((180 - 10, 180 + 10))
    results, _, _ = runhemisim(n, detectors, r, (0, 0, 0), ℓ; θ_range=θs, φ_range=φs)
    total_rate = (φs[2] - φs[1]) * abs(1 / 3 * (cos(θs[2])^3 - cos(θs[1])^3)) * 70 * ℓ^2 # 1 second normalisation    

    coin_hits_eff = sum(results[:, 1] .& results[:, 2])
    coin_hits_eff_err = √coin_hits_eff
    coin_hits_eff *= total_rate / n
    coin_hits_eff_err *= total_rate / n
    
    # println(coin_hits_eff)
    
    @test abs(coin_hits - coin_hits_eff) ≈ 0 atol = √(coin_hits_err^2 + coin_hits_eff_err^2)
    return nothing
end

testefficientcoverage()

end
