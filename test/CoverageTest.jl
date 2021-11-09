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
Function that generates the planar sampling and report the counts, errors, and expected values.
"""
function test_planar_coincidence()
    sim_num = Int(3e5)
    # size of the square plane from which rays are generated [m]
    x_bound = 0.05 / 2
    
    # Set up the detectors
    z_pos = range(-1.0, -5.0, step=-0.25)
    
    total_rate = pi / 2 * I0 * (2x_bound)^2 # 1 second normalisation
    
    detectors = [RectBox(0.05, 0.05, 0.001, position=(0.0, 0.0, z), efficiency=1.0) for z in z_pos]
    results = Array{Bool,2}(undef, sim_num, length(detectors))
    
    println("\n--- Planar simulation started ---")
    if Threads.nthreads() == 1
        @time collect(runsim!(i, detectors, results, x_bound, (0, 0, 0)) for i in 1:sim_num)
    else
        @time Folds.collect(runsim!(i, detectors, results, x_bound, (0, 0, 0)) for i in 1:sim_num)
    end

    hits = sum(results, dims=1)
    single_hits = hits * total_rate / sim_num
    single_hits_err = @. single_hits / √hits

    single_hits = vec(single_hits)
    single_hits_err = vec(single_hits_err)
    r = 0.0282
    expected_rate = @. pi^2 * r^2 * I0 * (r^2 / z_pos^2)
    return (abs.(z_pos), single_hits, single_hits_err, expected_rate)
end


function test_hemi_coincidence()
    # number of rays to be simulated
    sim_num = Int(3e5)
    # Radius is r (should have 25 cm^2 projected flat surface)
    r = 0.0282
    ℓ = 0.05
    
    # Set up the detectors
    z_pos = range(-1.0, -5.0, step=-0.25)
        
    total_rate = 2pi / 3 * I0 * ℓ^2 # 1 second normalisation    
    
    detectors = [RectBox(0.05, 0.05, 0.001, position=(0.0, 0.0, z), efficiency=1.0) for z in z_pos]
    results = Array{Bool,2}(undef, sim_num, length(detectors))
    
    println("\n--- Hemispherical simulation started ---")
    if Threads.nthreads() == 1
        @time collect(runhemisim!(i, detectors, results, r, (0, 0, 0), ℓ) for i in 1:sim_num)
    else
        @time Folds.collect(runhemisim!(i, detectors, results, r, (0, 0, 0), ℓ) for i in 1:sim_num)
    end

    hits = sum(results, dims=1)
    single_hits = hits * total_rate / sim_num
    single_hits_err = @. single_hits / √hits

    single_hits = vec(single_hits)
    single_hits_err = vec(single_hits_err)
    expected_rate = @. pi^2 * r^2 * I0 * (r^2 / z_pos^2)
    return (abs.(z_pos), single_hits, single_hits_err, expected_rate)
end

# # Test 1/R² rule
# (_, rates, errs, expected) = test_planar_coincidence()
# @test all(rates - 3.5errs <= expected <= rates + 3.5errs)
# (_, rates, errs, expected) = test_hemi_coincidence()
# @test all(rates - 3.5errs <= expected <= rates + 3.5errs)


end
