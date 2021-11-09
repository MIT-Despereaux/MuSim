using Test
using CoincidenceSim
using Folds

using StatsBase

# Submodules
# Note: Use explicit import when needed
import CoincidenceSim: _rand_vec!, _rand_cos3

# %%
"""
Helper function that returns the bin_centers, counts, errs and expected values in each bins from the desired distribution.
"""
function test_rand_dist(dist::Function, expected::Function)
    n = Int(3e5)
    vals = zeros(n)
    _rand_vec!(dist, vals)
    bin_edges = collect(range(0, pi / 2, length=30))
    bin_edges .+= pi / 2
    bin_centers = (bin_edges .+ ((bin_edges[2] - bin_edges[1]) / 2))[1:end - 1]
    counts = fit(Histogram, vals, bin_edges).weights
    errs = map(x -> max(x, 1), counts ./ .√counts)
    counts /= (bin_edges[2] - bin_edges[1]) * n
    errs /= (bin_edges[2] - bin_edges[1]) * n
    expected_counts = expected(bin_centers)

    return (bin_centers, counts, errs, expected_counts)
end
    

I0 = 102
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

# %%
# using Plots; gr();
# default(size=(800, 800))

# %%
function plot_hist()
    # Plotting the theoretical and data histogram
    (bin_centers, counts, errs, expected) = test_rand_dist(_rand_cos3, x -> @. abs(4cos(x)^3 * sin(x)))

    p = plot(bin_centers, counts, yerr=errs, label="Sim.")
    plot!(bin_centers, expected, label="Th.")
    return p
end

# %%
function plot_hist2()
    # Plotting uniform distribution test
    (pos, single_hits, single_hits_err, predicted_hits) = test_planar_coincidence()
    
    p = plot(pos, single_hits, yerr=single_hits_err, label="Sim.", xlabel="Distance from source [m]", ylabel="μ Rate [Hz]")
    plot!(pos, predicted_hits, label="Th.")
    title!("1/R² Rule Test for Planar Generation Surface")
    return p
end
    
# %%
function plot_hist3()
    (pos, single_hits, single_hits_err, predicted_hits) = test_hemi_coincidence()
    
    plot(pos, single_hits, yerr=single_hits_err, label="Sim.", xlabel="Distance from source [m]", ylabel="μ Rate [Hz]")
    plot!(pos, predicted_hits, label="Th.")
    title!("1/R² Rule Test for Hemispherical Generation Surface")
end

# %%
# All tests are here
@time @testset "Coincidence" begin
    # Test drawing from a custom distribution
    (_, counts, errs, expected) = test_rand_dist(_rand_cos3, x -> @. abs(4cos(x)^3 * sin(x)))
    @test all(counts - 3.5errs <= expected <= counts + 3.5errs)

    # Test 1/R² rule
    (_, rates, errs, expected) = test_planar_coincidence()
    @test all(rates - 3.5errs <= expected <= rates + 3.5errs)
    (_, rates, errs, expected) = test_hemi_coincidence()
    @test all(rates - 3.5errs <= expected <= rates + 3.5errs)
end
