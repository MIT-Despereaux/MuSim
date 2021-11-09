### Scratches and plots

using Plots; gr();
using Plots.PlotMeasures
default(size=(600, 600), dpi=100, margin=3.0mm)
println("Default style set.")
using MuSim
import MuSim:_randcos3

include("DistributionTest.jl")

# %%
using Random
Random.seed!(42)
function plot_hist()
    # Plotting the theoretical and data histogram
    expected = x -> @. abs(4cos(x)^3 * sin(x))
    (bin_centers, counts, errs, expected) = DistributionTest.test_zenangle_dist(_randcos3, expected, π / 2, π)

    p1 = plot(bin_centers, counts, yerr=errs, label="Sim.")
    plot!(bin_centers, expected, label="Th.")
    return p1
end

p1 = plot_hist()

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
