### Scratches and plots


using MuSim
import MuSim: _randcos2

# %%
include("DistributionTest.jl")

using Random
Random.seed!(42)
# Plotting the theoretical and data histogram
expected_f = x -> @. abs(3cos(x)^2 * sin(x))
(bin_centers, counts, errs, expected) = DistributionTest.testzenangledist(_randcos3, expected_f, π / 2, π)

p1 = plot(bin_centers, counts, yerr = errs, label = "Sim.")
plot!(bin_centers, expected, label = "Th.")

# %%
include("CoverageTest.jl")
# Plotting uniform distribution test
(pos, single_hits, single_hits_err, predicted_hits) = CoverageTest.testplanarcoverage()

p1 = plot(pos, single_hits, yerr = single_hits_err, label = "Sim.", xlabel = "Distance from source [m]", ylabel = "μ Rate [s⁻¹]")
plot!(pos, predicted_hits, label = "Th.")
title!("1/R² Rule Test for Planar Generation Surface")

# %%
(pos, single_hits, single_hits_err, predicted_hits) = CoverageTest.testhemisphericalcoverage()

p1 = plot(pos, single_hits, yerr = single_hits_err, label = "Sim.", xlabel = "Distance from source [m]", ylabel = "μ Rate [Hz]")
plot!(pos, predicted_hits, label = "Th.")
title!("1/R² Rule Test for Hemispherical Generation Surface")

# %%

