### Tests for sampling from a distribution

module DistributionTest

include(joinpath(dirname(@__FILE__), "testutils.jl"))

using Test
using MuSim
using StatsBase

# Submodules
# Note: Use explicit import when needed
import MuSim:_randvec!, _randcos2, _randcos3

# %%
"""
Helper function that returns the bin_centers, counts, errs and expected values in each bins from the desired zenith angle distribution.
"""
function test_zenangle_dist(dist::Function, expected::Function, l_bound::Real, u_bound::Real, n::Int=35000)
    vals = zeros(n)
    _randvec!(dist, vals)
    bin_edges = collect(range(l_bound, u_bound, length=30))
    bin_centers = (bin_edges .+ ((bin_edges[2] - bin_edges[1]) / 2))[1:end - 1]
    counts = fit(Histogram, vals, bin_edges).weights
    errs = map(x -> max(x, 1), .√counts)
    counts /= (bin_edges[2] - bin_edges[1]) * n
    errs /= (bin_edges[2] - bin_edges[1]) * n
    expected_counts = expected(bin_centers)

    return (bin_centers, counts, errs, expected_counts)
end

# %%
# Test drawing from a custom distribution
expected_f = x -> @. abs(3cos(x)^2 * sin(x))
(_, counts, errs, expected) = test_zenangle_dist(_randcos2, expected_f, π / 2, π)
println(counts)
@test all(counts - 3.5errs <= expected <= counts + 3.5errs)
expected_f = x -> @. abs(4cos(x)^3 * sin(x))
(_, counts, errs, expected) = test_zenangle_dist(_randcos3, expected_f, π / 2, π)
@test all(counts - 3.5errs <= expected <= counts + 3.5errs)

end # module
