### Tests for sampling from a distribution.

module DistributionTest

# %%
include("testutils.jl")

using Test
# using Infiltrator

# The following code is necessary to fix VSCode julia local module linting
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MuSim.jl")
    using .MuSim
else
    using MuSim
end

# Submodules
# Note: Use explicit import when needed
import MuSim: _randcos2, _randcos3

using StatsBase


# %%
"""
Helper function that returns the bin_centers, counts, errs and expected values in each bins from the desired θ angle distribution.
"""
function test_θdist(dist::Function, expected::Function, l_bound::Real, u_bound::Real, n::Int=35000)
    vals = [dist(l_bound, u_bound) for i in 1:n]
    bin_edges = collect(range(l_bound, u_bound, length=30))
    bin_centers = (bin_edges.+((bin_edges[2]-bin_edges[1])/2))[1:end-1]
    counts = fit(Histogram, vals, bin_edges).weights
    errs = map(x -> max(x, 1), .√counts)
    counts /= (bin_edges[2] - bin_edges[1]) * n
    errs /= (bin_edges[2] - bin_edges[1]) * n
    expected_counts = expected(bin_centers)

    return (bin_centers, counts, errs, expected_counts)
end

function main()
    initrand()
    # Test drawing from a custom distribution
    expected_f = x -> @. abs(3cos(x)^2 * sin(x))
    (_, counts, errs, expected) = test_θdist(_randcos2, expected_f, π / 2, π)
    @test all(counts - 3.5errs <= expected <= counts + 3.5errs)
    expected_f = x -> @. abs(4cos(x)^3 * sin(x))
    (_, counts, errs, expected) = test_θdist(_randcos3, expected_f, π / 2, π)
    @test all(counts - 3.5errs <= expected <= counts + 3.5errs)
end

!isinteractive() && main()

end # module
