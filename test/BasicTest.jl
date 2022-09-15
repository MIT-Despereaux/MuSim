### Tests for basic structures and functions.

module BasicTest

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

import MuSim: _cart2unitsph

# %%
"""
Compare hashes.
"""
function testrectboxhash()
    o1 = RectBox("o1", 1, 1, 1, position=(0, 0, 0), orientation=deg2rad.((90, 0)), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test hash(o1) == hash(o2)
end


"""
Compare two LabObjects.
"""
function testrectboxequal()
    o1 = RectBox("o1", 1, 1, 1, position=(0, 0, 0), orientation=deg2rad.((90, 0)), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


function testsphereequal()
    o1 = Sphere("o1", 1, position=(0, 0, 0), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


function testcylinderequal()
    o1 = Cylinder("o1", 1, 1, position=(0, 0, 0), orientation=deg2rad.((90, 0)), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end

function testcart2unitsph()
    @test _cart2unitsph(1, 0, 0) |> Tuple == (π / 2, 0)
    @test _cart2unitsph(0, 1, 0) |> Tuple == (π / 2, π / 2)
    @test _cart2unitsph(0, 0, 1) |> Tuple == (0, 0)
    @test _cart2unitsph(1 / √2, 1 / √2, 1) |> Tuple == (π / 4, π / 4)
end


# %%
testrectboxhash()
testrectboxequal()
testsphereequal()
testcylinderequal()
testcart2unitsph()

# %%
end # module
