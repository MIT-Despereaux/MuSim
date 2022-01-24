### Tests for basic structures and functions.

module BasicTest

include(joinpath(dirname(@__FILE__), "testutils.jl"))

using Test
# The following code is necessary to fix VSCode julia local module linting
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MuSim.jl")
    using .MuSim
else
    using MuSim
end


# %%
"""
Compare hashes.
"""
function testrectboxhash()
    o1 = RectBox("o1", 1, 1, 1, position = (0, 0, 0), orientation = deg2rad.((90, 0)), efficiency = 1.0, material = "Unknown")
    o2 = deepcopy(o1)
    @test hash(o1) == hash(o2)
end


"""
Compare two LabObjects.
"""
function testrectboxequal()
    o1 = RectBox("o1", 1, 1, 1, position = (0, 0, 0), orientation = deg2rad.((90, 0)), efficiency = 1.0, material = "Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


function testsphereequal()
    o1 = Sphere("o1", 1, position = (0, 0, 0), efficiency = 1.0, material = "Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


function testcylinderequal()
    o1 = Cylinder("o1", 1, 1, position = (0, 0, 0), orientation = deg2rad.((90, 0)), efficiency = 1.0, material = "Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


# %%
testrectboxequal()
testsphereequal()
testcylinderequal()

end # module
