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

import MuSim: _cart2unitsph, _relativedir

# %%
"""
Compare hashes.
"""
function test_rectboxhash()
    o1 = RectBox("o1", 1, 1, 1, position=(0, 0, 0), orientation=deg2rad.((90, 0)), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test hash(o1) == hash(o2)
end


"""
Compare two LabObjects.
"""
function test_rectboxequal()
    o1 = RectBox("o1", 1, 1, 1, position=(0, 0, 0), orientation=deg2rad.((90, 0)), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


function test_sphereequal()
    o1 = Sphere("o1", 1, position=(0, 0, 0), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end


function test_cylinderequal()
    o1 = Cylinder("o1", 1, 1, position=(0, 0, 0), orientation=deg2rad.((90, 0)), efficiency=1.0, material="Unknown")
    o2 = deepcopy(o1)
    @test o1 == o2
end

function test_cart2unitsph()
    @test _cart2unitsph(1, 0, 0) |> Tuple == (π / 2, 0)
    @test _cart2unitsph(0, 1, 0) |> Tuple == (π / 2, π / 2)
    @test _cart2unitsph(0, 0, 1) |> Tuple == (0, 0)
    @test _cart2unitsph(1 / √2, 1 / √2, 1) |> Tuple == (π / 4, π / 4)
end

function test_relativedir()
    box1 = RectBox("A", 0.001, 0.001, 0.001, position=(0, 0, 0), orientation=deg2rad.((0, 0)))
    box2 = RectBox("B", 2.0, 2.0, 0.001, position=(0, 0, √2), orientation=deg2rad.((0, 0)))
    @test all(.≈(_relativedir(box1, box2), (0, 0, π / 4, 6π / 4), atol=1e-3))
    box1 = RectBox("A", 2.0, 2.0, 0.001, position=(0, 0, 0), orientation=deg2rad.((0, 0)))
    box2 = RectBox("B", 0.001, 0.001, 0.001, position=(0, 0, √2), orientation=deg2rad.((0, 0)))
    @test all(.≈(_relativedir(box1, box2), (0, 0, π / 4, 6π / 4), atol=1e-3))
    box1 = RectBox("A", 0.001, 0.001, 0.001, position=(0, 0, 0), orientation=deg2rad.((0, 0)))
    box2 = RectBox("B", 2.0, 1.0, 0.001, position=(0, 0, √2), orientation=deg2rad.((30, 0)))
    @test all(.≈(_relativedir(box1, box2), (0, 0, atan(2 / (2√2 - 1)), 5π / 3), atol=2e-3))
end


# %%
test_rectboxhash()
test_rectboxequal()
test_sphereequal()
test_cylinderequal()
test_cart2unitsph()
test_relativedir()

# %%
end # module
