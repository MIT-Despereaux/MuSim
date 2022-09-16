### Tests for simulation utilities.

module SimUtilsTest

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

# %%
OUT_DIR = abspath(joinpath(dirname(@__FILE__), "test_cache"))
println(OUT_DIR)
mkpath(OUT_DIR)

# %%
"""
Integration test for the essential functions.
"""
function testrunexp(n::Int=1000000)
    # Construct two planar surfaces oriented at θ
    w = 0.5
    d = 1.0 # r distance between centers (on a sphere)
    θ = π / 5
    θ2 = θ + deg2rad(4)
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    box3 = RectBox("C", w, w, 0.001; position=(sin(θ2) * d * 1.5, 0, cos(θ2) * d * 1.5), orientation=(θ2, 0))
    dets = [box1, box2, box3]

    sim_configs = Dict{String,Any}[]
    # The list of tangent plane length for ray generation to scan over
    # Total number of simulations
    config = Dict{String,Any}("sim_num" => n)
    # The radius of the hemisphere in [m]
    config["center"] = nothing
    config["detectors"] = dets
    push!(sim_configs, config)

    # Total number of simulations
    config = Dict{String,Any}("sim_num" => n)
    # The radius of the hemisphere in [m]
    config["ℓ"] = 1.0
    config["r"] = 2.0
    config["center"] = (0, 0, 0)
    config["detectors"] = dets
    push!(sim_configs, config)

    res, sim_configs = runexp(OUT_DIR, dets, sim_configs; batch_size=100000)

    # Check the two results are the same
    det_order = Dict{String,Int}()
    for (i, d) in enumerate(dets)
        det_order[d.name] = i
    end

    res1 = calcgeo(det_order, res[1], ["A", "B"], true)
    res2 = calcgeo(det_order, res[2], ["A", "B"], false)
    # println("res1 = $res1")
    # println("res2 = $res2")
    @test all(.≈(res1, res2, rtol=0.1))

    res1 = calcgeo(det_order, res[1], ["A", "C"], true)
    res2 = calcgeo(det_order, res[2], ["A", "C"], false)
    # println("res1 = $res1")
    # println("res2 = $res2")
    @test all(.≈(res1, res2, rtol=0.1))

    βs1 = geometricio(OUT_DIR, res[1], sim_configs[1])
    βs2 = geometricio(OUT_DIR, res[2], sim_configs[2])
    println("βs1 = $βs1")
    println("βs2 = $βs2")
    @test all(map(k -> ≈(βs1[k], βs2[k], rtol=0.1), collect(keys(βs1))))
end

# %%
testrunexp()

# %%
end
