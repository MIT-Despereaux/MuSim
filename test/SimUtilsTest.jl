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
using JLD2

# %%
OUT_DIR = abspath(joinpath(dirname(@__FILE__), "test_cache"))
println(OUT_DIR)
mkpath(OUT_DIR)

# %%
"""
Integration test for the essential functions.
"""
function test_runexp(n::Int=1000000)
    # Construct two planar surfaces oriented at θ
    w = 0.75
    d = 1.0 # r distance between centers (on a sphere)
    θ1 = 29.75π / 60
    θ2 = θ1 + deg2rad(4)
    box1 = RectBox("A", w, w, 0.1; orientation=(θ1, 0))
    box2 = RectBox("B", w, w, 0.1; position=(sin(θ1) * d, 0, cos(θ1) * d), orientation=(θ1, 0))
    box3 = RectBox("C", w, w, 0.1; position=(sin(θ2) * d * 1.2, 0, cos(θ2) * d * 1.2), orientation=(θ2, 0))
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
    config = Dict{String,Any}("sim_num" => n * 10)
    # The radius of the hemisphere in [m]
    config["ℓ"] = 1.0
    config["r"] = 100
    config["center"] = (0, 0, 0)
    config["detectors"] = dets
    push!(sim_configs, config)

    res, sim_configs = runexp(OUT_DIR, sim_configs; overwrite=true)

    # Check the two results are the same
    det_order = Dict{String,Int}(d.name => i for (i, d) in enumerate(dets))

    res1 = calculateβ(det_order, res[1], ["A", "B"], true) .* sim_configs[1]["ℓ"]^2
    res2 = calculateβ(det_order, res[2], ["A", "B"], false) .* sim_configs[2]["ℓ"]^2
    println("res1 = $res1")
    println("res2 = $res2")
    @test all(.≈(res1, res2, rtol=0.1))

    res1 = calculateβ(det_order, res[1], ["A", "C"], true) .* sim_configs[1]["ℓ"]^2
    res2 = calculateβ(det_order, res[2], ["A", "C"], false) .* sim_configs[2]["ℓ"]^2
    println("res1 = $res1")
    println("res2 = $res2")
    @test all(.≈(res1, res2, rtol=0.1))
end

# %%
test_runexp()

# %%
"""
Integration test for geometric factor errors.
"""
function test_compose(n::Int=1000000)
    w = 0.2
    d = 0.5
    θ = deg2rad(43)
    γ = deg2rad(138)
    box1 = RectBox("A", w, w, 0.1; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.1; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    box3 = RectBox("C", w, w, 0.1; position=(sin(γ) * d, 0, cos(γ) * d), orientation=(γ, 0))
    dets = [box1, box2, box3]

    config = Dict{String,Any}("sim_num" => n * 10)
    config["detectors"] = dets
    config["center"] = (0, 0, 0)
    config["ℓ"] = 1.5
    config["r"] = 2.0
    res, sim_configs = runexp(OUT_DIR, [config])
    βs1 = βio(OUT_DIR, res[1], sim_configs[1]; savefile=true)
    βs1 = Dict(βs1)
    println("βs1 = $βs1")

    βs2 = composeβ(OUT_DIR, dets, n)
    βs2 = Dict(βs2)
    println("βs2 = $βs2")

    @test all(map(k -> ≈(βs1[k], βs2[k], rtol=0.05), collect(keys(βs1))))
end

# %%
test_compose()

# %%
end
