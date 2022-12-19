### Tests for simulation utilities.

module SimUtilsTest

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

using JLD2

# %%
"""
Setup function for test geometry.
"""
function setup_geometry()
    hL = 70e-3
    dL = 70e-3
    lenL = 600e-3

    hS = 20e-3
    dS = 72e-3
    lenS = 510e-3

    hCW = 10e-3
    dCW = 50e-3
    lenCW = 50e-3
    hhCW = 30e-3 # height of the scintillator from the bottom

    gCB = 8e-3
    gBD = 11e-3
    gDA = 6e-3

    (xC, yC, zC) = (0.625, 0.220, 0.8682)

    zB = zC + hL / 2 + gCB + hS / 2
    zD = zB + hS / 2 + gBD + hL / 2
    zA = zD + hL / 2 + gDA + hS / 2
    zCW = zC + hL / 2 + hhCW

    xL = xC
    xS = xC - lenL / 2 + lenS / 2
    xCW = xS + lenS / 2 + lenCW / 2 + 0.002


    loc_det_C = (xL, yC, zC)
    loc_det_B = (xS, yC, zB)
    loc_det_D = (xL, yC, zD)
    loc_det_A = (xS, yC, zA)
    loc_det_CWD = (xCW, yC, zCW)

    loc_Chip = (0.338, 0.229, 1.4148)

    det_CWD = RectBox("CWD", lenCW, dCW, hCW, position=loc_det_CWD, material="POP Doped Polystyrene")
    det_A = RectBox("A", lenS, dS, hS, position=loc_det_A, material="EJ-200")
    det_D = RectBox("D", lenL, dL, hL, position=loc_det_D, material="EJ-200")
    det_B = RectBox("B", lenS, dS, hS, position=loc_det_B, material="EJ-200")
    det_C = RectBox("C", lenL, dL, hL, position=loc_det_C, material="EJ-200")
    det_Chip = RectBox("Chip", 350e-6, 0.005, 0.005, position=loc_Chip, material="Unknown")

    detectors = [det_A, det_B, det_C, det_D, det_Chip]
    return detectors
end


# %%
"""
Comparison test for the essential methods.
"""
function test_comparison1(output_dir; n_sim::Int=Int(1e6))
    # Construct two planar surfaces oriented at θ
    w = 0.8
    d = 1.0 # r distance between centers (on a sphere)
    θ1 = 11π / 60
    θ2 = θ1 + deg2rad(2)
    box1 = RectBox("A", w, w, 0.1; orientation=(θ1, 0))
    box2 = RectBox("B", w, w, 0.1; position=(sin(θ1) * d, 0, cos(θ1) * d), orientation=(θ1, 0))
    box3 = RectBox("C", w, w, 0.1; position=(sin(θ2) * d * 1.2, 0, cos(θ2) * d * 1.2), orientation=(θ2, 0))
    dets = [box1, box2, box3]

    sim_configs = Dict{String,Any}[]
    # Total number of simulations
    config = Dict{String,Any}("sim_num" => n_sim)
    # The radius of the hemisphere in [m]
    config["ℓ"] = 2.5
    config["R"] = 10.0
    config["center"] = (0, 0, 0)
    config["detectors"] = dets
    push!(sim_configs, config)

    res, sim_configs = runexp(output_dir, sim_configs; overwrite=true)

    # Check the two results are the same
    det_order = Dict{String,Int}(d.name => i for (i, d) in enumerate(dets))

    res1 = calculateβ(sim_configs[1], det_order, res[1], ["A", "B"])
    res2 = calculateβ_MC(sim_configs[1], ["A", "B"])
    println("res_AB_1 = $res1")
    println("res_AB_2 = $res2")
    @test all(.≈(res1, res2, rtol=0.1))

    res1 = calculateβ(sim_configs[1], det_order, res[1], ["A", "C"])
    res2 = calculateβ_MC(sim_configs[1], ["A", "C"])
    println("res_AC_1 = $res1")
    println("res_AC_2 = $res2")
    @test all(.≈(res1, res2, rtol=0.1))
end


"""
Integration test for geometric factor errors.
"""
function test_comparison2(output_dir; n_sim::Int=Int(1e6))
    dets = setup_geometry()

    config = Dict{String,Any}("sim_num" => n_sim)
    config["detectors"] = dets
    config["center"] = (0.338, 0.229, 1.4148)
    config["ℓ"] = 3.0
    config["R"] = 100.0
    res, sim_configs = runexp(output_dir, [config])
    βs1 = βio(output_dir, res[1], sim_configs[1]; savefile=true, overwrite=false)
    βs1 = Dict(βs1)

    βs2 = βio_MC(output_dir, sim_configs[1]; savefile=true, overwrite=true)
    βs2 = Dict(βs2)

    for k in keys(βs1)
        println("βs1[$k] = $(βs1[k])")
        println("βs2[$k] = $(βs2[k])")
    end
end


# """
# Test for detector position sampling.
# """
# function test_detpos(n_sim::Int=Int(1e6))
#     w = 0.3
#     d = 0.25
#     θ = deg2rad(43)
#     γ = deg2rad(125)
#     box1 = RectBox("A", w, w, 0.1; orientation=(θ, 0))
#     box2 = RectBox("B", w, w, 0.1; position=(-sin(θ) * d, 0, cos(θ) * d), orientation=(θ, π))
#     box3 = RectBox("C", w, w, 0.1; position=(sin(γ) * d, 0, cos(γ) * d), orientation=(γ, 0))
#     dets = [box1, box2, box3]

#     gen_det_list = gendetectorpos(dets, 0.0001)
#     for det_list in gen_det_list
#         (cfg_hash, βs) = composeβ(OUT_DIR, det_list, n_sim)
#         println("config_hash: $cfg_hash")
#         println("βs = $βs")
#     end
# end


function main()
    # initrand()
    OUT_DIR = abspath(joinpath(dirname(@__FILE__), "test_cache"))
    println(OUT_DIR)
    mkpath(OUT_DIR)

    # test_comparison1(OUT_DIR)
    # test_comparison2(OUT_DIR, n_sim=Int(1e6))
    # test_detpos()
end

end
