using StatsBase, Distributions, Interpolations, Plots, MuSim, LinearAlgebra

# %%
function setup_geometry()
    # This is the reference location (where the sim will be centered)

    hL = 70e-3
    wL = 70e-3
    lenL = 600e-3

    hS = 20e-3
    wS = 72e-3
    lenS = 510e-3

    hCW = 10e-3
    wCW = 50e-3
    lenCW = 50e-3
    hhCW = 30e-3 # height of the scintillator from the bottom

    xA = 2.54e-2 + hS / 2
    zA = 1.788 + wS / 2
    yA = lenL / 2

    xD = 2.54e-2 + hL / 2
    zD = 1.029 + wL / 2
    yD = lenL / 2

    xB = 2.54e-2 + hL + 0.505 + wS / 2
    zB = 1.029 + hS / 2
    yB = lenS / 2

    xC = 2.54e-2 + hL + 0.6 + wL / 2
    zC = 1.029 + hL / 2
    yC = lenL / 2


    loc_det_C = (xC, yC, zC)
    loc_det_B = (xB, yB, zB)
    loc_det_D = (xD, yD, zD)
    loc_det_A = (xA, yA, zA)

    loc_chip = (0.338, 0.229, 1.4148)
    # loc_chip = ref_loc = (0.0, 0.0, 1.0)

    # det_CWD = RectBox("CWD", lenCW, dCW, hCW, position=loc_det_CWD, material="POP Doped Polystyrene")
    det_A = RectBox("A", hS, lenS, wS, position=loc_det_A, material="EJ-200")
    det_D = RectBox("D", wL, lenL, hL, position=loc_det_D, material="EJ-200")
    det_B = RectBox("B", wS, lenS, hS, position=loc_det_B, material="EJ-200")
    det_C = RectBox("C", wL, lenL, hL, position=loc_det_C, material="EJ-200")
    det_Chip = RectBox("Chip", 350e-6, 0.005, 0.005, position=loc_chip, material="Unknown")


    # remove detectors here as needed
    detectors = [det_A, det_B, det_C, det_D, det_Chip]
    return detectors
end

# %%
sim_configs = []
config = Dict{String,Any}("sim_num" => Int(3e7))
# The radius of the hemisphere in [m]
config["ℓ"] = 3.0
config["R"] = 100.0
# config["center"] = (0.338, 0.3, 1.4148)
config["center"] = (0.6, 0.3, 1.0)
dets = setup_geometry()
config["detectors"] = dets
push!(sim_configs, config)

OUT_DIR = abspath(joinpath(dirname(@__FILE__), "test_cache"))
mkpath(OUT_DIR)
results_raw, sim_configs_raw = runexp(OUT_DIR, sim_configs; lite=false, overwrite=false)

# %%
included_dets = [x for x in dets if (x.name == "A" || x.name == "D")]
res = analytic_R(included_dets)

# %%
sim_configs = []
config = Dict{String,Any}("sim_num" => Int(2e7))
# The radius of the hemisphere in [m]
config["ℓ"] = 3.0
config["R"] = 101.0
# config["center"] = (0.338, 0.3, 1.4148)
config["center"] = (0.6, 0.3, 1.0)
dets = setup_geometry()
config["detectors"] = dets
push!(sim_configs, config)
results_new, sim_configs_new = runexp(OUT_DIR, sim_configs; lite=false, overwrite=true, mc_config=res.config)

# %%
det = 4
θs_raw = [v[1] for i in range(1, 20) for (k, v) in results_raw[1][i][3][det] if (results_raw[1][i][1][k, det] && results_raw[1][i][1][k, 1])]
θs_new = [v[1] for i in range(1, 15) for (k, v) in results_new[1][i][3][det] if (results_new[1][i][1][k, det] && results_new[1][i][1][k, 1])]
histogram(θs_raw, bins=250, alpha=0.5, normalize=:pdf, label="θ_hemi")
histogram!(θs_new, bins=250, alpha=0.5, normalize=:pdf, label="θ_mc")

# %%
φ_raw = [v[2] for i in range(1, 20) for (k, v) in results_raw[1][i][3][det] if (results_raw[1][i][1][k, det] && results_raw[1][i][1][k, 1])]
φ_new = [v[2] for i in range(1, 15) for (k, v) in results_new[1][i][3][det] if (results_new[1][i][1][k, det] && results_new[1][i][1][k, 1])]
histogram(φ_raw, bins=250, alpha=0.5, normalize=:pdf, label="φ_hemi")
histogram!(φ_new, bins=250, alpha=0.5, normalize=:pdf, label="φ_mc")

# %%
det = 4
Ds_raw = [norm(diff(v, dims=1)) for i in range(1, 20) for (k, v) in results_raw[1][i][2][det] if (results_raw[1][i][1][k, det] && results_raw[1][i][1][k, 1])]
Ds_new = [norm(diff(v, dims=1)) for i in range(1, 10) for (k, v) in results_new[1][i][2][det] if (results_new[1][i][1][k, det] && results_new[1][i][1][k, 1])]
# weights_new = Float64[]
# for i in range(1, 10)
#     for j in range(1, length(results_new[1][i][end]))
#         if (results_new[1][i][1][j, det] && results_new[1][i][1][j, 1])
#             push!(weights_new, results_new[1][i][end][j])
#         end
#     end
# end
# println(weights_raw)
# println(weights_new)
histogram(Ds_raw, bins=250, alpha=0.2, normalize=:pdf, label="D_hemi")
histogram!(Ds_new, bins=250, alpha=0.2, normalize=:pdf, label="D_mc")
# histogram!(Ds_new, bins=500, weights=weights_new, alpha=0.2, normalize=:pdf, legend=:topleft)
