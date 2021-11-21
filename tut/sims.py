# This demo script shows how to call julia functions from python.

from julia import Main

print("The following lines are run from julia with code block:")

results = Main.eval("""
    using MuSim
    sim_num = Int(100)
    r = 0.80
    ℓ = 0.75
    det1 = RectBox("Name", 0.01, 0.05, 0.05, position=(0, 0, 0), orientation=deg2rad.((0, 0)), efficiency=0.98, material="POP Doped Polystyrene")
    detectors = [det1]
    runhemisim(sim_num, detectors, r, (0, 0, 0), ℓ)
    """)
# Note: results[0] is NOT a sparse matrix so it will take up a lot of space.
# One could get around this by getting as much work done in julia, therefore
# outputting only high level analysis variables.
print(results)

# Running a julia script directly and get the script output.
print("The following lines are run from julia with include:")
results = Main.eval("""
    include("simulations.jl")
    """)
print("Result sparse array density: {:.5f}".format(
    results['sim_res'].sparse.density))
