# This is a generic script for simulating and saving results of a bunch of detectors.
# Actual analysis (such as coincidence rates, etc) are done by reading the resultant pkl files.

include(joinpath(dirname(@__FILE__), "../test/testutils.jl"))
initplots()

# Set up output directory for the simulation
OUT_DIR = abspath(joinpath(dirname(@__FILE__), "./demo"));
println("Output directory: $OUT_DIR")
mkpath(OUT_DIR)

# The following code is necessary to fix VSCode julia local module linting
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MuSim.jl")
    using .MuSim
else
    using MuSim
end

# Using pandas from python to save DataFrames to pkl files
using PythonCall
pd = pyimport("pandas")
using Printf

# %%
# Constants
I₀ = 57.0 # in Hz m⁻² Sr⁻¹
b₀ = 0.635 # in Hz

# %%
# Construct a simulation configuration vector, each entry being a dictionary [variable name -> value]
sim_configs = Dict{String,Real}[]
# Use hemispherical generation for efficiency
# The list of tangent plane length for ray generation to scan over
for ℓ = 0.01:0.02:0.10
    # Total number of simulations
    config = Dict{String,Real}("sim_num" => Int(1e4))
    # The radius of the hemisphere in [m]
    config["r"] = 0.8
    config["ℓ"] = ℓ
    push!(sim_configs, config)
end

"""
Runs the simulation in given configurations. It scans over the ℓ parameter
in the given list (provided in sim_configs). 
The output are vectors of simulation results.
"""
function runexp(detectors, sim_configs; θs=deg2rad.((0, 90)), φs=deg2rad.((0, 360)))
    # Run the simulation using the setup and return raw results
    # Print number of threads
    @printf("Total threads: %d\n", Threads.nthreads())

    # Initilize the result vectors
    # Each entry in the vector is a full simulation output
    results = Vector{Any}(undef, length(sim_configs))
    crx_pts = similar(results)
    ray_dirs = similar(results)

    for (i, config) in enumerate(sim_configs)
        ℓ = config["ℓ"]
        r = config["r"]
        sim_num = config["sim_num"]
        println("\n--- Simulation events $sim_num, ℓ#$i, ℓ=$ℓ r=$r started ---")
        # results: size(sim_num, detectors) sparse matrix with trues and falses.
        # crx_pts: size(detectors) vector of dictionary with [event → two cross points of this event] as [key → value] pairs.
        # ray_dirs: size(detectors) vector of dictionary with [event → (θ, φ) of the muon direction].
        results[i], crx_pts[i], ray_dirs[i] = runhemisim(sim_num, detectors, r, (0, 0, 0), ℓ; θ_range=θs, φ_range=φs)
    end

    return results, crx_pts, ray_dirs
end


# %%
# Set up the detectors
dz = 0.05
dy = 0.0
# Adjusts the coverage to increase simulation efficiency
# Note this will make single detector rates invalid and only coincidence rates will be correct.
θs = deg2rad.((5, 15)) # from 5° to 15°
φs = deg2rad.((-25, 25)) # ±25°

det1 = RectBox("Detector 1", 0.01, 0.05, 0.05, position=(0.3, dy, dz),
    orientation=deg2rad.((0, 0)), efficiency=0.98,
    material="POP Doped Polystyrene")
chip = RectBox("Chip", 350e-6, 0.005, 0.005, position=(0, dy, 0),
    orientation=deg2rad.((0, 0)), efficiency=1.0,
    material="Unknown")
det2 = RectBox("Detector 2", 0.01, 0.05, 0.05, position=(-0.3, dy, -dz),
    orientation=deg2rad.((0, 0)), efficiency=0.98,
    material="POP Doped Polystyrene")

# Put the detectors into a list
detectors = [det1, chip, det2]
for d in detectors
    println(d)
end

results, crx_pts, ray_dirs = runexp(detectors, sim_configs; θs=θs, φs=φs);

# %%
# Visualization
"""
This function plots the cross points for the simulation as a scatter plot.
It assumes the points are coming from a single simulation.
Returns the plot for further manipulation.
"""
function plotcrxpts(res, crx, detectors)
    all_pts = [] # the n_hit x n_det array that contains all the hits for the detectors
    # Note: this array is NOT regular since n_hit will differ by different detectors
    # Loop through all the detectors to push the cross points
    for (j, d) in enumerate(detectors)
        # The column for this particular detector
        r = res[:, j]
        println("Total events through detector $j: $(sum(r))")
        # Select the cross points for this detector
        crosses = crx[j]
        pts = vcat([val for val in values(crosses)]...)
        push!(all_pts, pts)
    end

    # Configure palette for the plot
    pal = cgrad(:matter, length(detectors), categorical=true)
    p1 = Plots.plot(dpi=300, size=(1000, 1000), palette=pal)
    for (j, d) in enumerate(detectors)
        pts = all_pts[j]
        # Set the maximum number of points to plot for each detector
        max_pts = min(size(pts)[1], 250)
        scatter!(pts[1:max_pts, 1], pts[1:max_pts, 2], pts[1:max_pts, 3],
            aspect_ratio=:equal, markersize=1.5, markeralpha=0.7,
            markerstrokewidth=0.2, label="Det. $j")
    end
    # Adjust the axes limits
    xaxis!("x", (-0.5, 0.5))
    yaxis!("y", (-0.5, 0.5))
    plot!(zaxis=("z", (-0.5, 0.5)))
    savefig(p1, joinpath(OUT_DIR, "Hitpoint_Scatter.png"))
    return p1
end

"""
This function plots the rates for some selected detector combinations 
with scanned ℓ_list as the x-axis. If the length of the list is 1 then abort and return nothing.
The combinations should be in a list in ["1", "12", "14", ...] format, 
with the numbers denoting the detector columns.
Assumes that the given sim_configs contains only scan over ℓ.
Returns the plot for further manipulation.
"""
function plotratestab(results, sim_configs, combinations, θs, φs)
    # Plot the detector rate against the simulation configs
    if length(sim_configs) == 1
        @warn "sim_configs length is 1, abort plotting..."
        return nothing
    end
    rates = zeros(length(sim_configs), length(combinations))
    for (i, config) in enumerate(sim_configs)
        ℓ = config["ℓ"]
        # Normalize the flux
        total_rate = (φs[2] - φs[1]) * abs(1 / 3 * (cos(θs[2])^3 - cos(θs[1])^3)) * I₀ * ℓ^2 # 1 second normalisation
        res = results[i]
        n_sim = size(res)[1]
        # Loop through the combinations
        for (j, c) in enumerate(combinations)
            # Initilize a temporary result vector
            res_tmp = trues(n_sim)
            for char in c
                col = parse(Int, char)
                # Check the column number is in bound
                @assert col <= size(res)[2] "Column $col exceeds the number of detector in simulation."
                res_tmp = res_tmp .& view(res, :, col)
            end
            # Fill the rates matrix
            rates[i, j] = sum(res_tmp)
        end
        rates[i, :] *= (total_rate / n_sim)
    end

    pal = palette([:purple, :green], length(combinations))
    p1 = Plots.plot(palette=pal)
    plot!([c["ℓ"] for c in sim_configs], rates, label=hcat(["Comb. $c" for c in combinations]...), legend=:topleft)
    xlabel!("Tangent generation plane side length [m]")
    ylabel!("Muon rate [s⁻¹]")
    title!("Vertical Muon Rate = $I₀ [m⁻² s⁻¹ Sr⁻¹]")
    savefig(p1, joinpath(OUT_DIR, "Rates_Stability.pdf"))
    return p1
end

# Picks out the desired simulation result from the scanned list
i = length(sim_configs)
res = results[i]
crx = crx_pts[i]
dirs = ray_dirs[i]
config = sim_configs[i]
# Plot the combinations
combinations = ["12", "13", "23"]
p1 = plotcrxpts(res, crx, detectors)
p2 = plotratestab(results, sim_configs, combinations, θs, φs)

# %%
# Save simulation results as .pkl files (pandas DataFrame)
"""
This function converts the simulation results to DataFrame and saves it as pkl files.
Names of the files are sim_res, sim_crx, sim_dirs, etc.
Overwrite controls whether to overwrite existing data.
The function reads the files on disk if they are already there and returns
a dictionary containing the DataFrames.
"""
function expio(res, crx, dirs, config, detectors; overwrite=false)
    # Construct the metadata
    metadata = "Total number of simulations: $(size(res)[1])\n"
    detector_names = String[]
    for d in detectors
        metadata *= sprint(show, d)
        metadata *= repeat("-", 50)
        metadata *= "\n"
        push!(detector_names, d.name)
    end
    metadata *= "Simulation Configurations:\n"
    for (k, v) in config
        metadata *= "$k -> $v \n"
    end
    println(metadata)

    # Construct the DataFrames and sort in place
    df_res = pd.DataFrame(res, columns=detector_names, dtype=pd.SparseDtype("bool", false))

    df_crx = pd.DataFrame(crx, index=detector_names, dtype=pd.SparseDtype("object")).transpose()
    df_crx.sort_index(inplace=true)

    df_dirs = pd.DataFrame(dirs, index=detector_names, dtype=pd.SparseDtype("object")).transpose()
    df_dirs.sort_index(inplace=true)

    file_names = ["sim_res", "sim_crx", "sim_dirs"]
    files = [df_res, df_crx, df_dirs]


    # Save the metadata and DFs to disk
    io = open(joinpath(OUT_DIR, "metadata.txt"), "w")
    write(io, metadata)
    close(io)

    sim_res = Dict{String,PyObject}()
    for (n, d) in zip(file_names, files)
        f_name = joinpath(OUT_DIR, n * ".pkl")
        if !isfile(f_name) || overwrite
            d.to_pickle(f_name)
        else
            d = read_pickle(f_name)
        end
        sim_res[n] = d
    end
    return sim_res
end

# %%
sim_res = expio(res, crx, dirs, config, detectors; overwrite=true)

# %%
# This demo shows what saved data looks like
# Plotting the histogram of simulated θ angles that hit both the chip and detector 1
df_dirs = sim_res["sim_dirs"]
df_raw = df_dirs["Detector 1"]
df_mask = df_dirs["Detector 1"].notna().values .& df_dirs["Detector 2"].notna().values
println("Total selected events: $(sum(df_mask))")
df = df_raw.sparse.to_dense().values[df_mask]
vals = [df[i][1] for i = 1:length(df)]
p3 = Plots.histogram(rad2deg.(vals .- π / 2), bin=20)
xlabel!("θ of the incident ray")
ylabel!("Counts")
savefig(p3, joinpath(OUT_DIR, "Demo_Histogram.pdf"))

return sim_res
