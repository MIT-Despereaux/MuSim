### This part contains functions used for simulations

using FLoops, BangBang, MicroCollections
using SparseArrays
using DataStructures
using Distributions
using Random
using JLD2

"""
Returns the relative direction and error going from center of T1 to center of T2, in spherical coordinates.
"""
function _relativedir(T1::RectBox{T}, T2::RectBox{T}) where {T<:Real}
    center_1 = T1.position
    center_2 = T2.position
    dir = center_2 - center_1
    (θ, φ) = _cart2unitsph(dir...)
    θ_max = θ
    θ_min = θ
    φ_max = φ
    φ_min = φ
    for x_sign in (-1, 1)
        for y_sign in (-1, 1)
            for z_sign in (-1, 1)
                v = (x_sign * T2.delta_x / 2, y_sign * T2.delta_y / 2, z_sign * T2.delta_z / 2) |> SVector{3,Float64}
                v = _rotate(v, T2.orientation..., 0)
                dir2 = center_2 - center_1 + v
                (θ2, φ2) = _cart2unitsph(dir2...)
                θ_max = max(θ_max, θ2)
                θ_min = min(θ_min, θ2)
                φ_max = max(φ_max, φ2)
                φ_min = min(φ_min, φ2)
            end
        end
    end
    Δθ = θ_max - θ_min
    Δφ = φ_max - φ_min
    return (θ, φ, Δθ, Δφ)
end


"""
Determines the simulation parameters by centering on the first detector.
"""
function _get_ℓ_r(detectors::Vector{RectBox{T}}) where {T<:Real}
    center = detectors[1].position
    (x, y, z) = (detectors[1].delta_x, detectors[1].delta_y, detectors[1].delta_z)
    ℓ = √(max(x * y, x * z, y * z)) * 2
    r = max(x, y, z) * 2
    return ℓ, r
end


"""
Runs the simulation with a hemispherical generating surface. 
Optional parameters θ_range and φ_range denotes the range of solid angles
for the Rays to spawn. The solid angle vector extends from the center of
the hemisphere to its surface. 
By default, the simulation will center on the first detector.
"""
function runhemisim(n_sim::Int, detectors::Vector{T},
    R::Real, ℓ::Real;
    exec=ThreadedEx(), center::Union{NTuple{3,Real},Nothing}=nothing) where {T<:LabObject{<:Real}}
    if center === nothing
        center = detectors[1].position |> Tuple
        optimize = true
    else
        optimize = false
    end
    if optimize
        # Find the relative direction between other detectors and the first detector
        rel_dirs = []
        for i in 2:length(detectors)
            push!(rel_dirs, _relativedir(detectors[1], detectors[i]))
        end
        mix_dist_θ = MixtureModel([Normal(θ > π / 2 ? π - θ : θ, Δθ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_θ = truncated(mix_dist_θ, 0, π / 2)
        mix_dist_φ = MixtureModel([Normal(φ, Δφ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_φ = truncated(mix_dist_φ, 0, 2π)
    else
        mix_dist_θ = nothing
        mix_dist_φ = nothing
    end
    let optimize = optimize, mix_dist_θ = mix_dist_θ, mix_dist_φ = mix_dist_φ, center = center
        @floop exec for i = 1:n_sim
            # Private mutable variables
            @init begin
                crosses = SortedDict{Float64,SVector{3,Float64}}()
                hit_vec = falses(length(detectors))
                i_vec = zeros(length(detectors))
                j_vec = zeros(length(detectors))
                crx = Vector{Union{Missing,Dict}}(missing, length(detectors))
                dir = Vector{Union{Missing,Dict}}(missing, length(detectors))
                angles = zeros(2)
            end
            # Set the hit_vec to false to prepare for a new ray
            hit_vec .*= false
            # Generate a random ray
            if optimize
                θ_sample = rand(mix_dist_θ)
                φ_sample = rand(mix_dist_φ)
                angles = (θ_sample, φ_sample)
                θs = π - θ_sample
                φs = π + φ_sample
            else
                θs = nothing
                φs = nothing
            end
            ray = Ray(R, center, ℓ; θ=θs, φ=φs)
            # Loop through each dectector to fill the pre-allocated vectors
            for (j, d) in enumerate(detectors)
                hit = isthrough!(ray, d, crosses)
                if hit
                    i_vec[j] = i
                    j_vec[j] = j
                    hit_vec[j] = hit
                    points = hcat(values(crosses)...)'
                    # Note the points are vertically concatenated
                    crx[j] = Dict(i => points)
                    dir[j] = Dict(i => Float64[ray.azimuth_agl, ray.polar_agl])
                    empty!(crosses)
                end
            end
            # Hit indices to be filled
            ii = i_vec[hit_vec]
            jj = j_vec[hit_vec]
            # Reduce the accumulators
            # Note that for threaded executions, the initial values will be
            # assigned to the variables (e.g. ii) at the end of the loop
            # This means e.g. ii will be Int[] at some time
            @reduce() do (i_list = Float64[]; ii),
            (j_list = Float64[]; jj),
            (crx_pts = [Dict{Int,Matrix{Float64}}() for i = 1:length(detectors)]; crx),
            (ray_dirs = [Dict{Int,Vector{Float64}}() for i = 1:length(detectors)]; dir),
            (ang_tuple = Float64[]; angles)
                append!(i_list, ii)
                append!(j_list, jj)
                append!!(ang_tuple, angles)
                for (j, c) in enumerate(crx)
                    if !ismissing(c) && !isempty(c)
                        merge!(crx_pts[j], crx[j])
                        merge!(ray_dirs[j], dir[j])
                    end
                end
            end
        end
        results = sparse(i_list, j_list, trues(length(i_list)), n_sim, length(detectors))
        if optimize
            return (results, crx_pts, ray_dirs, mix_dist_θ, mix_dist_φ, ang_tuple)
        else
            return (results, crx_pts, ray_dirs)
        end
    end
end


"""
Runs the simulation with a hemispherical generating surface, outputting only the sparse matrix.
    """
function runhemisimlite(n_sim::Int, detectors::Vector{T},
    R::Real, ℓ::Real;
    exec=ThreadedEx(), center::Union{NTuple{3,Real},Nothing}=nothing) where {T<:LabObject{<:Real}}
    if center === nothing
        center = detectors[1].position |> Tuple
        optimize = true
    else
        optimize = false
    end
    if optimize
        # Find the relative direction between other detectors and the first detector
        rel_dirs = []
        for i in 2:length(detectors)
            push!(rel_dirs, _relativedir(detectors[1], detectors[i]))
        end
        mix_dist_θ = MixtureModel([Normal(θ > π / 2 ? π - θ : θ, Δθ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_θ = truncated(mix_dist_θ, 0, π / 2)
        mix_dist_φ = MixtureModel([Normal(φ, Δφ) for (θ, φ, Δθ, Δφ) in rel_dirs])
        mix_dist_φ = truncated(mix_dist_φ, -π, π)
    else
        mix_dist_θ = nothing
        mix_dist_φ = nothing
    end
    # See https://juliafolds.github.io/FLoops.jl/dev/howto/avoid-box/#avoid-box for the "let" phrase.
    let optimize = optimize, mix_dist_θ = mix_dist_θ, mix_dist_φ = mix_dist_φ, center = center
        @floop exec for i = 1:n_sim
            # Private mutable variables
            @init begin
                crosses = SortedDict{Float64,SVector{3,Float64}}()
                hit_vec = falses(length(detectors))
                i_vec = zeros(length(detectors))
                j_vec = zeros(length(detectors))
                angles = zeros(2)
            end
            # Set the hit_vec to false to prepare for a new ray
            hit_vec .*= false
            # Generate a random ray
            if optimize
                θ_sample = rand(mix_dist_θ)
                φ_sample = rand(mix_dist_φ)
                angles = (θ_sample, φ_sample)
                θs = π - θ_sample
                φs = π + φ_sample
            else
                θs = nothing
                φs = nothing
            end
            ray = Ray(R, center, ℓ; θ=θs, φ=φs)
            # Loop through each dectector to fill the pre-allocated vectors
            for (j, d) in enumerate(detectors)
                hit = isthrough!(ray, d, crosses)
                if hit
                    i_vec[j] = i
                    j_vec[j] = j
                    hit_vec[j] = hit
                    empty!(crosses)
                end
            end
            # Hit indices to be filled
            ii = i_vec[hit_vec]
            jj = j_vec[hit_vec]
            # Reduce the accumulators
            @reduce() do (i_list = Float64[]; ii),
            (j_list = Float64[]; jj),
            (ang_tuple = Float64[]; angles)
                append!!(i_list, ii)
                append!!(j_list, jj)
                append!!(ang_tuple, angles)
            end
        end
        results = sparse(i_list, j_list, trues(length(i_list)), n_sim, length(detectors))
        if optimize
            return (results, mix_dist_θ, mix_dist_φ, ang_tuple)
        else
            return (results)
        end
    end
end


"""
Runs the simulation in given configurations. 
It scans over all the sim_configs. 
The batch_size parameter denotes the number of simulations between caching,
defaults to 1e6.
The cached file blocks contain relevant metadata for the simulation 
and adopt the naming convention "sim_cache_H#.jld2":
H -- Hash of the config dict;
The outputs are vectors of simulation results.
"""
function runexp(output_dir, detectors, sim_configs; batch_size=Int(1e6), lite=true, overwrite=false)
    # Run the simulation using the setup and return raw results
    # Print number of threads
    @printf("Total threads: %d\n", Threads.nthreads())

    # Initilize the result vectors
    # Each entry in the vector is a full simulation output
    # results: size(sim_num, detectors) sparse matrix with trues and falses.
    # crx_pts: size(detectors) vector of dictionary with [event → two cross points of this event] as [key → value] pairs.
    # ray_dirs: size(detectors) vector of dictionary with [event → (θ, φ) of the muon direction].
    results = Vector{Any}(undef, length(sim_configs))

    for (i, config) in enumerate(sim_configs)
        config_hash = hash(config)
        println("Config hash: $(config_hash)")
        f_name = joinpath(output_dir, "sim_cache_H$(config_hash).jld2")
        if !overwrite && isfile(f_name)
            println("Cache found, loading...")
            @load f_name config detectors results_tmp
        else
            results_tmp = []
        end

        total_batch = Int(ceil(config["sim_num"] / batch_size))
        for b in 1:total_batch
            if size(results_tmp)[1] >= b
                continue
            end
            batch_sim_num = min(config["sim_num"] - (b - 1) * batch_size, batch_size)
            center = config["center"]
            if center === nothing
                # Automatically determine ℓ, r if center is not provided
                ℓ, r = _get_ℓ_r(detectors)
            else
                ℓ = config["ℓ"]
                r = config["r"]
            end
            config["ℓ"] = ℓ
            config["r"] = r
            println("\n--- Simulation events $batch_sim_num, batch#$b, config#$i, ℓ=$ℓ r=$r started ---")
            if lite
                push!(results_tmp, runhemisimlite(batch_sim_num, detectors, r, ℓ; center=center))
            else
                push!(results_tmp, runhemisim(batch_sim_num, detectors, r, ℓ; center=center))
            end
            @save f_name config detectors results_tmp
        end
        results[i] = results_tmp
    end

    return results
end


# """
# Helper function that copies the detector setup.
# """
# function copydets(detectors::Vector{T}, n_copies::Int) where {T<:LabObject{<:Real}}
#     dets = []
#     for i in 1:n_copies
#         for (j, d) in enumerate(detectors)
#             idx = length(detectors) * (i - 1) + j
#             new_det = deepcopy(d)
#             # Rotate the detector by the correct angle
#             angle = 2pi / n_copies * (i - 1)
#             pos = new_det.position
#             ori_end = pos + objorient(new_det)
#             new_pos = _rotate(pos, 0, angle, 0)
#             new_ori_end = _rotate(ori_end, 0, angle, 0)
#             new_ori = new_ori_end - new_pos
#             new_det.position = new_pos
#             new_det.orientation = _cart2unitsph(new_ori...)
#             push!(dets, new_det)
#         end
#     end
#     dets = Vector{LabObject{<:Real}}(dets)
#     return dets
# end


# """
# This function converts the simulation results to DataFrame and saves it as pkl files.
# Names of the files are sim_res, sim_crx, sim_dirs, etc.
# Overwrite controls whether to overwrite existing data.
# The function reads the files on disk if they are already there and returns
# a dictionary containing the DataFrames.
# """
# function expio(output_dir, detectors, config; overwrite=false, res=nothing, crx=nothing, dirs=nothing)
#     config_hash = hash(config)
#     println("Config hash: $(config_hash)")
#     # Construct the metadata
#     metadata = "Total number of simulations: $(size(res)[1])\n"
#     n_copies = size(res)[2] ÷ length(detectors)
#     metadata *= "# of copies: $(n_copies)\n"
#     metadata *= "Simulation Configurations:\n"
#     detector_names = String[]
#     for (k, v) in config
#         if k == "detectors"
#             for grp in 1:n_copies
#                 for d in v
#                     metadata *= sprint(show, d)
#                     metadata *= repeat("-", 50)
#                     metadata *= "\n"
#                     push!(detector_names, "$(d.name)_copy$(grp)")
#                 end
#             end
#         else
#             metadata *= "$k -> $v \n"
#         end
#     end
#     println(metadata)

#     sim_res = Dict{String,Any}()
#     # Construct the DataFrames and sort in place
#     if res !== nothing
#         f_name = joinpath(output_dir, "sim_res_H$(config_hash).jld2")
#         if isfile(f_name) && !overwrite
#             println("$(f_name) found, loading...")
#             @load f_name res
#         else
#             @save f_name res
#         end
#         sim_res["sim_res"] = res
#     end

#     if crx !== nothing
#         f_name = joinpath(output_dir, "sim_crx_H$(config_hash).jld2")
#         if isfile(f_name) && !overwrite
#             println("$(f_name) found, loading...")
#             @load f_name crx
#         else
#             @save f_name crx
#         end
#         sim_res["sim_crx"] = crx
#     end

#     if dirs !== nothing
#         f_name = joinpath(output_dir, "sim_dirs_H$(config_hash).jld2")
#         if isfile(f_name) && !overwrite
#             println("$(f_name) found, loading...")
#             @load f_name dirs
#         else
#             @save f_name dirs
#         end
#         sim_res["sim_dirs"] = dirs
#     end

#     # Save the metadata and config to disk
#     # Convert the dict to python dict
#     py_config = Dict{String,Any}(config)
#     # Deletes unsupported PyObject entry
#     pop!(py_config, "detectors")
#     # Open a python IO and save
#     f_name = joinpath(OUT_DIR, "config_H$(config_hash).pkl")
#     @py pickle.dump(py_config, open(f_name, "wb"))
#     open(joinpath(OUT_DIR, "metadata_H$(config_hash).txt"), "w") do io
#         write(io, metadata)
#     end
#     return sim_res
# end


"""
Helper function that takes in the result table and the detector "pairs"
and outputs the inclusive and exclusive geometric factor. 
The pair should be a list of string, and can be of any order.
Requires a dict containing the "detect_order": det name -> column#.
Detector name comes from the uncopied detector list, with convention "Det_?", where ? is the detector name.
"""
function calcgeo(det_order, res, pair, optimize)
    exc_geo_factor = 0
    inc_geo_factor = 0
    sort!(pair)
    total_n = 0
    for res_tup in res
        if typeof(res_tup) <: Tuple
            res_spmat = res_tup[1]
        else
            res_spmat = res_tup
        end
        total_n += size(res_spmat, 1)
        if optimize
            angles = reshape(res_tup[end], 2, :)
            dist_θ = res_tup[end-2]
            dist_φ = res_tup[end-1]
        end
        init_idx = det_order[pair[1]]
        exc_pair_res = res_spmat[:, init_idx]
        inc_pair_res = res_spmat[:, init_idx]
        for (det, col) in det_order
            if det == pair[1]
                continue
            end
            bit_array = view(res_spmat, :, col) |> BitArray
            if det in pair
                inc_pair_res .&= bit_array
                exc_pair_res .&= bit_array
            else
                exc_pair_res = exc_pair_res .> bit_array
            end
        end
        if optimize
            for i in eachindex(inc_pair_res)
                if inc_pair_res[i]
                    local θ = angles[1, i]
                    local φ = angles[2, i]
                    inc_geo_factor += 3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ))
                end
            end
            for i in eachindex(exc_pair_res)
                if exc_pair_res[i]
                    local θ = angles[1, i]
                    local φ = angles[2, i]
                    exc_geo_factor += 3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ))
                end
            end
        else
            exc_geo_factor += sum(exc_pair_res)
            inc_geo_factor += sum(inc_pair_res)
        end
    end
    println("Combination: $(pair)")

    println("Total number of events: $(total_n)")
    exc_res = exc_geo_factor / total_n
    inc_res = inc_geo_factor / total_n
    println("Exclusive beta: $(exc_res)")
    println("Inclusive beta: $(inc_res)")
    return inc_res, exc_res
end


"""
This function processes the input result table and outputs the naive geometric
factors as a dict. It also saves the dict as a pkl file.
"""
function geometricio(output_dir, detectors, linked_dets, res, config; overwrite=false)
    config_hash = hash(config)
    println("Config hash: $(config_hash)")
    # Check the output table
    f_name = joinpath(output_dir, "comb_table_H$(config_hash).csv")
    if isfile(f_name) && !overwrite
        println("$(f_name) found, loading...")
        df_gfac = pd.read_csv(f_name, index_col=0)
        df_gfac = df_gfac.transpose().iloc[0]
        geo_factors = df_gfac.to_dict()
    else
        geo_factors = Dict{String,Float64}()
        det_order = Dict{String,Int}()
        all_dets = []
        linked_pairs = []
        for (i, det) in enumerate(detectors)
            name = split(det.name, "_", limit=2)[end]
            det_order[name] = i
            push!(all_dets, name)
        end
        for p in linked_dets
            linked_tuple = Tuple(split(d.name, "_", limit=2)[end] for d in p)
            push!(linked_pairs, linked_tuple)
        end
        # List all combinations
        betas = combinations(all_dets)
        for β in betas
            # join(β, delim) for concatenation
            (inc_res, exc_res) = calcgeo(det_order, res, β)
            geo_factors["inc_beta_$(join(β, "_"))"] = inc_res * config["ℓ"]^2
            geo_factors["exc_beta_$(join(β, "_"))"] = exc_res * config["ℓ"]^2
        end
        df_gfac = pd.Series(geo_factors)
        df_gfac.to_csv(f_name)
    end
    return geo_factors
end

