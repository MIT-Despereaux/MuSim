### This part contains functions used for simulations
using HCubature.QuadGK

@kwdef struct SimOutput{T<:Real}
    col_names::Vector{String}
    mat::Matrix{T}
    sim_num::Int
    μPDFSettings::Any = nothing
end


"""
Runs the simulation with a hemispherical generating surface.
The output contains the following sparse matrices: traversed distances; Muon energies; cross points; ray directions;
This is the "vanilla" version without any optimizations.
"""
function runhemisim(N::Int, detectors::Vector{T}, R::Real, ℓ::Real;
    exec=ThreadedEx(),
    center::Union{NTuple{3,Real},SVector{3}}=(0, 0, 0),
    s::μPDFSettings=μPDFSettings(),
    cos2::Bool=false
) where {T<:LabObject{<:Real}}
    # Generate the HMC Eμ and θ samples
    p = μPDF()
    samples = drawsamples(p, s, N=N)
    # See https://juliafolds.github.io/FLoops.jl/dev/howto/avoid-box/#avoid-box for the "let" phrase.
    let center = center, cos2 = cos2
        @floop exec for i = 1:N
            # Private mutable variables
            @init begin
                ray = Ray(0.0, 0.0)
                crosses = SortedDict{Float64,SVector{3,Float64}}()
                hits = Vector{Vector{Float64}}(undef, length(detectors))
                hits_selection = falses(length(detectors))
            end
            hits_selection .*= false
            # Generate a random ray
            if cos2
                θ̃ = nothing
            else
                θ̃ = samples[2, i]
            end
            modifyray!(ray, R, center, ℓ; θ=θ̃)
            # Loop through each dectector to fill the pre-allocated vectors
            for (j, d) in enumerate(detectors)
                through = isthrough!(ray, d, crosses)
                if through
                    hits_selection[j] = true
                    single_hit = Float64[]
                    push!(single_hit, i)
                    push!(single_hit, j)
                    d = norm(last(crosses)[2] - first(crosses)[2])
                    push!(single_hit, d)
                    push!(single_hit, samples[1, i])
                    append!(single_hit, first(crosses)[2])
                    append!(single_hit, last(crosses)[2])
                    push!(single_hit, ray.zenith, ray.polar)
                    empty!(crosses)
                    hits[j] = single_hit
                end
            end
            # Reduce the accumulators
            selected_hits = @view hits[hits_selection]
            @reduce() do (total_hits = EmptyVector{Float64}(); selected_hits)
                append!!(total_hits, selected_hits)
            end
        end
        col_names = ["row",
            "col",
            "dist",
            "energy",
            "x0",
            "y0",
            "z0",
            "x1",
            "y1",
            "z1",
            "zenith",
            "polar"]
        res = SimOutput(col_names, vcat(total_hits'...), N, s)
        return res
    end
end


"""
Runs the simulation with a hemispherical generating surface.
The output contains the following sparse matrices: traversed distances; Muon energies;
Note the T/F table can be inferred from the sparse matrices.
This is the "vanilla" version without any optimizations.
"""
function runhemisimlite(N::Int, detectors::Vector{T}, R::Real, ℓ::Real;
    exec=ThreadedEx(),
    center::Union{NTuple{3,Real},SVector{3}}=(0, 0, 0),
    s::μPDFSettings=μPDFSettings(),
    cos2::Bool=false
) where {T<:LabObject{<:Real}}
    # Generate the HMC Eμ and θ samples
    p = μPDF()
    samples = drawsamples(p, s, N=N)
    # See https://juliafolds.github.io/FLoops.jl/dev/howto/avoid-box/#avoid-box for the "let" phrase.
    let center = center, cos2 = cos2
        @floop exec for i = 1:N
            # Private mutable variables
            @init begin
                ray = Ray(0.0, 0.0)
                crosses = SortedDict{Float64,SVector{3,Float64}}()
                hit_vec = falses(length(detectors))
                i_vec = zeros(length(detectors))
                j_vec = zeros(length(detectors))
                dist = zeros(length(detectors))
                energy = zeros(length(detectors))
            end
            # Set the hit_vec to false to prepare for a new ray
            hit_vec .*= false
            # Generate a random ray
            if cos2
                θ̃ = nothing
            else
                θ̃ = samples[2, i]
            end
            modifyray!(ray, R, center, ℓ; θ=θ̃)
            # Loop through each dectector to fill the pre-allocated vectors
            for (j, d) in enumerate(detectors)
                hit = isthrough!(ray, d, crosses)
                if hit
                    i_vec[j] = i
                    j_vec[j] = j
                    hit_vec[j] = true
                    d = norm(last(crosses)[2] - first(crosses)[2])
                    dist[j] = d
                    energy[j] = samples[1, i]
                    empty!(crosses)
                end
            end
            # Hit indices to be filled
            ii = i_vec[hit_vec]
            jj = j_vec[hit_vec]
            dd = dist[hit_vec]
            ee = energy[hit_vec]
            # Reduce the accumulators
            @reduce() do (i_list = Float64[]; ii),
            (j_list = Float64[]; jj),
            (dist_list = Float64[]; dd),
            (energy_list = Float64[]; ee)
                append!!(i_list, ii)
                append!!(j_list, jj)
                append!!(dist_list, dd)
                append!!(energy_list, ee)
            end
        end
        dist_table = sparse(i_list, j_list, dist_list, N, length(detectors))
        # energy_table = sparse(i_list, j_list, energy_list, N, length(detectors))
        # res = SimOutput(dist_table, energy_table, nothing, nothing, s)
        return dist_table
    end
end


"""
Runs one simulation in a given configuration. 
The cached file adopts the naming convention "sim_cache_H#.jld2":
H -- Hash of the config dict;
"""
function runexp(output_dir, config;
    lite=true,
    overwrite=false,
    mc_config::Union{MCIntegration.Configuration,Nothing}=nothing
)
    # Run the simulation using the setup and return raw results
    # Print number of threads
    @printf("Total threads: %d\n", Threads.nthreads())

    detectors = config["detectors"]
    center = config["center"]
    ℓ = config["ℓ"]
    R = config["R"]
    cos2 = config["cos2"]
    config_hash = hash(config)
    sim_num = config["sim_num"]
    println("Config hash: $(config_hash)")
    f_name = joinpath(output_dir, "sim_cache_H$(config_hash).jld2")
    if !overwrite && isfile(f_name)
        println("Cache found, loading...")
        @load f_name result
    else
        println("--- Simulation events $sim_num, ℓ=$ℓ, R=$R started ---")
        if lite
            result = runhemisimlite(batch_sim_num, detectors, R, ℓ; center=center, cos2=cos2)
        else
            result = runhemisim(batch_sim_num, detectors, R, ℓ; center=center, cos2=cos2)
        end
        @save f_name result
    end
    return result
end


"""
Helper function that takes in the result table corresponding
to a member of sim_config and the detector "combinations",
and outputs the inclusive and exclusive geometry factor (β). 
Output: the β factors are in # of hits, and the multiplication factor is provided at the end.
The combination should be a list of string, and can be of any order.
Requires a dict containing the "detect_order": det name -> column#.
"""
function calculateβ(config, det_order, res_vec, comb)
    inclusive_β = 0
    exclusive_β = 0
    sort!(comb)
    total_n = 0
    for sim_out in res_vec
        # Construct the T/F table
        I = sim_out.mat[:, 1]
        J = sim_out.mat[:, 2]
        res_spmat = sparse(I, J, trues(length(I)), sim_out.sim_num, length(det_order))
        total_n += sim_out.sim_num
        init_idx = det_order[comb[1]]
        inc_comb_res = res_spmat[:, init_idx]
        exc_comb_res = res_spmat[:, init_idx]
        for (det, col) in det_order
            if det == comb[1]
                continue
            end
            bit_array = view(res_spmat, :, col) |> BitArray
            if det in comb
                exc_comb_res .&= bit_array
                inc_comb_res .&= bit_array
            else
                exc_comb_res = exc_comb_res .> bit_array
            end
        end
        inclusive_β += sum(inc_comb_res)
        exclusive_β += sum(exc_comb_res)
    end
    println("Combination: $(comb)")

    println("Total number of events: $(total_n)")
    inclusive_β_err = sqrt(inclusive_β)
    exclusive_β_err = sqrt(exclusive_β)
    mult = 1 / total_n * config["ℓ"]^2

    return inclusive_β * mult, inclusive_β_err * mult, exclusive_β * mult, exclusive_β_err * mult, inclusive_β, exclusive_β
end


# """
# The MC counterpart of calculating the β factors.
# """
# function calculateβ_MC(config, comb)
#     println("Combination: $(comb)")
#     dets = config["detectors"]
#     included_dets = Vector{RectBox{Float64}}()
#     for d in dets
#         if d.name in comb
#             push!(included_dets, d)
#         end
#     end
#     missing_dets = setdiff(dets, included_dets)
#     inclusive_β_res = analytic_R(included_dets)
#     exclusive_β_res = analytic_R(included_dets, excluded_detectors=missing_dets, config=inclusive_β_res.config)
#     println("Inclusive beta: $(inclusive_β_res.mean) ± $(inclusive_β_res.stdev)")
#     println("Exclusive beta: $(exclusive_β_res.mean) ± $(exclusive_β_res.stdev)")
#     return inclusive_β_res.mean, inclusive_β_res.stdev, exclusive_β_res.mean, exclusive_β_res.stdev
# end


"""
This function processes the input result table and outputs the geometric
factors as a dict. It also saves the dict as a csv file.
"""
function βio(output_dir, res_vec, config; savefile=false, overwrite=false)
    config_hash = hash(config)
    println("βio config hash: $(config_hash)")
    # Check the output table
    f_name = joinpath(output_dir, "comb_table_H$(config_hash).csv")
    if !overwrite && isfile(f_name)
        println("$(f_name) found, loading...")
        geo_factors = CSV.File(f_name) |> DataFrame
    else
        gf_name = String[]
        gf_β = Float64[]
        gf_β_err = Float64[]
        gf_n_hits = Int64[]
        gf_analytic = Float64[]

        detectors = config["detectors"]
        det_order = Dict{String,Int}(d.name => i for (i, d) in enumerate(detectors))
        # List all combinations
        combs = combinations(collect(keys(det_order)))
        settings = res_vec[1].μPDFSettings
        p = μPDF()
        f = x -> (p)(x; return_log=false)
        f_vert = x -> (p)([x, 0]; return_log=false, return_jac=false)
        totalI = hcubature(f, (settings.E_range[1], settings.θ_range[1]), (settings.E_range[2], settings.θ_range[2]))[1]
        vertI = quadgk(f_vert, settings.E_range[1], settings.E_range[2])[1]
        factor = totalI / vertI * 3
        println("Factor: $(factor)")
        for c in combs
            sort!(c)
            # join(β, delim) for concatenation
            (inc_β, inc_β_err, exc_β, exc_β_err, inc_n_hits, exc_n_hits) = calculateβ(config, det_order, res_vec, c)
            inc_β *= factor
            inc_β_err *= factor
            exc_β *= factor
            exc_β_err *= factor
            println("Inclusive beta: $(inc_β) ± $(inc_β_err)")
            println("Exclusive beta: $(exc_β) ± $(exc_β_err)")
            println("Inclusive relative err (%): $(100.0 / sqrt(inc_n_hits))")
            if length(c) == 1
                d = detectors[det_order[c[1]]]
                analytic = analytic_R(d)
                analytic /= (2π / 3)
                println("Analytic rate: $(analytic)")
                push!(gf_analytic, analytic)
                push!(gf_analytic, 0.0)
            else
                append!(gf_analytic, [0.0, 0.0])
            end
            push!(gf_name, "inc_beta_$(join(c, "_"))")
            push!(gf_name, "exc_beta_$(join(c, "_"))")
            push!(gf_β, inc_β)
            push!(gf_β, exc_β)
            push!(gf_β_err, inc_β_err)
            push!(gf_β_err, exc_β_err)
            push!(gf_n_hits, inc_n_hits)
            push!(gf_n_hits, exc_n_hits)
        end
        geo_factors_dict = OrderedDict("combination" => gf_name,
            "beta" => gf_β,
            "beta_err" => gf_β_err,
            "n_hits" => gf_n_hits,
            "n_total" => repeat([config["sim_num"]], length(gf_name)),
            "factor" => repeat([factor], length(gf_name)),
            "analytic" => gf_analytic)
        geo_factors = DataFrame(geo_factors_dict)
        println(geo_factors)
        if savefile
            sort!(geo_factors, "combination")
            CSV.write(f_name, geo_factors)
        end
    end
    return geo_factors
end


# """
# MC Counterpart of βio.
# """
# function βio_MC(output_dir, config; savefile=false, overwrite=false)
#     config_hash = hash(config)
#     println("βio config hash: $(config_hash)")
#     # Check the output table
#     f_name = joinpath(output_dir, "comb_table_MC_H$(config_hash).csv")
#     if !overwrite && isfile(f_name)
#         println("$(f_name) found, loading...")
#         geo_factors = CSV.File(f_name) |> Dict{String,Float64}
#     else
#         geo_factors = Dict{String,Float64}()
#         detectors = config["detectors"]
#         det_order = Dict{String,Int}(d.name => i for (i, d) in enumerate(detectors))
#         # List all combinations
#         combs = combinations(collect(keys(det_order)))
#         for c in combs
#             sort!(c)
#             # join(β, delim) for concatenation
#             (inc_β, inc_β_err, exc_β, exc_β_err) = calculateβ_MC(config, c)
#             println("Inclusive beta: $(inc_β) ± $(inc_β_err)")
#             println("Exclusive beta: $(exc_β) ± $(exc_β_err)")
#             println("Inclusive relative err (%): $(100.0 * inc_β_err/ inc_β)")
#             geo_factors["inc_beta_$(join(c, "_"))"] = inclusive_β
#             geo_factors["exc_beta_$(join(c, "_"))"] = exclusive_β
#         end
#         if savefile
#             geo_factors = sort(collect(geo_factors), by=x -> x[1])
#             CSV.write(f_name, geo_factors)
#         end
#     end
#     return geo_factors
# end


# """
# This function saves simulation outputs to jld2 files.
# Names of the files are sim_res, sim_crx, sim_dirs, etc.
# Overwrite controls whether to overwrite existing data.
# """
# function expio(output_dir, config; overwrite=false, res=nothing, crx=nothing, dirs=nothing)
#     config_hash = hash(config)
#     println("Config hash: $(config_hash)")
#     # Construct the metadata
#     metadata = ""
#     metadata *= "Simulation Configurations:\n"
#     detector_names = String[]
#     for (k, v) in config
#         if k == "detectors"
#             for d in v
#                 metadata *= sprint(show, d)
#                 metadata *= repeat("-", 50)
#                 metadata *= "\n"
#                 push!(detector_names, "$(d.name)")
#             end
#         else
#             metadata *= "$k -> $v \n"
#         end
#     end
#     println(metadata)

#     sim_res = Dict{String,Any}()
#     if res !== nothing
#         println(res)
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

#     # Save the metadata to disk
#     f_name = joinpath(output_dir, "metadata_H$(config_hash).txt")
#     if !isfile(f_name) || overwrite
#         open(f_name, "w") do io
#             write(io, metadata)
#         end
#     end

#     return sim_res
# end


# --- Scratch ---


# """
# Compose the simulation assuming the first "detector" is the chip.
# This function aims to output the "correct" geometric factor.
# """
# function composeβ(output_dir, detectors::Vector{RectBox{T}}, n_sim::Int; overwrite=false, include_chip=true) where {T<:Real}
#     config = Dict{String,Any}("sim_num" => n_sim)
#     config["detectors"] = detectors
#     config_hash = hash(config)
#     println("composeβ config hash: $(config_hash)")
#     # Check the output table
#     f_name = joinpath(output_dir, "comb_table_H$(config_hash).csv")
#     if isfile(f_name) && !overwrite
#         println("$(f_name) found, loading...")
#         βs = CSV.File(f_name) |> Dict{String,Float64}
#     else
#         βs = Dict{String,Float64}()
#         for i in eachindex(detectors)
#             dets = circshift(detectors, -i)
#             config["detectors"] = dets
#             config["center"] = nothing
#             if include_chip
#                 if i == length(detectors)
#                     config["sim_num"] = n_sim * 10
#                 end
#             end
#             res, sim_configs = runexp(output_dir, [config])
#             merge!(βs, βio(output_dir, res[1], sim_configs[1]))
#         end
#         config["sim_num"] = n_sim
#         for i in eachindex(detectors)
#             dets = circshift(detectors, -i)
#             config["detectors"] = dets
#             config["center"] = dets[1].position |> Tuple
#             ℓ, r = _get_ℓ_r(dets)
#             config["ℓ"] = ℓ
#             config["r"] = r
#             if include_chip
#                 if i == length(detectors)
#                     config["sim_num"] = n_sim * 10
#                 end
#             end
#             res, sim_configs = runexp(output_dir, [config])
#             merge!(βs, βio(output_dir, res[1], sim_configs[1]; first_only=true))
#         end
#     end
#     βs = sort(collect(βs), by=x -> x[1])
#     println("Saving βs to $(f_name)...")
#     println("βs: $(βs)")
#     CSV.write(f_name, βs)
#     return (config_hash, βs)
# end


# """
# Returns the relative direction and error going from center of T1 to center of T2, in spherical coordinates.
# !!! NEEDS FIX: FROM EDGE OF T1 to EDGE OF T2
# """
# function _relativedir(T1::RectBox{T}, T2::RectBox{T}) where {T<:Real}
#     center_1 = T1.position
#     center_2 = T2.position
#     dir = center_2 - center_1
#     (θ, φ) = _cart2unitsph(dir...)
#     θ_max = θ
#     θ_min = θ
#     φ_max = φ
#     φ_min = φ
#     for x_sign_2 in (-1, 1)
#         for y_sign_2 in (-1, 1)
#             for z_sign_2 in (-1, 1)
#                 for x_sign_1 in (-1, 1)
#                     for y_sign_1 in (-1, 1)
#                         for z_sign_1 in (-1, 1)
#                             v2 = (x_sign_2 * T2.delta_x / 2, y_sign_2 * T2.delta_y / 2, z_sign_2 * T2.delta_z / 2) |> SVector{3,Float64}
#                             v2 = _rotate(v2, T2.orientation..., 0)
#                             v1 = (x_sign_1 * T1.delta_x / 2, y_sign_1 * T1.delta_y / 2, z_sign_1 * T1.delta_z / 2) |> SVector{3,Float64}
#                             v1 = _rotate(v1, T1.orientation..., 0)
#                             dir2 = center_2 + v2 - (center_1 + v1)
#                             (θ2, φ2) = _cart2unitsph(dir2...)
#                             θ_max = max(θ_max, θ2)
#                             θ_min = min(θ_min, θ2)
#                             φ_max = max(φ_max, φ2)
#                             φ_min = min(φ_min, φ2)
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     Δθ = θ_max - θ_min
#     Δφ = φ_max - φ_min
#     return (θ, φ, Δθ, Δφ)
# end


# """
# Determines the simulation parameters by centering on the first detector.
# """
# function _get_ℓ_r(detectors::Vector{RectBox{T}}) where {T<:Real}
#     (x, y, z) = (detectors[1].delta_x, detectors[1].delta_y, detectors[1].delta_z)
#     ℓ = √(max(x * y, x * z, y * z)) * 5
#     r = 100
#     return ℓ, r
# end


# """
# Helper function that samples the detector positions in accordance of "fuzzy" detectors.
# "δr" is the 1σ position error.
# """
# function gendetectorpos(detectors::Vector{RectBox{T}}, δr::Real; n::Int=10, seed::Union{Nothing,Int}=nothing) where {T<:Real}
#     generated_detectors = []
#     if seed !== nothing
#         Random.seed!(seed)
#     end
#     for i in 1:n
#         dets = deepcopy(detectors)
#         for d in dets
#             d.position += rand(Normal(0, δr), 3)
#         end
#         push!(generated_detectors, dets)
#     end
#     return generated_detectors
# end
