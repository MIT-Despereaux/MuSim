### This part stores properties of the generated particles and their interactions.


"""
A mutable ray object (just a line) that could go through multiple LabObjects.
(θ, φ) denotes the unit direction vector, so (0, π) goes up.
"""
@kwdef mutable struct Ray{T<:Real}
    zenith::T # the angle θ in spherical coordinates in radian
    polar::T # angle φ in spherical coordinates in radian
    start_position::MVector{3,Float64} = @MVector [0, 0, 0]
end

"""
Convenience constructor.
"""
function Ray(θ::T, φ::U; start_pos::NTuple{3,<:Real}=(0, 0, 0)) where {T,U}
    args = promote(θ, φ)
    return Ray(zenith=args[1], polar=args[2], start_position=MVector{3,Float64}(start_pos))
end


"""
Returns a random number with pdf: Normalized(cos²(x)sin(x)) ∈ [x_min, x_max). We use inverse sampling for speed.
!!!Note: sampling on a sphere requires taking into account of the Jacobian sin(θ)dθdφ = dΩ,
The constraint is within any solid angle dΩ, the random point counts are constant: so ∫f(θ,φ)dΩ = 1 => pdf = f(θ,φ)sin(θ).
"""
function _randcos2(x_min::Real=π / 2, x_max::Real=π)
    @assert 0 <= x_min <= x_max <= π
    # Note (-rand)^(1/3) doesn't work; it gives a complex root (not sure why Julia dev does this...)
    N = 1 / 3 * cos(x_min)^3 - 1 / 3 * cos(x_max)^3
    u = rand() * N
    res = (1 / 3 * cos(x_min)^3 - u)
    return acos(cbrt(3 * res))
end


"""
Returns the distribution f(θ) = cos³(θ)sin(θ). The extra cos factor is required for sampling on a horizontal plane with rays coming
non-vertically (for the angle between the flux vector and norm vector of the surface).
"""
function _randcos3(x_min::Real=π / 2, x_max::Real=π)
    @assert π / 2 <= x_min <= x_max <= π
    # For x_min smaller than π/2 cos(x)^4 will not be monotonic which causes trouble
    N = 1 / 4 * cos(x_min)^4 - 1 / 4 * cos(x_max)^4
    u = rand() * N
    res = (1 / 4 * cos(x_min)^4 - u)
    return acos(-(4 * res)^(1 / 4))
end


"""
Generates the uniform distribution on the sphere for variable θ. f(θ) = sin(θ). Note the range is [0, π).
"""
_randsin() = acos(1 - 2rand())


"""
Randomly modifies a ray with start uniform x, y ∈ (-max_x, max_x) U (-max_y, max_y), with max_bounds = (max_x, max_y)
The center is given as a tuple.
If not given, the θ and φ are drawn from cos³sin and uniform distributions respectively.
"""
function modifyray!(r::Ray{<:Real}, max_bounds::NTuple{2,Real}, center::Union{NTuple{3,Real},SVector{3}};
    θ::Union{Nothing,Real}=nothing,
    φ::Union{Nothing,Real}=nothing)::Nothing
    x = rand() * 2max_bounds[1] - max_bounds[1] + center[1]
    y = rand() * 2max_bounds[2] - max_bounds[2] + center[2]
    if θ === nothing
        # Generate by the distribution
        θ = _randcos3()
    end
    if φ === nothing
        φ = rand() * 2π
    end
    r.zenith = θ
    r.polar = φ
    (r.start_position.x, r.start_position.y, r.start_position.z) = (x, y, center[3])
    return nothing
end


"""
Randomly modifies a ray on the upper hemisphere with radius R. The direction is along the normal vector.
If not given, the θ and φ are drawn from cos²sin and uniform distributions respectively.
If not given, the x and y are drawn from uniform distribution within [-ℓ/2, ℓ/2].
See http://arxiv.org/abs/1912.05462.
"""
function modifyray!(r::Ray{<:Real}, R::Real, center::Union{NTuple{3,Real},SVector{3}}, ℓ::Real;
    θ::Union{Nothing,Real}=nothing,
    φ::Union{Nothing,Real}=nothing,
    x::Union{Nothing,Real}=nothing,
    y::Union{Nothing,Real}=nothing)::Nothing
    if θ === nothing
        # Generate by the distribution
        θ = _randcos2()
    end
    if φ === nothing
        φ = rand() * 2π
    end
    start_θ = π - θ
    start_φ = φ + π
    x_vec = SA_F64[sin(start_θ)*cos(start_φ), sin(start_θ)*sin(start_φ), cos(start_θ)] * R
    if x === nothing
        x = rand() * ℓ - ℓ / 2
    end
    if y === nothing
        y = rand() * ℓ - ℓ / 2
    end
    @assert -ℓ / 2 <= x <= ℓ / 2 && -ℓ / 2 <= y <= ℓ / 2
    pos = @SVector [x, y]
    θ_hat = SA_F64[cos(start_θ)*cos(start_φ), cos(start_θ)*sin(start_φ), -sin(start_θ)]
    φ_hat = SA_F64[-sin(start_φ), cos(start_φ), 0]
    dx_vec = pos[1] * θ_hat + pos[2] * φ_hat
    x_vec += (dx_vec + (center |> SVector))
    r.zenith = θ
    r.polar = φ
    (r.start_position.x, r.start_position.y, r.start_position.z) = (x_vec.x, x_vec.y, x_vec.z)
    return nothing
end


"""
Returns the direction vector in Cartesian coordinates.
"""
function raydir(r::Ray{<:Real})::SVector{3,Float64}
    return _unitsph2cart(r.zenith, r.polar)
end


"""
Returns the intersection point a ray and an inifinite plane. 
It's a tuple of intersection time and intersection position vector.
"""
function rayplaneint(dir::SVector{3,Float64}, r_pos::SVector{3,Float64}, plane_norm::SVector{3,Float64}, plane_pos::SVector{3,Float64})::Union{Tuple,Nothing}
    β = dir ⋅ plane_norm
    if β == 0
        return nothing
    else
        t = (plane_pos - r_pos) ⋅ plane_norm / β
        tup = (t, r_pos + t * dir)
        return tup
    end
end


"""
Sample from the Continuous variable from MCIntegration.
"""
function randContinuous(x::MCIntegration.Continuous{G}; n::Int=1, inv_itp=nothing) where {G<:AbstractVector}
    if inv_itp === nothing
        xs = range(0, stop=1; length=length(x.grid))
        inv_itp = Interpolations.scale(interpolate(x.grid, BSpline(Linear())), xs)
    end
    if n > 1
        return inv_itp.(rand(n))
    else
        return inv_itp(rand())
    end
end


"""
Interpolate the Continuous variable from MCIntegration.
"""
function interpContinuous(x::MCIntegration.Continuous{G}) where {G<:AbstractVector}
    xs = collect(range(0, stop=1; length=length(x.grid)))
    itp = interpolate(x.grid, xs, FritschButlandMonotonicInterpolation())
    return itp
end


"""
Calculate the PDF of the Continuous variable from MCIntegration.
"""
function probContinuous(x::MCIntegration.Continuous{G}, x0::Union{Real,AbstractArray}) where {G<:AbstractVector}
    itp = interpContinuous(x)
    return map(x -> Interpolations.gradient(itp, x)[1], x0)
end


"""
Define the Muon flux pdf problem, in the form of a struct.
"""
@kwdef struct μPDF
    # Parameter units: E₀ [GeV], ϵ [GeV]
    ρ::Vector{Float64} = [0.102573, -0.068287, 0.958633, 0.0407253, 0.817285]
    ϵ::Vector{Float64} = [115, 850]
    E₀::Float64 = 3.64
    γ::Float64 = 2.7
end


"""
Define the loglikelihood of μPDF, given E and θ.
See https://arxiv.org/pdf/1509.06176.pdf.
Note the final pdf needs to be multiplied by the Jacobian sin(θ).
"""
function (p::μPDF)(params::Union{NamedTuple,AbstractVector}; return_log=true, return_jac=true)
    if params isa AbstractVector
        E, θ = params
    elseif params isa NamedTuple
        E = params[:E]
        θ = params[:θ]
    end
    z = @. √((cos(θ)^2 + (p.ρ[1])^2 + p.ρ[2] * cos(θ)^(p.ρ[3]) + p.ρ[4] * cos(θ)^(p.ρ[5])) / (1 + (p.ρ[1]^2) + p.ρ[2] + p.ρ[4]))
    q = @. (1 + p.E₀ / (E * z^1.29))
    flux_pdf = @. (E * q)^(-p.γ) * (1 / (1 + 1.1 * E * z / p.ϵ[1]) + 0.054 / (1 + 1.1 * E * z / p.ϵ[2]))
    if return_log
        return @. log(flux_pdf * sin(θ))
    else
        if return_jac
            return @. flux_pdf * sin(θ)
        else
            return flux_pdf
        end
    end
end


struct μPDFSettings
    # E and θ range
    # Parameter units: E [GeV], θ [rad]
    E_range::Tuple{Float64,Float64}
    θ_range::Tuple{Float64,Float64}
    trans::Any
    rng::AbstractRNG
end


"""
Constructor for the μPDFSettings. Use this for initializing the μPDFSettings.
"""
function μPDFSettings(; E_range=(1.0, 1000.0), θ_range=(0.0, π / 2), seed::Union{Int,Nothing}=nothing)
    t = as((E=as(Real, E_range[1], E_range[2]), θ=as(Real, θ_range[1], θ_range[2])))
    if seed === nothing
        rng = Random.GLOBAL_RNG
    else
        rng = Random.MersenneTwister(seed)
    end
    return μPDFSettings(E_range, θ_range, t, rng)
end


"""
This function provides a batch sampler of μPDF by using HMC.
The "μPDFSettings" input is a struct containing the HMCSettings.
The output will be a 2 by N matrix.
"""
function drawsamples(p::μPDF, s::μPDFSettings; N::Int=1000)
    trans = s.trans
    trans_logℓ = TransformedLogDensity(trans, p)
    ∇ℓ = ADgradient(:ForwardDiff, trans_logℓ)
    results = mcmc_with_warmup(s.rng, ∇ℓ, N; reporter=NoProgressReport())
    res = TransformVariables.transform.(trans, eachcol(results.posterior_matrix))
    res = hcat([getindex.(res, i) for i in 1:length(res[1])]...)'
    return res
end


# --- Scratch ---


# """
# A sampler for an arbitrary univariate distribution.
# It takes in a function with possible lower and upper bounds and performs inverse sampling.
# If provided, it also saves the results into the cache vector.
# If the total number of provided cache exceeds (by default) 1e5 the function will draw from the saved cache.
# """
# function unisampler!(f::Function; low_bound=-Inf, upp_bound=Inf, cache=nothing, cache_size=Int(1e5))
#     if cache !== nothing && length(cache) > cache_size
#         return rand(cache)
#     end
#     (total, err) = quadgk(f, low_bound, upp_bound)
#     x = rand()
#     eqn = y -> unicdf(f, y, low_bound=low_bound, upp_bound=upp_bound, total=total) - x
#     res = find_zero(eqn, (low_bound, upp_bound), A42())
#     if cache !== nothing
#         push!(cache, res)
#     end
#     return res
# end


# """
# A sampler for an arbitrary univariate distribution using its cdf.
# NOTE: THIS METHOD IS NOT STABLE. 
# """
# function unisamplercdf!(ecdf; low_bound=-Inf, upp_bound=Inf, cache=nothing, cache_size=Int(1e5))
#     if cache !== nothing && length(cache) > cache_size
#         return rand(cache)
#     end
#     x = rand()
#     eqn = y -> ecdf(y) - x
#     res = find_zero(eqn, (low_bound, upp_bound), Bisection())
#     if cache !== nothing
#         push!(cache, res)
#     end
#     return res
# end


# """
# Return the cdf of a univariate distribution. Note this is not normalized.
# """
# function unicdf(f::Function, x::Real; low_bound=-Inf, upp_bound=Inf, total=1.0)
#     @assert low_bound <= x <= upp_bound "$x exceeding function support [$low_bound, $upp_bound]."
#     (res, err) = quadgk(f, low_bound, x)
#     res /= total
#     return res
# end
