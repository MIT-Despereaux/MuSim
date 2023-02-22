using StatsBase, Distributions, Interpolations, Plots, LinearAlgebra, Random

# %%
# Creating original samples...
C = Cauchy(0, 0.01)
mean_V = [0, 0]
diag_V = [0.1, 0.5]
N = MvNormal(mean_V, diag_V)
x_range = 0:0.005:1
y_range = -1:0.01:1


# %%
# Define the log likelihood function
using LogDensityProblems, LogDensityProblemsAD, TransformVariables, TransformedLogDensities, DynamicHMC

function (D::ContinuousUnivariateDistribution)(θ::NamedTuple)
    return logpdf(D, θ.x)
end

function (MvD::MultivariateDistribution)(θ::NamedTuple)
    return logpdf(MvD, [θ.x, θ.y])
end


# %%
# Define the constraint in the input space
trans_C = as((x=asℝ₊,))
trans_N = as((x=asℝ₊, y=asℝ))

trans_logℓ_C = TransformedLogDensity(trans_C, C)
trans_logℓ_N = TransformedLogDensity(trans_N, N)

∇ℓ_C = ADgradient(:ForwardDiff, trans_logℓ_C)
∇ℓ_N = ADgradient(:ForwardDiff, trans_logℓ_N)

# %%
results_C = mcmc_with_warmup(Random.default_rng(), ∇ℓ_C, 1000)
results_N = mcmc_with_warmup(Random.default_rng(), ∇ℓ_N, 1000)


# %%
# Plot the pdf of Cauchy
plot(x_range, pdf.(C, x_range) * 2, label="Cauchy")
# Plot the posterior after the transformation
post_x = first.(transform.(trans_C, eachcol(results_C.posterior_matrix)))
# Plot the histogram in log
p1 = histogram!(post_x, label="Posterior", bins=x_range, normed=true, alpha=0.3)
display(p1)

# %%
x_range = 0:0.05:1
y_range = -1:0.05:1

# Plot the posterior after the transformation
post_x = [i.x for i in (transform.(trans_N, eachcol(results_N.posterior_matrix)))]
post_y = [i.y for i in (transform.(trans_N, eachcol(results_N.posterior_matrix)))]
# Plot the histogram in log
p1 = histogram2d(post_x, post_y, label="Posterior", bins=(x_range, y_range), normed=true, color=:plasma)

pdf_MvNormal = [pdf(N, [x, y]) for y in y_range, x in x_range]
# Create the heatmap plot
heatmap!(x_range, y_range, pdf_MvNormal * 2,
    aspect_ratio=:equal,
    xlabel="x",
    ylabel="y",
    color=:plasma,
    alpha=0.3,
    title="2D Distribution")


# %%
using MuSim, BenchmarkTools

p = μPDF()
s = μPDFSettings()
drawsamples(p, s, N=10000)

# %%
# Plot the theoretical and histogram
E_range = LinRange(s.E_range..., 101)
θ_range = LinRange(s.θ_range..., 101)

# %%
# Plot the histogram in log
p1 = histogram2d(res[1, :], res[2, :], label="Posterior", bins=(E_range, θ_range), normed=true, color=:plasma, alpha=0.1)

pdf_theo = exp.(p.([(E=e, θ=t) for t in θ_range, e in E_range])) / 0.01398009970294201
# Create the heatmap plot
contour!(E_range, θ_range, pdf_theo,
    xlabel="x",
    ylabel="y",
    levels=100,
    color=:plasma,
    title="2D Distribution")

# %%
# Plot the traces
# plot(res[1, :], label="E", yscale=:log10)
function E_pdf(p, E)
    return (p.E₀ + E)^(-p.n) * (1 + E / p.ϵ)^(-1) * 50
end

function θ_pdf(p, θ)
    D = √(p.Rd_ratio^2 * cos(θ)^2 + 2p.Rd_ratio + 1) - p.Rd_ratio * cos(θ)
    return D^(-p.n + 1)
end

# histogram(res[1, :], label="E", yscale=:log10, legend=false, normalize=:pdf)
# plot!(E_range, [E_pdf(p, e) for e in E_range], label="E theoretical", yscale=:log10, legend=false)
histogram(res[2, :], label="θ", yscale=:log10, legend=false, normalize=:pdf)
plot!(θ_range, [θ_pdf(p, t) for t in θ_range], label="θ theoretical", yscale=:log10, legend=false)

# %%
itp = linear_interpolation(edges, (eCDF(edges)), extrapolation_bc=Flat())

# %%

result = []
for n in 1:10000
    push!(result, unisamplercdf!(itp; low_bound=floor(minimum(eCDF)), upp_bound=ceil(maximum(eCDF))))
end

grad = []
for e in edges
    push!(grad, gradient(itp, e)[1])
end

plot(edges, grad)
histogram!(obs, bins=100, normalize=:pdf, alpha=0.3)
histogram!(result, bins=100, normalize=:pdf, alpha=0.3)

# %%
X = MCIntegration.Continuous(0.0, 2.0)
res = integrate((X, c) -> (X[1]^2 + X[2]^2 < 1.5); var=X, dof=2)
println(res.mean)
plot(res.config.var[1].grid)

# Cdf --> Pmf..
# pmf_obs = Float64[]
# for i in 1:length(obs)
#     if i == 1
#         append!(pmf_obs, cdf_obs[1])
#     else
#         append!(pmf_obs, cdf_obs[i] - cdf_obs[i-1])
#     end
# end
# # Creating a full distribution object..
# d = DiscreteNonParametric(obs, pmf_obs)
# # Sampling from the distribution object..
# out = rand(d, 100)
# # Test it...
# plot(sort(out))
# # quantile(d, 0.5)
# # quantile(d, 1 - 0.025)

# Current problem:
# 1. Interpolation of eCDF and sampling: interpolation is not monotonic and stable, and sampling also could not use CDF directly.
# 2. MC integration does not output explicit optimized variables. Have a hard time to calculate derived quantities.
