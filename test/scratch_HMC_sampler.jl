# %%
using StatsBase, Distributions, Interpolations, Plots, LinearAlgebra, Random
using LogDensityProblems, LogDensityProblemsAD, TransformVariables, TransformedLogDensities, DynamicHMC

# %%
# Creating original samples...
C = Cauchy(0, 0.01)
mean_V = [0, 0]
diag_V = [0.1, 0.5]
N = MvNormal(mean_V, diag_V)
x_range = 0:0.005:1
y_range = -1:0.01:1

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
using MuSim, HCubature, Plots

p = μPDF()
s = μPDFSettings()
res = drawsamples(p, s, N=10000)

# %%
# Plot the theoretical and histogram
E_range = LinRange(s.E_range..., 101)
θ_range = LinRange(s.θ_range..., 101)

# %%
# Plot the histogram in log
p1 = histogram2d(res[1, :], res[2, :], label="Posterior", bins=(E_range, θ_range), normed=true, color=:plasma, alpha=0.2)

normalization = hcubature(x -> exp(p(x)), (s.E_range[1], s.θ_range[1]), (s.E_range[2], s.θ_range[2]))[1]
pdf_theo = exp.(p.([(E=e, θ=t) for t in θ_range, e in E_range])) / normalization
# Create the heatmap plot
p1 = contour!(E_range, θ_range, pdf_theo,
    xlabel="E",
    ylabel="θ",
    levels=100,
    color=:plasma,
    title="2D Distribution")
display(p1)

# %%
# Plot the traces
# plot(res[1, :], label="E", yscale=:log10)
normalization = hcubature(x -> exp(p(x)), (s.E_range[1], s.θ_range[1]), (s.E_range[2], s.θ_range[2]))

# %%
function E_pdf(p, E)
    f = y -> hquadrature(x -> exp.(p((E=y, θ=x))), s.θ_range[1], s.θ_range[2])[1]
    return f(E) / normalization[1]
end

histogram(res[1, :], label="E", normalize=:pdf)
plot!(E_range, [E_pdf(p, e) for e in E_range], label="E theoretical")

# %%
function θ_pdf(p, θ)
    f = y -> hquadrature(x -> exp.(p((E=x, θ=y))), s.E_range[1], s.E_range[2])[1]
    return f(θ) / normalization[1]
end

histogram(res[2, :], label="θ", normalize=:pdf)
plot!(θ_range, [θ_pdf(p, t) for t in θ_range], label="θ theoretical")

