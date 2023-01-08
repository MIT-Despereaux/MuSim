using StatsBase, Distributions, Interpolations, Plots, MuSim, LinearAlgebra, MCIntegration

# %%
θ = res.config.var[1]
itp = interpContinuous(θ)
plot(0:0.005:dist.grid[end], itp.(0:0.005:dist.grid[end]))
plot(0:0.005:dist.grid[end], probContinuous(θ, 0:0.005:dist.grid[end]))
obs = randContinuous(θ; n=1000000)
histogram!(obs, bins=500, normalize=:pdf, alpha=0.5)

# %%
φ = res.config.var[2]
itp = interpContinuous(φ)
plot(0:0.005:dist.grid[end], itp.(0:0.005:dist.grid[end]))
plot(0:0.005:dist.grid[end], probContinuous(φ, 0:0.005:dist.grid[end]))
obs = randContinuous(φ; n=1000000)
histogram!(obs, bins=500, normalize=:pdf, alpha=0.5)

# %%
# Creating original samples...
C = Cauchy(0, 0.01)
obs = rand(C, 100000)
eCDF = ecdf(obs)

# %%
hist = fit(Histogram, obs, nbins=100)
edge_diff = hist.edges[1][2] - hist.edges[1][1]
edges = hist.edges[1][1]:edge_diff:hist.edges[1][end]

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
