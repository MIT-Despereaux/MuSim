using StatsBase, Distributions, MuSim
# Creating original samples...
N = Normal(0, 1)
obs = rand(N, 10000)
sort!(obs)
# Cdf..
myCdf = ecdf(obs)

cdf_obs = map(x -> myCdf(x), obs)
result = []
for n in 1:10000
    push!(result, unisamplercdf!(myCdf; low_bound=round(minimum(myCdf)), upp_bound=round(maximum(myCdf))))
end

histogram(obs, bins=100, normalize=:pdf)
histogram!(result, bins=100, normalize=:pdf)

# %%
# using MCIntegration

# res = integrate((x, c) -> log(x[1]) / sqrt(x[1]))
# plot(res.config.var[1].histogram)

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
