### Run tests

println("Testing distributions...")
t = @elapsed include("DistributionTest.jl")
println("Time elapse: $t seconds")
