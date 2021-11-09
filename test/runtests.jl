### Run tests

println("Testing distributions...")
t = @elapsed include("DistributionTest.jl")
println("Time elapse: $t seconds")

println("Testing coverages...")
t = @elapsed include("CoverageTest.jl")
println("Time elapse: $t seconds")
