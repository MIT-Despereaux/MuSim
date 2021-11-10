### Run tests

println("Testing distributions...")
t = @elapsed include("DistributionTest.jl")
println("Time elapsed: $t seconds")

println("Testing coverages...")
t = @elapsed include("CoverageTest.jl")
println("Time elapsed: $t seconds")
