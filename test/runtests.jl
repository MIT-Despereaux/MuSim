### Run tests

println("Testing basic structures...")
t = @elapsed include("BasicTest.jl")
println("Time elapsed: $t seconds")

println("Testing distributions...")
t = @elapsed include("DistributionTest.jl")
println("Time elapsed: $t seconds")

println("Testing coverages...")
t = @elapsed include("CoverageTest.jl")
println("Time elapsed: $t seconds")

println("Testing simulation utilities...")
t = @elapsed include("SimUtilsTest.jl")
println("Time elapsed: $t seconds")
