### Test utensils


"""
Initialize plotting.
"""
function initplots()
    # Use @eval to execute the block as macro under global environment.
    @eval begin
        using Plots
        ENV["GKSwstype"] = "100" # Special fix for running without graphical output (only saving files).
        gr()
        using Plots.PlotMeasures
        default(size=(600, 600), dpi=100, margin=3.0mm)
    end
    println("Default style set.")
end

"""
Initialize random seed.
"""
function initrand()
    @eval begin
        using Random
        println("Random seed set to 1234.")
        Random.seed!(1234)
    end
end
