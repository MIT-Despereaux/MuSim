### Test utensils

using Random
Random.seed!(42)

"""
Initialize plotting.
"""
function initplots()
    # Use @eval to execute the block as macro under global environment.
    @eval begin
        using Plots
        gr()
        using Plots.PlotMeasures
        default(size = (600, 600), dpi = 100, margin = 3.0mm)
    end
    println("Default style set.")
end
