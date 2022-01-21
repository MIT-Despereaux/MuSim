### Test utensils

using Random
Random.seed!(42)

"""
Initialize plotting.
"""
function initplots()
    using Plots
    gr()
    using Plots.PlotMeasures
    default(size = (600, 600), dpi = 100, margin = 3.0mm)
    println("Default style set.")
end
