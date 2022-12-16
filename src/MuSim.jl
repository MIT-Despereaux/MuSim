module MuSim

"""
This small module is for simulating muon interactions with the detectors.
"""


include("Particles.jl")
export Ray, modifyray!, raydir, rayplaneint, μenergy!, unisampler!, unisamplercdf!

include("RayTracing.jl")
export LabObject, RectBox, Sphere, Cylinder, isthrough!, objorient, analytic_R

include("SimUtils.jl")
export gendetectorpos, runhemisim, runhemisimlite, runexp, calculateβ, βio, composeβ, expio

end
