module MuSim

"""
This small module is for simulating muon interactions with the detectors.
"""


include("Particles.jl")
export Ray, raydir, rayplaneint, μenergy!

include("RayTracing.jl")
export LabObject, RectBox, Sphere, Cylinder, isthrough!

include("SimUtils.jl")
export runsim, runhemisim, runhemisimlite

end
