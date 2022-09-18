module MuSim

"""
This small module is for simulating muon interactions with the detectors.
"""


include("Particles.jl")
export Ray, raydir, rayplaneint, μenergy!

include("RayTracing.jl")
export LabObject, RectBox, Sphere, Cylinder, isthrough!, objorient

include("SimUtils.jl")
export gensimconfigs, runhemisim, runhemisimlite, runexp, calculateβ, βio, composeβ

end
