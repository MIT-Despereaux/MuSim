module MuSim

"""
This small module is for simulating muon interactions with the detectors.
"""


include("Particles.jl")
export Ray, raydir, rayplaneint, μenergy!

include("RayTracing.jl")
export LabObject, RectBox, Sphere, Cylinder, isthrough!, objorient

include("SimUtils.jl")
export runhemisim, runhemisimlite, runexp, calcgeo, geometricio

end
