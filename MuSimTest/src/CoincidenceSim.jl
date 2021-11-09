module CoincidenceSim

"""
This small module is for simulating cosmic ray coincidences.
"""


include("Particles.jl")
export Ray, raydir, rayplaneint, μenergy!

include("RayTracing.jl")
export LabObject, RectBox, Sphere, Cylinder, isthrough!

include("SimUtils.jl")
export runsim, runhemisim

end
