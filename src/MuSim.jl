module MuSim

"""
This small module is for simulating muon interactions with the detectors.
"""


include("Particles.jl")
export Ray, modifyray!, raydir, rayplaneint, randContinuous, interpContinuous, probContinuous, μPDF, μPDFSettings, drawsamples

include("RayTracing.jl")
export LabObject, RectBox, Sphere, Cylinder, isthrough!, objorient, analytic_R

include("SimUtils.jl")
export gendetectorpos, runhemisim, runhemisimlite, runexp, calculateβ, calculateβ_MC, βio, βio_MC, expio

end
