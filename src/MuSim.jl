module MuSim

"""
This small module is for simulating muon interactions with the detectors.
"""

import Base: +, -, ==, hash, show, @kwdef
using Printf
using LinearAlgebra
using Random
using StaticArrays
using DataStructures
using SparseArrays
using Roots
using JLD2
using StatsBase
using HCubature
using MCIntegration
using Interpolations
using LogDensityProblems
using LogDensityProblemsAD
using TransformVariables
using TransformedLogDensities
using DynamicHMC
using FLoops
using BangBang
using MicroCollections
using Distributions
using Combinatorics
using CSV
using DataFrames


export Ray, modifyray!, raydir, rayplaneint, randContinuous, interpContinuous,
    probContinuous, μPDF, μPDFSettings, drawsamples, LabObject, RectBox, Sphere,
    Cylinder, isthrough!, objorient, analytic_R, gendetectorpos, runhemisim,
    runhemisimlite, runexp, calculateβ, calculateβ_MC, βio, βio_MC, expio


include("Particles.jl")
include("RayTracing.jl")
include("SimUtils.jl")


# end module
end
