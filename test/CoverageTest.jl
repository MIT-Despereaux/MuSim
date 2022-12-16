### Tests for sampling geometry validations.

module CoverageTest

# %%
include("testutils.jl")

using Test
# using Infiltrator

# The following code is necessary to fix VSCode julia local module linting
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MuSim.jl")
    using .MuSim
else
    using MuSim
end

using Distributions, FLoops
# using Plots

# %%
"""
Test for basic coverage of one box. Returns beta and beta_err multiplied by simulation ℓ^2.
"""
function test_coverage1_hemisimlite(n_sim::Int=500000)
    I₀ = 58
    (x, y, z) = (1.2, 1.9, 2.2)
    ori = deg2rad.((45.0, 112.5))

    # Construct a X x Y x Z box (unit m)
    box = [RectBox("A", x, y, z, orientation=ori)]
    r = max(x, y, z) * 3
    ℓ = max(x * y, x * z, y * z) * 1.75
    # println("r = $r, ℓ = $ℓ")
    results = runhemisimlite(n_sim, box, r, ℓ; center=(0, 0, 0))

    expected_rate = analytic_R(box[1], I₀=I₀)

    β = sum(results) / n_sim
    β_err = √(β * (1 - β) / n_sim)
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    # println("rate = $rate ± $rate_err")
    # println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end

# %%
"""
Test for coincidence between multiple detectors. 
"""
function test_coverage2_hemisimlite(n_sim::Int=1000000)
    I₀ = 58

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ = 27π / 60
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    dets = [box1, box2]

    r = d * 2
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    @time results = runhemisimlite(n_sim, dets, r, ℓ; center=(0, 0, 0))
    res = view(results, :, 1) .& view(results, :, 2)

    β = sum(res) / n_sim
    β_err = √(β * (1 - β) / n_sim)
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    # println("rate = $rate ± $rate_err")
    # println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end


# %%
"""
Test for optimization.
"""
function test_coverage3_hemisimlite(n_sim::Int=1000000)
    I₀ = 1

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ = 20π / 60
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    dets = [box1, box2]

    r = d * 2
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    @time (results, dist_θ, dist_φ, angles) = runhemisimlite(n_sim, dets, r, ℓ)
    # println("dist_θ = $dist_θ, dist_φ = $dist_φ")
    angles = reshape(angles, 2, :)
    res = view(results, :, 1) .& view(results, :, 2)
    # p1 = histogram(angles[1, :], bins=500, label="θ")
    # p2 = histogram(angles[2, :], bins=500, label="φ")
    # display(p1)
    # display(p2)
    # Importance sampling
    β = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            β += 3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ))
        end
    end
    β /= n_sim
    # Calculate the variance
    σβ = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            σβ += (3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ)) - β)^2
        else
            σβ += (β)^2
        end
    end
    β_err = (√σβ) / n_sim
    β *= ℓ^2
    β_err *= ℓ^2
    # println("β = $β ± $β_err")
    rate = β * 2π / 3 * I₀
    rate_err = β_err * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end


# %%
"""
Test for optimization for the full simulation.
"""
function test_coverage4_hemisim(n_sim::Int=100000)
    I₀ = 58

    # Construct two planar surfaces oriented at θ
    w = 0.05
    d = 1.0 # r distance between centers (on a sphere)
    θ = 0.0
    θ2 = θ + deg2rad(3)
    box1 = RectBox("A", w, w, 0.001; orientation=(θ, 0))
    box2 = RectBox("B", w, w, 0.001; position=(sin(θ) * d, 0, cos(θ) * d), orientation=(θ, 0))
    box3 = RectBox("C", w, w, 0.001; position=(sin(θ2) * d * 1.5, 0, cos(θ2) * d * 1.5), orientation=(θ2, 0))
    dets = [box1, box2, box3]

    r = d * 2
    ℓ = w * 2
    # println("r = $r, ℓ = $ℓ")
    (results, _, _, dist_θ, dist_φ, angles) = runhemisim(n_sim, dets, r, ℓ)
    angles = reshape(angles, 2, :)
    res = view(results, :, 1) .& view(results, :, 2)
    # Importance sampling
    β = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            β += 3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ))
        end
    end
    β /= n_sim
    # Calculate the variance
    σβ = 0
    for i in eachindex(res)
        if res[i]
            local θ = angles[1, i]
            local φ = angles[2, i]
            σβ += (3 / (2π) * cos(θ)^2 * sin(θ) / (pdf(dist_θ, θ) * pdf(dist_φ, φ)) - β)^2
        else
            σβ += (β)^2
        end
    end
    β_err = (√σβ) / n_sim
    # println("β = $β ± $β_err")
    rate = β * ℓ^2 * 2π / 3 * I₀
    rate_err = β_err * ℓ^2 * 2π / 3 * I₀
    expected_rate = w^2 * I₀ * cos(θ)^2 * w^2 / d^2
    println("rate = $rate ± $rate_err")
    println("expected rate = $expected_rate")
    @test rate ≈ expected_rate atol = 2rate_err
end


function setup_geometry()
    hL = 70e-3
    dL = 70e-3
    lenL = 600e-3

    hS = 20e-3
    dS = 72e-3
    lenS = 510e-3

    hCW = 10e-3
    dCW = 50e-3
    lenCW = 50e-3
    hhCW = 30e-3 # height of the scintillator from the bottom

    gCB = 8e-3
    gBD = 11e-3
    gDA = 6e-3

    (xC, yC, zC) = (0.625, 0.220, 0.8682)

    zB = zC + hL / 2 + gCB + hS / 2
    zD = zB + hS / 2 + gBD + hL / 2
    zA = zD + hL / 2 + gDA + hS / 2
    zCW = zC + hL / 2 + hhCW

    xL = xC
    xS = xC - lenL / 2 + lenS / 2
    xCW = xS + lenS / 2 + lenCW / 2 + 0.002


    loc_det_C = (xL, yC, zC)
    loc_det_B = (xS, yC, zB)
    loc_det_D = (xL, yC, zD)
    loc_det_A = (xS, yC, zA)
    loc_det_CWD = (xCW, yC, zCW)

    loc_Chip = (0.338, 0.229, 1.4148)

    det_CWD = RectBox("CWD", lenCW, dCW, hCW, position=loc_det_CWD, material="POP Doped Polystyrene")
    det_A = RectBox("A", lenS, dS, hS, position=loc_det_A, material="EJ-200")
    det_D = RectBox("D", lenL, dL, hL, position=loc_det_D, material="EJ-200")
    det_B = RectBox("B", lenS, dS, hS, position=loc_det_B, material="EJ-200")
    det_C = RectBox("C", lenL, dL, hL, position=loc_det_C, material="EJ-200")
    det_Chip = RectBox("Chip", 350e-6, 0.005, 0.005, position=loc_Chip, material="Unknown")

    detectors = [det_A, det_B, det_C, det_D, det_CWD, det_Chip]
    return detectors
end


function main()
    initrand()
    test_coverage1_hemisimlite()
    test_coverage2_hemisimlite()
    # test_coverage3_hemisimlite()
    # test_coverage4_hemisim()
end

!isinteractive() && main()

end
