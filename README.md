# Ray Tracing
This folder contains a package for quick muon-detector simulations. 

# Installation for development
Add the current package to Julia general envronment for dev, it will be recorded in the general Manifest file:
```julia
]dev /path/to/package/
```
To remove (since this package is local and not registered under the official Julia package repo):
```julia
]rm packagename
```

Note: the required packages include `PyCall` which deals with python bindings.
Its default python executable and environment can be set when installed. 
I just bind it to my conda production environment which has `jupyter` installed.
To bind a specific python environment, run the following commands in julia:
```julia
using Pkg;
Pkg.activate()

if Base.find_package("PyCall") === nothing
    # Check if the production environment has been created
    prod_env = "prod_linux" # <-- this is the environment where you want PyCall to bind to
    file = "eval \"\$(conda shell.bash hook)\" \nconda activate $prod_env && which python"
    file_name = joinpath(ENV["HOME"], ".julia/config/tmp.sh")
    open(file_name, "w") do f
        println(f, file)
        flush(f)
    end
    if !success(Cmd(`bash $file_name`))
        @warn "Conda $prod_env cannot be activated. Exit building PyCall..."
    else
        Pkg.add("PyCall")
        pypath = ENV["PYTHON"] = readchomp(`bash $file_name`)
        println("Using $pypath as python path.")
        Pkg.build("PyCall")
    end
    rm(file_name)
else
```

You would then have the `PyCall` package installed under your root environment.
All subsequent python wrapper package would use the same python binding as 
the `PyCall` package.

Reference: https://github.com/JuliaPy/PyCall.jl

# Usage

To run a simulation, simply add the `MuSim` package to your Julia root environment (see above) and use it as a normal package.
For example:
```julia
using MuSim
sim_num = Int(1e3)
r = 0.80
ℓ = 0.75
det1 = RectBox(0.01, 0.05, 0.05, position=(0, 0, 0), orientation=deg2rad.((0, 0)), efficiency=0.98, material="POP Doped Polystyrene")
detectors = [det1]
runhemisim(sim_num, detectors, r, (0, 0, 0), ℓ)
```

# TODOs:
- [High Priority] Complete the tests.
- [Medium Priority] Explore the interface with python.
- [Low Priority] Construct a multip-dimensional sampler for arbitrary distributions
    - Motivation: how to efficiently sample an arbitrary multivariate distribution given it's analytical form?
    - We also needs to sample its integrated probability density, for example: sample the energy spectrum of a muon after we know its solid angle.