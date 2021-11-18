# MuSim
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

# If there is no "PyCall" in the root environment, then proceed to install it
if Base.find_package("PyCall") === nothing
    # Check if the production environment has been created
    conda_env = "prod_linux" # <-- this is the environment where you want PyCall to bind to
    # Create a file with the following content:
    # 1. Activate conda environment specified by "conda_env"
    # 2. Get the python executable path by running "which python"
    file = "eval \"\$(conda shell.bash hook)\" \nconda activate $conda_env && which python"
    # Save the file to $HOME/.julia/config directory
    file_name = joinpath(ENV["HOME"], ".julia/config/tmp.sh")
    open(file_name, "w") do f
        println(f, file)
        flush(f)
    end
    # Run the file and get the python executable path
    if !success(Cmd(`bash $file_name`))
        @warn "Conda $prod_env cannot be activated. Exit building PyCall..."
    else
        # Install the PyCall package in the root environment
        Pkg.add("PyCall")
        pypath = ENV["PYTHON"] = readchomp(`bash $file_name`)
        println("Using $pypath as python path.")
        Pkg.build("PyCall")
    end
    # Delete the tmp file
    rm(file_name)
end
```

You would then have the `PyCall` package installed under your root environment.
All subsequent python wrapper package would use the same python binding as 
the `PyCall` package.

Reference: https://github.com/JuliaPy/PyCall.jl

# Usage
To run a simulation, simply add the `MuSim` package to your Julia root environment (see above) and use it as a normal package.
Tutorial files are in the `./tut` directory, but the following gives a very quick example.
```julia
using MuSim
sim_num = Int(1e3)
r = 0.80
ℓ = 0.75
det1 = RectBox(0.01, 0.05, 0.05, position=(0, 0, 0), orientation=deg2rad.((0, 0)), efficiency=0.98, material="POP Doped Polystyrene")
detectors = [det1]
runhemisim(sim_num, detectors, r, (0, 0, 0), ℓ)
```

# Python Interface
Since most analyses will be run in python, calling this package from python is documented below.


# TODOs:
- [x] Complete writing the tests.
- [ ] [Medium Priority] Explore the interface with python.
- [ ] [Medium Priority] Construct a multip-dimensional sampler for arbitrary distributions
    - Motivation: how to efficiently sample an arbitrary multivariate distribution given it's analytical form?
    - We also needs to sample its integrated probability density, for example: sample the energy spectrum of a muon after we know its solid angle.