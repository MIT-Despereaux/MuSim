# MuSim
This folder contains a package for quick muon-detector simulations. 


# Installation for development
Add the current package to Julia general environment for dev, it will be recorded in the general Manifest file:
```julia
using Pkg; Pkg.develop("/path/to/package/")
```
To remove (since this package is local and not registered under the official Julia package repo):
```julia
using Pkg; Pkg.rm("PackageName")
```

# Python Binding
**Currently this section is under re-assessment for choosing between `PythonCall` and `PyCall`. Note that at the current state calling julia function from python should be avoided.**

Using `PyCall`: 
- Its default python executable and environment can be set when installed. I just bind it to my conda production environment which has `jupyter` installed. To bind a specific python environment, run the following commands in julia root environment:
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
All subsequent python wrapper package would use the same python binding as the `PyCall` package.

Reference: https://github.com/JuliaPy/PyCall.jl

Using `PythonCall`:
- This is much simpler: just install the package under the julia root environment.
- After that, follow this section from the official document to choose the python executable you want: https://cjdoris.github.io/PythonCall.jl/stable/pythoncall/#pythoncall-config.

**Please note that saving data in python-compatible formats current is not in the scope of this project. Please use `PythonCall` or `PyCall` in your own environment.**


# Usage
To run a simulation, simply add the `MuSim` package to your Julia root environment (see above) and use it as a normal package.
Examplary tests are in the `./test` directory, but the following gives a very quick example.
```julia
using MuSim
sim_num = Int(1e3)
r = 0.80
ℓ = 0.75
det1 = RectBox("Name", 0.01, 0.05, 0.05, position=(0, 0, 0), orientation=deg2rad.((0, 0)), efficiency=0.98, material="POP Doped Polystyrene")
detectors = [det1]
@time runhemisim(sim_num, detectors, r, (0, 0, 0), ℓ)
```
