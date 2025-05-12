# AquaSkyLES.jl

### Instructions

First [install Julia](https://julialang.org/downloads/); suggested version 1.10. See [juliaup](https://github.com/JuliaLang/juliaup) README for how to install 1.10 and make that version the default.

Then clone this repository

```bash
git clone git@github.com:navidcy/AquaSkyLES.jl.git
```

Open Julia from within the local directory of the repo via:

```bash
julia --project
```

The first time, you need to install any dependencies:

```julia
julia> using Pkg; Pkg.instantiate()
```

Now you are ready to run any of the examples!

For instance,

```julia
julia> include("examples/free_convection.jl")
```

produces

https://github.com/user-attachments/assets/dc45d188-6c61-4eb5-95fb-9a51c6f99013
