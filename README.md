# AquaSkyLES.jl

Oceananigans.jl comes above the ocean's surface...

### Using AquaSkyLES

First [install Julia](https://julialang.org/downloads/).

Then clone this repository

```bash
$ git clone git@github.com:navidcy/AquaSkyLES.jl.git
```

Open Julia from within the local directory of the repo via:

```bash
$ julia --project
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

https://github.com/user-attachments/assets/cc52de82-3d06-4742-adc2-01129fa198d7
