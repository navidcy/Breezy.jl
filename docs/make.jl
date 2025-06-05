push!(LOAD_PATH,"../src/")

using AquaSkyLES
using Documenter

DocMeta.setdocmeta!(AquaSkyLES, :DocTestSetup, :(using AquaSkyLES); recursive=true)

makedocs(sitename="AquaSkyLES",
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo="github.com/gregorywagner/AquaSkyLES.jl",
    devbranch="main",
    push_preview=true,
) 