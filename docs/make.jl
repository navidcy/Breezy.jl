push!(LOAD_PATH,"../src/")

using Breeze
using Documenter

DocMeta.setdocmeta!(Breeze, :DocTestSetup, :(using Breeze); recursive=true)

makedocs(sitename="Breeze",
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo="github.com/gregorywagner/Breeze.jl",
    devbranch="main",
    push_preview=true,
) 