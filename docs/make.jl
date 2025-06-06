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
    repo="github.com/NumericalEarth/Breeze.jl",
    devbranch="main",
    push_preview=true,
)
