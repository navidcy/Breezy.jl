using AquaSkyLES
using Documenter

DocMeta.setdocmeta!(AquaSkyLES, :DocTestSetup, :(using AquaSkyLES); recursive=true)

makedocs(;
    modules=[AquaSkyLES],
    authors="Gregory Wagner",
    repo="https://github.com/gregorywagner/AquaSkyLES.jl/blob/{commit}{path}#{line}",
    sitename="AquaSkyLES.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gregorywagner.github.io/AquaSkyLES.jl",
        edit_link="main",
        assets=String[],
    ),
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