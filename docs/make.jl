using TCCFunctions
using Documenter

DocMeta.setdocmeta!(TCCFunctions, :DocTestSetup, :(using TCCFunctions); recursive=true)

makedocs(;
    modules=[TCCFunctions],
    authors="renilsonlisboa <renilsonlisboajunior@gmail.com> and contributors",
    repo="https://github.com/renilsonlisboa/TCCFunctions.jl/blob/{commit}{path}#{line}",
    sitename="TCCFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://renilsonlisboa.github.io/TCCFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/renilsonlisboa/TCCFunctions.jl",
)
