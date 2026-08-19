using ForestInventory
using Documenter

DocMeta.setdocmeta!(ForestInventory, :DocTestSetup, :(using ForestInventory); recursive=true)

makedocs(;
    modules=[ForestInventory],
    authors="Marcos Daniel da Silva <marcosdasilva@5a.tec.br> and contributors",
    sitename="ForestInventory.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaForests.github.io/ForestInventory.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaForests/ForestInventory.jl",
    devbranch="master",
)
