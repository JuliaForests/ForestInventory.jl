using ForestInventory
using Documenter

DocMeta.setdocmeta!(ForestInventory, :DocTestSetup, :(using ForestInventory); recursive=true)

makedocs(;
  modules=[ForestInventory],
  authors="Marcos Daniel da Silva <marcosdasilva@5a.tec.br> and contributors",
  sitename="ForestInventory.jl",
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://JuliaForests.github.io/ForestInventory.jl",
    edit_link="master",
    assets=String[],
  ),
  checkdocs=:exports,
  pages=[
    "ForestInventory Package" => "forestinventory.md",
    "Getting Started" => "tutorial.md",
    "API Reference" => "reference.md",
    "Index" => "index.md",
  ],
)

deploydocs(;
  repo="github.com/JuliaForests/ForestInventory.jl",
  devbranch="master",
)
