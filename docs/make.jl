using Documenter
using TightBindingLattice

makedocs(
    modules=[TightBindingLattice],
    doctest=true,
    sitename="TightBindingLattice.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
  )

deploydocs(
    repo = "github.com/kyungminlee/TightBindingLattice.jl.git",
  )
