using Documenter
using Literate
using TightBindingLattice

example_directory = joinpath(@__DIR__, "..", "examples")
output_directory = joinpath(@__DIR__, "src/generated")
for filename in ["example_group_isomorphism.jl",
                 "example_honeycomb_symmetry.jl",
                 "example_orthocube.jl",
                 "example_kagome_symmetry.jl",
                 "example_little_group_2d.jl",
                 "example_little_group_3d.jl",
                 "example_little_symmetry_kspace_honeycomb.jl",
                 "example_point_group_4mm.jl",
                 "example_point_group_elements.jl",]
    example = joinpath(example_directory, filename)
    Literate.markdown(example, output_directory)
end

makedocs(
    modules=[TightBindingLattice],
    doctest=true,
    sitename="TightBindingLattice.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Group isomorphism" => "generated/example_group_isomorphism.md",
            "Orthocube" => "generated/example_orthocube.md",
            "Point group elements" => "generated/example_point_group_elements.md",
            "Point group 4mm (C₄ᵥ)" => "generated/example_point_group_4mm.md",
            "Little group in 2d" => "generated/example_little_group_2d.md",
            "Little group in 3d" => "generated/example_little_group_3d.md",
            "Honeycomb lattice" => "generated/example_honeycomb_symmetry.md",
            "Honeycomb lattice in k-space" => "generated/example_little_symmetry_kspace_honeycomb.md",
            "Kagome lattice" => "generated/example_kagome_symmetry.md",
        ]
    ]
)

deploydocs(
        repo = "github.com/kyungminlee/TightBindingLattice.jl.git",
        devbranch = "dev",
)
