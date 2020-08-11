using Documenter
using Literate
using TightBindingLattice

example_directory = joinpath(@__DIR__, "..", "examples")
output_directory = joinpath(@__DIR__, "src/generated")
for filename in [
    "example_group_isomorphism.jl",
    "example_honeycomb_symmetry.jl",
    "example_orthocube.jl",
    "example_kagome_symmetry.jl",
    "example_little_group_2d.jl",
    "example_little_group_3d.jl",
    "example_little_symmetry_kspace_honeycomb.jl",
    "example_point_group_4mm.jl",
    "example_point_group_elements.jl",
]
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
        "lattice.md",
        "Symmetry" => [
            "Symmetry/operation-space.md",
            "Symmetry/symmetry-space.md",
            "Symmetry/irrep.md",
        ],
        "Examples" => [
            "Group Isomorphism" => "generated/example_group_isomorphism.md",
            "OrthoCube" => "generated/example_orthocube.md",
            "Point Group Elements" => "generated/example_point_group_elements.md",
            "Point Group 4mm (C₄ᵥ)" => "generated/example_point_group_4mm.md",
            "Little Group in 2D" => "generated/example_little_group_2d.md",
            "Little Group in 3D" => "generated/example_little_group_3d.md",
            "Honeycomb Lattice" => "generated/example_honeycomb_symmetry.md",
            "Honeycomb Lattice in k-space" => "generated/example_little_symmetry_kspace_honeycomb.md",
            "Kagome Lattice" => "generated/example_kagome_symmetry.md",
        ],
        "API Reference" => [
            "Basic" => "API/basic.md",
            "Group" => "API/group.md",
            "Symmetry Operation" => "API/symmetry-operation.md",
            "Symmetry" => "API/symmetry.md",
            "Symmetry Embedding" => "API/symmetry-embedding.md",
            "Symmetry Irrep Component" => "API/symmetry-irrep-component.md",
        ]
    ]
)

deploydocs(
    repo = "github.com/kyungminlee/TightBindingLattice.jl.git",
    devbranch = "dev",
)
