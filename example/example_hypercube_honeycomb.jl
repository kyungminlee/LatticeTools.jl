using TightBindingLattice
using YAML
using LinearAlgebra

scale_matrix = [4 -2; 2 2]
@assert det(scale_matrix) ≈ 12

unitcell = make_unitcell([1 -0.5; 0 0.5*sqrt(3.0)])
addorbital!(unitcell, "A", carte2fract(unitcell, [0.5, 0.5/sqrt(3)]))
addorbital!(unitcell, "B", carte2fract(unitcell, [0.5,-0.5/sqrt(3)]))

lattice = make_lattice(unitcell, scale_matrix)
hypercube = lattice.hypercube

#supercell = make_supercell(unitcell, scale_matrix)

#hypercube = HypercubicLattice(scale_matrix)
translation_symmetry = TranslationSymmetry(hypercube)

data_directory = joinpath("..", "PointGroupData")
data_yaml = YAML.load_file(joinpath(data_directory, "PointGroup3D-19.yaml"))
point_symmetry = project(PointSymmetry(data_yaml), [1 0 0; 0 1 0])

#@show translation_symmetry
#@show point_symmetry

transperm = get_orbital_permutations(lattice)
pointperm = get_orbital_permutations(lattice, point_symmetry)

println("# Translation group period lengths")
println(translation_symmetry.group.period_lengths)

println("# Translation permutations")
for t in transperm
    println(t)
end
println("* orders are different. (hypercubic order vs. orthogonal order)")
println()

println("# Point group period lengths")
println(point_symmetry.group.period_lengths)

println("# Point permutations")
for p in pointperm
    println(p)
end
println("* orders are the same.")
println()


hypercube = orthogonalize(hypercube)
translation_symmetry = TranslationSymmetry(hypercube)

data_directory = joinpath("..", "PointGroupData")
data_yaml = YAML.load_file(joinpath(data_directory, "PointGroup3D-19.yaml"))
point_symmetry = project(PointSymmetry(data_yaml), [1 0 0; 0 1 0])

transperm = get_orbital_permutations(lattice)
pointperm = get_orbital_permutations(lattice, point_symmetry)

println("# Translation group period lengths")
println(translation_symmetry.group.period_lengths)

println("# Translation permutations")
for t in transperm
    println(t)
end
println("* orders are same now")
println()

println("# Point group period lengths")
println(point_symmetry.group.period_lengths)

println("# Point permutations")
for p in pointperm
    println(p)
end
println("* orders are the same.")
println()
