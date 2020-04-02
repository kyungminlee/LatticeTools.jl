using TightBindingLattice
using YAML
using LinearAlgebra

scale_matrix = [4 -2; 2 2]
@assert det(scale_matrix) ≈ 12

unitcell = make_unitcell([1 -0.5; 0 0.5*sqrt(3.0)])
addorbital!(unitcell, "A", carte2fract(unitcell, [0.5, 0.5/sqrt(3)]))
addorbital!(unitcell, "B", carte2fract(unitcell, [0.5,-0.5/sqrt(3)]))

lattice = make_lattice(unitcell, scale_matrix)

#hypercube = lattice.hypercube
#hypercube = orthogonalize(hypercube)
translation_symmetry = TranslationSymmetry(lattice)

#data_directory = joinpath("..", "PointGroupData")
#data_yaml = YAML.load_file(joinpath(data_directory, "PointGroup3D-19.yaml"))
point_symmetry = project(
    PointSymmetryDatabase.get(19),
    [1 0 0; 0 1 0])

transperm = get_orbital_permutations(lattice, translation_symmetry)
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

@show translation_symmetry.hypercube.coordinates[translation_symmetry.generators]
@show translation_symmetry.orthogonal_coordinates
@show translation_symmetry.hypercube.coordinates



for momentum_index in 1:num_irreps(translation_symmetry)
    @show momentum_index
    @show translation_symmetry.hypercube.coordinates[momentum_index]
    @show translation_symmetry.orthogonal_coordinates[momentum_index]
    @show little_group(translation_symmetry, momentum_index, point_symmetry)
end


tsym = translation_symmetry
psym = point_symmetry

for momentum_index in 1:12
    @show momentum_index
    for o in get_irrep_iterator(lattice, translation_symmetry, momentum_index, 1)
        @show o
    end
    println()
end

for tsym_irrep_index in 1:12

    psym_little = little_symmetry(tsym, tsym_irrep_index, psym)
    @show tsym_irrep_index
    @show psym_little.hermann_mauguinn
    @show num_irreps(psym_little)

    for psym_irrep_index in 1:num_irreps(psym_little)
        println(tsym.irreps[tsym_irrep_index].name,"    ", psym_little.irreps[psym_irrep_index].name)
        for o in get_irrep_iterator(lattice,
                                    tsym, psym_little,
                                    tsym_irrep_index, psym_irrep_index,
                                    1, 1)
            #@show o
        end
        #println()
    end
    println()
end
