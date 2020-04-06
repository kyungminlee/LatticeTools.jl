using TightBindingLattice

# Little group

unitcell = make_unitcell([1.0 0.0 0.0; 0.0 1.0 0.0; 0 0 1]; OrbitalType=String)
addorbital!(unitcell, "Ox", FractCoord([0,0,0], [0.5, 0.0, 0.0]))
addorbital!(unitcell, "Oy", FractCoord([0,0,0], [0.0, 0.5, 0.0]))

psym = project(PointSymmetryDatabase.get(15), [1 0 0; 0 1 0; 0 0 1])
#idx_C4 = 3
#@show psym
#@show typeof(psym)

lattice = make_lattice(unitcell, [4 0 0; 0 4 0; 0 0 3])
perms = get_orbital_permutations(lattice, psym)

tsym = TranslationSymmetry(lattice)

@show num_irreps(tsym)

@show tsym.orthogonal_coordinates


for idx in 1:num_irreps(tsym)
    #@show little_symmetry(tsym, idx, psym)
    psym_little = little_symmetry(tsym, idx, psym)
    psym_little2 = TightBindingLattice.little_symmetry_iso(tsym, idx, psym)
    @show idx
    @show tsym.hypercube.coordinates[idx]

    @show psym_little.hermann_mauguinn, group_order(psym_little)
    @show psym_little2.hermann_mauguinn, group_order(psym_little2)
    @show iscompatible(tsym, idx, psym)
    @show iscompatible(tsym, idx, psym_little)
    @show iscompatible(tsym, idx, psym_little2)
end
