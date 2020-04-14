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
tsym = TranslationSymmetry(lattice)

perms = get_orbital_permutations(lattice, psym)



@show num_irreps(tsym)

@show tsym.orthogonal_coordinates

# x = collect(get_irrep_iterator(lattice, tsym, 1, :))
# y = collect(get_irrep_iterator(lattice, tsym, 1, 1))
# @show x == y
# exit()

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


count = 1
for ssic in get_irrep_components(tsym, psym)
    global count
    count += 1
end
@show count
