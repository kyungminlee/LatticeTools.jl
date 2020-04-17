using TightBindingLattice

# Little group

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
#idx_C4 = 3

lattice = make_lattice(unitcell, [4 0; 0 4])
perms = get_orbital_permutations(lattice, psym)

tsym = TranslationSymmetry(lattice)

lge = little_group_elements(tsym, psym)
