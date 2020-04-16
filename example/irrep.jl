using TightBindingLattice


unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])

lattice = make_lattice(unitcell, [4 0; 0 4])
tsym = TranslationSymmetry(lattice)

tsic = TranslationSymmetryIrrepComponent(tsym, 2, 1)
@show collect( get_irrep_iterator(lattice, tsic) )

#psic = PointSymmetryIrrepComponent(psym, 1, 1)
#@show collect( get_irrep_iterator(lattice, psic) )

psym_little = little_symmetry(tsic, psym)
psic = PointSymmetryIrrepComponent(psym_little, 1, 1)
ssic = SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)
@show collect( get_irrep_iterator(lattice, ssic) )


psic = PointSymmetryIrrepComponent(psym, 1, 1)
ssic = SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)

#=
for tsym_irrep in 1:num_irreps(tsym)
  psym_little = little_symmetry(tsym, tsym_irrep, psym)
  k = tsym.hypercube.coordinates[tsym_irrep]
  @test iscompatible(tsym, tsym_irrep, psym) == (k in [[0,0], [2,2]])
  @test iscompatible(tsym, tsym_irrep, psym_little)
  lg_matrep = psym.matrix_representations[little_group_element_indices(tsym, tsym_irrep, psym)]
  @test !isnothing(group_isomorphism(little_group(tsym, tsym_irrep, psym),
                                     FiniteGroup(group_multiplication_table(lg_matrep))))

end

=#
