# # Little Group in 3D

using LatticeTools

unitcell = makeunitcell([1.0 0.0 0.0; 0.0 1.0 0.0; 0 0 1]; SiteType=String)
addsite!(unitcell, "Ox", FractCoord([0,0,0], [0.5, 0.0, 0.0]))
addsite!(unitcell, "Oy", FractCoord([0,0,0], [0.0, 0.5, 0.0]))

# Currently, makelattice in three-dimension is not well supported

lattice = makelattice(unitcell, [4 0 0; 0 4 0; 0 0 3], [1 0 0; 0 1 0; 0 0 1])
tsym = FiniteTranslationSymmetry(lattice.hypercube, [1 0 0; 0 1 0; 0 0 1])
psym = project(PointSymmetryDatabase.get(15), [1 0 0; 0 1 0; 0 0 1])

for idx in 1:num_irreps(tsym)
    kf = tsym.fractional_momenta[idx]
    k = lattice.unitcell.reducedreciprocallatticevectors * kf
    psym_little1 = little_symmetry(tsym, idx, psym)
    psym_little2 = LatticeTools.little_symmetry_iso(tsym, idx, psym)
    println("- irrep_index: $(idx)")
    println("  momentum: $(k)")
    println("  little_point_group1: { name: \"$(psym_little1.hermann_mauguin)\", order: $(group_order(psym_little1)) }")
    println("  little_point_group2: { name: \"$(psym_little2.hermann_mauguin)\", order: $(group_order(psym_little2)) }")
    println("  is_psym_compatible: $(iscompatible(tsym, idx, psym))")
    println("  is_psym_little1_compatible: $(iscompatible(tsym, idx, psym_little1))")
    println("  is_psym_little2_compatible: $(iscompatible(tsym, idx, psym_little2))")
end
