
using TightBindingLattice
tsym = TranslationSymmetry([4 0; 0 4])
psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])
ssym = tsym ⋊ psym

# @show ssym

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
addorbital!(unitcell, "A", FractCoord([0, 0], [0.5, 0.0]))
addorbital!(unitcell, "B", FractCoord([0, 0], [0.0, 0.5]))
lattice = make_lattice(unitcell, [4 0; 0 4])

ssymbed = embed(lattice, ssym)

for ssic in get_irrep_components(ssymbed)
    L = length(collect(get_irrep_iterator(ssic)))
    println(symmetry_name(ssic.normal_symmetry), "\t", symmetry_name(ssic.rest_symmetry), "\t", L)
end