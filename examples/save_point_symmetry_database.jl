using TightBindingLattice
using JLD2


POINT_SYMMETRY_DATABASE = [PointSymmetryDatabase.get(gn) for gn in 1:32]
@save "PointGroup3D.jld2" POINT_SYMMETRY_DATABASE
