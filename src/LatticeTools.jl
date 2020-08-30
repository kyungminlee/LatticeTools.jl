module LatticeTools

include("Basic/basic.jl")
include("Basic/coord.jl")
include("Basic/orthocube.jl")  # replacement for hypercube
include("Basic/unitcell.jl")
include("Basic/lattice.jl")

include("Group/abstractgroup.jl")
include("Group/finitegroup.jl")
include("Group/permutation.jl")

include("SymmetryOperation/abstractoperation.jl")
include("SymmetryOperation/identityoperation.jl")
include("SymmetryOperation/translationoperation.jl")
include("SymmetryOperation/pointoperation.jl")
include("SymmetryOperation/spaceoperation.jl")
# include("SymmetryOperation/productoperation.jl")

include("Symmetry/abstractsymmetry.jl")
include("Symmetry/translationsymmetry.jl")
include("Symmetry/pointsymmetry.jl")
include("Symmetry/compatibility.jl")
include("Symmetry/littlesymmetry.jl")

include("SymmetryEmbedding/sitepermutation.jl")
include("SymmetryEmbedding/sitemap.jl")
include("SymmetryEmbedding/localunitary.jl")
include("SymmetryEmbedding/symmetryembedding.jl")

include("Irrep/irrep.jl")

export PointSymmetryDatabase
include("Symmetry/pointsymmetrydatabase.jl")
export IrrepDatabase
include("Irrep/irrepdatabase.jl")

include("Symmorphic/symmorphicsymmetry.jl")
include("Symmorphic/symmorphicsymmetryembedding.jl")
include("Symmorphic/symmorphicirrep.jl")


include("momentumpath.jl")

end # module
