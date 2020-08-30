using Test
using LatticeTools

include("test_basic.jl")
include("test_parser.jl")
include("test_coord.jl")
include("test_unitcell.jl")
#include("test_unitcell_deprecated.jl")
include("test_orthocube.jl")
include("test_lattice.jl")
include("test_permutation.jl")

include("test_group.jl")

include("test_symmetryoperation.jl")

include("test_translationsymmetry.jl")
include("test_pointsymmetry.jl")

include("test_irrep.jl")
include("test_irrepdatabase.jl")

include("test_sitemap.jl")
include("test_sitepermutation.jl")
include("test_symmetryembedding.jl")

include("test_compatibility.jl")

include("test_symmorphic.jl")

include("test_momentumpath.jl")
