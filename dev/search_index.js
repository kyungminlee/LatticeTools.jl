var documenterSearchIndex = {"docs":
[{"location":"#TightBindingLattice-1","page":"TightBindingLattice","title":"TightBindingLattice","text":"","category":"section"},{"location":"#","page":"TightBindingLattice","title":"TightBindingLattice","text":"Implements tight binding","category":"page"},{"location":"#","page":"TightBindingLattice","title":"TightBindingLattice","text":"Modules = [TightBindingLattice]","category":"page"},{"location":"#TightBindingLattice.CarteCoord","page":"TightBindingLattice","title":"TightBindingLattice.CarteCoord","text":"CarteCoord\n\nCartesian coordinates. Vector{Float64}.\n\n\n\n\n\n","category":"type"},{"location":"#TightBindingLattice.FractCoord","page":"TightBindingLattice","title":"TightBindingLattice.FractCoord","text":"FractCoord\n\nFractional coordinates.\n\nMembers\n\nwhole ::Vector{Int}: Integer part of fractional coordinates\nfraction ::Vector{Float64}: [0,1) part of fractional coordinates\n\n\n\n\n\n","category":"type"},{"location":"#TightBindingLattice.UnitCell","page":"TightBindingLattice","title":"TightBindingLattice.UnitCell","text":"UnitCell{O}\n\nParameters\n\nO: type of \"orbital\". Any type can be used, but we recommend using String or tuple of String and Int      for compatibility with JSON.\n\nMembers\n\nlatticevectors ::Array{Float64, 2}: Lattice vectors\nreducedreciprocallatticevectors ::Array{Float64, 2}: Reduced reciprocal lattice vectors (transpose of inverse of latticevectors)\nreciprocallatticevectors ::Array{Float64, 2}: Reciprocal lattice vectors. 2π * reducedreciprocallatticevectors\norbitals ::Vector{Tuple{T, FractCoord}}: List of orbitals within unit cell\norbitalindices ::Dict{T, Int}: Indices of orbitals\n\n\n\n\n\n","category":"type"},{"location":"#TightBindingLattice.addorbital!-Union{Tuple{O}, Tuple{UnitCell{O},O,FractCoord}} where O","page":"TightBindingLattice","title":"TightBindingLattice.addorbital!","text":"addorbital!\n\nAdd an orbital to the unit cell.\n\nArguments\n\nuc ::UnitCell{T}\norbitalname ::{T}\norbitalcoord ::FractCoord\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.carte2fract-Tuple{AbstractArray{#s23,2} where #s23<:AbstractFloat,Array{Float64,1}}","page":"TightBindingLattice","title":"TightBindingLattice.carte2fract","text":"carte2fract\n\nArguments\n\nlatticevectors ::AbstractArray{<:AbstractFloat, 2}: square matrix whose columns are lattice vectors.\ncc ::CarteCoord: cartesian coordinates\ntol ::Real=sqrt(eps(Float64)): tolerance\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.carte2fract-Tuple{UnitCell,Array{Float64,1}}","page":"TightBindingLattice","title":"TightBindingLattice.carte2fract","text":"carte2fract\n\nArguments\n\nlatticevectors ::Array{Float64, 2}\ncc ::CarteCoord\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.dimension-Tuple{FractCoord}","page":"TightBindingLattice","title":"TightBindingLattice.dimension","text":"dimension\n\nDimension of the fractional coordinates\n\nArguments\n\nfc ::FractCoord: Fractional coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.dimension-Tuple{UnitCell}","page":"TightBindingLattice","title":"TightBindingLattice.dimension","text":"dimension\n\nSpatial dimension of the unit cell.\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.fract2carte-Tuple{AbstractArray{#s31,2} where #s31<:AbstractFloat,FractCoord}","page":"TightBindingLattice","title":"TightBindingLattice.fract2carte","text":"fract2carte\n\nArguments\n\nlatticevectors ::AbstractArray{<:AbstractFloat, 2}: square matrix whose columns are lattice vectors.\nfc ::FractCoord: fractional coordinates\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.fract2carte-Tuple{UnitCell,FractCoord}","page":"TightBindingLattice","title":"TightBindingLattice.fract2carte","text":"fract2carte\n\nArguments\n\nlatticevectors ::Array{Float64, 2}\nfc ::FractCoord\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbital-Tuple{UnitCell,Integer}","page":"TightBindingLattice","title":"TightBindingLattice.getorbital","text":"getorbital\n\nArguments\n\nuc ::UnitCell{T}\nindex ::Integer\n\nReturn\n\n(orbitalname, fractcoord)\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbital-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O","page":"TightBindingLattice","title":"TightBindingLattice.getorbital","text":"getorbital\n\nGet the orbital (its orbital name and its fractional coordinates) with the given name.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\nReturn\n\n(orbitalname, fractcoord)\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbitalcoord-Tuple{UnitCell,Integer}","page":"TightBindingLattice","title":"TightBindingLattice.getorbitalcoord","text":"getorbitalcoord\n\nArguments\n\nuc ::UnitCell\nidx ::Integer\n\nReturn\n\nfractcoord\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbitalcoord-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O","page":"TightBindingLattice","title":"TightBindingLattice.getorbitalcoord","text":"getorbitalcoord\n\nGet the fractional coordinates of the orbital with the given name.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\nReturn\n\nfractcoord\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbitalindex-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O","page":"TightBindingLattice","title":"TightBindingLattice.getorbitalindex","text":"getorbitalindex\n\nGet index of the given orbital.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbitalindexcoord-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O","page":"TightBindingLattice","title":"TightBindingLattice.getorbitalindexcoord","text":"getorbitalindexcoord\n\nArguments\n\nuc ::UnitCell{T}\nname ::T\n\nReturn\n\n(index, fractcoord)\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.getorbitalname-Tuple{UnitCell,Integer}","page":"TightBindingLattice","title":"TightBindingLattice.getorbitalname","text":"getorbitalname\n\nArguments\n\nuc ::UnitCell\nindex ::Integer\n\nReturn\n\norbitalname\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.hasorbital-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O","page":"TightBindingLattice","title":"TightBindingLattice.hasorbital","text":"hasorbital{T}\n\nTest whether the unit cell contains the orbital of given name.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.linpath-Tuple{Array{Array{Float64,1},1}}","page":"TightBindingLattice","title":"TightBindingLattice.linpath","text":"momentumpath\n\nGenerate a list of momenta\n\nArguments\n\nanchorpoints\n(Optional) nseg - number of points in each segment\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.make_supercell-Union{Tuple{O}, Tuple{UnitCell{O},Array{#s28,2} where #s28<:Integer}} where O","page":"TightBindingLattice","title":"TightBindingLattice.make_supercell","text":"Returns\n\nnew unit cell\nembedding function\nwhich takes the unitcell displacement in the orignal displacement\nwhich returns a tuple of (supercell lattice displacement, sublattice subscript)\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.make_unitcell-Tuple{AbstractArray{#s28,2} where #s28<:AbstractFloat}","page":"TightBindingLattice","title":"TightBindingLattice.make_unitcell","text":"UnitCell\n\nConstruct an n-dimensional lattice.\n\nArguments\n\nlatticevectors ::AbstractArray{<:AbstractFloat, 2}: Lattice vectors\nOrbitalType::DataType\n\nOptional Arguments\n\ntol=sqrt(eps(Float64)): Epsilon\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.make_unitcell-Tuple{AbstractFloat}","page":"TightBindingLattice","title":"TightBindingLattice.make_unitcell","text":"UnitCell\n\nConstruct a one-dimensional lattice.\n\nArguments\n\nlatticeconstant ::Float64: Lattice constant\nOrbitalType: List of orbitals\n\nOptional Arguments\n\ntol=sqrt(eps(Float64)): Tolerance\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.momentumgrid-Tuple{UnitCell,AbstractArray{#s33,1} where #s33<:Integer}","page":"TightBindingLattice","title":"TightBindingLattice.momentumgrid","text":"momentumgrid\n\nGenerate an n-dimensional grid of momenta of given shape\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.momentumpath-Tuple{UnitCell,Array{Array{Float64,1},1}}","page":"TightBindingLattice","title":"TightBindingLattice.momentumpath","text":"momentumpath\n\nThe anchorpoints are given in units of the reciprocal lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.numorbital-Tuple{UnitCell}","page":"TightBindingLattice","title":"TightBindingLattice.numorbital","text":"numorbital\n\nNumber of orbitals of the unit cell.\n\nArguments\n\nuc ::UnitCell\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.whichunitcell-Union{Tuple{O}, Tuple{UnitCell{O},O,Array{Float64,1}}} where O","page":"TightBindingLattice","title":"TightBindingLattice.whichunitcell","text":"whichunitcell\n\nReturn\n\nR ::Vector{Int}: which unit cell the specificied orbital/cartesian coordinates belongs to.\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.whichunitcell-Union{Tuple{O}, Tuple{UnitCell{O},O,FractCoord}} where O","page":"TightBindingLattice","title":"TightBindingLattice.whichunitcell","text":"whichunitcell\n\nReturn\n\nR ::Vector{Int}: which unit cell the specificied orbital/cartesian coordinates belongs to.\n\n\n\n\n\n","category":"method"},{"location":"#Base.isapprox-Tuple{FractCoord,FractCoord}","page":"TightBindingLattice","title":"Base.isapprox","text":"isapprox(x, y; rtol::Real=atol>0 ? 0 : √eps, atol::Real=0, nans::Bool=false, norm::Function)\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingLattice.hypercubic_cluster-Tuple{Array{#s24,2} where #s24<:Integer}","page":"TightBindingLattice","title":"TightBindingLattice.hypercubic_cluster","text":"hypercubic_cluster\n\nGenerate a hypercubic cluster\n\n. . . . . .   . . . o . .  . o . . . .  . . . . o .  . . o . . .  . . . . . .\n\n\n\n\n\n","category":"method"}]
}