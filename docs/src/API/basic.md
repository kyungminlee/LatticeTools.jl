# Basics

```@meta
CurrentModule = LatticeTools
```

## Coordinate Types

```@docs
CarteCoord
FractCoord
```

### Functions
```@docs
FractCoord(::Vector{Int}, ::Vector{Float64})
FractCoord(::Vector{Float64})
FractCoord(::Int)
fract2carte(::Matrix{Float64}, ::FractCoord)
carte2fract(::Matrix{Float64}, ::CarteCoord)
```

## UnitCell

```@docs
LatticeTools.UnitCell
LatticeTools.makeunitcell
```

### Functions
```@docs
LatticeTools.dimension
LatticeTools.numsite
LatticeTools.sitecount
LatticeTools.addsite!
LatticeTools.hassite
LatticeTools.getsite
LatticeTools.getsiteindex
LatticeTools.getsitecoord
LatticeTools.getsiteindexcoord
LatticeTools.getsitename
LatticeTools.carte2fract(::UnitCell, ::CarteCoord)
LatticeTools.fract2carte(::UnitCell, ::FractCoord)
LatticeTools.whichunitcell
LatticeTools.momentumgrid(::UnitCell, ::Vector{Int})
LatticeTools.findsiteindex
```

## Hypercube

```@docs
LatticeTools.Hypercube
```

```@docs
LatticeTools.dimension(::Hypercube)
LatticeTools.volume(::Hypercube)
LatticeTools.isequiv(::Hypercube, ::Hypercube)
LatticeTools.find_generators(::Hypercube)
LatticeTools.find_generators_2d(::Hypercube)
LatticeTools.generate_coordinates(::Hypercube, ::Matrix{Int})
```

## Lattice

```@docs
LatticeTools.Lattice
LatticeTools.makelattice
```

```@docs
LatticeTools.dimension(::Lattice)
```
