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

## OrthoCube

```@docs
LatticeTools.OrthoCube
```

```@docs
LatticeTools.dimension(::OrthoCube)
LatticeTools.volume(::OrthoCube)
LatticeTools.isequiv(::OrthoCube, ::OrthoCube)
LatticeTools.find_generators(::OrthoCube)
LatticeTools.find_generators_2d(::OrthoCube)
LatticeTools.generate_coordinates(::OrthoCube, ::Matrix{Int})
```

## Lattice

```@docs
LatticeTools.Lattice
LatticeTools.makelattice
```

```@docs
LatticeTools.dimension(::Lattice)
```
