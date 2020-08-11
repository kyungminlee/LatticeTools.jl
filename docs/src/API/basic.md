# Basics

```@meta
CurrentModule = TightBindingLattice
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
TightBindingLattice.UnitCell
TightBindingLattice.makeunitcell
```

### Functions
```@docs
TightBindingLattice.dimension
TightBindingLattice.numsite
TightBindingLattice.sitecount
TightBindingLattice.addsite!
TightBindingLattice.hassite
TightBindingLattice.getsite
TightBindingLattice.getsiteindex
TightBindingLattice.getsitecoord
TightBindingLattice.getsiteindexcoord
TightBindingLattice.getsitename
TightBindingLattice.carte2fract(::UnitCell, ::CarteCoord)
TightBindingLattice.fract2carte(::UnitCell, ::FractCoord)
TightBindingLattice.whichunitcell
TightBindingLattice.momentumgrid(::UnitCell, ::Vector{Int})
TightBindingLattice.findsiteindex
```

## OrthoCube

```@docs
TightBindingLattice.OrthoCube
```

```@docs
TightBindingLattice.dimension(::OrthoCube)
TightBindingLattice.volume(::OrthoCube)
TightBindingLattice.isequiv(::OrthoCube, ::OrthoCube)
TightBindingLattice.find_generators(::OrthoCube)
TightBindingLattice.find_generators_2d(::OrthoCube)
TightBindingLattice.generate_coordinates(::OrthoCube, ::Matrix{Int})
```

## Lattice

```@docs
TightBindingLattice.Lattice
TightBindingLattice.makelattice
```

```@docs
TightBindingLattice.dimension(::Lattice)
```
