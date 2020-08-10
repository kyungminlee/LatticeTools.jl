# Basic

```@meta
CurrentModule = TightBindingLattice
```

## Coordinate Types

```@docs
TightBindingLattice.FractCoord
```

### Functions
```@docs
TightBindingLattice.fract2carte
TightBindingLattice.carte2fract
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
TightBindingLattice.momentumgrid(::UnitCell, ::AbstractVector{<:Integer})
TightBindingLattice.findsiteindex
```

## OrthoCube

```@docs
TightBindingLattice.OrthoCube
```

```@docs
TightBindingLattice.dimension
TightBindingLattice.volume
TightBindingLattice.isequiv
TightBindingLattice.find_generators
TightBindingLattice.find_generators_2d
TightBindingLattice.generate_coordinates
```

## Lattice

```@docs
TightBindingLattice.Lattice
TightBindingLattice.make_lattice
TightBindingLattice.makelattice
```

```@docs
TightBindingLattice.dimension
```