# Lattice Definition

The most basic feature **LatticeTools.jl** provides is defining a lattice system, which consist of:
1. a unit cell, defined by the lattice vectors and the basis sites
2. Bravais lattice

[`UnitCell`](@ref) represents the first component, while the second component, its cluster shape and periodic boundary condition for finite a finite size lattice, makes use of [`Hypercube`](@ref).


## UnitCell

The unit cell of a lattice is represented by [`UnitCell`](@ref).
The lattice vectors, which defines the size and shape of the unit cell, is stored in a matrix as column vectors.
The basis of the unit cell is represented as a list of sites, with ther names and locations within the unit cell.
The "name" of a site can be of any arbitrary type (except integer in order to avoid confusion with the site "index").
Their locations are stored in fractional coordinates in units of the lattice vectors as objects of type [`FractCoord`](@ref).

Creating a [`UnitCell`](@ref) is straight forward, using [`makeunitcell`](@ref) and [`addsite!`](@ref).
The following example defines a two-dimensional square unit cell with two sites, site A at location [0.1, 0],
and site B at location [0.0, 0.1].
```@example example-unitcell
using LatticeTools
unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "A", FractCoord([0,0], [0.1, 0.0]))
addsite!(unitcell, "B", FractCoord([0,0], [0.0, 0.1]))
unitcell
```
Here [`FractCoord`](@ref) represents location in the fractional coordinates,
i.e. in units of the lattice vectors, with `whole` ($\in \mathbb{Z}^{D}$) and
`fraction` ($\in [0, 1)^D$) parts.


## Hypercube

A finite size lattice with periodic boundary condition can be constructed as a quotient set of an infinite lattice, with equivalence relation defined by a "supercell translation".
In such a case, the equivalence relation within the basis under translation is identity; only the equivalence relation within the Bravais lattice needs to be worked out.
[`Hypercube`](@ref) provides methods for calculating equivalence relation.

For example, when the super cell is defined by lattice vectors `[3, -1]` and `[1, 3]` (in units of the lattice vectors of the original unit cell), 
```@repl
using LatticeTools # hide
cube = Hypercube([3 1; -1 3])
cube.wrap([4, 0])
```
Here `cube.wrap` is similar to `divrem`:
`[1, 0]` is the super cell translation of `[4, 0]`, and `[1, 1]` is the remaining translation.


## Lattice

Now the two can be combined into a [`Lattice`](@ref) using [`makelattice`](@ref).

```@repl
using LatticeTools # hide
unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String) # hide
addsite!(unitcell, "A", FractCoord([0,0], [0.1, 0.0])) # hide
addsite!(unitcell, "B", FractCoord([0,0], [0.0, 0.1])) # hide
lattice = makelattice(unitcell, [3 1; -1 3])
```