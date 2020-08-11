# Defining a Lattice

A **lattice** is defined by the following components:
1. lattice vectors
2. sites
3. Bravais lattice

The components 1 and 2 are stored in [`UnitCell`](@ref), while the component 3
makes use of [`OrthoCube`](@ref).


## UnitCell

A [`UnitCell`](@ref) contains the lattice vectors, stored as a matrix whose columns are
unit translation vectors, and list of sites, with their names and locations.
The "name" can be of any arbitrary type (except integer, to avoid confusion
with the site index), and the locations are stored as type [`FractCoord`](@ref).

Creating a [`UnitCell`](@ref) is straight forward: The following example defines
a two-dimensional square unit cell with two sites, site A at location [0.1, 0],
and site B at location [0.0, 0.1].

```@example example-unitcell
using TightBindingLattice
unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "A", FractCoord([0,0], [0.1, 0.0]))
addsite!(unitcell, "B", FractCoord([0,0], [0.0, 0.1]))
unitcell
```

Here [`FractCoord`](@ref) represents location in the fractional coordinates,
i.e. in units of the lattice vectors, with `whole` ($\in \mathbb{Z}^{D}$) and
`fraction` ($\in [0, 1)^D$) parts.


## OrthoCube

The Bravais lattice is represented by [`OrthoCube`](@ref).


## Lattice

Now the two can be combined into a [`Lattice`](@ref).