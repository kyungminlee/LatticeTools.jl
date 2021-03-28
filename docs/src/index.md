# Welcome to LatticeTools.jl


## Overview

**LatticeTools.jl** is a Julia package for defining lattices and performing space symmetry analysis on them, useful for creating and solving quantum Hamiltonians.

## Features
- **Unit cell definition**: Supports definition of unit cells with arbitrary shape in arbitrary dimension and with basis sites. Orbitals within a site is not yet implemented.
- **Lattice definition**: Supports definition of a lattice with a unit cell with periodic boundary condition.
- **Momentum space**: Provides methods for transformations between real space and momentum space.
- **Symmetry analysis**: Supports group representation theoretical analysis on translation and point symmetry. Currently, only symmorphic space groups are supported.


## Installation

`LatticeTools` is currently not included in Julia's default package registry.
To install, add the package registry `KyugminLeeRegistry` and then install the package:
```julia-repl
(@v1.5) pkg> registry add https://github.com/kyungminlee/KyungminLeeRegistry.git
(@v1.5) pkg> add LatticeTools
```


## Quick Example

```@example example-unitcell
using LatticeTools

unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))
lattice = makelattice(unitcell, [4 0; 0 4]);
tsym = FiniteTranslationSymmetry(lattice);
psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0]);

println("Irreps and Little Groups")
println("------------------------")
for tsic in get_irrep_components(tsym)
    idx = tsic.irrep_index
    kf = tsym.fractional_momenta[idx]
    k = lattice.unitcell.reducedreciprocallatticevectors * kf
    psym_little = little_symmetry(tsym, idx, psym)
    println("- irrep_index: $(idx)")
    println("  momentum: $(k)")
    println("  little_point_group: { name: \"$(psym_little.hermann_mauguin)\", order: $(group_order(psym_little)) }")
    println("  is_psym_compatible: $(iscompatible(tsym, idx, psym))")
    println("  is_psym_little_compatible: $(iscompatible(tsym, idx, psym_little))")
end
```