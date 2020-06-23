# Lattice

Definition of lattice involves
1. definition of unitcell with lattice vectors
2. definition of sites within a unitcell.
3. definition of Bravais lattice

## UnitCell and Sites

```@example example-unitcell
using TightBindingLattice
unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "A", FractCoord([0,0], [0.1, 0.0]))
addsite!(unitcell, "B", FractCoord([0,0], [0.0, 0.1]))
unitcell
```

## OrthoCube


