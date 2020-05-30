# Lattice

Definition of lattice involves
1. definition of unitcell with lattice vectors
2. definition of orbitals within a unitcell.
3. definition of Bravais lattice

## UnitCell and Orbitals

```@example example-unitcell
using TightBindingLattice
unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
addorbital!(unitcell, "A", FractCoord([0,0], [0.1, 0.0]))
addorbital!(unitcell, "B", FractCoord([0,0], [0.0, 0.1]))
unitcell
```

## OrthoCube


