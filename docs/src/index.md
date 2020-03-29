# TightBindingLattice

Lattice.


`Group` has group structure. Examples are `FiniteAbelianGroup`, `GenericGroup`
- Multiplication table
- order (period length) of every element


Symmetry is the "presentation" of group.


`Symmetry` on the other hand, is a representation of the group structure in some space.
Examples are: `TranslationSymmetry`, `PointSymmetry`. Members are
- group
- generators
- conjugacy_classes
- character_table
- irreps
- element names

`PointSymmetry` has additional info
- matrix_representations (i.e. representation in units of lattice vectors)
- Schoenflies
- Hermann Mauguinn







```@autodocs
Modules = [TightBindingLattice]
```
