# Symmetry Operations

[`AbstractSpaceSymmetryOperation{S}`](@ref) represent spatial symmetry,
including translation, and point operation, and the combination of the two.
Here `S` is the type for the coordinates.

```@setup tbl
using LatticeTools
```

## Translation Operation

The type [`TranslationOperation{S}`](@ref) represents translation operation
on the space of coordinates of type `S`.

A translation operation can be of integer type or non-integer type.
```@repl tbl
TranslationOperation([1, 2])
TranslationOperation([0.5, 0.0, 0.5])
TranslationOperation([1//2, 3//4])
```
On a lattice system with discrete translation symmetry, it is convenient to use
integer translations of type [`TranslationOperation{Int}`](@ref),
which represent lattice translations.

A translation operation can be applied to coordinates
```@repl tbl
t = TranslationOperation([1, 2])
t([3, 4])
apply_operation(t, [5, 6])
```

## Point Operation

Point operation is represented by the type [`PointOperation{S}`](@ref),
which includes rotation, inversion, and mirror operations.
A lattice system where lattice vectors are not orthogonal,
the matrix of the point operation need not be an orthogonal matrix.
For example, a sixfold rotation on a triangular Bravais lattice
may be written as the following

![Sixfold rotation](rotation.png)

```@repl tbl
p = PointOperation([1 -1; 1 0])
p([1, 0])
p([0, 1])
```

```@repl tbl
p = PointOperation([0 1; 1 0])
p([2, 5])
```

## Space Operation

```@repl tbl
s = SpaceOperation([1 -1; 1 0], [1, 0])
s([1, 0])
```