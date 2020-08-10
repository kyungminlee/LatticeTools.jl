# Symmetry Operations

`AbstractSpaceSymmetryOperation{S}` (where `S` is the type for the coordinates)
represent spatial symmetry, including translation, and point operation, and the
combination of the two.

```@setup tbl
using TightBindingLattice
```

## Translation Operation

The type `TranslationOperation{S}` represents translation operation on the space
of coordinates of type `S`.

A translation operation can be of integer type.
```@example tbl
TranslationOperation([1, 2])
```

It can also represent translation of non-integer type
```@example tbl
TranslationOperation([0.5, 0.0, 0.5])
```

```@example tbl
TranslationOperation([1//2, 3//4])
```

A translation operation can be applied to a coordinate
```@example tbl
t = TranslationOperation([1, 2])
@show t([3, 4])
@show apply_operation(t, [3, 4])
```

## Point Operation

```@example tbl
p = PointOperation([0 1; 1 0])
@show p([3, 4])
```