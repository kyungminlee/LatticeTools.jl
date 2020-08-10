# Symmetry Operations

`AbstractSpaceSymmetryOperation{S}` (where `S` is the type for the coordinates)
represent spatial symmetry, including translation, and point operation, and the
combination of the two.

## Translation Operation

The type `TranslationOperation{S}` represents translation operation on the space
of coordinates of type `S`.

```
julia> using TightBindingLattice

julia> TranslationOperation([0.5, 0.0, 0.5])
TranslationOperation{Float64}([0.5, 0.0, 0.5])

julia> TranslationOperation([1, 2])
TranslationOperation{Int64}([1, 2])
```

```
julia> using TightBindingLattice

julia> t = TranslationOperation([1, 2])
TranslationOperation{Int64}([1, 2])

julia> t([3, 4])
2-element Array{Int64,1}:
 4
 6

julia> apply_operation(t, [3, 4])
2-element Array{Int64,1}:
 4
 6
```

## Point Operation

```@example
julia> p = PointOperation([0 1; 1 0])
PointOperation{Int64}([0 1; 1 0])

julia> p([3, 4])
2-element Array{Int64,1}:
 4
 3
```