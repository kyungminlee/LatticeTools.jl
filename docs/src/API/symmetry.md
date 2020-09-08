# Symmetry

```@meta
CurrentModule = LatticeTools
```


## AbstractSymmetry

```@docs
AbstractSymmetry
```

## TranslationSymmetry

```@docs
TranslationSymmetry
```

### Constructors
```@docs
TranslationSymmetry(::Matrix{Int})
TranslationSymmetry(::Lattice)
TranslationSymmetry(::Hypercube)
TranslationSymmetry(::Hypercube, ::Matrix{Int})
```

### Common Functions for Symmetry
```@docs
dimension(::TranslationSymmetry)
elements(::TranslationSymmetry)
element(::TranslationSymmetry, ::Any)
element_names(::TranslationSymmetry)
element_name(::TranslationSymmetry, ::Any)
group(::TranslationSymmetry)
group_order(::TranslationSymmetry, ::Any)
group_multiplication_table(::TranslationSymmetry)
generator_indices(::TranslationSymmetry)
generator_elements(::TranslationSymmetry)
symmetry_name(::TranslationSymmetry)
symmetry_product(::TranslationSymmetry)
```

### Irreducible Representations
```@docs
character_table(::TranslationSymmetry)
irreps(::TranslationSymmetry)
irrep(::TranslationSymmetry, ::Any)
num_irreps(::TranslationSymmetry)
numirreps(::TranslationSymmetry)
irrepcount(::TranslationSymmetry)
irrep_dimension(::TranslationSymmetry, ::Int)
```

### Functions Specific to Translation Symmetry
```@docs
fractional_momentum(::TranslationSymmetry)
isbragg(::Vector{Int}, ::Vector{Int}, ::Vector{Int})
isbragg(::Vector{Int}, ::Vector{Int}, ::Vector{Vector{Int}})
isbragg(::Vector{Rational{Int}}, ::Vector{Int})
isbragg(::Vector{Rational{Int}}, ::Vector{Vector{Int}})
isbragg(::TranslationSymmetry, ::Int, ::TranslationOperation{Int})
isbragg(::TranslationSymmetry, ::Int, ::Vector{TranslationOperation{Int}})
```

## PointSymmetry

```@docs
PointSymmetry
```

### Common Functions for Symmetry
```@docs
dimension(::PointSymmetry)
elements(::PointSymmetry)
element(::PointSymmetry, ::Any)
element_names(::PointSymmetry)
element_name(::PointSymmetry, ::Any)
group(::PointSymmetry)
group_order(::PointSymmetry)
group_multiplication_table(::PointSymmetry)
generator_indices(::PointSymmetry)
generator_elements(::PointSymmetry)
symmetry_name(::PointSymmetry)
symmetry_product(::PointSymmetry)
```

### Irreducible Representations
```@docs
character_table(::PointSymmetry)
irreps(::PointSymmetry)
irrep(::PointSymmetry, ::Int)
num_irreps(::PointSymmetry)
numirreps(::PointSymmetry)
irrepcount(::PointSymmetry)
irrep_dimension(::PointSymmetry, ::Int)
```

### Projection to Subspace
```@docs
project(::PointSymmetry, ::Matrix{Int})
```

## Little Symmetry

```@docs
little_group_elements
little_group
little_symmetry
```


## Compatibility

```@docs
iscompatible
```


## PointSymmetryDatabase

```@docs
LatticeTools.PointSymmetryDatabase
LatticeTools.PointSymmetryDatabase.get2d
LatticeTools.PointSymmetryDatabase.get3d
LatticeTools.PointSymmetryDatabase.find2d
LatticeTools.PointSymmetryDatabase.find3d
```


## IrrepDatabase

```@docs
LatticeTools.IrrepDatabase
```


## Symmorphic Symmetry

```@docs
SymmorphicSymmetry
```

## Operators

```@docs
⋊
⋉
```

```@docs
symmetry_product(::SymmorphicSymmetry{TranslationSymmetry,PointSymmetry,SpaceOperation{Int}})
```

```@docs
group(::SymmorphicSymmetry)
group_order(::SymmorphicSymmetry)
group_multiplication_table(::SymmorphicSymmetry)
elements(::SymmorphicSymmetry)
element(::SymmorphicSymmetry, ::Any)
element_names(::SymmorphicSymmetry)
element_name(::SymmorphicSymmetry, ::Any)
fractional_momentum(::SymmorphicSymmetry{TranslationSymmetry, PointSymmetry, SpaceOperation})
generator_indices(::SymmorphicSymmetry)
generator_elements(::SymmorphicSymmetry)
```
