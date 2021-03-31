# Symmetry

```@meta
CurrentModule = LatticeTools
```


## AbstractSymmetry

```@docs
AbstractSymmetry
```

## FiniteTranslationSymmetry

```@docs
FiniteTranslationSymmetry
```

### Constructors
```@docs
FiniteTranslationSymmetry(::Matrix{Int})
FiniteTranslationSymmetry(::Lattice)
FiniteTranslationSymmetry(::Hypercube)
FiniteTranslationSymmetry(::Hypercube, ::Matrix{Int})
```

### Common Functions for Symmetry
```@docs
dimension(::FiniteTranslationSymmetry)
elements(::FiniteTranslationSymmetry)
element(::FiniteTranslationSymmetry, ::Any)
element_names(::FiniteTranslationSymmetry)
element_name(::FiniteTranslationSymmetry, ::Any)
group(::FiniteTranslationSymmetry)
group_order(::FiniteTranslationSymmetry, ::Any)
group_multiplication_table(::FiniteTranslationSymmetry)
generator_indices(::FiniteTranslationSymmetry)
generator_elements(::FiniteTranslationSymmetry)
symmetry_name(::FiniteTranslationSymmetry)
symmetry_product(::FiniteTranslationSymmetry)
```

### Irreducible Representations
```@docs
character_table(::FiniteTranslationSymmetry)
irreps(::FiniteTranslationSymmetry)
irrep(::FiniteTranslationSymmetry, ::Any)
num_irreps(::FiniteTranslationSymmetry)
numirreps(::FiniteTranslationSymmetry)
irrepcount(::FiniteTranslationSymmetry)
irrep_dimension(::FiniteTranslationSymmetry, ::Int)
```

### Functions Specific to Translation Symmetry
```@docs
fractional_momentum(::FiniteTranslationSymmetry)
isbragg(::Vector{Int}, ::Vector{Int}, ::Vector{Int})
isbragg(::Vector{Int}, ::Vector{Int}, ::Vector{Vector{Int}})
isbragg(::Vector{Rational{Int}}, ::Vector{Int})
isbragg(::Vector{Rational{Int}}, ::Vector{Vector{Int}})
isbragg(::FiniteTranslationSymmetry, ::Int, ::TranslationOperation{Int})
isbragg(::FiniteTranslationSymmetry, ::Int, ::Vector{TranslationOperation{Int}})
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
symmetry_product(::SymmorphicSymmetry{FiniteTranslationSymmetry,PointSymmetry,SpaceOperation{Int}})
```

```@docs
group(::SymmorphicSymmetry)
group_order(::SymmorphicSymmetry)
group_multiplication_table(::SymmorphicSymmetry)
elements(::SymmorphicSymmetry)
element(::SymmorphicSymmetry, ::Any)
element_names(::SymmorphicSymmetry)
element_name(::SymmorphicSymmetry, ::Any)
fractional_momentum(::SymmorphicSymmetry{FiniteTranslationSymmetry, PointSymmetry, SpaceOperation})
generator_indices(::SymmorphicSymmetry)
generator_elements(::SymmorphicSymmetry)
```
