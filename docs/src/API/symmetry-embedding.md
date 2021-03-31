# Symmetry Embedding

## SitePermutation

```@docs
SitePermutation
```

```@docs
embed(::Lattice, ::TranslationOperation{Int})
embed(::Lattice, ::PointOperation{Int})
embed(::Lattice, ::SpaceOperation{Int, Int})
isidentity(::SitePermutation)
```


## SymmetryEmbedding

```@docs
SymmetryEmbedding
```

```@docs
elements(::SymmetryEmbedding)
element(::SymmetryEmbedding)
symmetry(::SymmetryEmbedding)
generator_elements(::SymmetryEmbedding)
```

```@docs
group_order(::SymmetryEmbedding, ::Any)
group_multiplication_table(::SymmetryEmbedding, ::Any)
generator_indices(::SymmetryEmbedding, ::Any)
element_names(::SymmetryEmbedding, ::Any)
element_name(::SymmetryEmbedding, ::Any)
character_table(::SymmetryEmbedding, ::Any)
irreps(::SymmetryEmbedding, ::Any)
irrep(::SymmetryEmbedding, ::Any)
irrep_dimension(::SymmetryEmbedding, ::Any)
num_irreps(::SymmetryEmbedding, ::Any)
numirreps(::SymmetryEmbedding, ::Any)
irrepcount(::SymmetryEmbedding, ::Any)
```

```@docs
fractional_momentum(::SymmetryEmbedding{FiniteTranslationSymmetry})
iscompatible(::SymmetryEmbedding{FiniteTranslationSymmetry}, ::SymmetryEmbedding{PointSymmetry})
iscompatible(::SymmetryEmbedding{FiniteTranslationSymmetry}, ::Int, ::SymmetryEmbedding{PointSymmetry})
little_group_elements(::SymmetryEmbedding{FiniteTranslationSymmetry}, ::SymmetryEmbedding{PointSymmetry})
little_group_elements(::SymmetryEmbedding{FiniteTranslationSymmetry}, ::Int, ::SymmetryEmbedding{PointSymmetry})
little_symmetry(::SymmetryEmbedding{FiniteTranslationSymmetry}, ::SymmetryEmbedding{PointSymmetry})
little_symmetry(::SymmetryEmbedding{FiniteTranslationSymmetry}, ::Int, ::SymmetryEmbedding{PointSymmetry})
symmetry_name(::SymmetryEmbedding)
```

## Embedding a Symmetry onto a Lattice

```@docs
embed(::Lattice, ::FiniteTranslationSymmetry) 
embed(::Lattice, ::PointSymmetry)
translation_symmetry_embedding(::Lattice)
```

## findsitemap

```@docs
findsitemap
```

