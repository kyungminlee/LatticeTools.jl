# Symmetry Irrep Components

## IrrepComponent

```@docs
IrrepComponent
```

```@docs
group_order(::IrrepComponent)
get_irrep_components(::AbstractSymmetry)
get_irrep_iterator(::IrrepComponent)
```

## Little Group

```@docs
little_group_elements(::IrrepComponent{FiniteTranslationSymmetry}, ::PointSymmetry)
little_group(::IrrepComponent{FiniteTranslationSymmetry}, ::PointSymmetry)
little_group(::IrrepComponent{SymmetryEmbedding{FiniteTranslationSymmetry}}, ::SymmetryEmbedding{PointSymmetry})
little_symmetry(::IrrepComponent{FiniteTranslationSymmetry}, ::PointSymmetry)
iscompatible(::IrrepComponent{FiniteTranslationSymmetry}, ::PointSymmetry)
```


## SymmorphicIrrepComponent

```@docs
SymmorphicIrrepComponent
```

```@docs
group_order(::SymmorphicIrrepComponent)
get_irrep_components(::SymmorphicSymmetry)
get_irrep_components(::SymmorphicSymmetryEmbedding)
get_irrep_iterator(::SymmorphicIrrepComponent)
```