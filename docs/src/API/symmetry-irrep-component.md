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
little_group_elements(::IrrepComponent{TranslationSymmetry}, ::PointSymmetry)
little_group(::IrrepComponent{TranslationSymmetry}, ::PointSymmetry)
little_group(::IrrepComponent{SymmetryEmbedding{TranslationSymmetry}}, ::SymmetryEmbedding{PointSymmetry})
little_symmetry(::IrrepComponent{TranslationSymmetry}, ::PointSymmetry)
iscompatible(::IrrepComponent{TranslationSymmetry}, ::PointSymmetry)
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