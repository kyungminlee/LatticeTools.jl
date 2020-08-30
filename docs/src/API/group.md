# Group

```@meta
CurrentModule = LatticeTools
```

## AbstractGroup

```@docs
AbstractGroup
```


## FiniteGroup

```@docs
FiniteGroup
```

```@docs
element(::FiniteGroup, ::Any)
elements(::FiniteGroup)
element_name(::FiniteGroup, ::Any)
element_names(::FiniteGroup)
group_order(::FiniteGroup)
group_order(::FiniteGroup, ::Any)
period_length(::FiniteGroup, ::Any)
group_multiplication_table(::FiniteGroup)
isabelian(::FiniteGroup)
group_product(::FiniteGroup)
group_inverse(::FiniteGroup)
LatticeTools.group_inverse(::FiniteGroup, ::Int)
conjugacy_class(::FiniteGroup, ::Int)
generate_subgroup(::FiniteGroup, ::Int)
issubgroup(::FiniteGroup, ::Set{Int})
minimal_generating_set(::FiniteGroup)
group_isomorphism(::FiniteGroup, ::FiniteGroup)
group_multiplication_table(::Vector, ::Function)
ishomomorphic(::FiniteGroup, ::Vector)
```

## Permutation

```@docs
Permutation
Base.:(*)(::Permutation, ::Permutation)
generate_group(::Permutation...)
isidentity(::Permutation)
```