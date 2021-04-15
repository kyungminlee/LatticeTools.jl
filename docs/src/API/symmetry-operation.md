# Symmetry Operation


## AbstractSpaceSymmetryOperation
```@docs
AbstractSpaceSymmetryOperation
domaintype(::AbstractSpaceSymmetryOperation)
```

## TranslationOperation
```@docs
TranslationOperation
```

### Constructors
```@docs
TranslationOperation{Int}(::Vector{Int})
TranslationOperation(::Vector{Int})
```

### Properties
```@docs
isidentity(::TranslationOperation)
istranslation(::TranslationOperation)
ispoint(::TranslationOperation)
dimension(::TranslationOperation)
```

### Apply
```@docs
apply_operation(::TranslationOperation{Int}, ::Vector{Int})
```


## PointOperation
```@docs
PointOperation
```

### Constructors
```@docs
PointOperation(::Matrix{Int})
```

### Properties
```@docs
isidentity(::PointOperation)
istranslation(::PointOperation)
ispoint(::PointOperation)
dimension(::PointOperation)
```

### Apply
```@docs
apply_operation(::PointOperation{Int}, ::Vector{Int})
```


## SpaceOperation
```@docs
SpaceOperation
```

### Constructors
```@docs
SpaceOperation(::Matrix{Int}, ::Vector{Int})
SpaceOperation(::PointOperation{Int}, ::TranslationOperation{Int})
```

### Properties
```@docs
isidentity(::SpaceOperation)
istranslation(::SpaceOperation)
ispoint(::SpaceOperation)
dimension(::SpaceOperation)
```

### Apply
```@docs
apply_operation(::SpaceOperation{Int}, ::Vector{Int})
```

