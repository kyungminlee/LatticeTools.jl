export iscompatible
## 1. Bragg condition

"""
    iscompatible(tsym, tsym_irrep_index, identity_translation)

Test whether the identity translation is compatible with the irreducible representation
of the translation symmetry, i.e. (1) lattice is compatible with translation symmetry, 
and (2) 
"""
function iscompatible(tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      identity_translation::TranslationOperation{<:Integer})
    orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
    orthogonal_shape = tsym.orthogonal_shape
    return isbragg(orthogonal_shape, orthogonal_momentum, identity_translation.displacement)
end

function iscompatible(tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      identity_translations::AbstractVector{<:TranslationOperation{<:Integer}})
    orthogonal_shape = tsym.orthogonal_shape
    orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
    return all(isbragg(orthogonal_shape, orthogonal_momentum, t.displacement)
                   for t in identity_translations)
end


## 2. Hypercube and Operation

function iscompatible(hypercube::HypercubicLattice, op::TranslationOperation{<:Integer})::Bool
    return true
end

"""
    iscompatible(hypercube, point_operation)

    (R ⋅ L) mod L = 0   and   (R ⋅ L) / L = 1
"""
function iscompatible(hypercube::HypercubicLattice, op::PointOperation{<:Integer})::Bool
    R, r = hypercube.wrap(op.matrix * hypercube.scale_matrix)
    return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
end


## 3. Hypercube and Symmetry

function iscompatible(hypercube::HypercubicLattice, tsym::TranslationSymmetry)::Bool
    R, r = iszero(hypercube.wrap(tsym.hypercube.scale_matrix))
    return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
end

function iscompatible(hypercube::HypercubicLattice, psym::PointSymmetry)::Bool
    return all(iscompatible(hypercube, pop) for pop in elements(psym))
end


## 4. TranslationSymmetry and PointOperation/PointSymmetry

function iscompatible(tsym::TranslationSymmetry, pop::PointOperation{<:Integer})::Bool
    # sm = tsym.hypercube.scale_matrix
    # smi = tsym.hypercube.inverse_scale_matrix
    # return all(mod(x,1) == 0 for x in smi * pop.matrix * sm)
    return iscompatible(tsym.hypercube, pop)
end

function iscompatible(tsym::TranslationSymmetry, psym::PointSymmetry)::Bool
    return iscompatible(tsym.hypercube, psym)
end


## 5. Translation Irreps and PointOperation/PointSymmetry

function iscompatible(tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      pop::PointOperation{<:Integer})::Bool
    !iscompatible(tsym, pop) && return false
    reciprocal_matrix = ExactLinearAlgebra.inverse(transpose(pop.matrix))
    k1 = tsym.fractional_momenta[tsym_irrep_index]
    k2 = (x -> mod(x, 1)).(reciprocal_matrix * k1)
    return k1 == k2
end

function iscompatible(tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      psym::PointSymmetry)::Bool
    return all(iscompatible(tsym, tsym_irrep_index, pop) for pop in generator_elements(psym))
end


## 6. Lattice and Operation

function iscompatible(lattice::Lattice, op::TranslationOperation{<:Integer})
    return true
end

function iscompatible(lattice::Lattice, op::PointOperation{<:Integer})
    return iscompatible(lattice.hypercube, op) &&
           !isnothing(findorbitalmap(lattice.unitcell, op))
end


## 7. Lattice and Symmetry

"""
    iscompatible(lattice, translation_symmetry)

Test whether lattice and the symmetry are compatible.
For translation symmetry, this means that the hypercubic lattice for both are the same.
"""
function iscompatible(lattice::Lattice, tsym::TranslationSymmetry)
    return lattice.hypercube == tsym.hypercube
end


function iscompatible(lattice::Lattice, psym::PointSymmetry)::Bool
    return iscompatible(lattice.hypercube, psym) &&
           all(!isnothing(findorbitalmap(lattice.unitcell, pop)) for pop in generator_elements(psym))
end





















# function iscompatible(# lattice::Lattice,
#                       tsym::TranslationSymmetry,
#                       tsym_irrep_index::Integer,
#                       identity_translation::AbstractVector{<:Integer})
#     # !iscompatible(lattice, tsym) && return false
#     orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
#     orthogonal_shape = tsym.orthogonal_shape
#     return iscompatible(orthogonal_momentum, orthogonal_shape, identity_translation)
# end


# function iscompatible(# lattice::Lattice,
#                       tsym::TranslationSymmetry,
#                       tsym_irrep_index::Integer,
#                       identity_translations::AbstractVector{<:AbstractVector{<:Integer}})
#     orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
#     orthogonal_shape = tsym.orthogonal_shape
#     return all(iscompatible(orthogonal_momentum, orthogonal_shape, t) for t in identity_translations)
# end


# function iscompatible(orthogonal_shape::AbstractVector{<:Integer},
#                       orthogonal_momentum::AbstractVector{<:Integer},
#                       identity_translation::TranslationOperation{<:Integer})
#     value = Rational{Int}(0)
#     for (n, i, j) in zip(orthogonal_shape, orthogonal_momentum, identity_translation.displacement)
#         value += i * j // n
#     end
#     return mod(value, 1) == 0
# end


# """
#     iscompatible(orthogonal_momentum, orthogonal_shape, identity_translations)

# Test whether the given momentum is compatible with *all* the identity translations.
# """
# function iscompatible(orthogonal_shape::AbstractVector{<:Integer},
#                       orthogonal_momentum::AbstractVector{<:Integer},
#                       identity_translations::AbstractVector{<:TranslationOperation{<:Integer}})
#     return all(iscompatible(orthogonal_shape, orthogonal_momentum, t) for t in identity_translations)
# end


## Between Something and PointSymmetry

# function iscompatible(hypercube::HypercubicLattice, matrix::AbstractMatrix{<:Integer})::Bool
#     _, elems = hypercube.wrap(matrix * hypercube.scale_matrix)
#     all(iszero(elems)) # all, since scale_matrix
# end

# function iscompatible(tsym::TranslationSymmetry, matrix::AbstractMatrix{<:Integer})::Bool
#     sm = tsym.hypercube.scale_matrix
#     smi = tsym.hypercube.inverse_scale_matrix
#     return all(mod(x,1) == 0 for x in smi * matrix * sm)
# end

