export iscompatible

import LinearAlgebraX

## 1. Hypercube and Operation

"""
    iscompatible(hypercube, op::TranslationOperation{<:Integer})

Check whether the translation operation is compatible with the hypercube,
i.e. leaves the superlattice invariant. If the dimensions do not match, throw exception.
Otherwise always return true.
    
"""
function iscompatible(hypercube::Hypercube, op::TranslationOperation{<:Integer})::Bool
    if dimension(hypercube) != dimension(op)
        throw(DimensionMismatch("hypercube and operation have different dimension"))
    end
    return true
end

"""
    iscompatible(hypercube, op::PointOperation{<:Integer})

Check whether the point operation is compatible with the hypercube,
i.e. leaves the superlattice invariant.
"""
function iscompatible(hypercube::Hypercube, op::PointOperation{<:Integer})::Bool
    R, r = hypercube.wrap(op.matrix * hypercube.shape_matrix)
    return iszero(r) && isone(abs(LinearAlgebraX.detx(R)))
end


## 2. Hypercube and Symmetry

"""
    iscompatible(hypercube, symmetry::FiniteTranslationSymmetry)

Check whether the translation symmetry is compatible with the hypercube,
i.e. has the same superlattice.
"""
function iscompatible(hypercube::Hypercube, tsym::FiniteTranslationSymmetry)::Bool
    R, r = hypercube.wrap(tsym.hypercube.shape_matrix)
    return iszero(r) && isone(abs(LinearAlgebraX.detx(R)))
end

"""
    iscompatible(hypercube, symmetry::PointSymmetry)

Check whether the point symmetry is compatible with the hypercube,
i.e. all of its elements leaves the superlattice invariant.
"""
function iscompatible(hypercube::Hypercube, psym::PointSymmetry)::Bool
    return all(iscompatible(hypercube, pop) for pop in elements(psym))
end


## 3. FiniteTranslationSymmetry and PointOperation/PointSymmetry

"""
    iscompatible(tsym::FiniteTranslationSymmetry, pop::PointOperation{<:Integer})

Check whether the point operation is compatible with the translation symmetry,
i.e. leaves the superlattice invariant.
"""
function iscompatible(tsym::FiniteTranslationSymmetry, pop::PointOperation{<:Integer})::Bool
    return iscompatible(tsym.hypercube, pop)
end

"""
    iscompatible(tsym::FiniteTranslationSymmetry, psym::PointSymmetry)

Check whether the point symmetry is compatible with the translation symmetry,
i.e. all of its elements leaves the superlattice invariant.
"""
function iscompatible(tsym::FiniteTranslationSymmetry, psym::PointSymmetry)::Bool
    return iscompatible(tsym.hypercube, psym)
end

"""
    iscompatible(psym::PointSymmetry, tsym::FiniteTranslationSymmetry)

Check whether the point symmetry is compatible with the translation symmetry,
i.e. all of its elements leaves the superlattice invariant.
"""
function iscompatible(psym::PointSymmetry, tsym::FiniteTranslationSymmetry)::Bool
    return iscompatible(tsym, psym)
end


## 4. Translation Irreps and PointOperation/PointSymmetry

"""
    iscompatible(tsym, tsym_irrep_index, point_operation)

Check whether the point operation is compatible with the given irrep of translation (a.k.a. momentum).
In specific, check whether
(1) the point operation is compatible with the translation symmetry, and
(2) the momentum is invariant under the point operation.
"""
function iscompatible(
    tsym::FiniteTranslationSymmetry,
    tsym_irrep_index::Integer,
    pop::PointOperation{<:Integer},
)::Bool
    !iscompatible(tsym, pop) && return false
    reciprocal_matrix = LinearAlgebraX.invx(transpose(pop.matrix))
    k1 = tsym.fractional_momenta[tsym_irrep_index]
    k2 = (x -> mod(x, 1)).(reciprocal_matrix * k1)
    return k1 == k2
end


"""
    iscompatible(tsym, tsym_irrep_index, psym)

Check whether the point symmetry is compatible with the given irrep of translation (a.k.a. momentum).
In specific, check whether
(1) the point operation is compatible with the translation symmetry, and
(2) the momentum is invariant under all point operations.
"""
function iscompatible(
    tsym::FiniteTranslationSymmetry,
    tsym_irrep_index::Integer,
    psym::PointSymmetry,
)::Bool
    return all(iscompatible(tsym, tsym_irrep_index, pop) for pop in generator_elements(psym))
end


## 5. Lattice and Operation

"""
    iscompatible(lattice, translation_operation)
"""
function iscompatible(lattice::Lattice, op::TranslationOperation{<:Integer})
    if dimension(lattice) != dimension(op)
        throw(DimensionMismatch("lattice and op have different dimension"))
    end
    return true
end


"""
    iscompatible(lattice, point_operation)
"""
function iscompatible(lattice::Lattice, op::PointOperation{<:Integer})
    return iscompatible(lattice.hypercube, op) &&
           !isnothing(findsitemap(lattice.unitcell, op))
end


## 6. Lattice and Symmetry

"""
    iscompatible(lattice, translation_symmetry)

Test whether lattice and the symmetry are compatible.
For translation symmetry, this means that the hypercubic lattice for both are the same.
"""
function iscompatible(lattice::Lattice, tsym::FiniteTranslationSymmetry)
    return iscompatible(lattice.hypercube, tsym)
end

"""
    iscompatible(lattice, point_symmetry)
"""
function iscompatible(lattice::Lattice, psym::PointSymmetry)::Bool
    return (
        iscompatible(lattice.hypercube, psym) &&
        all(pop -> !isnothing(findsitemap(lattice.unitcell, pop)), generator_elements(psym))
    )
end
