export iscompatible


## 1. OrthoCube and Operation

"""
    iscompatible(orthocube, operator::TranslationOperation{<:Integer})
"""
function iscompatible(orthocube::OrthoCube, op::TranslationOperation{<:Integer})::Bool
    if dimension(orthocube) != dimension(op)
        throw(DimensionMismatch("hypercube and operation have different dimension"))
    end
    return true
end

"""
    iscompatible(orthocube, operator::PointOperation{<:Integer})
"""
function iscompatible(orthocube::OrthoCube, op::PointOperation{<:Integer})::Bool
    R, r = orthocube.wrap(op.matrix * orthocube.shape_matrix)
    return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
end


## 2. OrthoCube and Symmetry

"""
    iscompatible(orthocube, symmetry::TranslationSymmetry)
"""
function iscompatible(orthocube::OrthoCube, tsym::TranslationSymmetry)::Bool
    R, r = orthocube.wrap(tsym.orthocube.shape_matrix)
    return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
end

"""
    iscompatible(orthocube, symmetry::PointSymmetry)
"""
function iscompatible(orthocube::OrthoCube, psym::PointSymmetry)::Bool
    return all(iscompatible(orthocube, pop) for pop in elements(psym))
end


## 3. TranslationSymmetry and PointOperation/PointSymmetry

"""
    iscompatible(tsym::TranslationSymmetry, pop::PointOperation{<:Integer})
"""
function iscompatible(tsym::TranslationSymmetry, pop::PointOperation{<:Integer})::Bool
    return iscompatible(tsym.orthocube, pop)
end

"""
    iscompatible(tsym::TranslationSymmetry, psym::PointSymmetry)
"""
function iscompatible(tsym::TranslationSymmetry, psym::PointSymmetry)::Bool
    return iscompatible(tsym.orthocube, psym)
end

"""
    iscompatible(psym::PointSymmetry, tsym::TranslationSymmetry)
"""
function iscompatible(psym::PointSymmetry, tsym::TranslationSymmetry)::Bool
    return iscompatible(tsym, psym)
end


## 4. Translation Irreps and PointOperation/PointSymmetry

"""
    iscompatible(tsym, tsym_irrep_index, point_operation)
"""
function iscompatible(
    tsym::TranslationSymmetry,
    tsym_irrep_index::Integer,
    pop::PointOperation{<:Integer},
)::Bool
    !iscompatible(tsym, pop) && return false
    reciprocal_matrix = ExactLinearAlgebra.inverse(transpose(pop.matrix))
    k1 = tsym.fractional_momenta[tsym_irrep_index]
    k2 = (x -> mod(x, 1)).(reciprocal_matrix * k1)
    return k1 == k2
end

"""
    iscompatible(tsym, tsym_irrep_index, psym)
"""
function iscompatible(
    tsym::TranslationSymmetry,
    tsym_irrep_index::Integer,
    psym::PointSymmetry,
)::Bool
    return all(
        iscompatible(tsym, tsym_irrep_index, pop)
        for pop in generator_elements(psym)
    )
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
    return iscompatible(lattice.orthocube, op) &&
           !isnothing(findsitemap(lattice.unitcell, op))
end


## 6. Lattice and Symmetry

"""
    iscompatible(lattice, translation_symmetry)

Test whether lattice and the symmetry are compatible.
For translation symmetry, this means that the hypercubic lattice for both are the same.
"""
function iscompatible(lattice::Lattice, tsym::TranslationSymmetry)
    # return lattice.hypercube == tsym.hypercube
    return iscompatible(lattice.orthocube, tsym)
end

"""
    iscompatible(lattice, point_symmetry)
"""
function iscompatible(lattice::Lattice, psym::PointSymmetry)::Bool
    return (
        iscompatible(lattice.orthocube, psym) &&
        all(pop -> !isnothing(findsitemap(lattice.unitcell, pop)), generator_elements(psym))
    )
end



# function iscompatible(hypercube::HypercubicLattice, op::TranslationOperation{<:Integer})::Bool
#     if dimension(hypercube) != dimension(op)
#         throw(DimensionMismatch("hypercube and operation have different dimension"))
#     end
#     return true
# end

# """
#     iscompatible(hypercube, point_operation)

#     (R ⋅ L) mod L = 0   and   (R ⋅ L) / L = 1
# """
# function iscompatible(hypercube::HypercubicLattice, op::PointOperation{<:Integer})::Bool
#     R, r = hypercube.wrap(op.matrix * hypercube.scale_matrix)
#     return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
# end


# function iscompatible(hypercube::HypercubicLattice, tsym::TranslationSymmetry)::Bool
#     R, r = hypercube.wrap(tsym.hypercube.scale_matrix)
#     return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
# end

# function iscompatible(hypercube::HypercubicLattice, psym::PointSymmetry)::Bool
#     return all(iscompatible(hypercube, pop) for pop in elements(psym))
# end
