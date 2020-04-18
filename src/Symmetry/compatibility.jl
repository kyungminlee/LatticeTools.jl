export iscompatible


## 1. Hypercube and Operation

function iscompatible(hypercube::HypercubicLattice, op::TranslationOperation{<:Integer})::Bool
    if dimension(hypercube) != dimension(op)
        throw(DimensionMismatch("hypercube and operation have different dimension"))
    end
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


## 2. Hypercube and Symmetry

function iscompatible(hypercube::HypercubicLattice, tsym::TranslationSymmetry)::Bool
    R, r = hypercube.wrap(tsym.hypercube.scale_matrix)
    return iszero(r) && abs(ExactLinearAlgebra.determinant(R)) == 1
end

function iscompatible(hypercube::HypercubicLattice, psym::PointSymmetry)::Bool
    return all(iscompatible(hypercube, pop) for pop in elements(psym))
end


## 3. TranslationSymmetry and PointOperation/PointSymmetry

function iscompatible(tsym::TranslationSymmetry, pop::PointOperation{<:Integer})::Bool
    # sm = tsym.hypercube.scale_matrix
    # smi = tsym.hypercube.inverse_scale_matrix
    # return all(mod(x,1) == 0 for x in smi * pop.matrix * sm)
    return iscompatible(tsym.hypercube, pop)
end

function iscompatible(tsym::TranslationSymmetry, psym::PointSymmetry)::Bool
    return iscompatible(tsym.hypercube, psym)
end


## 4. Translation Irreps and PointOperation/PointSymmetry

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


## 5. Lattice and Operation

function iscompatible(lattice::Lattice, op::TranslationOperation{<:Integer})
    if dimension(lattice) != dimension(op)
        throw(DimensionMismatch("lattice and op have different dimension"))
    end
    return true
end

function iscompatible(lattice::Lattice, op::PointOperation{<:Integer})
    return iscompatible(lattice.hypercube, op) &&
           !isnothing(findorbitalmap(lattice.unitcell, op))
end


## 6. Lattice and Symmetry

"""
    iscompatible(lattice, translation_symmetry)

Test whether lattice and the symmetry are compatible.
For translation symmetry, this means that the hypercubic lattice for both are the same.
"""
function iscompatible(lattice::Lattice, tsym::TranslationSymmetry)
    # return lattice.hypercube == tsym.hypercube
    return iscompatible(lattice.hypercube, tsym)
end


function iscompatible(lattice::Lattice, psym::PointSymmetry)::Bool
    return iscompatible(lattice.hypercube, psym) &&
           all(!isnothing(findorbitalmap(lattice.unitcell, pop))
                   for pop in generator_elements(psym))
end
