export AbstractSymmetryOperation
export combinable

abstract type AbstractSymmetryOperation end

combinable(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation) = false
