export AbstractSymmetry
# export AbstractUnitSymmetry

"""
    AbstractSymmetry{OperationType<:AbstractSpaceSymmetryOperation}

Abstract symmetry whose elements are of type `OperationType`.
"""
abstract type AbstractSymmetry{OperationType<:AbstractSpaceSymmetryOperation} end
# abstract type AbstractUnitSymmetry{OperationType<:AbstractSpaceSymmetryOperation} <: AbstractSymmetry{OperationType} end
