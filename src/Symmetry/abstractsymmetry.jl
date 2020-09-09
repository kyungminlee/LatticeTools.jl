export AbstractSymmetry

"""
    AbstractSymmetry{OperationType<:AbstractSpaceSymmetryOperation}

Abstract symmetry whose elements are of type `OperationType`.
"""
abstract type AbstractSymmetry{OperationType<:AbstractSpaceSymmetryOperation} end
