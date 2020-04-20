export AbstractSymmetry
# export AbstractUnitSymmetry
# export AbstractProductSymmetry

abstract type AbstractSymmetry{OperationType<:AbstractSymmetryOperation} end


# abstract type AbstractUnitSymmetry{OperationType<:AbstractSymmetryOperation} <: AbstractSymmetry{OperationType} end
# abstract type AbstractProductSymmetry{OperationType<:AbstractSymmetryOperation} <: AbstractSymmetry{OperationType} end
