export AbstractSymmetry

abstract type AbstractSymmetry end

ConjugacyClassType = NamedTuple{(:name, :elements), Tuple{String, Vector{Int}}}
IrrepType = NamedTuple{(:name, :matrices), Tuple{String, Vector{Matrix{ComplexF64}}}}
