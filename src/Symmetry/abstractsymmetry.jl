
abstract type AbstractSymmetry end

ConjugacyClassType = NamedTuple{(:name,), Tuple{String}}
RepresentationType = NamedTuple{(:name, :matrices), Tuple{String, Vector{Matrix{ComplexF64}}}}
