# Local degrees of freedom

module LDOF

export LocalDof
export ProductDof


struct LocalDof
    name::String
    names::Vector{String}
end

import Base.length
length(ldof::LocalDof) = length(ldof.names)

import Base.iterate
iterate(ldof::LocalDof) = iterate(ldof.names)


struct ProductDof
    localdof::Vector{LocalDof}
end

import Base.length
length(pdof::ProductDof) = prod(length(ldof) for ldof in pdof.localdof)

import Base.iterate
iterate(pdof::ProductDof) = iterate(Iterators.product(pdof.localdof...))


# import Base.×
export ×
(×)(lhs::LocalDof, rhs::LocalDof) = ProductDof([lhs, rhs])
(×)(lhs::ProductDof, rhs::LocalDof) = ProductDof([lhs..., rhs])
(×)(lhs::LocalDof, rhs::ProductDof) = ProductDof([lhs, rhs...])
(×)(lhs::ProductDof, rhs::ProductDof) = ProductDof([lhs..., rhs...])

end

# spin = LocalDof("spin", ["up", "dn"])
# orbital = LocalDof("orbital", ["px", "py", "pz"])

#=
module LDOF

export LocalDof
export ProductDof

struct LocalDof{names} end
struct ProductDof{local_degrees<:Tuple} end

# import Base.×
export ×
(×)(lhs::Type{LocalDof{S1}}, rhs::Type{LocalDof{S2}}) where {S1, S2} = ProductDof{Tuple{lhs, rhs}}
(×)(lhs::Type{ProductDof{T1}}, rhs::Type{LocalDof{S2}}) where {T1, S2} = ProductDof{Tuple{T..., rhs}}
(×)(lhs::Type{LocalDof{S1}}, rhs::Type{ProductDof{T2}}) where {S1, T2} = ProductDof{Tuple{lhs, T...}}
(×)(lhs::Type{ProductDof{T1}}, rhs::Type{ProductDof{T2}}) where {T1, T2} = ProductDof{Tuple{T1..., T2...}}

end

using .LDOF

spin = LocalDof{(:up, :dn)}
orbital = LocalDof{(:px, :py)}
=#