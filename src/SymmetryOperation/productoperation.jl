export DirectProductOperation
export isidentity
import LinearAlgebra.×
export directproduct

"""
    DirectProductOperation{Ops}

Direct product of symmetry operations.
math```
  g = (g_1, g_2, \\ldots)
```
"""
struct DirectProductOperation{Ops<:Tuple{Vararg{AbstractSymmetryOperation}}}<:AbstractSymmetryOperation
    operations::Ops
    function DirectProductOperation(ops::Vararg{AbstractSymmetryOperation})
        T = typeof(ops)
        return new{T}(ops)
    end
    function DirectProductOperation(ops::T) where {T<:Tuple{Vararg{AbstractSymmetryOperation}}}
        return new{T}(ops)
    end
end

×(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation) = DirectProductOperation(lhs, rhs)
×(lhs::DirectProductOperation, rhs::AbstractSymmetryOperation) = DirectProductOperation(lhs.operations..., rhs)
×(lhs::AbstractSymmetryOperation, rhs::DirectProductOperation) = DirectProductOperation(lhs, rhs.operations...)
×(lhs::DirectProductOperation, rhs::DirectProductOperation) = DirectProductOperation(lhs.operations..., rhs.operations...)

Base.hash(x::P, h::UInt) where {P<:DirectProductOperation} = Base.hash(P, Base.hash(x.operations, h))

Base.:(*)(lhs::P, rhs::P) where {P<:DirectProductOperation} = DirectProductOperation(lhs.operations .* rhs.operations)
Base.:(^)(obj::DirectProductOperation, n::Integer) = DirectProductOperation(obj.operations.^n)
Base.:(==)(lhs::P, rhs::P) where {P<:DirectProductOperation} = all(lhs.operations .== rhs.operations)

Base.inv(obj::DirectProductOperation) = DirectProductOperation(Base.inv.(obj.operations))

function Base.isapprox(lhs::P, rhs::P; atol::Real=0, rtol::Real=Base.rtoldefault(Float64)) where {P<:DirectProductOperation}
    return all(isapprox.(lhs.operations, rhs.operations; atol=atol, rtol=rtol))
end

isidentity(obj::DirectProductOperation) = all(isidentity, obj.operations)

function directproduct(::Type{E}, products::Function...) where {E<:DirectProductOperation}
    function product(lhs::E, rhs::E)
        return DirectProductOperation([p(l, r) for (p, l, r) in zip(products, lhs.operations, rhs.operations)]...)
    end
    return product
end


# export ProductOperation

# export apply_operation
# export canonize
# export iscanonical
# export domaintype
# export isidentity

# struct ProductOperation{S<:Real, F<:Tuple}
#     dimension::Int
#     factors::F
#     function ProductOperation(factors::AbstractSpaceSymmetryOperation...)
#         F = typeof(factors)
#         dim = 0
#         S = Bool
#         for f in factors
#             S = promote_type(domaintype(f), S)
#             df = dimension(f)
#             if df != 0
#                 if dim == 0
#                     dim = df
#                 else
#                     dim != df && throw(ArgumentError("All factors must have the same dimension"))
#                 end
#             else
#                 # this should not happen since `identity * something` is always `something`
#             end
#         end
#         new{S, F}(dim, factors)
#     end

#     function ProductOperation{S}(factors::AbstractSpaceSymmetryOperation...) where {S<:Real}
#         F = typeof(factors)
#         dim = 0
#         for f in factors
#             df = dimension(f)
#             if df != 0
#                 if dim == 0
#                     dim = df
#                 else
#                     dim != df && throw(ArgumentError("All factors must have the same dimension"))
#                 end
#             else
#                 # this should not happen since `identity * something` is always `something`
#             end
#         end
#         new{S, F}(dim, factors)
#     end
# end


# ## properties
# dimension(arg::ProductOperation) = arg.dimension
# domaintype(arg::ProductOperation{S,F}) where {S,F} = S
# combinable(arg::ProductOperation) = false
# isidentity(arg::ProductOperation) = isa(canonize(arg), IdentityOperation) # <- TODO: Think about whether this is a good idea


# ## operators
# import Base.==
# """
#     ==(lhs::ProductOperation, rhs::ProductOperation)

# Test for LITERAL equality.
# """
# function ==(lhs::ProductOperation, rhs::ProductOperation)
#     dimension(lhs) != dimension(rhs) && return false
#     length(lhs.factors) != length(rhs.factors) && return false
#     return all(l == r for (l, r) in zip(lhs.factors, rhs.factors))
# end

# import Base.*
# function (*)(lhs::AbstractSpaceSymmetryOperation, rhs::AbstractSpaceSymmetryOperation)
#     return ProductOperation(lhs, rhs)
# end

# function (*)(lhs::ProductOperation, rhs::AbstractSpaceSymmetryOperation)
#     return ProductOperation(lhs.factors..., rhs)
# end

# function (*)(lhs::AbstractSpaceSymmetryOperation, rhs::ProductOperation)
#     return ProductOperation(lhs, rhs.factors...)
# end

# function (*)(lhs::ProductOperation, rhs::ProductOperation)
#     return ProductOperation(lhs.factors..., rhs.factors...)
# end

# import Base.^
# function (^)(lhs::ProductOperation{S,F}, rhs::Integer)::ProductOperation where {S,F<:Tuple}
#     if rhs == 0
#         return ProductOperation()
#     elseif rhs > 0
#         return ProductOperation(vcat([collect(lhs.factors) for i in 1:rhs]...)...)
#     else
#         lhs_inv = inv(lhs)
#         return ProductOperation(vcat([collect(lhs_inv.factors) for i in 1:(-rhs)]...)...)
#     end
# end

# import Base.inv
# function inv(arg::ProductOperation)
#     return ProductOperation(reverse(inv.(arg.factors))...)
# end


# ## apply
# function apply_operation(symop::ProductOperation{S,F}, coord::AbstractVector{S}) where {S<:Real, F}
#     return foldr(apply_operation, symop.factors; init=coord)
# end

# function (symop::ProductOperation{S,F})(coord::AbstractVector{S}) where {S, F}
#     return foldr(apply_operation, symop.factors; init=coord)
# end


# ## canonical
# iscanonical(pop::ProductOperation) = false
# function iscanonical(
#             pop::ProductOperation{<:Real,
#                                   <:Tuple{<:PointOperation, <:TranslationOperation}})
#     return iscanonical(pop.factors[1]) && iscanonical(pop.factors[2])
# end

# function canonize(arg::ProductOperation)
#     _reorder(lhs::IdentityOperation, rhs::AbstractSpaceSymmetryOperation) = (rhs, lhs)
#     _reorder(lhs::AbstractSpaceSymmetryOperation, rhs::IdentityOperation) = (lhs, rhs)
#     _reorder(lhs::TranslationOperation, rhs::TranslationOperation) = (lhs, rhs)
#     _reorder(lhs::PointOperation, rhs::PointOperation) = (lhs, rhs)
#     _reorder(lhs::PointOperation, rhs::TranslationOperation) = (lhs, rhs)
#     function _reorder(lhs::TranslationOperation, rhs::PointOperation)
#         rhs_inv = inv(rhs)
#         (rhs, TranslationOperation(rhs_inv.matrix * lhs.displacement))
#     end

#     factors = collect(arg.factors)
#     # first, bubble sort
#     # order: Point, Translation, Identity (implemented in _reorder)
#     n = length(factors)
#     for m in n:-1:2, i in 1:(m-1)
#         factors[i], factors[i+1] = _reorder(factors[i], factors[i+1])
#     end

#     new_factors = AbstractSpaceSymmetryOperation[]
#     op = IdentityOperation()
#     for f in factors
#         if combinable(op, f)
#             op *= f
#         else
#             op = canonize(op)
#             op != IdentityOperation() && push!(new_factors, op)
#             op = f
#         end
#     end
#     op = canonize(op)
#     op != IdentityOperation() && push!(new_factors, op)

#     if length(new_factors) == 0
#         return IdentityOperation()
#     elseif length(new_factors) == 1
#         return first(new_factors)
#     else
#         return ProductOperation(new_factors...)
#     end
# end
