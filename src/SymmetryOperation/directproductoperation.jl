export DirectProductOperation
export isidentity
export directproduct
export ×ˢ

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

Base.hash(x::P, h::UInt) where {P<:DirectProductOperation} = Base.hash(P, Base.hash(x.operations, h))

Base.:(*)(lhs::P, rhs::P) where {P<:DirectProductOperation} = DirectProductOperation(lhs.operations .* rhs.operations)
Base.:(^)(obj::DirectProductOperation, n::Integer) = DirectProductOperation(obj.operations.^n)
Base.:(==)(lhs::P, rhs::P) where {P<:DirectProductOperation} = all(lhs.operations .== rhs.operations)

Base.inv(obj::DirectProductOperation) = DirectProductOperation(Base.inv.(obj.operations))

function Base.isapprox(lhs::P, rhs::P; atol::Real=0, rtol::Real=Base.rtoldefault(Float64)) where {P<:DirectProductOperation}
    return all(isapprox.(lhs.operations, rhs.operations; atol=atol, rtol=rtol))
end

isidentity(obj::DirectProductOperation) = all(isidentity, obj.operations)


×ˢ(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation) = DirectProductOperation(lhs, rhs)
×ˢ(lhs::DirectProductOperation, rhs::AbstractSymmetryOperation) = DirectProductOperation(lhs.operations..., rhs)
×ˢ(lhs::AbstractSymmetryOperation, rhs::DirectProductOperation) = DirectProductOperation(lhs, rhs.operations...)
×ˢ(lhs::DirectProductOperation, rhs::DirectProductOperation) = DirectProductOperation(lhs.operations..., rhs.operations...)
function ×ˢ(lhs::AbstractVector{<:AbstractSymmetryOperation}, rhs::AbstractVector{<:AbstractSymmetryOperation})
    return [×ˢ(l, r) for l in lhs, r in rhs]
end

directproduct(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation) = DirectProductOperation(lhs, rhs)
directproduct(lhs::DirectProductOperation, rhs::AbstractSymmetryOperation) = DirectProductOperation(lhs.operations..., rhs)
directproduct(lhs::AbstractSymmetryOperation, rhs::DirectProductOperation) = DirectProductOperation(lhs, rhs.operations...)
directproduct(lhs::DirectProductOperation, rhs::DirectProductOperation) = DirectProductOperation(lhs.operations..., rhs.operations...)
function directproduct(lhs::AbstractVector{<:AbstractSymmetryOperation}, rhs::AbstractVector{<:AbstractSymmetryOperation})
    return [directproduct(l, r) for l in lhs, r in rhs]
end

function directproduct(::Type{E}, products::Function...) where {E<:DirectProductOperation}
    function product(lhs::E, rhs::E)
        return DirectProductOperation([p(l, r) for (p, l, r) in zip(products, lhs.operations, rhs.operations)]...)
    end
    return product
end
