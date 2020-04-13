export ProductOperation

export apply_operation
export canonize
export iscanonical


struct ProductOperation{F<:Tuple}
    factors::F
    function ProductOperation(factors::AbstractSymmetryOperation...)
        F = typeof(factors)
        new{F}(factors)
    end
end

import Base.==

"""
    ==(lhs::ProductOperation, rhs::ProductOperation)

Test for LITERAL equality.
"""
function ==(lhs::ProductOperation, rhs::ProductOperation)
    if length(lhs.factors) != length(rhs.factors)
        return false
    end
    return all(l == r for (l, r) in zip(lhs.factors, rhs.factors))
end


import Base.*
function (*)(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation)
    return ProductOperation(lhs, rhs)
end

function (*)(lhs::ProductOperation, rhs::AbstractSymmetryOperation)
    ProductOperation(lhs.factors..., rhs)
end

function (*)(lhs::AbstractSymmetryOperation, rhs::ProductOperation)
    ProductOperation(lhs, rhs.factors...)
end

function (*)(lhs::ProductOperation, rhs::ProductOperation)
    return ProductOperation(lhs.factors..., rhs.factors...)
end

import Base.inv
function inv(arg::ProductOperation)
    return ProductOperation(reverse(inv.(arg.factors))...)
end

function (^)(lhs::ProductOperation{T}, rhs::Integer)::ProductOperation where {T<:Tuple}
    if rhs == 0
        return ProductOperation()
    elseif rhs > 0
        return ProductOperation(vcat([collect(lhs.factors) for i in 1:rhs]...)...)
    else
        lhs_inv = inv(lhs)
        return ProductOperation(vcat([collect(lhs_inv.factors) for i in 1:(-rhs)]...)...)
    end
end



function apply_operation(symop::ProductOperation, coord::AbstractVector{<:Real})
    return foldr(apply_operation, symop.factors; init=coord)
end


function (symop::ProductOperation)(coord::AbstractVector{<:Real})
    return foldr(apply_operation, symop.factors; init=coord)
end


function canonize(arg::ProductOperation)
    _reorder(lhs::IdentityOperation, rhs::AbstractSymmetryOperation) = (rhs, lhs)
    _reorder(lhs::AbstractSymmetryOperation, rhs::IdentityOperation) = (lhs, rhs)
    _reorder(lhs::TranslationOperation, rhs::TranslationOperation) = (lhs, rhs)
    _reorder(lhs::PointOperation, rhs::PointOperation) = (lhs, rhs)
    _reorder(lhs::PointOperation, rhs::TranslationOperation) = (lhs, rhs)
    function _reorder(lhs::TranslationOperation, rhs::PointOperation)
        rhs_inv = inv(rhs)
        (rhs, TranslationOperation(rhs_inv.matrix * lhs.displacement))
    end

    factors = collect(arg.factors)
    # first, bubble sort
    # order: Point, Translation, Identity (implemented in _reorder)
    n = length(factors)
    for m in n:-1:2, i in 1:(m-1)
        factors[i], factors[i+1] = _reorder(factors[i], factors[i+1])
    end

    new_factors = AbstractSymmetryOperation[]
    op = IdentityOperation()
    for f in factors
        if combinable(op, f)
            op *= f
        else
            op = canonize(op)
            if op != IdentityOperation()
                push!(new_factors, op)
            end
            op = f
        end
    end
    op = canonize(op)
    if op != IdentityOperation()
        push!(new_factors, op)
    end
    
    if length(new_factors) == 0
        return IdentityOperation()
    elseif length(new_factors) == 1
        return first(new_factors)
    else
        return ProductOperation(new_factors...)
    end
end

iscanonical(pop::ProductOperation) = false
function iscanonical(pop::ProductOperation{<:Tuple{<:PointOperation, <:TranslationOperation}})
    return iscanonical(pop.factors[1]) && iscanonical(pop.factors[2])
end