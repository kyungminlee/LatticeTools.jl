export ProductOperation

export inverse
export apply_symmetry
export canonize
export iscanonical


struct ProductOperation{F<:Tuple}
    factors::F
    function ProductOperation(factors::AbstractSymmetryOperation...)
        F = typeof(factors)
        new{F}(factors)
    end
end


import Base.*
function (*)(lhs::AbstractSymmetryOperation, rhs::AbstractSymmetryOperation)
    return ProductOperation(lhs, rhs)
end

function (*)(lhs::ProductOperation, rhs::AbstractSymmetryOperation)
    isempty(lhs.factors) && return rhs
    if isa(last(lhs.factors), typeof(rhs))
        ProductOperation(lhs.factors[1:end-1]..., lhs.factors[end]*rhs)
    else
        ProductOperation(lhs.factors..., rhs)
    end
end

function (*)(lhs::AbstractSymmetryOperation, rhs::ProductOperation)
    isempty(rhs.factors) && return lhs
    if isa(first(rhs.factors), typeof(lhs))
        ProductOperation(lhs*rhs.factors[1], rhs.factors[2:end]...)
    else
        ProductOperation(lhs, rhs.factors...)
    end
end

function (*)(lhs::ProductOperation, rhs::ProductOperation)
    return ProductOperation(lhs.factors..., rhs.factors...)
end


function inverse(arg::ProductOperation)
    return ProductOperation(reverse(inverse.(arg.factors))...)
end


function apply_symmetry(symop::ProductOperation, coord::AbstractVector{<:Real})
    return foldr(apply_symmetry, symop.factors; init=coord)
end


function canonize(arg::ProductOperation)
    _reorder(lhs::IdentityOperation, rhs::AbstractSymmetryOperation) = (rhs, lhs)
    _reorder(lhs::AbstractSymmetryOperation, rhs::IdentityOperation) = (lhs, rhs)
    _reorder(lhs::TranslationOperation, rhs::TranslationOperation) = (lhs, rhs)
    _reorder(lhs::PointOperation, rhs::PointOperation) = (lhs, rhs)
    _reorder(lhs::PointOperation, rhs::TranslationOperation) = (lhs, rhs)
    function _reorder(lhs::TranslationOperation, rhs::PointOperation)
        rhs_inv = inverse(rhs)
        (rhs, TranslationOperation(rhs_inv.matrix * lhs.displacement))
    end

    factors = collect(arg.factors)
    # first, bubble sort
    # order: Identity, Point, Translation (implemented in _reorder)
    n = length(factors)
    for m in n:-1:2
        for i in 1:(m-1)
            factors[i], factors[i+1] = _reorder(factors[i], factors[i+1])
        end
    end

    out = prod(factors)
    if !isa(out, ProductOperation)
        return out
    end

    factors = [canonize(f) for f in out.factors]
    
    filter!(x->!isa(x, IdentityOperation), factors)
    pop = filter(x->isa(x, PointOperation), factors)
    top = filter(x->isa(x, TranslationOperation), factors)

    if isempty(pop)
        if isempty(top)
            return IdentityOperation()
        else
            @assert length(top) == 1
            return first(top)
        end
    else
        @assert length(pop) == 1
        if isempty(top)
            return first(pop)
        else
            @assert length(top) == 1
            return ProductOperation(first(pop), first(top))
        end
    end
end

iscanonical(pop::ProductOperation) = false
function iscanonical(pop::ProductOperation{<:Tuple{<:PointOperation, <:TranslationOperation}})
    return iscanonical(pop.factors[1]) && iscanonical(pop.factors[2])
end