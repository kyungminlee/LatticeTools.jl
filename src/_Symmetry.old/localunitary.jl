
export LocalSpecialUnitaryOperation

import LinearAlgebra

struct LocalSpecialUnitaryOperation{M<:AbstractMatrix} <: AbstractSymmetryOperation
    operators ::Dict{Int, M}
    function LocalSpecialUnitaryOperation(p::AbstractDict{Int, M})
        for (k, v) in p
            if size(v,1) != size(v, 2)
                throw(ArgumentError("element at index $p is not a square matrix")
            elseif !isapprox(LinearAlgebra.det(v), 1)
                throw(ArgumentError("element at index $p is not a special unitary matrix (has determinant $(LinearAlgebra.det(v)))")
            elseif !isapprox(v * adjoint(v), LinearAlgebra.I)
                throw(ArgumentError("element at index $p is not a special unitary matrix")
            end
        end
        return new(p)
    end

    function LocalSpecialUnitaryOperation(p::Pair{<:Integer, <:AbstractMatrix}...)
        return LocalSpecialUnitaryOperation(Dict(p))
    end

    function LocalSpecialUnitaryOperation(kv)
        return LocalSpecialUnitaryOperation(Dict(kv))
    end

end

import Base.*

function Base.*(lhs ::LocalSpecialUnitaryOperation{M1}, rhs ::LocalSpecialUnitaryOperation{M2}) where {M1, M2}
    keys_lhs = setdiff(keys(lhs), keys(rhs))
    keys_rhs = setdiff(keys(rhs), keys(lhs))
    keys_intersect = intersect(keys(lhs), keys(rhs))
    M3 = promote_type(M1, M2)
    out = Dict{Int, M3}()
    for k in keys_lhs
        out[k] = lhs.operators[k]
    end
    for k in keys_rhs
        out[k] = rhs.operator[k]
    end
    for k in keys_intersect
        out[k] = lhs.operators[k] * rhs.operators[k]
    end
    return LocalSpecialUnitaryOperation(out)
end

export ProductSymmetryOperation
struct ProductSymmetryOperation <: AbstractSymmetryOperation
    operators ::Vector{AbstractSymmetryOperation}
end

function Base.*(lhs ::LocalSpecialUnitaryOperation{M}, rhs ::Permutation)
    return ProductSymmetryOperation([lhs, rhs])
end

function Base.*(lhs ::Permutation, rhs ::LocalSpecialUnitaryOperation{M})
    new_lhs = LocalSpecialUnitaryOperation(Dict{Int, M}(lhs.permutation[i] -> u for (i, u) in rhs.operators))
    new_rhs = lhs
    return ProductSymmetryOperation([new_lhs, new_rhs])
end

function simplify(arg::Permutation; tol::Real=0)
    return arg
end

function simplify(arg::LocalSpecialUnitaryOperation{M}; tol=Base.rtoldefault(Float64))
    out = Dict{Int, M}()
    for (i, u) in arg.operators
        if !isappox(u, LinearAlgebra.I; atol=tol, rtol=tol)
            out[i] = u
        end
    end
    return LocalSpecialUnitaryOperation(out)
end

function simplify(arg::ProductSymmetryOperation; tol=Base.rtoldefault(Float64)) where M
    new_operators = [simplify(op; tol=tol) for op in arg.operators]
end

function simplify(arg::ProductSymmetryOperation{M1, M2}; tol=Base.rtoldefault(Float64)) where {M1, M2}

end
