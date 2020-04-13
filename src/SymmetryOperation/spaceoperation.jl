export SpaceOperation

export apply_operation
export canonize
export iscanonical

struct SpaceOperation{S<:Real}
    matrix::Matrix{S}
    displacement::Vector{S}

    function SpaceOperation(matrix::AbstractMatrix{S},
                            displacement::AbstractVector{S}) where {S<:Real}
        return new{S}(matrix, displacement)
    end

    function SpaceOperation(point::PointOperation{S},
                            translation::TranslationOperation{S}) where {S<:Real}
        return new{S}(point.matrix, translation.displacement)
    end
end


import Base.*
function (*)(lhs::ProductOperation, rhs::TranslationOperation)
    SpaceOperation(lhs.matrix, rhs.displacement)
end


function (*)(lhs::TranslationOperation, rhs::ProductOperation)
    rhs_inv = inv(rhs)
    SpaceOperation(rhs, TranslationOperation(rhs_inv.matrix * lhs.displacement))
end


import Base.inv
function inv(arg::SpaceOperation)
    matrix_inv = ExactLinearAlgebra.inverse(arg.matrix)
    return SpaceOperation(matrix_inv, -arg.matrix * arg.displacement)
end