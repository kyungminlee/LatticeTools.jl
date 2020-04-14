export SpaceOperation
export apply_operation
export canonize
export iscanonical

using LinearAlgebra

struct SpaceOperation{S<:Real}
    matrix::Matrix{S}
    displacement::Vector{S}

    function SpaceOperation(matrix::AbstractMatrix{S},
                            displacement::AbstractVector{S}) where {S<:Real}
        dim = size(matrix, 1)
        if size(matrix, 2) != dim
            throw(DimensionMismatch("matrix is not square: dimensions are $(size(matrix))"))
        elseif length(displacement) != dim
            throw(DimensionMismatch("dimensions of matrix and displacement do not match: matrix dimensions are $(size(matrix)), and vector dimension is $(length(displacement))"))
        end
        return new{S}(matrix, displacement)
    end

    function SpaceOperation(point::PointOperation{S},
                            translation::TranslationOperation{S}) where {S<:Real}
        return new{S}(point.matrix, translation.displacement)
    end

    function SpaceOperation(arg::ProductOperation{S, F}) where {S, F}
        dim = dimension(arg)
        if dim == 0
            throw(ArgumentError("Cannot create space operation out of an identity (dimension undetermined)"))
        end
        
        canon = canonize(arg)
        if isa(canon, IdentityOperation)
            return SpaceOperation(Matrix(I, dim, dim), zeros(S, dim))
        elseif isa(canon, TranslationOperation)
            return SpaceOperation(Matrix(I, dim, dim), canon.displacement)
        elseif isa(canon, PointOperation)
            return SpaceOperation(canon.matrix, zeros(S, dim))
        else
            @assert length(canon.factors) == 2
            return SpaceOperation(canon.factors[1], canon.factors[2])
        end
    end
end


import Base.inv
function inv(arg::SpaceOperation)
    matrix_inv = ExactLinearAlgebra.inverse(arg.matrix)
    return SpaceOperation(matrix_inv, -arg.matrix * arg.displacement)
end