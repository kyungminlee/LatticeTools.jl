export SpaceOperation
export apply_operation
export dimension
export isidentity, istranslation, ispoint

import LinearAlgebra
import LinearAlgebraX

"""
    SpaceOperation{Tp<:Real, Tt<:Real}

Represent a spatial symmetry operation of the following form:
```math
S_{[M, \\mathbf{R}]}: \\mathbf{r} \\mapsto M \\cdot ( \\mathbf{r} + \\mathbf{R} )
```

# Fields
* `matrix`: Rotation matrix
* `displacement`: Displacement vector
"""
struct SpaceOperation{Tp<:Real, Tt<:Real} <:AbstractSpaceSymmetryOperation{Tt}
    matrix::Matrix{Tp}
    displacement::Vector{Tt}

    function SpaceOperation{Tp, Tt}(dim::Integer) where {Tp<:Real, Tt<:Real}
        return new{Tp, Tt}(Matrix{Tp}(LinearAlgebra.I, dim, dim), zeros(Tt, dim))
    end

    function SpaceOperation(::Type{Tp}, ::Type{Tt}, dim::Integer) where {Tp<:Real, Tt<:Real}
        return new{Tp, Tt}(Matrix{Tp}(LinearAlgebra.I, dim, dim), zeros(Tt, dim))
    end

    function SpaceOperation{Tp, Tt}(
        matrix::AbstractMatrix,
        displacement::AbstractVector,
    ) where {Tp<:Real, Tt<:Real}
        dim = size(matrix, 1)
        if size(matrix, 2) != dim
            throw(DimensionMismatch("matrix is not square: dimensions are $(size(matrix))"))
        elseif length(displacement) != dim
            throw(DimensionMismatch(
                "dimensions of matrix and displacement do not match:"*
                " matrix dimensions are $(size(matrix)),"*
                " and vector dimension is $(length(displacement))"
            ))
        end
        return new{Tp, Tt}(matrix, displacement)
    end

    @doc """
        SpaceOperation(matrix::AbstractMatrix{Tp}, displacement::AbstractVector{Tt}) where {Tp<:Real, Tt<:Real}

    Construct a space operation by `matrix` and `displacement`.
    """
    function SpaceOperation(
        matrix::AbstractMatrix{Tp},
        displacement::AbstractVector{Tt},
    ) where {Tp<:Real, Tt<:Real}
        return SpaceOperation{Tp, Tt}(matrix, displacement)
    end


    @doc """
        SpaceOperation([point], [translation])

    Construct a space operation by `point` and `translation`.

    # Arguments
    * `point::PointOperation{Tp}`
    * `translation::TranslationOperation{Tt}`
    """
    function SpaceOperation(
        point::PointOperation{Tp},
        translation::TranslationOperation{Tt},
    ) where {Tp<:Real, Tt<:Real}
        return SpaceOperation{Tp, Tt}(point.matrix, translation.displacement)
    end

    function SpaceOperation(point::PointOperation{Tp}) where {Tp}
        dim = dimension(point)
        return SpaceOperation{Tp, Tp}(point.matrix, zeros(Tp, dim))
    end

    function SpaceOperation(translation::TranslationOperation{Tt}) where {Tt}
        dim = dimension(translation)
        return SpaceOperation{Int, Tt}(
            Matrix{Int}(LinearAlgebra.I, dim, dim),
            translation.displacement,
        )
    end
end

## Conversion

function Base.convert(::Type{SpaceOperation{Tp1, Tt1}}, arg::SpaceOperation{Tp2, Tt2}) where {Tp1, Tt1, Tp2, Tt2}
    return SpaceOperation(convert(Matrix{Tp1}, arg.matrix), convert(Vector{Tt1}, arg.displacement))
end

function Base.convert(::Type{SpaceOperation{Tp1, Tt1}}, op::IdentityOperation{T}) where {Tp1, Tt1, T}
    dim = dimension(op)
    mat = Matrix{Tp1}(LinearAlgebra.I, (dim, dim))
    dis = zeros(Tt1, size(mat, 1))
    return SpaceOperation(mat, dis)
end

function Base.convert(::Type{SpaceOperation{Tp1, Tt1}}, op::PointOperation{Tp2}) where {Tp1, Tt1, Tp2}
    dim = dimension(op)
    mat = convert(Matrix{Tp1}, op.matrix)
    dis = zeros(Tt1, dim)
    return SpaceOperation(mat, dis)
end

function Base.convert(::Type{SpaceOperation{Tp1, Tt1}}, op::TranslationOperation{Tt2}) where {Tp1, Tt1, Tt2}
    dim = dimension(op)
    mat = Matrix{Tp1}(LinearAlgebra.I, (dim, dim))
    dis = convert(Vector{Tt1}, op.displacement)
    return SpaceOperation(mat, dis)
end

function Base.convert(::Type{SpaceOperation{Tp1, Tt1}}, arg::AbstractVector{Tt2}) where {Tp1, Tt1, Tt2}
    dim = length(arg)
    mat = Matrix{Tp1}(LinearAlgebra.I, (dim, dim))
    dis = convert(Vector{Tt1}, arg)
    return SpaceOperation(mat, dis)
end

function Base.convert(::Type{SpaceOperation{Tp1, Tt1}}, arg::AbstractMatrix{Tp2}) where {Tp1, Tt1, Tp2}
    dim = size(arg, 1)
    mat = convert(Matrix{Tp1}, arg)
    dis = zeros(Tt1, dim)
    return SpaceOperation(mat, dis)
end


function Base.promote_rule(::Type{SpaceOperation{Tp1, Tt1}}, ::Type{SpaceOperation{Tp2, Tt2}}) where {Tp1, Tt1, Tp2, Tt2}
    Tp = promote_type(Tp1, Tp2)
    Tt = promote_type(Tt1, Tt2)
    return SpaceOperation{Tp, Tt}
end

function Base.promote_rule(::Type{SpaceOperation{Tp, Tt}}, ::Type{IdentityOperation{T}}) where {Tp, Tt, T}
    return SpaceOperation{Tp, Tt}
end

function Base.promote_rule(::Type{SpaceOperation{Tp1, Tt1}}, ::Type{TranslationOperation{Tt2}}) where {Tp1, Tt1, Tt2}
    Tt = promote_type(Tt1, Tt2)
    return SpaceOperation{Tp1, Tt}
end

function Base.promote_rule(::Type{SpaceOperation{Tp1, Tt1}}, ::Type{PointOperation{Tp2}}) where {Tp1, Tt1, Tp2}
    Tp = promote_type(Tp1, Tp2)
    return SpaceOperation{Tp, Tt1}
end

function Base.promote_rule(::Type{TranslationOperation{Tt}}, ::Type{PointOperation{Tp}}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}
end


## Properties

"""
    isidentity(arg::SpaceOperation)

Check whether `arg` is identity (i.e. identity point and identity translation)
"""
isidentity(arg::SpaceOperation) = isone(arg.matrix) && iszero(arg.displacement)

"""
    istranslation(arg::SpaceOperation)

Check whether `arg` is a translation operation (i.e., point operation is identity)
"""
istranslation(arg::SpaceOperation) = isone(arg.matrix)

"""
    ispoint(arg::SpaceOperation)

Check whether `arg` is a point operation (i.e. translation is identity)
"""
ispoint(arg::SpaceOperation) = iszero(arg.displacement)

"""
    dimension(arg::SpaceOperation)

Return the spatial dimension of `arg`.
"""
dimension(arg::SpaceOperation) = length(arg.displacement)


function Base.hash(op::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    h = hash(SpaceOperation{Tp, Tt})
    h = hash(op.matrix, h)
    h = hash(op.displacement, h)
    return h
end


## operators
function Base.:(==)(lhs::SpaceOperation{Tp, Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return lhs.matrix == rhs.matrix && lhs.displacement == rhs.displacement
end

function Base.:(==)(lhs::SpaceOperation{Tp, Tt}, rhs::TranslationOperation{Tt}) where {Tp, Tt}
    return isone(lhs.matrix) && lhs.displacement == rhs.displacement
end

function Base.:(==)(lhs::SpaceOperation{Tp, Tt}, rhs::PointOperation{Tp}) where {Tp, Tt}
    return iszero(lhs.displacement) && lhs.matrix == rhs.matrix
end

function Base.:(==)(lhs::TranslationOperation{Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return isone(rhs.matrix) && lhs.displacement == rhs.displacement
end

function Base.:(==)(lhs::PointOperation{Tp}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return iszero(rhs.displacement) && lhs.matrix == rhs.matrix
end

function Base.:(==)(sop::SpaceOperation{Tp, Tt}, iden::IdentityOperation{<:Union{Tp, Tt}}) where {Tp, Tt}
    return dimension(sop) == dimension(iden) && iszero(sop.displacement) && isone(sop.matrix)
end

function Base.:(==)(iden::IdentityOperation{<:Union{Tp, Tt}}, sop::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    return dimension(sop) == dimension(iden) && iszero(sop.displacement) && isone(sop.matrix)
end

function Base.:(*)(lhs::PointOperation{Tp}, rhs::TranslationOperation{Tt}) where {Tp, Tt}
    return SpaceOperation{Tp, Tt}(lhs.matrix, rhs.displacement)
end

function Base.:(*)(lhs::TranslationOperation{Tt}, rhs::PointOperation{Tp}) where {Tp, Tt}
    rhs_matrix_inv = Matrix{Tp}(LinearAlgebraX.invx(rhs.matrix))
    return SpaceOperation{Tp, Tt}(rhs.matrix, rhs_matrix_inv * lhs.displacement)
end

#   ML ⋅ ( MR ⋅ ( x + ρR ) + ρL )
# = ML MR x + ML MR ρR + ML ρL
# = ML MR (x + ρR + inv(MR) ⋅ ρL )
function Base.:(*)(lhs::SpaceOperation{Tp, Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix = lhs.matrix * rhs.matrix
    rhs_matrix_inv = Matrix{Tp}(LinearAlgebraX.invx(rhs.matrix))
    displacement = rhs.displacement + rhs_matrix_inv * lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function Base.:(*)(lhs::SpaceOperation{Tp, Tt}, rhs::TranslationOperation{Tt}) where {Tp, Tt}
    matrix = lhs.matrix
    displacement = rhs.displacement + lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function Base.:(*)(lhs::TranslationOperation{Tt}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix = rhs.matrix
    rhs_matrix_inv = Matrix{Tp}(LinearAlgebraX.invx(rhs.matrix))
    displacement = rhs.displacement + rhs_matrix_inv * lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function Base.:(*)(lhs::SpaceOperation{Tp, Tt}, rhs::PointOperation{Tp}) where {Tp, Tt}
    matrix = lhs.matrix * rhs.matrix
    rhs_matrix_inv = Matrix{Tp}(LinearAlgebraX.invx(rhs.matrix))
    displacement = rhs_matrix_inv * lhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function Base.:(*)(lhs::PointOperation{Tp}, rhs::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix = lhs.matrix * rhs.matrix
    displacement = rhs.displacement
    return SpaceOperation{Tp, Tt}(matrix, displacement)
end

function Base.inv(arg::SpaceOperation{Tp, Tt}) where {Tp, Tt}
    matrix_inv = Matrix{Tp}(LinearAlgebraX.invx(arg.matrix))
    return SpaceOperation{Tp, Tt}(matrix_inv, -arg.matrix * arg.displacement)
end

function Base.:(^)(op::SpaceOperation{Tp, Tt}, power::Integer) where {Tp, Tt}
    if power == 0
        return SpaceOperation(Tp, Tt, dimension(op))
    elseif power < 0
        op_inv = inv(op)
        return op_inv^(-power)
    else
        return op * op^(power-1)
    end
end


"""
    apply_operation(op::SpaceOperation{Tp, Tt}, coord::AbstractArray{<:Union{Tp, Tt}}) where {Tp, Tt}
"""
function apply_operation(op::SpaceOperation{Tp, Tt}, coord::AbstractArray{<:Union{Tp, Tt}}) where {Tp, Tt}
    return op.matrix * (coord .+ op.displacement)
end

function (op::SpaceOperation{Tp, Tt})(coord::AbstractArray{<:Union{Tp, Tt}}) where {Tp, Tt}
    return op.matrix * (coord .+ op.displacement)
end
