# COV_EXCL_START
struct LocalUnitary{S<:Number}<:AbstractSymmetryOperationEmbedding
    operations::Dict{Int, Matrix{S}}
    function LocalUnitary(ops::AbstractDict{<:Integer, <:AbstractMatrix{S}}) where {S<:Number}
        error("Local unitary not implemented")
        return new{S}(ops)
    end

    function LocalUnitary(pairs::Pair{<:Integer, <:AbstractMatrix{<:Number}}...)
        error("Local unitary not implemented")
        LocalUnitary(Dict(pairs...))
    end
end
# COV_EXCL_STOP
