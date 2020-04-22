export linpath
export momentumpath

"""
    momentumpath

Generate a list of momenta

# Arguments
* `anchorpoints`
* (Optional) `nseg` - number of points in each segment
"""
function linpath(anchors::AbstractVector...; nseg::Integer=100)
    n = length(anchors)
    if n <= 1
        throw(ArgumentError("Number of anchor points should be more than one"))
    end
    d = length(anchors[1])
    for i in 2:n
        if length(anchors[i]) != d
            throw(DimensionMismatch("All anchor points should have the same dimension"))
        end
    end
    out = Matrix{Float64}(undef, (d, (n-1) * nseg + 1))
    for i in 1:n-1
        from = anchors[i]
        to = anchors[i+1]
        dvec = (to - from) ./ nseg
        for j in 0:(nseg-1)
            out[:, (i-1)*nseg + j + 1] = from + dvec * j
        end
    end
    out[:, end] = anchors[end]
    return out
end

# """
#     momentumpath
#
# The anchorpoints are given in units of the reciprocal lattice vectors.
# """
# function momentumpath(unitcell::UnitCell,
#                       anchors::AbstractVector{<:AbstractVector{<:Number}})
#     real_anchorpoints = [unitcell.reciprocallatticevectors * ap for ap in anchorpoints]
#     return linpath(real_anchorpoints...)
# end
