export OrthoCube
export find_generators
export generate_coordinates
export volume
export isequiv


"""
    OrthoCube(shape)

Represent a hypercubic (Bravais) lattice.

# Fields
* `shape_matrix`: a matrix whose columns are the lattice vectors of the
* `inverse_shape_matrix`: `(shape_matrix)`
* `wrap`: periodic wrapping function which takes an integer array, and maps it onto a
  site in the Bravais lattice. Analogous to `fldmod`.
"""
struct OrthoCube
    shape_matrix::Matrix{Int}
    inverse_shape_matrix::Matrix{Rational{Int}}
    wrap::Function

    function OrthoCube(shape_matrix::AbstractMatrix{<:Integer})
        dim, dim2 = size(shape_matrix)
        if dim != dim2
            msg = "scale_matrix is not square: dimensions are ($dim, $dim2)"
            throw(DimensionMismatch(msg))
        end
        det = ExactLinearAlgebra.determinant(shape_matrix)
        det == 0 && throw(ArgumentError("scale matrix null"))

        inverse_shape_matrix = ExactLinearAlgebra.inverse(shape_matrix)
        function wrap(r::AbstractArray{<:Integer}, mode::RoundingMode=RoundDown)
            rnd = (x) -> round(Int, x, mode)
            R = rnd.(inverse_shape_matrix * r)
            r2 = r - shape_matrix * R
            return R, r2
        end
        return new(shape_matrix, inverse_shape_matrix, wrap)
    end
end

"""
    dimension(ortho::OrthoCube)

Return the spatial dimension of the orthocube.
"""
dimension(ortho::OrthoCube) = size(ortho.shape_matrix, 1)

"""
    volume(ortho::OrthoCube)

Return the signed volume of the orthocube, defined by the determinant of the shape.
"""
volume(ortho::OrthoCube) = ExactLinearAlgebra.determinant(ortho.shape_matrix)


function Base.:(==)(lhs::OrthoCube, rhs::OrthoCube)
    return (lhs.shape_matrix == rhs.shape_matrix)
end


"""
    isequiv(lhs::OrthoCube, rhs::OrthoCube)

Check whether the two orthocubes are equivalent.
"""
function isequiv(lhs::OrthoCube, rhs::OrthoCube)
    R, r = lhs.wrap(rhs.shape_matrix)
    if abs(ExactLinearAlgebra.determinant(R)) != 1 || !iszero(r)
        return false
    end
    return true
end

# function isequiv_old(lhs::OrthoCube, rhs::OrthoCube)
#     det_lhs = LatticeTools.ExactLinearAlgebra.determinant(lhs.shape_matrix)
#     det_rhs = LatticeTools.ExactLinearAlgebra.determinant(rhs.shape_matrix)
#     det_lhs != det_rhs && return false
#     inv_lhs = LatticeTools.ExactLinearAlgebra.inverse(lhs.shape_matrix)
#     inv_rhs = LatticeTools.ExactLinearAlgebra.inverse(rhs.shape_matrix)
#     return all(isinteger.(inv_lhs * rhs.shape_matrix)) && all(isinteger.(inv_rhs * lhs.shape_matrix))
# end


"""
    find_generators(ortho::OrthoCube)

Find generators of an OrthoCube. For one-dimensional lattice, which is isomorphic to Zₙ, the
generator is +1. For two-dimensional lattices, this function invokes
[`find_generators_2d`](@ref). Higher dimensions are currently not supported.
"""
function find_generators(ortho::OrthoCube)
    if dimension(ortho) == 1
        return ones(Int, 1, 1)
    elseif dimension(ortho) == 2
        return find_generators_2d(ortho)
    else
        error("Currently, only 1D and 2D lattices are supported")
    end
end


"""
    find_generators_2d(ortho::OrthoCube)

Find translation generators of an OrthoCube.
"""
function find_generators_2d(ortho::OrthoCube)
    dimension(ortho) != 2 && throw(ArgumentError("shape matrix should be 2x2"))

    orthovolume = abs(volume(ortho))
    allowed_pairs = [[1,0]]
    append!(
        allowed_pairs,
        [x, y] for y in 1:orthovolume, x in 0:orthovolume if gcd(x,y) == 1
    )

    function make_loop(t::Vector{Int})
        axis = Vector{Int}[zero(t)]
        r = ortho.wrap(t)[2]
        while !iszero(r)
            push!(axis, r)
            r = ortho.wrap(r .+ t)[2]
        end
        @assert allunique(axis)
        return axis
    end

    for r1p in allowed_pairs
        r1p[1] == 0 && continue
        r1 = [r1p[1], -r1p[2]]
        axis1 = make_loop(r1)
        mod(orthovolume, length(axis1)) != 0 && continue
        for r2 in allowed_pairs
            ExactLinearAlgebra.determinant(hcat(r1, r2)) != 1 && continue
            axis2 = make_loop(r2)
            orthovolume != length(axis1) * length(axis2) && continue
            coordinates = vec([ortho.wrap(sum(x))[2]
                                   for x in Iterators.product(axis1, axis2)])
            allunique(coordinates) && return hcat(r1, r2)
        end
    end
    error("Failed to find decomposition") # COV_EXCL_LINE
end


"""
    generate_coordinates(ortho, generators)

Generate a list of coordinates of the orthocube

# Arguments
* `ortho::OrthoCube`
* `generator_translations::AbstractMatrix{<:Integer}`
"""
function generate_coordinates(
    ortho::OrthoCube,
    generator_translations::AbstractMatrix{<:Integer},
)
    let dim = dimension(ortho)
        if size(generator_translations) != (dim, dim)
            throw(DimensionMismatch("OrthoCube and generator_translations have different dimensions"))
        elseif abs(ExactLinearAlgebra.determinant(generator_translations)) != 1
            throw(ArgumentError("generator_translations is not unimodular"))
        end
    end
    function make_loop(t::AbstractVector{<:Integer})
        axis = Vector{Int}[zero(t)]
        r = ortho.wrap(t)[2]
        while !iszero(r)
            push!(axis, r)
            r = ortho.wrap(r .+ t)[2]
        end
        return axis
    end
    axes = [make_loop(t) for t in eachcol(generator_translations)]
    coordinates = vec([ortho.wrap(sum(x))[2] for x in Iterators.product(axes...)])

    if length(coordinates) != abs(volume(ortho))
        msg = "translations $generator_translations generates $(length(coordinates)) " *
            "coordinates, while the volume is $(abs(volume(ortho)))"
        throw(ArgumentError(msg))
    end
    @assert allunique(coordinates)
    return coordinates
end
