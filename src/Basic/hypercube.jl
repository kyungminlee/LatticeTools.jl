export Hypercube
export find_generators
export generate_coordinates
export volume
export isequiv

import MathExpr


"""
    Hypercube(shape)

Represent a hypercubic (Bravais) lattice.

# Fields
* `shape_matrix`: a matrix whose columns are the lattice vectors of the
* `inverse_shape_matrix`: `(shape_matrix)`
* `wrap`: periodic wrapping function which takes an integer array, and maps it onto a
  site in the Bravais lattice. Analogous to `fldmod`.
"""
struct Hypercube
    shape_matrix::Matrix{Int}
    inverse_shape_matrix::Matrix{Rational{Int}}
    wrap::Function

    function Hypercube(shape_matrix::AbstractMatrix{<:Integer})
        dim, dim2 = size(shape_matrix)
        if dim != dim2
            msg = "scale_matrix is not square: dimensions are ($dim, $dim2)"
            throw(DimensionMismatch(msg))
        end
        det = MathExpr.determinant(shape_matrix)
        det == 0 && throw(ArgumentError("scale matrix null"))

        inverse_shape_matrix = MathExpr.inverse(shape_matrix)
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
    dimension(cube::Hypercube)

Return the spatial dimension of the hypercube.
"""
dimension(cube::Hypercube) = size(cube.shape_matrix, 1)


"""
    volume(cube::Hypercube)

Return the signed volume of the hypercube, defined by the determinant of the shape.
"""
volume(cube::Hypercube) = MathExpr.determinant(cube.shape_matrix)


function Base.:(==)(lhs::Hypercube, rhs::Hypercube)
    return (lhs.shape_matrix == rhs.shape_matrix)
end


"""
    isequiv(lhs::Hypercube, rhs::Hypercube)

Check whether the two hypercubes are equivalent.
"""
function isequiv(lhs::Hypercube, rhs::Hypercube)
    R, r = lhs.wrap(rhs.shape_matrix)
    if abs(MathExpr.determinant(R)) != 1 || !iszero(r)
        return false
    end
    return true
end


"""
    find_generators(cube::Hypercube)

Find generators of an Hypercube. For one-dimensional lattice, which is isomorphic to Zₙ, the
generator is +1. For two-dimensional lattices, this function invokes
[`find_generators_2d`](@ref). Higher dimensions are currently not supported.
"""
function find_generators(cube::Hypercube)
    if dimension(cube) == 1
        return ones(Int, 1, 1)
    elseif dimension(cube) == 2
        return find_generators_2d(cube)
    else
        error("Currently, only 1D and 2D lattices are supported")
    end
end


"""
    find_generators_2d(cube::Hypercube)

Find translation generators of an Hypercube.
"""
function find_generators_2d(cube::Hypercube)
    dimension(cube) != 2 && throw(ArgumentError("shape matrix should be 2x2"))

    cubevolume = abs(volume(cube))
    allowed_pairs = [[1,0]]
    append!(
        allowed_pairs,
        [x, y] for y in 1:cubevolume, x in 0:cubevolume if gcd(x,y) == 1
    )

    function make_loop(t::Vector{Int})
        axis = Vector{Int}[zero(t)]
        r = cube.wrap(t)[2]
        while !iszero(r)
            push!(axis, r)
            r = cube.wrap(r .+ t)[2]
        end
        @assert allunique(axis)
        return axis
    end

    for r1p in allowed_pairs
        r1p[1] == 0 && continue
        r1 = [r1p[1], -r1p[2]]
        axis1 = make_loop(r1)
        mod(cubevolume, length(axis1)) != 0 && continue
        for r2 in allowed_pairs
            MathExpr.determinant(hcat(r1, r2)) != 1 && continue
            axis2 = make_loop(r2)
            cubevolume != length(axis1) * length(axis2) && continue
            coordinates = vec([cube.wrap(sum(x))[2]
                                   for x in Iterators.product(axis1, axis2)])
            allunique(coordinates) && return hcat(r1, r2)
        end
    end
    error("Failed to find decomposition") # COV_EXCL_LINE
end


"""
    generate_coordinates(cube, generators)

Generate a list of coordinates of the hypercube

# Arguments
* `cube::Hypercube`
* `generator_translations::AbstractMatrix{<:Integer}`
"""
function generate_coordinates(
    cube::Hypercube,
    generator_translations::AbstractMatrix{<:Integer},
)
    let dim = dimension(cube)
        if size(generator_translations) != (dim, dim)
            throw(DimensionMismatch("Hypercube and generator_translations have different dimensions"))
        elseif abs(MathExpr.determinant(generator_translations)) != 1
            throw(ArgumentError("generator_translations is not unimodular"))
        end
    end
    function make_loop(t::AbstractVector{<:Integer})
        axis = Vector{Int}[zero(t)]
        r = cube.wrap(t)[2]
        while !iszero(r)
            push!(axis, r)
            r = cube.wrap(r .+ t)[2]
        end
        return axis
    end
    axes = [make_loop(t) for t in eachcol(generator_translations)]
    coordinates = vec([cube.wrap(sum(x))[2] for x in Iterators.product(axes...)])

    if length(coordinates) != abs(volume(cube))
        msg = "translations $generator_translations generates $(length(coordinates)) " *
            "coordinates, while the volume is $(abs(volume(cube)))"
        throw(ArgumentError(msg))
    end
    @assert allunique(coordinates)
    return coordinates
end
