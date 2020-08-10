export CarteCoord
export FractCoord
export dimension
export fract2carte, carte2fract


"""
    CarteCoord

Cartesian coordinates. `Vector{Float64}`.
"""
const CarteCoord = Vector{Float64}


"""
    FractCoord

Fractional coordinates.

# Members
* `whole ::Vector{Int}`: Integer part of fractional coordinates
* `fraction ::Vector{Float64}`: [0,1) part of fractional coordinates
"""
struct FractCoord  # TODO: Type parameter for integer and float
    whole::Vector{Int}
    fraction::Vector{Float64}

    function FractCoord(
        w::AbstractVector{<:Integer}, f::AbstractVector{<:Real};
        tol::Real=Base.rtoldefault(Float64)
    )
        if length(w) != length(f)
            throw(ArgumentError("w and f need to be of same the length"))
        elseif tol < 0
            throw(ArgumentError("tol must be non-negative"))
        elseif !all(x -> 0 <= x < 1+tol, f)
            throw(ArgumentError("f must be a list of floating points in [0, 1) (got $f)"))
        end

        neww = copy(w)
        newf = copy(f)
        for (idx, (wval, fval)) in enumerate(zip(neww, newf))
            if fval <= tol
                newf[idx] = 0.0
            elseif fval + tol >= 1.0
                newf[idx] = 0.0
                neww[idx] += 1
            end
        end
        return new(neww, newf)
    end

    function FractCoord(coord::AbstractVector{<:Real}; tol::Real=Base.rtoldefault(Float64))
        w = Int[fld(x,1) for x in coord]
        f = Float64[mod(x,1) for x in coord]
        return FractCoord(w, f; tol=tol)
    end

    function FractCoord(ndim::Integer)
        ndim <= 0 && throw(ArgumentError("ndim should be positive"))
        w = zeros(Int, ndim)
        f = zeros(Float64, ndim)
        return new(w, f)
    end
end


"""
    dimension(fc::FractCoord)

Dimension of the fractional coordinates

# Arguments
* `fc::FractCoord`: Fractional coordinates.
"""
dimension(fc::FractCoord) = length(fc.whole)

function Base.:(+)(fc::FractCoord)::FractCoord
    return fc
end

function Base.:(-)(fc::FractCoord)::FractCoord
    return FractCoord(-(fc.whole + fc.fraction))
end

function Base.:(+)(fractcoord::FractCoord, R::AbstractVector{<:Integer})::FractCoord
    return FractCoord(fractcoord.whole + R, fractcoord.fraction)
end

function Base.:(-)(fractcoord::FractCoord, R::AbstractVector{<:Integer})::FractCoord
    return FractCoord(fractcoord.whole - R, fractcoord.fraction)
end

function Base.:(+)(fc1::FractCoord, fc2::FractCoord)::FractCoord
    R = Int[fld(x+y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
    r = Float64[mod(x+y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
    return FractCoord(fc1.whole + fc2.whole + R, r)
end

function Base.:(-)(fc1::FractCoord, fc2::FractCoord) ::FractCoord
    R = Int[fld(x-y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
    r = Float64[mod(x-y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
    return FractCoord(fc1.whole - fc2.whole + R, r)
end

function Base.isapprox(
    fc1::FractCoord, fc2::FractCoord;
    atol::Real=Base.rtoldefault(Float64), rtol::Real=Base.rtoldefault(Float64)
)::Bool
    return (
        (fc1.whole == fc2.whole) &&
        isapprox(fc1.fraction, fc2.fraction; atol=atol, rtol=rtol)
    )
end

function Base.:(*)(mat::AbstractMatrix{<:Number}, fc::FractCoord)::FractCoord
    return FractCoord(mat * (fc.whole + fc.fraction))
end

function Base.:(==)(fc1::FractCoord, fc2::FractCoord)::Bool
    return (fc1.whole == fc2.whole) && (fc1.fraction == fc2.fraction)
end

function Base.show(io::IO, fc::FractCoord)
    print(io, "FractCoord(", fc.whole, " + ", fc.fraction,")")
end


"""
    fract2carte(latticevectors, fc)

# Arguments
* `latticevectors::AbstractArray{<:Real, 2}`: square matrix whose columns are lattice vectors.
* `fc::FractCoord`: fractional coordinates
"""
function fract2carte(latticevectors::AbstractMatrix{<:Real}, fc::FractCoord)::CarteCoord
    d1, d2 = size(latticevectors)
    d3 = dimension(fc)
    if d1 != d2
        throw(ArgumentError("latticevectors should be a square matrix "))
    elseif d1 != d3
        throw(ArgumentError("Dimensions of latticevectors and fc "))
    end
    mc = fc.whole + fc.fraction
    cc = latticevectors * mc
    return CarteCoord(cc)
end


"""
    carte2fract(latticevectors, cc; tol=√ϵ)

# Arguments
* `latticevectors::AbstractArray{<:Real, 2}`: square matrix whose columns are lattice vectors.
* `cc::CarteCoord`: cartesian coordinates
* `tol::Real=Base.rtoldefault(Float64)`: tolerance
"""
function carte2fract(
    latticevectors::AbstractMatrix{<:Real},
    cc::CarteCoord;
    tol::Real=Base.rtoldefault(Float64)
)::FractCoord
    fc = inv(latticevectors) * cc
    w = Int[fld(x, 1) for x in fc]
    f = Float64[mod(x, 1) for x in fc]
    return FractCoord(w, f; tol=tol)
end
