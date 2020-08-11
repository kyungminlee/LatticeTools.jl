export UnitCell
export make_unitcell, makeunitcell
export dimension

export numsite, sitecount,
       hassite,
       addsite!,
       getsite,
       getsiteindex,
       getsitecoord,
       getsiteindexcoord,
       getsitename

export carte2fract,
       fract2carte,
       whichunitcell,
       momentumgrid

export findsiteindex

import LinearAlgebra

export numorbital,
       hasorbital,
       addorbital!,
       getorbital,
       getorbitalindex,
       getorbitalcoord,
       getorbitalindexcoord,
       getorbitalname

@deprecate numorbital(args...) numsite(args...)
@deprecate hasorbital(args...) hassite(args...)
@deprecate addorbital!(args...) addsite!(args...)
@deprecate getorbital(args...) getsite(args...)
@deprecate getorbitalindex(args...) getsiteindex(args...)
@deprecate getorbitalcoord(args...) getsitecoord(args...)
@deprecate getorbitalindexcoord(args...) getsiteindexcoord(args...)
@deprecate getorbitalname(args...) getsitename(args...)
@deprecate findorbitalindex(args...; kwargs...) findsiteindex(args...; kwargs...)

"""
    UnitCell{O}

Represent a unitcell of a lattice, which contains sites at fixed locations
(does not yet implement multiple orbitals per site). It is recommended to use `makeunitcell`
rather than the constructor to make a UnitCell object.

# Parameters
- `O`: type of "site". Any type can be used, but we recommend using `String`
  or tuple of `String` and `Int` for compatibility with JSON.

# Fields
- `latticevectors::Array{Float64, 2}`: Lattice vectors
- `reducedreciprocallatticevectors::Array{Float64, 2}`: Reduced reciprocal lattice vectors
   (transpose of inverse of `latticevectors`)
- `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors.
  `2π * reducedreciprocallatticevectors`
- `sites::Vector{Tuple{T, FractCoord}}`: List of sites within unit cell
- `siteindices::Dict{T, Int}`: Indices of sites
"""
mutable struct UnitCell{S}
    latticevectors::Array{Float64, 2}
    sites::Vector{Tuple{S, FractCoord}}

    reducedreciprocallatticevectors::Array{Float64, 2}
    reciprocallatticevectors::Array{Float64, 2}
    siteindices::Dict{S, Int}

    function UnitCell{S}(
        latticevectors::AbstractArray{<:Real, 2},
        sites::AbstractVector{Tuple{S, FractCoord}},
        reducedreciprocallatticevectors::AbstractArray{<:Real, 2},
        reciprocallatticevectors::AbstractArray{<:Real, 2},
        siteindices::AbstractDict{S, Int}
    ) where {S}
        if (S <: Integer)
            throw(ArgumentError("Integer SiteType is disallowed to avoid confusion"))
        end
        new{S}(
            latticevectors,
            sites,
            reducedreciprocallatticevectors,
            reciprocallatticevectors,
            siteindices,
        )
    end
end


"""
    makeunitcell(latticevectors; SiteType=Any, tol=√ϵ)

Construct an n-dimensional lattice.

# Arguments
* `latticevectors`: Lattice vectors. Can be a nested array of lattice vectors,
  or a two-dimensional array whose columns are lattice vectors, or a real number.

# Optional Arguments
* `SiteType::DataType`
* `tol=√ϵ`: Epsilon
"""
function makeunitcell(
    latticevectors::AbstractArray{<:Real, 2};
    SiteType::DataType=Any,
    OrbitalType::DataType=Nothing,
    tol::Real=Base.rtoldefault(Float64)
)
    if OrbitalType != Nothing
        @warn "OrbitalType is deprecated. use SiteType instead"
        if SiteType == Any
            SiteType = OrbitalType
        else
            @error "Cannot specify OrbitalType and SiteType simultaneously"
        end
    end

    (ndim, ndim_) = size(latticevectors)
    if ndim != ndim_
        throw(ArgumentError("lattice vectors has dimension ($(ndim), $(ndim_))"))
    elseif abs(LinearAlgebra.det(latticevectors)) <= tol
        throw(ArgumentError("lattice vectors define zero volume $(latticevectors)"))
    elseif tol < 0
        throw(ArgumentError("tol must be non-negative"))
    end

    reduced_rlv = transpose(inv(latticevectors))
    sites = Tuple{SiteType, FractCoord}[]
    siteindices = Dict{SiteType, Int}()
    return UnitCell{SiteType}(latticevectors, sites, reduced_rlv, 2*π*reduced_rlv, siteindices)
end

function makeunitcell(
    latticevectors::AbstractVector{<:AbstractVector};
    SiteType::DataType=Any,
    tol::Real=Base.rtoldefault(Float64)
)
    lv = hcat(latticevectors...)
    return makeunitcell(lv; SiteType=SiteType, tol=tol)
end

function makeunitcell(
    latticeconstant::Real;
    SiteType::DataType=Any,
    OrbitalType::DataType=Nothing,
    tol::Real=Base.rtoldefault(Float64),
)
    return makeunitcell(
        reshape([latticeconstant], (1,1));
        SiteType=SiteType, OrbitalType=OrbitalType, tol=tol,
    )
end


make_unitcell = makeunitcell


"""
    dimension(uc)

Spatial dimension of the unit cell.
"""
dimension(uc::UnitCell) = size(uc.latticevectors, 1)


"""
    numsite(unitcell)

Number of sites of the unit cell.

# Arguments
* `uc ::UnitCell`
"""
numsite(uc::UnitCell) = length(uc.sites)



"""
    sitecount(unitcell)

Number of sites of the unit cell.

# Arguments
* `uc ::UnitCell`
"""
sitecount(uc::UnitCell) = length(uc.sites)


"""
    addsite!(unitcell, sitename, sitecoord)

Add an site to the unit cell.

# Arguments
* `uc ::UnitCell{T}`
* `sitename ::{T}`
* `sitecoord ::FractCoord`
"""
function addsite!(uc::UnitCell{O}, sitename::O, sitecoord::FractCoord) where {O}
    (ndim, ndim_) = size(uc.latticevectors)
    if dimension(sitecoord) != ndim
        throw(ArgumentError("sitecoord has wrong dimension"))
    elseif haskey(uc.siteindices, sitename)
        throw(ArgumentError( "duplicate site name"))
    end

    for site in uc.sites
        if isapprox(site[2].fraction, sitecoord.fraction)
            throw(ArgumentError("site $(site[1]) is already at $(sitecoord.fraction)"))
        end
    end

    push!(uc.sites, (sitename, sitecoord))
    index = length(uc.siteindices)+1
    uc.siteindices[sitename] = index
    return index
end


"""
    hassite(unitcell, name)

Test whether the unit cell contains the site of given name.

# Arguments
* `unitcell::UnitCell{O}`
* `name::O`
"""
function hassite(uc::UnitCell{O}, name::O) where {O}
    return haskey(uc.siteindices, name)
end


"""
    getsiteindex(unitcell, name)

Get index of the given site.

# Arguments
* `unitcell::UnitCell{O}`
* `name::O`
"""
function getsiteindex(uc::UnitCell{O}, name::O) where {O}
    return uc.siteindices[name]
end


"""
    getsite(unitcell, name)

Get the site (its site name and its fractional coordinates) with the given name.

# Arguments
* `unitcell::UnitCell{O}`
* `name::O`

# Return
* `(sitename, fractcoord)`
"""
function getsite(uc::UnitCell{O}, name::O) where {O}
    return uc.sites[uc.siteindices[name]]
end


"""
    getsitecoord(unitcell, name)

Get the fractional coordinates of the site with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `fractcoord`
"""
function getsitecoord(uc::UnitCell{O}, name::O) where {O}
    return getsite(uc, name)[2]
end


"""
    getsiteindexcoord(unitcell, name)

# Arguments
* `unitcell::UnitCell{T}`
* `name::T`

# Return
* `(index, fractcoord)`
"""
function getsiteindexcoord(uc::UnitCell{O}, name::O) where {O}
    index = getsiteindex(uc, name)
    coord = getsitecoord(uc, index)
    return (index, coord)
end


"""
    getsite(unitcell, index)

# Arguments
* `unitcell::UnitCell`
* `index::Integer`

# Return
* `(sitename, fractcoord)`
"""
function getsite(uc::UnitCell, index::Integer)
    return uc.sites[index]
end


"""
    getsitename(unitcell, index)

# Arguments
* `uc::UnitCell`
* `index::Integer`

# Return
* `sitename`
"""
function getsitename(uc::UnitCell{O}, index::Integer)::O where {O}
    return uc.sites[index][1]
end


"""
    getsitecoord(unitcell, index)

# Arguments
* `unitcell::UnitCell`
* `index::Integer`

# Return
* `FractCoord`
"""
function getsitecoord(uc::UnitCell, index::Integer)::FractCoord
    return uc.sites[index][2]
end


"""
    fract2carte(unitcell, fractcoord)

# Arguments
* `unitcell::UnitCell`
* `fractcoord::FractCoord`
"""
function fract2carte(unitcell::UnitCell, fc::FractCoord)::CarteCoord
    if dimension(unitcell) != dimension(fc)
        throw(ArgumentError("unitcell and fractcoord must have the same dimension"))
    end
    mc = fc.whole + fc.fraction
    cc = unitcell.latticevectors * mc
    return CarteCoord(cc)
end


"""
    carte2fract(unitcell, cartecoord; tol=√ϵ)

# Arguments
* `unitcell::UnitCell`
* `cc::CarteCoord`
* `tol::Real=Base.rtoldefault(Float64)`
"""
function carte2fract(
    unitcell::UnitCell, cc::CarteCoord;
    tol::Real=Base.rtoldefault(Float64)
)
    if dimension(unitcell) != length(cc)
        throw(ArgumentError("unitcell and cartecoord must have the same dimension"))
    end
    fc = transpose(unitcell.reducedreciprocallatticevectors) * cc
    w = Int[fld(x, 1) for x in fc]
    f = Float64[mod(x, 1) for x in fc]
    return FractCoord(w, f; tol=tol)
end


"""
whichunitcell(unitcell, name, cartecoord; tol=√ϵ)

# Return
- `Vector{Int}`: the integer coordinate of the unitcell that the specified site/cartesian
  coordinate belongs to.
"""
function whichunitcell(
    uc::UnitCell{O},
    name::O,
    cc::CarteCoord;
    tol::Real=Base.rtoldefault(Float64)
) where {O}
    fc1 = getsitecoord(uc, name)
    fc2 = carte2fract(uc, cc)
    if !isapprox(fc1.fraction, fc2.fraction; atol=tol)
        throw(ArgumentError("$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]"))
    end
    R = fc2.whole - fc1.whole
    return R
end


"""
    whichunitcell(unitcell, name, fractcoord; tol=√ϵ)

# Return
- `Vector{Int}`: the integer coordinate of the unitcell that the specified site/fractional
  coordinate belongs to.
"""
function whichunitcell(
    uc::UnitCell{O},
    name::O,
    fc::FractCoord;
    tol::Real=Base.rtoldefault(Float64)
) where {O}
    fc1 = getsitecoord(uc, name)
    fc2 = fc
    if !isapprox(fc1.fraction, fc2.fraction; atol=tol)
        throw(ArgumentError("$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]"))
    end
    R = fc2.whole - fc1.whole
    return R
end


#=
function zero(uc::UnitCell; dtype::DataType=ComplexF64)
    norb = numsite(uc)
    return Base.zeros(dtype, (norb, norb))
end
=#


"""
    momentumgrid(unitcell, shape)

Generate an n-dimensional grid of momenta of given shape.

Returns an n-dimensional array of the following form:
```math
   k[i_1, i_2, \\ldots ] = G \\cdot ( \\frac{i_1}{n_1}, \\frac{i_2}{n_2}, \\ldots )
````
where G is the reciprocal lattice vector
"""
function momentumgrid(uc::UnitCell, shape::AbstractVector{<:Integer})
    if length(shape) != dimension(uc)
        throw(ArgumentError("dimension mismatch"))
    elseif !all((x) -> x>0, shape)
        throw(ArgumentError("shape should be positive"))
    end
    ranges = [range(0,stop=1,length=n+1)[1:end-1] for n in shape]
    cubicgrid = map((x) -> [x...], Base.product(ranges...))
    momentumgrid = map((x) -> uc.reciprocallatticevectors * x, cubicgrid)
    return momentumgrid
end


function Base.:(==)(lhs::UnitCell{O}, rhs::UnitCell{O}) where O
    return (
        lhs.latticevectors == rhs.latticevectors &&
        lhs.sites == rhs.sites
    )
end


"""
    findsiteindex(unitcell::UnitCell, fc::FractCoord; tol=√ϵ)

Find the index of the site, and the unitcell coordinate at the specified fractional coordinate.

# Arguments
- `unitcell::UnitCell`
- `fc::FractCoord`

# Returns
- `(site_index, unitcell_vector)`, or `(-1, [])` if not found.
"""
function findsiteindex(
    unitcell::UnitCell,
    fc::FractCoord;
    tol::Real=Base.rtoldefault(Float64)
)
    i = findfirst(
        x -> isapprox(x, fc.fraction; atol=tol),
        [orbfc.fraction for (orbname, orbfc) in unitcell.sites]
    )
    if i !== nothing
        (orbname, orbfc) = unitcell.sites[i]
        return (i, fc.whole - orbfc.whole)
    else
        return (-1, Int[])
    end
end
