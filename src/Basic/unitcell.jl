export UnitCell
export make_unitcell
export dimension

export numsite,
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



"""
    UnitCell{O}

# Parameters
* `O`: type of "site". Any type can be used, but we recommend using `String` or tuple of `String` and `Int`
       for compatibility with JSON.

# Members
* `latticevectors ::Array{Float64, 2}`: Lattice vectors
* `reducedreciprocallatticevectors ::Array{Float64, 2}`: Reduced reciprocal lattice vectors (transpose of inverse of `latticevectors`)
* `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors. `2π * reducedreciprocallatticevectors`
* `sites ::Vector{Tuple{T, FractCoord}}`: List of sites within unit cell
* `siteindices ::Dict{T, Int}`: Indices of sites
"""
mutable struct UnitCell{O}
    latticevectors ::Array{Float64, 2}
    sites ::Vector{Tuple{O, FractCoord}}

    reducedreciprocallatticevectors ::Array{Float64, 2}
    reciprocallatticevectors ::Array{Float64, 2}
    siteindices ::Dict{O, Int}
    function UnitCell{O}(latticevectors ::AbstractArray{<:Real, 2},
                         sites ::AbstractVector{Tuple{O, FractCoord}},
                         reducedreciprocallatticevectors ::AbstractArray{<:Real, 2},
                         reciprocallatticevectors ::AbstractArray{<:Real, 2},
                         siteindices ::AbstractDict{O, Int}) where {O}
        if (O <: Integer)
            throw(ArgumentError("SiteType should not be integer to avoid confusion"))
        end
        new{O}(latticevectors,
               sites,
               reducedreciprocallatticevectors,
               reciprocallatticevectors,
               siteindices)
    end
end


"""
    UnitCell

Construct a one-dimensional lattice.

# Arguments
* `latticeconstant ::Float64`: Lattice constant
* `SiteType`: List of sites

# Optional Arguments
* `tol=√ϵ`: Tolerance
"""
function make_unitcell(latticeconstant ::Real;
                       SiteType::DataType=Any,
                       tol::Real=Base.rtoldefault(Float64))
    return make_unitcell(reshape([latticeconstant], (1,1));
                         SiteType=SiteType, tol=tol)
end


"""
    UnitCell

Construct an n-dimensional lattice.

# Arguments
* `latticevectors ::AbstractArray{<:AbstractFloat, 2}`: Lattice vectors
* `SiteType::DataType`

# Optional Arguments
* `tol=√ϵ`: Epsilon
"""
function make_unitcell(latticevectors ::AbstractArray{<:Real, 2};
                       SiteType::DataType=Any,
                       tol ::Real=Base.rtoldefault(Float64))
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
    return UnitCell{SiteType}(latticevectors, sites,
                                 reduced_rlv, 2*π*reduced_rlv, siteindices)
end


function make_unitcell(latticevectors::AbstractVector{<:AbstractVector};
                       SiteType::DataType=Any,
                       tol::Real=Base.rtoldefault(Float64))
    lv = hcat(latticevectors...)
    return make_unitcell(lv; SiteType=SiteType, tol=tol)
end


"""
    dimension

Spatial dimension of the unit cell.
"""
function dimension(uc ::UnitCell)
    return size(uc.latticevectors, 1)
end


"""
    numsite

Number of sites of the unit cell.

# Arguments
* `uc ::UnitCell`
"""
function numsite(uc ::UnitCell)
    return length(uc.sites)
end


"""
    addsite!

Add an site to the unit cell.

# Arguments
* `uc ::UnitCell{T}`
* `sitename ::{T}`
* `sitecoord ::FractCoord`
"""
function addsite!(uc ::UnitCell{O},
                     sitename ::O,
                     sitecoord ::FractCoord) where {O}
    (ndim, ndim_) = size(uc.latticevectors)
    if dimension(sitecoord) != ndim
        throw(ArgumentError("sitecoord has wrong dimension"))
    elseif haskey(uc.siteindices, sitename)
        throw(ArgumentError( "duplicate site name"))
    end
    push!(uc.sites, (sitename, sitecoord))
    index = length(uc.siteindices)+1
    uc.siteindices[sitename] = index
    return index
end


"""
    hassite{T}

Test whether the unit cell contains the site of given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`
"""
function hassite(uc ::UnitCell{O}, name ::O) where {O}
    return haskey(uc.siteindices, name)
end


"""
    getsiteindex

Get index of the given site.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`
"""
function getsiteindex(uc ::UnitCell{O}, name ::O) where {O}
    return uc.siteindices[name]
end


"""
    getsite

Get the site (its site name and its fractional coordinates) with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `(sitename, fractcoord)`
"""
function getsite(uc ::UnitCell{O}, name ::O) where {O}
    return uc.sites[ uc.siteindices[name] ]
end


"""
    getsitecoord

Get the fractional coordinates of the site with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `fractcoord`
"""
function getsitecoord(uc ::UnitCell{O}, name ::O) where {O}
    return getsite(uc, name)[2]
end


"""
    getsiteindexcoord

# Arguments
* `uc ::UnitCell{T}`
* `name ::T`

# Return
* `(index, fractcoord)`
"""
function getsiteindexcoord(uc ::UnitCell{O}, name::O) where {O}
    index = getsiteindex(uc, name)
    coord = getsitecoord(uc, index)
    return (index, coord)
end


"""
    getsite

# Arguments
* `uc ::UnitCell{T}`
* `index ::Integer`

# Return
* `(sitename, fractcoord)`
"""
function getsite(uc ::UnitCell, index::Integer)
    return uc.sites[index]
end


"""
    getsitename

# Arguments
* `uc ::UnitCell`
* `index ::Integer`

# Return
* `sitename`
"""
function getsitename(uc ::UnitCell, index ::Integer)
    return uc.sites[index][1]
end


"""
    getsitecoord

# Arguments
* `uc ::UnitCell`
* `idx ::Integer`

# Return
* `fractcoord`
"""
function getsitecoord(uc ::UnitCell, index ::Integer)
    return uc.sites[index][2]
end



"""
    fract2carte

# Arguments
* `latticevectors ::Array{Float64, 2}`
* `fc ::FractCoord`
"""
function fract2carte(unitcell ::UnitCell, fc ::FractCoord)
    if dimension(unitcell) != dimension(fc)
        throw(ArgumentError("unitcell and fractcoord must have the same dimension"))
    end
    mc = fc.whole + fc.fraction
    cc = unitcell.latticevectors * mc
    return CarteCoord(cc)
end


"""
    carte2fract

# Arguments
* `latticevectors ::Array{Float64, 2}`
* `cc ::CarteCoord`
"""
function carte2fract(unitcell::UnitCell, cc::CarteCoord;
                     tol::Real=Base.rtoldefault(Float64))
    if dimension(unitcell) != length(cc)
        throw(ArgumentError("unitcell and cartecoord must have the same dimension"))
    end
    fc = transpose(unitcell.reducedreciprocallatticevectors) * cc
    w = Int[fld(x, 1) for x in fc]
    f = Float64[mod(x, 1) for x in fc]
    return FractCoord(w, f; tol=tol)
end


"""
    whichunitcell

# Return
* `R ::Vector{Int}`: which unit cell the specificied site/cartesian coordinates belongs to.
"""
function whichunitcell(uc::UnitCell{O}, name::O, cc::CarteCoord;
                       tol::Real=Base.rtoldefault(Float64)) where {O}
    fc1 = getsitecoord(uc, name)
    fc2 = carte2fract(uc, cc)
    if !isapprox(fc1.fraction, fc2.fraction; atol=tol)
        throw(ArgumentError("$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]"))
    end
    R = fc2.whole - fc1.whole
    return R
end


"""
    whichunitcell

# Return
* `R ::Vector{Int}`: which unit cell the specificied site/cartesian coordinates belongs to.
"""
function whichunitcell(uc::UnitCell{O}, name::O, fc::FractCoord;
                       tol::Real=Base.rtoldefault(Float64)) where {O}
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
    momentumgrid

Generate an n-dimensional grid of momenta of given shape
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


import Base.==
function (==)(lhs::UnitCell{O}, rhs::UnitCell{O}) where O
  return (
    lhs.latticevectors == rhs.latticevectors &&
    lhs.sites == rhs.sites
  )
end


"""
    findsiteindex

Returns (site_index, unitcell_vector), or `(-1, [])` if not found.
"""
function findsiteindex(unitcell::UnitCell, fc::FractCoord; tol=Base.rtoldefault(Float64))
    i = findfirst(x -> isapprox(x, fc.fraction; atol=tol), [orbfc.fraction for (orbname, orbfc) in unitcell.sites])
    if i !== nothing
        (orbname, orbfc) = unitcell.sites[i]
        return (i, fc.whole - orbfc.whole)
    else
        return (-1, Int[])
    end
end
