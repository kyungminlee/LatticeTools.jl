export UnitCell
export make_unitcell
export dimension

export numsite,
       addsite!,
       hassite,
       getsite,
       getsiteindex,
       getsitecoord,
       getsiteindexcoord,
       getsitename

export numorbital,
       hasorbital,
       addorbital!,
       getorbital,
       getorbitalindex,
       getorbitalcoord,
       getorbitalindexcoord,
       getorbitalname

export carte2fract,
       fract2carte,
       whichunitcell,
       momentumgrid

export findsiteindex

import LinearAlgebra


"""
    UnitCell{O}

# Parameters
* `O`: type of "orbital". Any type can be used, but we recommend using `String` or tuple of `String` and `Int`
       for compatibility with JSON.

# Members
* `latticevectors ::Array{Float64, 2}`: Lattice vectors
* `reducedreciprocallatticevectors ::Array{Float64, 2}`: Reduced reciprocal lattice vectors (transpose of inverse of `latticevectors`)
* `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors. `2π * reducedreciprocallatticevectors`
* `orbitals ::Vector{Tuple{T, FractCoord}}`: List of orbitals within unit cell
* `orbitalindices ::Dict{T, Int}`: Indices of orbitals
"""
mutable struct UnitCell{S, O}
    latticevectors ::Array{Float64, 2}
    sites::Vector{Tuple{O, FractCoord}}
    orbitals ::Vector{Tuple{O, Int}}

    reducedreciprocallatticevectors ::Array{Float64, 2}
    reciprocallatticevectors ::Array{Float64, 2}
    siteindices::Dict{S, Int}
    orbitalindices::Dict{O, Int}
    function UnitCell{S, O}(latticevectors::AbstractArray{<:Real, 2},
                            sites::AbstractVector{Tuple{S, FractCoord}},
                            orbitals::AbstractVector{<:Tuple{O, <:Integer}},
                            reducedreciprocallatticevectors::AbstractArray{<:Real, 2},
                            reciprocallatticevectors::AbstractArray{<:Real, 2},
                            siteindices::AbstractDict{S, <:Integer},
                            orbitalindices::AbstractDict{O, <:Integer}) where {S, O}
        if (S <: Integer)
            throw(ArgumentError("SiteType should not be integer to avoid confusion"))
        end                    
        if (O <: Integer)
            throw(ArgumentError("OrbitalType should not be integer to avoid confusion"))
        end
        new{S, O}(latticevectors,
                  sites,
                  orbitals,
                  reducedreciprocallatticevectors,
                  reciprocallatticevectors,
                  siteindices,
                  orbitalindices)
    end
end


"""
    UnitCell

Construct a one-dimensional lattice.

# Arguments
* `latticeconstant ::Float64`: Lattice constant
* `OrbitalType`: List of orbitals

# Optional Arguments
* `tol=√ϵ`: Tolerance
"""
function make_unitcell(latticeconstant::Real;
                       sitetype::Type{SiteType}=Any,
                       orbitaltype::Type{OrbitalType}=Any,
                       tol::Real=Base.rtoldefault(Float64)) where {SiteType, OrbitalType}
    return make_unitcell(reshape([latticeconstant], (1,1));
                         sitetype=sitetype,
                         orbitaltype=orbitaltype,
                         tol=tol)
end


"""
    UnitCell

Construct an n-dimensional lattice.

# Arguments
* `latticevectors ::AbstractArray{<:AbstractFloat, 2}`: Lattice vectors
* `OrbitalType::DataType`

# Optional Arguments
* `tol=√ϵ`: Epsilon
"""
function make_unitcell(latticevectors::AbstractArray{<:Real, 2};
                       sitetype::Type{SiteType}=Any,
                       orbitaltype::Type{OrbitalType}=Any,
                       tol::Real=Base.rtoldefault(Float64)) where {SiteType, OrbitalType}
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

    orbitals = Tuple{OrbitalType, Int}[]
    orbitalindices = Dict{OrbitalType, Int}()
    return UnitCell{SiteType, OrbitalType}(
                latticevectors,
                sites, orbitals,
                reduced_rlv, 2*π*reduced_rlv,
                siteindices, orbitalindices)
end


function make_unitcell(latticevectors::AbstractVector{<:AbstractVector};
                       sitetype::Type{SiteType}=Any,
                       orbitaltype::Type{OrbitalType}=Any,
                       tol::Real=Base.rtoldefault(Float64)) where {SiteType, OrbitalType}
    lv = hcat(latticevectors...)
    return make_unitcell(lv; sitetype=sitetype, orbitaltype=orbitaltype, tol=tol)
end


"""
    dimension

Spatial dimension of the unit cell.
"""
function dimension(uc::UnitCell)
    return size(uc.latticevectors, 1)
end


function numsite(uc::UnitCell)
    return length(uc.sites)
end


"""
    numorbital

Number of orbitals of the unit cell.

# Arguments
* `uc ::UnitCell`
"""
function numorbital(uc::UnitCell)
    return length(uc.orbitals)
end


function addsite!(uc::UnitCell{S, O},
                 sitename::S,
                 sitecoord::FractCoord;
                 tol::Real=Base.rtoldefault(Float64)) where {S, O}
    (ndim, ndim_) = size(uc.latticevectors)

    if dimension(sitecoord) != ndim
        throw(ArgumentError("sitecoord has wrong dimension"))
    end

    for (n, c) in uc.sites
        n == sitename && throw(ArgumentError("duplicate site name"))
        isapprox(c, sitecoord; atol=tol) && throw(ArgumentError("sites with the same coordinates"))
    end
    push!(uc.sites, (sitename, sitecoord))
    index = length(uc.siteindices)+1
    uc.siteindices[sitename] = index
    return index
end


"""
    addorbital!

Add an orbital to the unit cell.

# Arguments
* `uc ::UnitCell{T}`
* `orbitalname ::{T}`
* `orbitalcoord ::FractCoord`
"""
function addorbital!(uc ::UnitCell{S, O},
                     sitename::S,
                     orbitalname::O,
                     ) where {S, O}
    haskey(uc.orbitalindices, orbitalname) && throw(ArgumentError("duplicate orbital name"))
    siteindex = get(uc.siteindices, sitename, 0)
    siteindex <= 0 && throw(ArgumentError("site $sitename not found"))
    push!(uc.orbitals, (orbitalname, siteindex))
    index = length(uc.orbitalindices)+1
    uc.orbitalindices[orbitalname] = index
    return index
end


"""
    hassite

Check whether the unit cell contains the site of given name.

# Arguments
* `uc::UnitCell{S, O}`
* `name::S`
"""
function hassite(uc ::UnitCell{S, O}, name::S) where {S, O}
    return haskey(uc.siteindices, name)
end


"""
    hasorbital

Check whether the unit cell contains the orbital of given name.

# Arguments
* `uc::UnitCell{S, O}`
* `name::O`
"""
function hasorbital(uc::UnitCell{S, O}, name::O) where {S, O}
    return haskey(uc.orbitalindices, name)
end


"""
    getsiteindex

Get index of the given site.

# Arguments
* `uc::UnitCell{S, O}`
* `name::S`
"""
function getsiteindex(uc::UnitCell{S, O}, name::S) where {S, O}
    return uc.siteindices[name]
end


"""
    getorbitalindex

Get index of the given orbital.

# Arguments
* `uc::UnitCell{O}`
* `name::O`
"""
function getorbitalindex(uc::UnitCell{S, O}, name::O) where {S, O}
    return uc.orbitalindices[name]
end


"""
    getsite

Get the site (its name and its fractional coordinates) with the given name.

# Arguments
* `uc::UnitCell{S, O}`
* `name::S`

# Return
* `(sitename, fractcoord)`
"""
function getsite(uc::UnitCell{S, O}, name::S) where {S, O}
    return uc.sites[ uc.siteindices[name] ]
end


"""
    getorbital

Get the orbital (its orbital name and its siteindex) with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `(orbitalname, siteindex)`
"""
function getorbital(uc ::UnitCell{S, O}, name::O) where {S, O}
    return uc.orbitals[ uc.orbitalindices[name] ]
end


"""
    getsitecoord

Get the fractional coordinate of the site with the given name
"""
function getsitecoord(uc::UnitCell{S, O}, name::S) where {S, O}
    return getsite(uc, name)[2]
end


"""
    getorbitalcoord

Get the fractional coordinates of the orbital with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `fractcoord`
"""
function getorbitalcoord(uc::UnitCell{S, O}, name::O) where {S, O}
    return uc.sites[ getorbitalsiteindex(uc, name) ][2]
end

function getsiteindexcoord(uc::UnitCell{S, O}, name::S) where {S, O}
    index = getsiteindex(uc, name)
    coord = getsitecoord(uc, index)
    return (index, coord)
end


"""
    getorbitalindexcoord

# Arguments
* `uc ::UnitCell{T}`
* `name ::T`

# Return
* `(index, fractcoord)`
"""
function getorbitalindexcoord(uc::UnitCell{S, O}, name::O) where {S, O}
    index = getorbitalindex(uc, name)
    coord = getorbitalcoord(uc, index)
    return (index, coord)
end


function getsite(uc::UnitCell, index::Integer)
    return uc.sites[index]
end

function getsitename(uc ::UnitCell, index ::Integer)
    return uc.sites[index][1]
end

function getsitecoord(uc ::UnitCell, index ::Integer)
    return uc.sites[index][2]
end


"""
    getorbital

# Arguments
* `uc ::UnitCell{T}`
* `index ::Integer`

# Return
* `(orbitalname, siteindex)`
"""
function getorbital(uc::UnitCell, index::Integer)
    return uc.orbitals[index]
end


"""
    getorbitalname

# Arguments
* `uc ::UnitCell`
* `index ::Integer`

# Return
* `orbitalname`
"""
function getorbitalname(uc::UnitCell, index::Integer)
    return uc.orbitals[index][1]
end


function getorbitalsiteindex(uc::UnitCell, index::Integer)
    return uc.orbitals[index][2]
end


"""
    getorbitalcoord

# Arguments
* `uc ::UnitCell`
* `idx ::Integer`

# Return
* `fractcoord`
"""
function getorbitalcoord(uc::UnitCell, index::Integer)
    return uc.sites[getorbitalsiteindex(uc, index)][2]
end


"""
    fract2carte

# Arguments
* `latticevectors ::Array{Float64, 2}`
* `fc ::FractCoord`
"""
function fract2carte(unitcell::UnitCell, fc::FractCoord)
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
* `R ::Vector{Int}`: which unit cell the specificied orbital/cartesian coordinates belongs to.
"""
function whichunitcell(uc::UnitCell{S, O}, name::S, cc::CarteCoord;
                       tol::Real=Base.rtoldefault(Float64)) where {S, O}
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
* `R ::Vector{Int}`: which unit cell the specificied orbital/cartesian coordinates belongs to.
"""
function whichunitcell(uc::UnitCell{S, O}, name::S, fc::FractCoord;
                       tol::Real=Base.rtoldefault(Float64)) where {S, O}
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
    norb = numorbital(uc)
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
function (==)(lhs::UnitCell{S, O}, rhs::UnitCell{S, O}) where {S, O}
    return (lhs.latticevectors == rhs.latticevectors &&
            lhs.sites == rhs.sites &&
            lhs.orbitals == rhs.orbitals)
end



"""
    findorbitalindex

Returns (orbital_index, unitcell_vector), or `(-1, [])` if not found.
"""
function findsiteindex(unitcell::UnitCell,
                       fc::FractCoord;
                       tol::Real=Base.rtoldefault(Float64))
    i = findfirst(x -> isapprox(x, fc.fraction; atol=tol),
                  [sitefc.fraction for (sitename, sitefc) in unitcell.sites])
    if i !== nothing
        (sitename, sitefc) = unitcell.sites[i]
        return (i, fc.whole - sitefc.whole)
    else
        return (-1, Int[])
    end
end
