export findsitemap
export findorbitalmap

@deprecate findorbitalmap(args...) findsitemap(args...)


"""
    findsitemap(unitcell, translation_operation)

Find which site gets mapped to which site with what lattice displacement.
Return a vector of tuples `(j, Rj)` at index `i`, representing that the site `i`
gets mapped to site `j` with lattice displacement `Rj`.

Since periodic boundary condition is always assumed, the return value always exists.

This function is used to find the mapping between sites under translation symmetry in a
given lattice. The overall transformation will be a combiation of (1) the result here, and
(2) Bravais lattice transformation. Since all the integer translations will be
taken care of by the Bravais lattice transformation, `Rj` here is all zero vectors.
"""
function findsitemap(
    unitcell::UnitCell,
    tsym_op::TranslationOperation{<:Integer},
)::Vector{Tuple{Int, Vector{Int}}}
    norb = numsite(unitcell)
    dim = dimension(unitcell)
    return [(i, zeros(Int, dim)) for i in 1:norb]
end


"""
    findsitemap(unitcell, translation_symmetry)

Find the sitemap for every element of the translation symmetry.
Return `[findsitemap(unitcell, m) for m in elements(translation_symmetry)]`.
"""
function findsitemap(unitcell::UnitCell, tsym::TranslationSymmetry)
    return [findsitemap(unitcell, m) for m in elements(tsym)]
end


"""
    findsitemap(unitcell, point_operation)

Find which site gets mapped to which site with what lattice displacement.
Return a vector of tuples `(j, Rj)` at index `i`, representing that the site `i`
gets mapped to site `j` with lattice displacement `Rj`.
If the unit cell is not compatible with the point operation, return `nothing`.
"""
function findsitemap(
    unitcell::UnitCell,
    psym_op::PointOperation{<:Integer},
)::Union{Nothing, Vector{Tuple{Int, Vector{Int}}}}
    norb = numsite(unitcell)
    map = Tuple{Int, Vector{Int}}[]
    for (orbname, orbfc) in unitcell.sites
        j, Rj = findsiteindex(unitcell, psym_op.matrix * orbfc)
        j <= 0 && return nothing  # throw(ArgumentError("site map not found with $unitcell and $psym_op"))
        push!(map, (j, Rj))
    end
    return map
end


"""
    findsitemap(unitcell, point_symmetry)

Find the sitemap for every element of the point symmetry.
Return `nothing` if the unit cell is incompatible with any element of the point symmetry.
"""
function findsitemap(
    unitcell::UnitCell,
    psym::PointSymmetry,
)::Union{Nothing, Vector{Vector{Tuple{Int, Vector{Int}}}}}
    out = Vector{Tuple{Int, Vector{Int}}}[]
    sizehint!(out, length(elements(psym)))
    for el in elements(psym)
        m = findsitemap(unitcell, el)
        isnothing(m) && return nothing
        push!(out, m)
    end
    return out
end
