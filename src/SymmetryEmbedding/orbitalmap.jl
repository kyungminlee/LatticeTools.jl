export findorbitalmap


"""
    findorbitalmap(unitcell, translation_operation)

Find which orbital gets mapped to which orbital with what lattice displacement.
Return a vector of tuples `(j, Rj)` at index `i`, representing that the orbital `i`
gets mapped to orbital `j` with lattice displacement `Rj`.

The overall transformation will be a combiation of (1) the result here, and
(2) Bravais lattice transformation. Since all the integer translations will be
taken care of by the Bravais lattice transformation, Rj here is all zero vectors.
"""
function findorbitalmap(unitcell::UnitCell,
                        tsym_op::TranslationOperation{<:Integer})::Vector{Tuple{Int, Vector{Int}}}
    norb = numorbital(unitcell)
    dim = dimension(unitcell)
    return [(i, zeros(Int, dim)) for i in 1:norb]
end


"""
    findorbitalmap(unitcell, translation_symmetry)
"""
function findorbitalmap(unitcell::UnitCell, tsym::TranslationSymmetry)
    return [findorbitalmap(unitcell, m) for m in elements(tsym)]
end


"""
    findorbitalmap(unitcell, point_operation)

Find which orbital gets mapped to which orbital with what lattice displacement.
Return a vector of tuples `(j, Rj)` at index `i`, representing that the orbital `i`
gets mapped to orbital `j` with lattice displacement `Rj`.
"""
function findorbitalmap(unitcell::UnitCell,
                        psym_op::PointOperation{<:Integer})::Union{Nothing, Vector{Tuple{Int, Vector{Int}}}}
    norb = numorbital(unitcell)
    map = Tuple{Int, Vector{Int}}[]
    for (orbname, orbfc) in unitcell.orbitals
        j, Rj = findorbitalindex(unitcell, psym_op.matrix * orbfc)
        j <= 0 && return nothing      # throw(ArgumentError("orbital map not found with $unitcell and $psym_op"))
        push!(map, (j, Rj))
    end
    return map
end


"""
    findorbitalmap(unitcell, point_symmetry)
"""
function findorbitalmap(unitcell::UnitCell, psym::PointSymmetry)::Union{Nothing, Vector{Vector{Tuple{Int, Vector{Int}}}}}
    out = Vector{Tuple{Int, Vector{Int}}}[]
    sizehint!(out, length(elements(psym)))
    for el in elements(psym)
        m = findorbitalmap(unitcell, el)
        isnothing(m) && return nothing
        push!(out, m)
    end
    return out
end