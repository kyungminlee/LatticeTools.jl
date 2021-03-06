export little_group_elements
export little_group
export little_symmetry
export little_symmetry_iso


"""
    little_group_elements(tsym, psym)

Return elements of `psym` compatible with the translation symmetry `tsym`.
"""
function little_group_elements(tsym::FiniteTranslationSymmetry, psym::PointSymmetry)::Vector{Int}
    return Int[
        i_elem
            for (i_elem, elem) in enumerate(elements(psym))
                if iscompatible(tsym, elem)
    ]
end


"""
    little_group_elements(tsym, tsym_irrep_index, psym)

Return elements of `psym` compatible with the translation symmetry `tsym` at
irrep `tsym_irrep_index`
"""
function little_group_elements(
    tsym::FiniteTranslationSymmetry,
    tsym_irrep_index::Integer,
    psym::PointSymmetry,
)::Vector{Int}
    return Int[
        i_elem
            for (i_elem, pop) in enumerate(elements(psym))
                if iscompatible(tsym, tsym_irrep_index, pop)
    ]
end


"""
    little_group(tsym, tsym_irrep_index, psym)
"""
function little_group(
    tsym::FiniteTranslationSymmetry,
    tsym_irrep_index::Integer,
    psym::PointSymmetry,
)::FiniteGroup
    lg_elems = little_group_elements(tsym, tsym_irrep_index, psym)
    return little_group(tsym, psym, lg_elems)
end


"""
    little_group(tsym, psym, elements)

Generate a little group with given elements.
The elements of the little group, which may be sparse, are compressed into consecutive integers.
"""
function little_group(
    tsym::FiniteTranslationSymmetry,
    psym::PointSymmetry,
    elements::AbstractVector{<:Integer},
)
    element_lookup = Dict(x=>i for (i, x) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[group_product(psym.group, elements[i], elements[j])]
    end
    return FiniteGroup(mtab)
end


"""
    little_symmetry(tsym, psym, lg_elements)

Generate a little symmetry.
First generate a little group with `lg_elements`, look up `PointSymmetryDatabase` using
the names of the elements, and then checks for group isomorphism.

# Arguments
* `tsym::FiniteTranslationSymmetry`: translation symmetry
* `psym::PointSymmetry`: point symmetry
* `lg_elements::AbstractVector{<:Integer}`: indices of the elements of the little group
"""
function little_symmetry(
    tsym::FiniteTranslationSymmetry,
    psym::PointSymmetry,
    lg_elements::AbstractVector{<:Integer},
)
    (lg_raw, lg_matrep_raw, lg_element_names_raw) = let
        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]
        (lg_raw, lg_matrep_raw, lg_element_names_raw)
    end

    psym2 = PointSymmetryDatabase.find(lg_element_names_raw)
    ϕ = group_isomorphism(lg_raw, psym2.group)
    isnothing(ϕ) && error("Group not isomorphic")

    lg_matrep = lg_matrep_raw[ϕ]
    lg_element_names = lg_element_names_raw[ϕ]

    return PointSymmetry(
        psym2.group,
        psym2.generators,
        psym2.conjugacy_classes,
        psym2.character_table,
        psym2.irreps,
        lg_element_names,
        lg_matrep,
        psym2.hermann_mauguin,
        psym2.schoenflies,
    )
end


"""
    little_symmetry(tsym, psym)
"""
function little_symmetry(tsym::FiniteTranslationSymmetry, psym::PointSymmetry)
    return little_symmetry(tsym, psym, little_group_elements(tsym, psym))
end


"""
    little_symmetry(tsym, tsym_irrep_index, psym)
"""
function little_symmetry(tsym::FiniteTranslationSymmetry, tsym_irrep::Integer, psym::PointSymmetry)
    if !iscompatible(tsym, psym)
        throw(ArgumentError("translation and point symmetries are not compatible"))
    end

    tsym_irrep == 1 && return little_symmetry(tsym, psym)

    (lg_raw, lg_matrep_raw, lg_element_names_raw) = let
        lg_elements = little_group_elements(tsym, tsym_irrep, psym)
        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]
        (lg_raw, lg_matrep_raw, lg_element_names_raw)
    end

    psym2 = PointSymmetryDatabase.find(lg_element_names_raw)
    ϕ = group_isomorphism(lg_raw, psym2.group)
    isnothing(ϕ) && error("Group not isomorphic")

    lg_matrep = lg_matrep_raw[ϕ]
    lg_element_names = lg_element_names_raw[ϕ]

    return PointSymmetry(
        psym2.group,
        psym2.generators,
        psym2.conjugacy_classes,
        psym2.character_table,
        psym2.irreps,
        lg_element_names,
        lg_matrep,
        psym2.hermann_mauguin,
        psym2.schoenflies,
    )
end


# COV_EXCL_START

"""
    little_symmetry_iso(tsym, tsym_irrep_index, psym)

Find little symmetry using group isomorphism
"""
function little_symmetry_iso(
    tsym::FiniteTranslationSymmetry,
    tsym_irrep_index::Integer,
    psym::PointSymmetry,
)
    @warn "Function little_symmetry_iso deprecated. Correctness not guaranteed." maxlog=1
    if !iscompatible(tsym, psym)
        throw(ArgumentError("translation and point symmetries are not compatible"))
    end
    tsym_irrep_index == 1 && return little_symmetry(tsym, psym)
    (lg_matrep, lg_element_names, lg_irrep) = let
        lg_elements = little_group_elements(tsym, tsym_irrep_index, psym)

        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]

        (lg_irrep, ϕ) = IrrepDatabase.find(lg_raw)

        lg_matrep = lg_matrep_raw[ϕ]
        lg_element_names = lg_element_names_raw[ϕ]
        (lg_matrep, lg_element_names, lg_irrep)
    end

    generators = minimal_generating_set(lg_irrep.group)
    hermann_mauguin = join(lg_element_names[generators])
    schoenflies = "unknown"

    simple_element_names = sort(simplify_name.(lg_element_names))
    for i in 1:32
        psym = PointSymmetryDatabase.get(i)
        if ( sort(simplify_name.(psym.element_names)) == simple_element_names
             && !isnothing(group_isomorphism(lg_irrep.group, psym.group)) )
            hermann_mauguin = psym.hermann_mauguin
            break
        end
    end

    return PointSymmetry(
        lg_irrep.group,
        generators,
        lg_irrep.conjugacy_classes,
        lg_irrep.character_table,
        lg_irrep.irreps,
        lg_element_names,
        lg_matrep,
        hermann_mauguin,
        schoenflies,
    )
end
# COV_EXCL_STOP