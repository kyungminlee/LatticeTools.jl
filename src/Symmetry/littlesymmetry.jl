export little_group_elements
export little_group
export little_symmetry
export little_symmetry_iso


## Symmetry reduction (little group etc.)

"""
    psym compatible with hypercube
"""
function little_group_elements(tsym::TranslationSymmetry, psym::PointSymmetry)
    lg_elements = [i for (i, elem) in enumerate(elements(psym))
                     if iscompatible(tsym, elem)]
    return lg_elements
end


function little_group_elements(tsym::TranslationSymmetry,
                               tsym_irrep_index::Integer,
                               psym::PointSymmetry)
    k1 = tsym.fractional_momenta[tsym_irrep_index]
    lg = Int[]
    for (i_elem, matrep) in enumerate(psym.matrix_representations)
        k2 = (x->mod(x,1)).(ExactLinearAlgebra.inverse(transpose(matrep)) * k1)
        k2 == k1 && push!(lg, i_elem)
    end
    return lg
end


function little_group(tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      psym::PointSymmetry)
    lg_elems = little_group_elements(tsym, tsym_irrep_index, psym)
    return little_group(tsym, psym, lg_elems)
end


"""
Generate a little group with given elements.
The elements of the little group, which may be sparse, are compressed into consecutive integers.
"""
function little_group(tsym::TranslationSymmetry,
                      psym::PointSymmetry,
                      elements::AbstractVector{<:Integer})
    element_lookup = Dict(x=>i for (i, x) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[ group_product(psym.group, elements[i], elements[j]) ]
    end
    return FiniteGroup(mtab)
end


function little_symmetry(tsym::TranslationSymmetry, psym::PointSymmetry)
    (lg_raw, lg_matrep_raw, lg_element_names_raw) = let
        lg_elements = little_group_elements(tsym, psym)
        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]
        (lg_raw, lg_matrep_raw, lg_element_names_raw)
    end

    group_index = PointSymmetryDatabase.find(lg_element_names_raw)
    psym2 = PointSymmetryDatabase.get(group_index)
    ϕ = group_isomorphism(lg_raw, psym2.group)
    isnothing(ϕ) && error("Group not isomorphic")

    lg_matrep = lg_matrep_raw[ϕ]
    lg_element_names = lg_element_names_raw[ϕ]

    PointSymmetry(psym2.group,
                  psym2.generators,
                  psym2.conjugacy_classes,
                  psym2.character_table,
                  psym2.irreps,
                  lg_element_names,
                  lg_matrep,
                  psym2.hermann_mauguinn,
                  psym2.schoenflies)
end


function little_symmetry(tsym::TranslationSymmetry, tsym_irrep::Integer, psym::PointSymmetry)
    tsym_irrep == 1 && return psym
    (lg_raw, lg_matrep_raw, lg_element_names_raw) = let
        lg_elements = little_group_elements(tsym, tsym_irrep, psym)
        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]
        (lg_raw, lg_matrep_raw, lg_element_names_raw)
    end

    group_index = PointSymmetryDatabase.find(lg_element_names_raw)
    psym2 = PointSymmetryDatabase.get(group_index)
    ϕ = group_isomorphism(lg_raw, psym2.group)
    isnothing(ϕ) && error("Group not isomorphic")

    lg_matrep = lg_matrep_raw[ϕ]
    lg_element_names = lg_element_names_raw[ϕ]

    PointSymmetry(psym2.group,
                  psym2.generators,
                  psym2.conjugacy_classes,
                  psym2.character_table,
                  psym2.irreps,
                  lg_element_names,
                  lg_matrep,
                  psym2.hermann_mauguinn,
                  psym2.schoenflies)
end


"""
    little_symmetry_iso(tsym, tsym_irrep_index, psym)

Find little symmetry using group isomorphism
"""
function little_symmetry_iso(tsym::TranslationSymmetry, tsym_irrep_index::Integer, psym::PointSymmetry)
    tsym_irrep_index == 1 && return psym
    (lg_irrep, lg_matrep, lg_element_names) = let
        lg_elements = little_group_elements(tsym, tsym_irrep_index, psym)

        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]

        (lg_irrep, ϕ) = IrrepDatabase.find(lg_raw)

        lg_matrep = lg_matrep_raw[ϕ]
        lg_element_names = lg_element_names_raw[ϕ]
        (lg_irrep, lg_matrep, lg_element_names)
    end

    generators = minimal_generating_set(lg_irrep.group)
    hermann_mauguinn = join(lg_element_names[generators])
    schoenflies = "unknown"

    simple_element_names = sort(simplify_name.(lg_element_names))
    for i in 1:32
        psym = PointSymmetryDatabase.get(i)
        if ( sort(simplify_name.(psym.element_names)) == simple_element_names )
             #&& !isnothing(group_isomorphism(lg_irrep.group, psym.group)) )
            hermann_mauguinn = psym.hermann_mauguinn
            break
        end
    end

    PointSymmetry(lg_irrep.group,
                  generators,
                  lg_irrep.conjugacy_classes,
                  lg_irrep.character_table,
                  lg_irrep.irreps,
                  lg_element_names,
                  lg_matrep,
                  hermann_mauguinn,
                  schoenflies)
end
