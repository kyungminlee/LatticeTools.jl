export SymmorphicIrrepComponent
export get_irrep_components
export get_irrep_iterator


"""
    SymmorphicIrrepComponent{S1<:SymmetryOrEmbedding, S2<:SymmetryOrEmbedding}

Irrep component of a symmorphic symmetry (embedding).

# Fields
* `normal::IrrepComponent{S1}`
* `rest::IrrepComponent{S2}`
"""
struct SymmorphicIrrepComponent{S1<:SymmetryOrEmbedding, S2<:SymmetryOrEmbedding}<:AbstractSymmetryIrrepComponent
    normal::IrrepComponent{S1}
    rest::IrrepComponent{S2}

    @doc """
        SymmorphicIrrepComponent(normal, rest)

    Construct a `SymmorphicIrrepComponent` of `normal` and `rest`.

    # Arguments
    * `normal::IrrepComponent{S1}`
    * `rest::IrrepComponent{S2}`
    )
    """
    function SymmorphicIrrepComponent(
        normal::IrrepComponent{S1},
        rest::IrrepComponent{S2},
    ) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        if !iscompatible(normal.symmetry, normal.irrep_index, rest.symmetry)
            throw(ArgumentError(
                "symmetry $(symmetry_name(rest.symmetry)) is not compatible with "*
                "symmetry $(symmetry_name(normal.symmetry)) at irrep $(normal.irrep_index)"
            ))
        end
        return new{S1, S2}(normal, rest)
    end

    """
        SymmorphicIrrepComponent(normal, rest)

    # Arguments
    * `normal::IrrepComponent{SymmetryEmbedding{S1}}`
    * `rest::IrrepComponent{SymmetryEmbedding{S2}}`
    """
    function SymmorphicIrrepComponent(
        normal::IrrepComponent{SymmetryEmbedding{S1}},
        rest::IrrepComponent{SymmetryEmbedding{S2}},
    ) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        if !iscompatible(normal.symmetry, normal.irrep_index, rest.symmetry)
            throw(ArgumentError(
                "symmetry $(symmetry_name(rest.symmetry)) is not compatible with "*
                "symmetry $(symmetry_name(normal.symmetry)) at irrep $(normal.irrep_index)"
            ))
        end
        return new{SymmetryEmbedding{S1}, SymmetryEmbedding{S2}}(normal, rest)
    end
end

"""
    group_order(arg::SymmorphicIrrepComponent)

Group order of the underlying symmorphic group of `arg`.
"""
GroupTools.group_order(arg::SymmorphicIrrepComponent) = group_order(arg.normal) * group_order(arg.rest)


"""
    get_irrep_components(sym::SymmorphicSymmetry)

Return an iterator which iterates over the irrep component of the symmorphic symmetry `sym`.
"""
function get_irrep_components(sym::SymmorphicSymmetry)
    return (
        SymmorphicIrrepComponent(normal_sic, rest_sic)
        for normal_sic in get_irrep_components(sym.normal)
        for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest))
    )
end


"""
    get_irrep_components(sym::SymmorphicSymmetryEmbedding)

Return an iterator which iterates over the irrep component of the symmorphic symmetry
embedding `sym`.
"""
function get_irrep_components(sym::SymmorphicSymmetryEmbedding)
    return (
        SymmorphicIrrepComponent(normal_sic, rest_sic)
        for normal_sic in get_irrep_components(sym.normal)
        for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest))
    )
end


"""
    get_irrep_iterator(ssic::SymmorphicIrrepComponent)

Return an iterator which iterates over the elements, together with their corresponding irrep
coefficient of the given irrep component.

# Return
* [(r⋅n, cᵣ⋅cₙ) for all r for all n]
"""
function get_irrep_iterator(ssic::SymmorphicIrrepComponent)
    normal_elements, normal_irrep_components = let
        nsym = ssic.normal.symmetry
        nii  = ssic.normal.irrep_index
        nic  = ssic.normal.irrep_component
        elements(nsym), [m[nic, nic] for m in irrep(nsym, nii)]
    end
    rest_elements, rest_irrep_components = let
        rsym = ssic.rest.symmetry
        rii  = ssic.rest.irrep_index
        ric  = ssic.rest.irrep_component
        elements(rsym), [m[ric, ric] for m in irrep(rsym, rii)]
    end
    return (
        (relem * nelem, rphase * nphase)
        for (relem, rphase) in zip(rest_elements, rest_irrep_components)
        for (nelem, nphase) in zip(normal_elements, normal_irrep_components)
    )
end
