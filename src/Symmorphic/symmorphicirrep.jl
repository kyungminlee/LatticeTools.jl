export SymmorphicIrrepComponent
export get_irrep_components
export get_irrep_iterator

struct SymmorphicIrrepComponent{S1<:SymmetryOrEmbedding, S2<:SymmetryOrEmbedding}<:AbstractSymmetryIrrepComponent
    normal::IrrepComponent{S1}
    rest::IrrepComponent{S2}

    function SymmorphicIrrepComponent(normal::IrrepComponent{S1}, rest::IrrepComponent{S2}
            ) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        if !iscompatible(normal.symmetry, normal.irrep_index, rest.symmetry)
            throw(ArgumentError("symmetry $(symmetry_name(rest.symmetry)) is not compatible with "*
                                "symmetry $(symmetry_name(normal.symmetry)) at irrep $(normal.irrep_index)"))
        end
        return new{S1, S2}(normal, rest)
    end

    function SymmorphicIrrepComponent(normal::IrrepComponent{SymmetryEmbedding{S1}},
                                      rest::IrrepComponent{SymmetryEmbedding{S2}}
            ) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        if !iscompatible(normal.symmetry, normal.irrep_index, rest.symmetry)
            throw(ArgumentError("symmetry $(symmetry_name(rest.symmetry)) is not compatible with "*
                                "symmetry $(symmetry_name(normal.symmetry)) at irrep $(normal.irrep_index)"))
        end
        return new{SymmetryEmbedding{S1}, SymmetryEmbedding{S2}}(normal, rest)
    end
end

group_order(arg::SymmorphicIrrepComponent) = group_order(arg.normal) * group_order(arg.rest)


function get_irrep_components(sym::Union{<:SymmorphicSymmetry, <:SymmorphicSymmetryEmbedding})
    (SymmorphicIrrepComponent(normal_sic, rest_sic)
        for normal_sic in get_irrep_components(sym.normal)
        for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest)))
end


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
    return ((relem * nelem, rphase * nphase)
                for (relem, rphase) in zip(rest_elements, rest_irrep_components)
                for (nelem, nphase) in zip(normal_elements, normal_irrep_components))
end
