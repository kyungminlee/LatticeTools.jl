export SymmorphicIrrepComponent
export get_irrep_components
export get_irrep_iterator

struct SymmorphicIrrepComponent{S1<:AbstractSymmetry, S2<:AbstractSymmetry}<:AbstractSymmetryIrrepComponent
    normal_symmetry::S1
    normal_irrep_index::Int
    normal_irrep_component::Int

    rest_symmetry::S2 # this needs to be little symmetry
    rest_irrep_index::Int
    rest_irrep_component::Int

    function SymmorphicIrrepComponent(
                normal_symmetry::S1, normal_irrep_index::Integer, normal_irrep_component::Integer,
                rest_symmetry::S2, rest_irrep_index::Integer, rest_irrep_component::Integer
            ) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        if !iscompatible(normal_symmetry, normal_irrep_index, rest_symmetry)
            throw(ArgumentError("symmetry $(symmetry_name(rest_symmetry)) is not compatible with "*
                                "symmetry $(symmetry_name(normal_symmetry)) at irrep $normal_irrep_index"))
        end
        return new{S1, S2}(normal_symmetry, normal_irrep_index, normal_irrep_component,
                           rest_symmetry, rest_irrep_index, rest_irrep_component)
    end
end


function get_irrep_components(sym::Union{<:SymmorphicSymmetry, <:SymmorphicSymmetryEmbedding})
    (SymmorphicIrrepComponent(normal_sic.symmetry, normal_sic.irrep_index, normal_sic.irrep_component,
                              rest_sic.symmetry, rest_sic.irrep_index, rest_sic.irrep_component)
        for normal_sic in get_irrep_components(sym.normal)
        for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest)))
end


function get_irrep_iterator(ssic::SymmorphicIrrepComponent)
    normal_elements, normal_irrep_components = let
        nsym = ssic.normal_symmetry
        nii  = ssic.normal_irrep_index
        nic  = ssic.normal_irrep_component
        elements(nsym), [m[nic, nic] for m in irrep(nsym, nii)]
    end
    rest_elements, rest_irrep_components = let
        rsym = ssic.rest_symmetry
        rii  = ssic.rest_irrep_index
        ric  = ssic.rest_irrep_component
        elements(rsym), [m[ric, ric] for m in irrep(rsym, rii)]
    end
    return ((relem * nelem, rphase * nphase)
                for (relem, rphase) in zip(rest_elements, rest_irrep_components)
                for (nelem, nphase) in zip(normal_elements, normal_irrep_components))
end
