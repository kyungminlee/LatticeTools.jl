export IrrepData

export AbstractSymmetryIrrepComponent
export IrrepComponent
export SymmorphicIrrepComponent

export group_order
export get_irrep_components
export get_irrep_iterator
export little_symmetry

struct IrrepData
    group::FiniteGroup
    conjugacy_classes::Vector{Vector{Int}}
    character_table::Matrix
    irreps::Vector{Vector{Matrix{ComplexF64}}}
end

abstract type AbstractSymmetryIrrepComponent end

struct IrrepComponent{SymmetryType}<:AbstractSymmetryIrrepComponent
    symmetry::SymmetryType
    irrep_index::Int
    irrep_component::Int
    function IrrepComponent(sym::S,
                            irrep_index::Integer,
                            irrep_compo::Integer=1) where {S<:AbstractSymmetry}
        if !(1 <= irrep_index <= num_irreps(sym))
            throw(ArgumentError("irrep index must be between 1 and $(num_irreps(sym))"))
        elseif !(1 <= irrep_compo <= irrep_dimension(sym, irrep_index))
            throw(ArgumentError("irrep component must be between 1 and $(irrep_dimension(sym, irrep_index))"))
        end
        return new{S}(sym, irrep_index, irrep_compo)
    end
end

import Base.==

function ==(lhs::IrrepComponent, rhs::IrrepComponent)
    lhs.symmetry == rhs.symmetry &&
    lhs.irrep_index == rhs.irrep_index &&
    lhs.irrep_component == rhs.irrep_component
end

group_order(sic::IrrepComponent) = group_order(sic.symmetry)

function get_irrep_components(sym::AbstractSymmetry)
    return (IrrepComponent(sym, irrep_index, 1)
                for irrep_index in 1:num_irreps(sym)
                for irrep_compo in 1:irrep_dimension(sym, irrep_index))
end

function get_irrep_iterator(sic::IrrepComponent)
    sym = sic.symmetry
    elems = elements(sic.symmetry)
    irrep_components = let irrep = irrep(sym, sic.irrep_index),
                           c = sic.irrep_component
                           [m[c, c] for m in irrep]
                       end
    @assert length(elems) == length(irrep_components)
    return zip(elems, irrep_components)
end


struct SymmorphicIrrepComponent{S1<:AbstractSymmetry, S2<:AbstractSymmetry} <:AbstractSymmetryIrrepComponent

    component1::IrrepComponent{S1} # e.g. Translation
    component2::IrrepComponent{S2} # e.g. Point

    function SymmorphicIrrepComponent(
                sic1::IrrepComponent{S1},
                sic2::IrrepComponent{S2}) where {
            S1<:Union{TranslationSymmetry, SymmetryEmbedding{TranslationSymmetry}},
            S2<:Union{PointSymmetry, SymmetryEmbedding{PointSymmetry}}}
        sym1 = sic1.symmetry
        sym_irrep_index1 = sic1.irrep_index
        sym2 = sic2.symmetry
        if !iscompatible(sym1, sym_irrep_index1, sym2)
            throw(ArgumentError("symmetry $(symmetry_name(sym2)) is not compatible with "*
                                "symmetry $(symmetry_name(sym1)) at irrep $sym_irrep_index1"))
        end
        return new{S1, S2}(sic1, sic2)
    end
end


group_order(sic::SymmorphicIrrepComponent) = group_order(sic.component1) * group_order(sic.component2)


function little_group_elements(tsic::IrrepComponent{TranslationSymmetry},
                               psym::PointSymmetry) ::Vector{Int}
    return little_group_elements(tsic.symmetry, tsic.irrep_index, psym)
end

function little_group_elements(tsic::IrrepComponent{SymmetryEmbedding{TranslationSymmetry}},
                               psymbed::SymmetryEmbedding{PointSymmetry}) ::Vector{Int}
    return little_group_elements(symmetry(tsic.symmetry),
                                 tsic.irrep_index,
                                 symmetry(psymbed))
end

function little_group(tsic::IrrepComponent{TranslationSymmetry},
                      psym::PointSymmetry)::FiniteGroup
    return little_group(tsic.symmetry, tsic.irrep_index, psym)
end

function little_group(tsic::IrrepComponent{SymmetryEmbedding{TranslationSymmetry}},
                      psymbed::SymmetryEmbedding{PointSymmetry})::FiniteGroup
    return little_group(symmetry(tsic.symmetry), tsic.irrep_index, symmetry(psymbed))
end


function little_symmetry(tsic::IrrepComponent{TranslationSymmetry},
                         psym::PointSymmetry)::PointSymmetry
    return little_symmetry(tsic.symmetry, tsic.irrep_index, psym)
end


function little_symmetry(tsic::IrrepComponent{SymmetryEmbedding{TranslationSymmetry}},
                         psymbed::SymmetryEmbedding{PointSymmetry})::SymmetryEmbedding{PointSymmetry}
    psym_little = little_symmetry(symmetry(tsic.symmetry),
                                  tsic.irrep_index,
                                  symmetry(psymbed))
    return embed(psymbed.lattice, psym_little)
end


function iscompatible(tsic::IrrepComponent{TranslationSymmetry},
                      psym::PointSymmetry)::Bool
    return iscompatible(tsic.symmetry, tsic.irrep_index, psym)
end


function iscompatible(tsic::IrrepComponent{SymmetryEmbedding{TranslationSymmetry}},
                      psymbed::SymmetryEmbedding{PointSymmetry})::Bool
    # TODO: Check lattice?
    return iscompatible(symmetry(tsic.symmetry),
                        tsic.irrep_index,
                        symmetry(psymbed))
end


function get_irrep_components(tsym::S1,
                              psym::S2) where {S1<:AbstractSymmetry,
                                               S2<:AbstractSymmetry}
    return (SymmorphicIrrepComponent(tsic, psic)
                for tsic in get_irrep_components(tsym)
                for psic in get_irrep_components(little_symmetry(tsic, psym)))
end


function get_irrep_iterator(lattice::Lattice,
                            ssic::SymmorphicIrrepComponent)

    tsym = ssic.component1.symmetry
    tsym_irrep_index = ssic.component1.irrep_index
    tsym_irrep_compo = ssic.component1.irrep_component

    psym = ssic.component2.symmetry
    psym_irrep_index = ssic.component2.irrep_index
    psym_irrep_compo = ssic.component2.irrep_component
    
    tsym_elements = elements(tsym)
    tsym_irrep = irrep(tsym, tsym_irrep_index)
    tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep]

    psym_elements = elements(psym)
    psym_irrep = irrep(psym, psym_irrep_index)
    psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep]

    return ((psym_elem * tsym_elem, psym_phase * tsym_phase)
                for (psym_elem, psym_phase) in zip(psym_elements, psym_irrep_components)
                for (tsym_elem, tsym_phase) in zip(tsym_elements, tsym_irrep_components))
end




#=

struct TranslationSymmetryIrrepComponent <:AbstractSymmetryIrrepComponent
    symmetry::TranslationSymmetry
    irrep_index::Int
    irrep_component::Int
    function TranslationSymmetryIrrepComponent(sym::TranslationSymmetry,
                                               irrep_index::Integer,
                                               irrep_compo::Integer=1)
        if !(1 <= irrep_index <= num_irreps(sym))
            throw(ArgumentError("translation symmetry irrep index should be between 1 and $(num_irreps(sym))"))
        elseif irrep_compo != 1
            throw(ArgumentError("irrep component of translation can only be 1 (since abelian)"))
        end
        return new(sym, irrep_index, irrep_compo)
    end
end


group_order(sic::TranslationSymmetryIrrepComponent) = group_order(sic.symmetry)


function get_irrep_components(tsym::TranslationSymmetry)
    return (TranslationSymmetryIrrepComponent(tsym, irrep_index, 1)
                for irrep_index in 1:num_irreps(tsym))
end

function get_irrep_components(tsym_embed::TranslationSymmetryEmbedding)
    get_irrep_components(tsym_embed.symmetry)
end


function get_irrep_iterator(lattice::Lattice,
                            sic::TranslationSymmetryIrrepComponent)
                            #tol::Real=Base.rtoldefault(Float64))
    sym = sic.symmetry
    permutations = get_orbital_permutations(lattice, sym)
    irrep_components = let sym_irrep = irrep(sym, sic.irrep_index),
                           c = sic.irrep_component
                           [m[c, c] for m in sym_irrep]
                       end
    return (
        (perm, phase) for (perm, phase) in zip(permutations, irrep_components)
        #if !isapprox(phase, 0; atol=tol)
    )
end


import Base.==

function ==(lhs::TranslationSymmetryIrrepComponent,
            rhs::TranslationSymmetryIrrepComponent)
    lhs.symmetry == rhs.symmetry && lhs.irrep_index == rhs.irrep_index
end


struct PointSymmetryIrrepComponent <:AbstractSymmetryIrrepComponent
    symmetry::PointSymmetry
    irrep_index::Int
    irrep_component::Int
    function PointSymmetryIrrepComponent(sym::PointSymmetry,
                                         irrep_index::Integer,
                                         irrep_compo::Integer)
        if !(1 <= irrep_index <= num_irreps(sym))
            throw(ArgumentError("point symmetry irrep index should be between 1 and $(num_irreps(sym))"))
        elseif !(1 <= irrep_compo <= irrep_dimension(sym, irrep_index))
            throw(ArgumentError("point symmetry irrep component wrong"))
        end
        new(sym, irrep_index, irrep_compo)
    end
end


group_order(sic::PointSymmetryIrrepComponent) = group_order(sic.symmetry)



function get_irrep_components(psym::PointSymmetry)
    return (PointSymmetryIrrepComponent(psym, irrep_index, 1)
                for irrep_index in 1:num_irreps(psym)
                for irrep_compo in 1:irrep_dimension(psym, irrep_index))
end

function get_irrep_components(psym_embed::PointSymmetryEmbedding)
    get_irrep_components(psym_embed.symmetry)
end


function get_irrep_iterator(lattice::Lattice,
                            sic::PointSymmetryIrrepComponent)
                            #tol::Real=Base.rtoldefault(Float64))
    sym = sic.symmetry
    permutations = get_orbital_permutations(lattice, sym)
    irrep_components = let irrep = irrep(sym, sic.irrep_index),
                           c = sic.irrep_component
                           [m[c, c] for m in irrep]
                       end
    return (
        (perm, phase) for (perm, phase) in zip(permutations, irrep_components)
        #if !isapprox(phase, 0; atol=tol)
    )
end



function get_irrep_components(tsym::TranslationSymmetry,
                              psym::PointSymmetry)
    return (SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)
                for tsic in get_irrep_components(tsym)
                for psic in get_irrep_components(little_symmetry(tsic, psym)))
end


function get_irrep_components(ssym_embed::SymmorphicSpaceSymmetryEmbedding)
    get_irrep_components(ssym_embed.translation_symmetry, ssym_embed.point_symmetry)
end




struct SymmorphicSpaceSymmetryIrrepComponent <:AbstractSymmetryIrrepComponent
    translation::TranslationSymmetryIrrepComponent
    point::PointSymmetryIrrepComponent

    function SymmorphicSpaceSymmetryIrrepComponent(
        tsic::TranslationSymmetryIrrepComponent,
        psic::PointSymmetryIrrepComponent)

        tsym = tsic.symmetry
        tsym_irrep_index = tsic.irrep_index
        psym = psic.symmetry
        if !iscompatible(tsym, tsym_irrep_index, psym)
            throw(ArgumentError("point symmetry $(psym.hermann_mauguinn) is not compatible with translation symmetry $(tsym.hypercube.scale_matrix) at irrep $tsym_irrep_index"))
        end
        return new(tsic, psic)
    end
end


group_order(sic::SymmorphicSpaceSymmetryIrrepComponent) = group_order(sic.translation) * group_order(sic.point)


function get_irrep_iterator(lattice::Lattice,
                            ssic::SymmorphicSpaceSymmetryIrrepComponent)
                            #tol::Real=Base.rtoldefault(Float64))

    tsym = ssic.translation.symmetry
    tsym_irrep_index = ssic.translation.irrep_index
    tsym_irrep_compo = ssic.translation.irrep_component

    psym = ssic.point.symmetry
    psym_irrep_index = ssic.point.irrep_index
    psym_irrep_compo = ssic.point.irrep_component

    tsym_permutations = get_orbital_permutations(lattice, tsym)
    tsym_irrep = irrep(tsym, tsym_irrep_index)
    tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep]

    psym_permutations = get_orbital_permutations(lattice, psym)
    psym_irrep = irrep(psym, psym_irrep_index)
    psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep]

    return (
        (psym_perm * tsym_perm ,  psym_phase * tsym_phase)
        for (psym_perm, psym_phase) in zip(psym_permutations,
                                           psym_irrep_components)
        for (tsym_perm, tsym_phase) in zip(tsym_permutations,
                                           tsym_irrep_components)
        #if !isapprox(psym_phase, 0; atol=tol) && !isapprox(psym_phase, 0; atol=tol)
    )
end


function little_group_elements(tsic::TranslationSymmetryIrrepComponent,
                               psym::PointSymmetry)
    return little_group_elements(tsic.symmetry, tsic.irrep_index, psym)
end


function little_group(tsic::TranslationSymmetryIrrepComponent,
                      psym::PointSymmetry)
    return little_group(tsic.symmetry, tsic.irrep_index, psym)
end


function little_symmetry(tsic::TranslationSymmetryIrrepComponent,
                         psym::PointSymmetry)
    return little_symmetry(tsic.symmetry, tsic.irrep_index, psym)
end


function iscompatible(tsic::TranslationSymmetryIrrepComponent,
                      psym::PointSymmetry)
    return iscompatible(tsic.symmetry, tsic.irrep_index, psym)
end

=#