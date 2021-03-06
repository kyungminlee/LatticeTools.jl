export IrrepData

export AbstractSymmetryIrrepComponent
export IrrepComponent

export group_order
export get_irrep_components
export get_irrep_iterator
export little_symmetry
export make_product_irrep

struct IrrepData
    group::FiniteGroup
    conjugacy_classes::Vector{Vector{Int}}
    character_table::Matrix
    irreps::Vector{Vector{Matrix{ComplexF64}}}
end


"""
    AbstractSymmetryIrrepComponent

Abstract type for symmetry irrep component.
"""
abstract type AbstractSymmetryIrrepComponent end


"""
    IrrepComponent{SymmetryType}

Represent an irrep component.
Consider a symmetry `S`, and  the i'th irrep Γᵢ, which has D dimensions.
The d-dimensional component (1 ≤ d ≤ D) is represented by `Irrep(S, i, D)`.
This type serves as a subscript for iterating over all components of all irreps.

# Parameters
* `SymmetryType<:SymmetryOrEmbedding`

# Fields
* `symmetry::SymmetryType`
* `irrep_index::Int`
* `irrep_component::Int`
"""
struct IrrepComponent{SymmetryType<:SymmetryOrEmbedding}<:AbstractSymmetryIrrepComponent
    symmetry::SymmetryType
    irrep_index::Int
    irrep_component::Int

    function IrrepComponent(
        sym::S,
        irrep_index::Integer,
        irrep_compo::Integer=1,
    ) where {S<:SymmetryOrEmbedding}
        if !(1 <= irrep_index <= num_irreps(sym))
            throw(ArgumentError("irrep index must be between 1 and $(num_irreps(sym))"))
        elseif !(1 <= irrep_compo <= irrep_dimension(sym, irrep_index))
            throw(ArgumentError("irrep component must be between 1 and $(irrep_dimension(sym, irrep_index))"))
        end
        return new{S}(sym, irrep_index, irrep_compo)
    end
end


function Base.:(==)(lhs::IrrepComponent, rhs::IrrepComponent)
    return (
        lhs.symmetry == rhs.symmetry &&
        lhs.irrep_index == rhs.irrep_index &&
        lhs.irrep_component == rhs.irrep_component
    )
end


"""
    group_order(sic::IrrepComponent)

Get order of the symmetry group of `sic`.
"""
GroupTools.group_order(sic::IrrepComponent) = group_order(sic.symmetry)


"""
    get_irrep_components(sym::SymmetryOrEmbedding)

Return a generator which gives `IrrepComponent(sym, irrep_index, irrep_component)`.
"""
function get_irrep_components(sym::SymmetryOrEmbedding)
    return (
        IrrepComponent(sym, irrep_index, irrep_compo)
        for irrep_index in 1:num_irreps(sym)
        for irrep_compo in 1:irrep_dimension(sym, irrep_index)
    )
end


"""
    get_irrep_iterator(sic::IrrepComponent)

Return a generator which gives (element, amplitude),
where amplitude is the irrep component of `element`.
"""
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


function get_irrep_iterator(sym::SymmetryOrEmbedding, index::Integer, component::Integer=1)
    elems = elements(sym)
    irrep_components = [m[component, component] for m in irrep(sym, index)]
    @assert length(elems) == length(irrep_components)
    return zip(elems, irrep_components)
end



"""
    little_group_elements(tsic, psym)

Return little group elements (i.e. indices) of `psym` corresponding to the irrep of translation symmetry specified by `tsic`.
`tsic` and `psym` are either
- `IrrepComponent{FiniteTranslationSymmetry}` and `PointSymmetry`, or
- `IrrepComponent{SymmetryEmbedding{FiniteTranslationSymmetry}}` and `SymmetryEmbedding{PointSymmetry}`
"""
function little_group_elements(
    tsic::IrrepComponent{FiniteTranslationSymmetry},
    psym::PointSymmetry,
) ::Vector{Int}
    return little_group_elements(tsic.symmetry, tsic.irrep_index, psym)
end


function little_group_elements(
    tsic::IrrepComponent{SymmetryEmbedding{FiniteTranslationSymmetry}},
    psymbed::SymmetryEmbedding{PointSymmetry}
)::Vector{Int}
    return little_group_elements(
        symmetry(tsic.symmetry),
        tsic.irrep_index,
        symmetry(psymbed)
    )
end


"""
    little_group(tsic, psym)

Return the `FiniteGroup` object that corresponds to the little group of `psym` at `tsic`.
"""
function little_group(
    tsic::IrrepComponent{FiniteTranslationSymmetry},
    psym::PointSymmetry,
)::FiniteGroup
    return little_group(tsic.symmetry, tsic.irrep_index, psym)
end


"""
    little_group(tsic, psymbed)

Return the `FiniteGroup` object that corresponds to the little group of `psymbed` at `tsic`.
"""
function little_group(
    tsic::IrrepComponent{SymmetryEmbedding{FiniteTranslationSymmetry}},
    psymbed::SymmetryEmbedding{PointSymmetry},
)::FiniteGroup
    return little_group(symmetry(tsic.symmetry), tsic.irrep_index, symmetry(psymbed))
end


"""
    little_symmetry(tsic, psym)

Return the `PointSymmetry` object that corresponds to the little group of psym at `tsic`.
"""
function little_symmetry(
    tsic::IrrepComponent{FiniteTranslationSymmetry},
    psym::PointSymmetry,
)::PointSymmetry
    return little_symmetry(tsic.symmetry, tsic.irrep_index, psym)
end


function little_symmetry(
    tsic::IrrepComponent{SymmetryEmbedding{FiniteTranslationSymmetry}},
    psymbed::SymmetryEmbedding{PointSymmetry},
)::SymmetryEmbedding{PointSymmetry}
    psym_little = little_symmetry(symmetry(tsic.symmetry), tsic.irrep_index, symmetry(psymbed))
    return embed(psymbed.lattice, psym_little)
end


"""
    iscompatible(tsic, psym)

Test whether the point symmetry `psym` is compatible with the irrep component `tsic`
(i.e. momentum) of the translation symmetry.
"""
function iscompatible(
    tsic::IrrepComponent{FiniteTranslationSymmetry},
    psym::PointSymmetry,
)::Bool
    return iscompatible(tsic.symmetry, tsic.irrep_index, psym)
end


function iscompatible(
    tsic::IrrepComponent{SymmetryEmbedding{FiniteTranslationSymmetry}},
    psymbed::SymmetryEmbedding{PointSymmetry},
)::Bool
    return (
        iscompatible(tsic.symmetry, psymbed) &&
        iscompatible(symmetry(tsic.symmetry), tsic.irrep_index, symmetry(psymbed))
    )
end



function make_product_irrep(irrep_iterators...; tol::Real=Base.rtoldefault(Float64))
    nested_symops_and_amplitude_list = [
        [(x, y) for (x, y) in irrep_iterator if !isapprox(y, zero(y); atol=tol)]
            for irrep_iterator in irrep_iterators
    ]
    symops_and_amplitudes = [
        (prod(reverse([op for (op, amp) in elems])), prod(amp for (op, amp) in elems))
            for elems in Iterators.product(nested_symops_and_amplitude_list...)
            # elems has the form of ((s, phis), (t, phit), (p, phit), ...)
            # we want the resulting element to be (p*t*s, ϕp*ϕt*ϕs). the phase commutes, while operations do not necessarily. 
    ]
    return symops_and_amplitudes
end



# struct SymmorphicSpaceIrrepComponent{S1<:SymmetryOrEmbedding, S2<:SymmetryOrEmbedding} <:AbstractSymmetryIrrepComponent

#     component1::IrrepComponent{S1} # e.g. Translation
#     component2::IrrepComponent{S2} # e.g. Point

#     function SymmorphicSpaceIrrepComponent(
#                 sic1::IrrepComponent{S1},
#                 sic2::IrrepComponent{S2}) where {
#             S1<:Union{FiniteTranslationSymmetry, SymmetryEmbedding{FiniteTranslationSymmetry}},
#             S2<:Union{PointSymmetry, SymmetryEmbedding{PointSymmetry}}}
#         sym1 = sic1.symmetry
#         sym_irrep_index1 = sic1.irrep_index
#         sym2 = sic2.symmetry
#         if !iscompatible(sym1, sym_irrep_index1, sym2)
#             throw(ArgumentError("symmetry $(symmetry_name(sym2)) is not compatible with "*
#                                 "symmetry $(symmetry_name(sym1)) at irrep $sym_irrep_index1"))
#         end
#         return new{S1, S2}(sic1, sic2)
#     end
# end


# group_order(sic::SymmorphicSpaceIrrepComponent) = group_order(sic.component1) * group_order(sic.component2)


# function get_irrep_components(tsym::S1,
#                               psym::S2) where {S1<:SymmetryOrEmbedding,
#                                                S2<:SymmetryOrEmbedding}
#     return (SymmorphicSpaceIrrepComponent(tsic, psic)
#                 for tsic in get_irrep_components(tsym)
#                 for psic in get_irrep_components(little_symmetry(tsic, psym)))
# end


# function get_irrep_iterator(ssic::SymmorphicSpaceIrrepComponent)

#     tsym = ssic.component1.symmetry
#     tsym_irrep_index = ssic.component1.irrep_index
#     tsym_irrep_compo = ssic.component1.irrep_component

#     psym = ssic.component2.symmetry
#     psym_irrep_index = ssic.component2.irrep_index
#     psym_irrep_compo = ssic.component2.irrep_component

#     tsym_elements = elements(tsym)
#     tsym_irrep = irrep(tsym, tsym_irrep_index)
#     tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep]

#     psym_elements = elements(psym)
#     psym_irrep = irrep(psym, psym_irrep_index)
#     psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep]

#     return ((psym_elem * tsym_elem, psym_phase * tsym_phase)
#                 for (psym_elem, psym_phase) in zip(psym_elements, psym_irrep_components)
#                 for (tsym_elem, tsym_phase) in zip(tsym_elements, tsym_irrep_components))
# end

# struct TranslationSymmetryIrrepComponent <:AbstractSymmetryIrrepComponent
#     symmetry::FiniteTranslationSymmetry
#     irrep_index::Int
#     irrep_component::Int
#     function TranslationSymmetryIrrepComponent(sym::FiniteTranslationSymmetry,
#                                                irrep_index::Integer,
#                                                irrep_compo::Integer=1)
#         if !(1 <= irrep_index <= num_irreps(sym))
#             throw(ArgumentError("translation symmetry irrep index should be between 1 and $(num_irreps(sym))"))
#         elseif irrep_compo != 1
#             throw(ArgumentError("irrep component of translation can only be 1 (since abelian)"))
#         end
#         return new(sym, irrep_index, irrep_compo)
#     end
# end


# group_order(sic::TranslationSymmetryIrrepComponent) = group_order(sic.symmetry)


# function get_irrep_components(tsym::FiniteTranslationSymmetry)
#     return (TranslationSymmetryIrrepComponent(tsym, irrep_index, 1)
#                 for irrep_index in 1:num_irreps(tsym))
# end

# function get_irrep_components(tsym_embed::TranslationSymmetryEmbedding)
#     get_irrep_components(tsym_embed.symmetry)
# end


# function get_irrep_iterator(lattice::Lattice,
#                             sic::TranslationSymmetryIrrepComponent)
#                             #tol::Real=Base.rtoldefault(Float64))
#     sym = sic.symmetry
#     permutations = get_site_permutations(lattice, sym)
#     irrep_components = let sym_irrep = irrep(sym, sic.irrep_index),
#                            c = sic.irrep_component
#                            [m[c, c] for m in sym_irrep]
#                        end
#     return (
#         (perm, phase) for (perm, phase) in zip(permutations, irrep_components)
#         #if !isapprox(phase, 0; atol=tol)
#     )
# end


# function Base.:(==)(lhs::TranslationSymmetryIrrepComponent,
#             rhs::TranslationSymmetryIrrepComponent)
#     lhs.symmetry == rhs.symmetry && lhs.irrep_index == rhs.irrep_index
# end


# struct PointSymmetryIrrepComponent <:AbstractSymmetryIrrepComponent
#     symmetry::PointSymmetry
#     irrep_index::Int
#     irrep_component::Int
#     function PointSymmetryIrrepComponent(sym::PointSymmetry,
#                                          irrep_index::Integer,
#                                          irrep_compo::Integer)
#         if !(1 <= irrep_index <= num_irreps(sym))
#             throw(ArgumentError("point symmetry irrep index should be between 1 and $(num_irreps(sym))"))
#         elseif !(1 <= irrep_compo <= irrep_dimension(sym, irrep_index))
#             throw(ArgumentError("point symmetry irrep component wrong"))
#         end
#         new(sym, irrep_index, irrep_compo)
#     end
# end


# group_order(sic::PointSymmetryIrrepComponent) = group_order(sic.symmetry)



# function get_irrep_components(psym::PointSymmetry)
#     return (PointSymmetryIrrepComponent(psym, irrep_index, 1)
#                 for irrep_index in 1:num_irreps(psym)
#                 for irrep_compo in 1:irrep_dimension(psym, irrep_index))
# end

# function get_irrep_components(psym_embed::PointSymmetryEmbedding)
#     get_irrep_components(psym_embed.symmetry)
# end


# function get_irrep_iterator(lattice::Lattice,
#                             sic::PointSymmetryIrrepComponent)
#                             #tol::Real=Base.rtoldefault(Float64))
#     sym = sic.symmetry
#     permutations = get_site_permutations(lattice, sym)
#     irrep_components = let irrep = irrep(sym, sic.irrep_index),
#                            c = sic.irrep_component
#                            [m[c, c] for m in irrep]
#                        end
#     return (
#         (perm, phase) for (perm, phase) in zip(permutations, irrep_components)
#         #if !isapprox(phase, 0; atol=tol)
#     )
# end



# function get_irrep_components(tsym::FiniteTranslationSymmetry,
#                               psym::PointSymmetry)
#     return (SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)
#                 for tsic in get_irrep_components(tsym)
#                 for psic in get_irrep_components(little_symmetry(tsic, psym)))
# end


# function get_irrep_components(ssym_embed::SymmorphicSpaceSymmetryEmbedding)
#     get_irrep_components(ssym_embed.translation_symmetry, ssym_embed.point_symmetry)
# end




# struct SymmorphicSpaceSymmetryIrrepComponent <:AbstractSymmetryIrrepComponent
#     translation::TranslationSymmetryIrrepComponent
#     point::PointSymmetryIrrepComponent

#     function SymmorphicSpaceSymmetryIrrepComponent(
#         tsic::TranslationSymmetryIrrepComponent,
#         psic::PointSymmetryIrrepComponent)

#         tsym = tsic.symmetry
#         tsym_irrep_index = tsic.irrep_index
#         psym = psic.symmetry
#         if !iscompatible(tsym, tsym_irrep_index, psym)
#             throw(ArgumentError("point symmetry $(psym.hermann_mauguin) is not compatible with translation symmetry $(tsym.hypercube.scale_matrix) at irrep $tsym_irrep_index"))
#         end
#         return new(tsic, psic)
#     end
# end


# group_order(sic::SymmorphicSpaceSymmetryIrrepComponent) = group_order(sic.translation) * group_order(sic.point)


# function get_irrep_iterator(lattice::Lattice,
#                             ssic::SymmorphicSpaceSymmetryIrrepComponent)
#                             #tol::Real=Base.rtoldefault(Float64))

#     tsym = ssic.translation.symmetry
#     tsym_irrep_index = ssic.translation.irrep_index
#     tsym_irrep_compo = ssic.translation.irrep_component

#     psym = ssic.point.symmetry
#     psym_irrep_index = ssic.point.irrep_index
#     psym_irrep_compo = ssic.point.irrep_component

#     tsym_permutations = get_site_permutations(lattice, tsym)
#     tsym_irrep = irrep(tsym, tsym_irrep_index)
#     tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep]

#     psym_permutations = get_site_permutations(lattice, psym)
#     psym_irrep = irrep(psym, psym_irrep_index)
#     psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep]

#     return (
#         (psym_perm * tsym_perm ,  psym_phase * tsym_phase)
#         for (psym_perm, psym_phase) in zip(psym_permutations,
#                                            psym_irrep_components)
#         for (tsym_perm, tsym_phase) in zip(tsym_permutations,
#                                            tsym_irrep_components)
#         #if !isapprox(psym_phase, 0; atol=tol) && !isapprox(psym_phase, 0; atol=tol)
#     )
# end


# function little_group_elements(tsic::TranslationSymmetryIrrepComponent,
#                                psym::PointSymmetry)
#     return little_group_elements(tsic.symmetry, tsic.irrep_index, psym)
# end


# function little_group(tsic::TranslationSymmetryIrrepComponent,
#                       psym::PointSymmetry)
#     return little_group(tsic.symmetry, tsic.irrep_index, psym)
# end


# function little_symmetry(tsic::TranslationSymmetryIrrepComponent,
#                          psym::PointSymmetry)
#     return little_symmetry(tsic.symmetry, tsic.irrep_index, psym)
# end


# function iscompatible(tsic::TranslationSymmetryIrrepComponent,
#                       psym::PointSymmetry)
#     return iscompatible(tsic.symmetry, tsic.irrep_index, psym)
# end
