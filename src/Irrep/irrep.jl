export IrrepData
export AbstractSymmetryIrrepComponent
export TranslationSymmetryIrrepComponent
export PointSymmetryIrrepComponent
export SymmorphicSpaceSymmetryIrrepComponent


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


function get_irrep_components(lattice::Lattice,
                              tsym::TranslationSymmetry)
    return (TranslationSymmetryIrrepComponent(tsym, irrep_index, 1)
                for irrep_index in 1:num_irreps(tsym))
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



function get_irrep_components(lattice::Lattice,
                              psym::PointSymmetry)
    return (PointSymmetryIrrepComponent(psym, irrep_index, 1)
                for irrep_index in 1:num_irreps(psym)
                for irrep_compo in 1:irrep_dimension(psym, irrep_index))
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



function get_irrep_components(lattice::Lattice,
                              tsym::TranslationSymmetry,
                              psym::PointSymmetry)
    return (SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)
                for tsic in get_irrep_components(lattice, tsym)
                for psic in get_irrep_components(lattice, little_symmetry(tsic, psym)) )
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
