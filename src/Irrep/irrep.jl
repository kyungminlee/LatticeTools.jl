export AbstractSymmetryIrrepComponent
export TranslationSymmetryIrrepComponent
export PointSymmetryIrrepComponent
export SymmorphicSpaceSymmetryIrrepComponent

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


function group_order(sic::TranslationSymmetryIrrepComponent)
    return group_order(sic.symmetry)
end


function get_irrep_iterator(lattice::Lattice,
                            sic::TranslationSymmetryIrrepComponent,
                            tol::Real=Base.rtoldefault(Float64))
    sym = sic.symmetry
    permutations = get_orbital_permutations(lattice, sym)
    irrep_components = let sym_irrep = irrep(sym, sic.irrep_index),
                           c = sic.irrep_component
                           [m[c, c] for m in sym_irrep.matrices]
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

function group_order(sic::PointSymmetryIrrepComponent)
    return group_order(sic.symmetry)
end


function get_irrep_iterator(lattice::Lattice,
                            sic::PointSymmetryIrrepComponent,
                            tol::Real=Base.rtoldefault(Float64))
    sym = sic.symmetry
    permutations = get_orbital_permutations(lattice, sym)
    irrep_components = let irrep = irrep(sym, sic.irrep_index),
                           c = sic.irrep_component
                           [m[c, c] for m in irrep.matrices]
                       end
    return (
        (perm, phase) for (perm, phase) in zip(permutations, irrep_components)
        #if !isapprox(phase, 0; atol=tol)
    )
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
            throw(ArgumentError("point symmetry not compatible with translation symmetry"))
        end
        return new(tsic, psic)
    end
end

function group_order(sic::SymmorphicSpaceSymmetryIrrepComponent)
    return group_order(sic.translation) * group_order(sic.point)
end


function get_irrep_iterator(lattice::Lattice,
    ssic::SymmorphicSpaceSymmetryIrrepComponent,
    tol::Real=Base.rtoldefault(Float64))

    tsym = ssic.translation.symmetry
    tsym_irrep_index = ssic.translation.irrep_index
    tsym_irrep_compo = ssic.translation.irrep_component

    psym = ssic.point.symmetry
    psym_irrep_index = ssic.point.irrep_index
    psym_irrep_compo = ssic.point.irrep_component

    if tsym_irrep_compo != 1 || psym_irrep_compo != 1
        @warn "Currently only supports Gamma point, trivial point irrep"
    end
    #@assert iscompatible(tsym, psym)

    tsym_permutations = get_orbital_permutations(lattice, tsym)
    tsym_irrep = irrep(tsym, tsym_irrep_index)
    tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep.matrices]

    psym_permutations = get_orbital_permutations(lattice, psym)
    psym_irrep = irrep(psym, psym_irrep_index)
    psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep.matrices]

    return (
        (psym_perm * tsym_perm ,  psym_phase * tsym_phase)
        for (psym_perm, psym_phase) in zip(psym_permutations,
                                           psym_irrep_components)
        for (tsym_perm, tsym_phase) in zip(tsym_permutations,
                                           tsym_irrep_components)
        #if !isapprox(psym_phase, 0; atol=tol) && !isapprox(psym_phase, 0; atol=tol)
    )
end


function little_symmetry(tsic::TranslationSymmetryIrrepComponent,
                         psym::PointSymmetry)
    return little_symmetry(tsic.symmetry, tsic.irrep_index, psym)
end
