export TranslationSymmetry

export group_order,
       group_multiplication_table,
       character_table,
       irreps,
       irrep,
       irrep_dimension,
       num_irreps,
       element_names,
       element_name


export get_orbital_permutations
export iscompatible

export get_irrep_iterator


struct TranslationSymmetry <:AbstractSymmetry
    hypercube::HypercubicLattice

    group::FiniteGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{ConjugacyClassType}
    character_table::Matrix{ComplexF64}
    irreps::Vector{IrrepType}
    element_names::Vector{String}

    # for quick
    orthogonal_shape::Vector{Int}
    orthogonal_coordinates::Vector{Vector{Int}}

    orthogonal_to_coordinate_map ::Dict{Vector{Int}, Vector{Int}}
    coordinate_to_orthogonal_map ::Dict{Vector{Int}, Vector{Int}}

    function TranslationSymmetry(shape::Matrix{<:Integer}; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(HypercubicLattice(shape))
    end

    function TranslationSymmetry(lattice::Lattice; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(lattice.hypercube)
    end

    function TranslationSymmetry(hypercube::HypercubicLattice; tol::Real=Base.rtoldefault(Float64))
        group = FiniteGroup(translation_group_multiplication_table(hypercube))
        ord_group = group_order(group)

        @assert isabelian(group)
        @assert ord_group == length(hypercube.coordinates)

        generators = minimal_generating_set(group)

        orthogonal_shape = [group.period_lengths[g] for g in generators] # in "generator" coordinates
        orthogonal_coordinates = vec([[x...] for x in Iterators.product([0:(d-1) for d in orthogonal_shape]...)])
        @assert prod(orthogonal_shape) == ord_group
        @assert length(orthogonal_coordinates) == ord_group

        orthogonal_to_coordinate_map = Dict{Vector{Int}, Vector{Int}}()
        coordinate_to_orthogonal_map = Dict{Vector{Int}, Vector{Int}}()

        let ortho_latvec = hcat(hypercube.coordinates[generators]...)
            for r_ortho in orthogonal_coordinates
                #_, i = hypercube.torus_wrap(ortho_latvec * r_ortho)
                #r = hypercube.coordinates[i]
                _, r = hypercube.wrap(ortho_latvec * r_ortho)
                orthogonal_to_coordinate_map[r_ortho] = r
                coordinate_to_orthogonal_map[r] = r_ortho
            end
        end

        # each element of an abelian group is a conjugacy class
        element_names = ["$(orthogonal_to_coordinate_map[t])" for t in orthogonal_coordinates]
        conjugacy_classes = [(name=x, elements=[i]) for (i,x) in enumerate(element_names)]

        momentum(oc::AbstractVector{<:Integer}) = [2π * x / d for (x, d) in zip(oc, orthogonal_shape)]
        character_table = ComplexF64[cis(-dot(momentum(kd), t))
                                     for kd in orthogonal_coordinates,
                                          t in orthogonal_coordinates]

        character_table = cleanup_number(character_table, Base.rtoldefault(Float64))

        # each element forms a conjugacy class
        irreps = IrrepType[]
        for (idx_rep, momentum) in enumerate(orthogonal_coordinates)
            matrices = Matrix{ComplexF64}[]
            for (idx_elem, orthogonal_translation) in enumerate(orthogonal_coordinates)
                push!(matrices, character_table[idx_rep, idx_elem] * ones(ComplexF64, 1, 1))
            end
            push!(irreps, (name="$momentum", matrices=matrices))
        end

        return new(hypercube, group, generators,
                   conjugacy_classes, character_table, irreps, element_names,
                   orthogonal_shape, orthogonal_coordinates,
                   orthogonal_to_coordinate_map, coordinate_to_orthogonal_map)
    end
end


group_order(sym::TranslationSymmetry) = group_order(sym.group)
group_multiplication_table(psym::TranslationSymmetry) = group_multiplication_table(psym.group)

element_names(sym::TranslationSymmetry) = sym.element_names
element_name(tsym::TranslationSymmetry, g) = tsym.element_names[g]

character_table(sym::TranslationSymmetry) = sym.character_table

irreps(sym::TranslationSymmetry) = sym.irreps
irrep(sym::TranslationSymmetry, idx) = sym.irreps[idx]
num_irreps(sym::TranslationSymmetry) = length(sym.irreps)
irrep_dimension(sym::TranslationSymmetry, idx::Integer) = size(sym.irreps[idx].matrices[1], 2)


function get_orbital_permutations(lattice::Lattice,
                                  translation_symmetry::TranslationSymmetry)
    if lattice.hypercube != translation_symmetry.hypercube
        throw(ArgumentError("lattice and translation symmetry not consistent"))
    end

    n_uc = length(lattice.hypercube.coordinates)
    n_orb = numorbital(lattice.unitcell)

    permutations = Permutation[]
    for trans_ortho in translation_symmetry.orthogonal_coordinates

        trans_coord = translation_symmetry.orthogonal_to_coordinate_map[trans_ortho]

        p = zeros(Int, n_uc * n_orb)
        for (orbital_index1, ((orbital_name1, uc_coord1), _)) in enumerate(lattice.supercell.orbitals)
            #_, uc_index2 = lattice.hypercube.torus_wrap(uc_coord1 + trans_coord)
            #uc_coord2 = lattice.hypercube.coordinates[uc_index2]
            _, uc_coord2 = lattice.hypercube.wrap(uc_coord1 + trans_coord)

            orbital_index1 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord1))
            orbital_index2 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord2))
            p[orbital_index1] = orbital_index2
        end
        push!(permutations, Permutation(p))
    end
    return permutations
end


function get_irrep_iterator(lattice::Lattice,
                            tsym::TranslationSymmetry,
                            tsym_irrep_index::Integer,
                            tsym_irrep_compo::Integer=1)

    tsym_permutations = get_orbital_permutations(lattice, tsym)
    tsym_irrep = irrep(tsym, tsym_irrep_index)
    tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep.matrices]
    return zip(tsym_permutations, tsym_irrep_components)
end


function iscompatible(orthogonal_momentum::AbstractVector{<:Integer},
                      orthogonal_shape::AbstractVector{<:Integer},
                      identity_translation::AbstractVector{<:Integer})
    value = Rational{Int}(0)
    for (i, j, k) in zip(orthogonal_momentum, identity_translation, orthogonal_shape)
        value += i * j // k
    end
    return mod(value, 1) == 0
end


function iscompatible(orthogonal_momentum::AbstractVector{<:Integer},
                      orthogonal_shape::AbstractVector{<:Integer},
                      identity_translations::AbstractVector{<:AbstractVector{<:Integer}})
    return all(iscompatible(orthogonal_momentum, orthogonal_shape, t) for t in identity_translations)
end


function iscompatible(lattice::Lattice,
                      tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      identity_translation::AbstractVector{<:Integer})
    orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
    orthogonal_shape = tsym.orthogonal_shape
    return iscompatible(orthogonal_momentum, orthogonal_shape, identity_translation)
end

function iscompatible(lattice::Lattice,
                      tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      identity_translations::AbstractVector{<:AbstractVector{<:Integer}})
    orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
    orthogonal_shape = tsym.orthogonal_shape
    return all(iscompatible(orthogonal_momentum, orthogonal_shape, t) for t in identity_translations)
end
