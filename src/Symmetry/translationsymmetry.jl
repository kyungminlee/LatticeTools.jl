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
export get_orbital_permutation
export iscompatible

export get_irrep_iterator


function decompose_lattice_2d(hypercube::HypercubicLattice)
    dimension(hypercube) != 2 && throw(ArgumentError("shape matrix should be 2x2"))
    group = FiniteGroup(translation_group_multiplication_table(hypercube))
    ord_group = group_order(group)
    allowed_pairs = Vector{Int}[[1,0]]
    sizehint!(allowed_pairs, ord_group^2*3÷4)
    for y in 1:ord_group, x in 0:ord_group
        if gcd(x, y) == 1
            push!(allowed_pairs, [x, y])
        end
    end
    for r1p in allowed_pairs
        if r1p[1] == 0
            continue
        end
        r1 = [r1p[1], -r1p[2]]
        j1 = hypercube.coordinate_indices[ hypercube.wrap(r1)[2] ]
        n1 = group_order(group, j1)
        for r2 in allowed_pairs
            j2 = hypercube.coordinate_indices[ hypercube.wrap(r2)[2] ]
            # condition 1/2: unimodular
            if ExactLinearAlgebra.determinant(hcat(r1, r2)) != 1
                continue
            end
            n2 = group_order(group, j2)
            nprod = length(generate_subgroup(group, [j1, j2]))
            # condition 2/2: orthogonal generators
            n1 * n2 == nprod && return hcat(r1, r2)
        end
    end
    error("Failed to find decomposition")
end

struct TranslationSymmetry <: AbstractSymmetry
    hypercube::HypercubicLattice
    elements::Vector{TranslationOperation{Int}}
    group::FiniteGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{Vector{Int}}
    character_table::Matrix{ComplexF64}
    irreps::Vector{Vector{Matrix{ComplexF64}}}
    element_names::Vector{String}

    # for quick
    generator_translations ::Matrix{Int}
    orthogonal_shape::Vector{Int}
    orthogonal_coordinates::Vector{Vector{Int}}

    orthogonal_to_coordinate_map ::Dict{Vector{Int}, Vector{Int}}
    coordinate_to_orthogonal_map ::Dict{Vector{Int}, Vector{Int}}

    orthogonal_scale_matrix::Matrix{Int}
    orthogonal_reduced_reciprocal_scale_matrix::Matrix{Rational{Int}}
    fractional_momenta::Vector{Vector{Rational{Int}}}

    function TranslationSymmetry(shape::Matrix{<:Integer}; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(HypercubicLattice(shape))
    end

    function TranslationSymmetry(lattice::Lattice; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(lattice.hypercube)
    end


    function TranslationSymmetry(hypercube::HypercubicLattice; tol::Real=Base.rtoldefault(Float64))
        if dimension(hypercube) == 1
            generator_translations = ones(Int, (1,1))
            return TranslationSymmetry(hypercube, generator_translations; tol=tol)
        elseif dimension(hypercube) == 2
            generator_translations = decompose_lattice_2d(hypercube)
            return TranslationSymmetry(hypercube, generator_translations; tol=tol)
        else
            error("Currenly only supports 1D and 2D")
        end
    end

    function TranslationSymmetry(hypercube::HypercubicLattice,
                                 generator_translations::AbstractMatrix{<:Integer};
                                 tol::Real=Base.rtoldefault(Float64))

        if ExactLinearAlgebra.determinant(generator_translations) != 1
            throw(ArgumentError("generator translation is not unimodular"))
        end

        group = FiniteGroup(translation_group_multiplication_table(hypercube))
        ord_group = group_order(group)

        @assert isabelian(group)
        @assert ord_group == length(hypercube.coordinates)

        elements = [TranslationOperation(v) for v in hypercube.coordinates]

        generators = Int[ hypercube.coordinate_indices[ hypercube.wrap(v)[2] ]
                             for v in eachcol(generator_translations) ]

        # if length(generators) < dimension(hypercube)
        #     @assert (length(generators) == 1 && dimension(hypercube) == 2) "Currently only supports 2D"
        #     ρ1 = hypercube.coordinates[first(generators)]
        #     R, r = hypercube.wrap(ρ1 * group_order(group, first(generators)))
        #     @assert iszero(r)
        #     gcd_val, R2p = extended_gcd(R...)
        #     @assert gcd_val == 1
        #     R2 = [-R2p[2], R2p[1]]
        #     ρ2 = hypercube.scale_matrix * R2
        #     # a x + b y = 1   =>  | x -b | = 1
        #     #                     | y  a |
        #     push!(generators, 1)
        #     generator_translations = hcat(ρ1, ρ2)
        # else
        #     generator_translations = hcat(hypercube.coordinates[generators]...)
        # end

        orthogonal_shape = [group.period_lengths[g] for g in generators] # in "generator" coordinates
        orthogonal_coordinates = vec([[x...] for x in Iterators.product([0:(d-1) for d in orthogonal_shape]...)])

        @assert prod(orthogonal_shape) == ord_group
        @assert length(orthogonal_coordinates) == ord_group

        orthogonal_to_coordinate_map = Dict{Vector{Int}, Vector{Int}}()
        coordinate_to_orthogonal_map = Dict{Vector{Int}, Vector{Int}}()

        let ortho_latvec = generator_translations #hcat(hypercube.coordinates[generators]...)
            for r_ortho in orthogonal_coordinates
                _, r = hypercube.wrap(ortho_latvec * r_ortho)
                orthogonal_to_coordinate_map[r_ortho] = r
                coordinate_to_orthogonal_map[r] = r_ortho
            end
        end
        orthogonal_scale_matrix = hcat(
                    [group_order(group, g) * v
                        for (g, v) in zip(generators, eachcol(generator_translations))]...
                )
        orthogonal_reduced_reciprocal_scale_matrix = ExactLinearAlgebra.inverse(transpose(orthogonal_scale_matrix))

        fractional_momenta = let mo = x -> mod(x, 1)
            [mo.( orthogonal_reduced_reciprocal_scale_matrix * orthogonal_integer_momentum )
                 for orthogonal_integer_momentum in orthogonal_coordinates]
        end

        # each element of an abelian group is a conjugacy class
        element_names = ["$(orthogonal_to_coordinate_map[t])" for t in orthogonal_coordinates]
        conjugacy_classes = [[i] for (i,x) in enumerate(element_names)]

        momentum(oc::AbstractVector{<:Integer}) = [2π * x / d for (x, d) in zip(oc, orthogonal_shape)]
        character_table = ComplexF64[cis(-dot(momentum(kd), t))
                                     for kd in orthogonal_coordinates,
                                          t in orthogonal_coordinates]

        character_table = cleanup_number(character_table, tol)

        # each element forms a conjugacy class
        irreps = Vector{Matrix{ComplexF64}}[]
        for (idx_rep, momentum) in enumerate(orthogonal_coordinates)
            matrices = Matrix{ComplexF64}[]
            for (idx_elem, orthogonal_translation) in enumerate(orthogonal_coordinates)
                push!(matrices, character_table[idx_rep, idx_elem] * ones(ComplexF64, 1, 1))
            end
            push!(irreps, matrices)
        end

        return new(hypercube, elements, group, generators,
                   conjugacy_classes, character_table, irreps, element_names,
                   generator_translations,
                   orthogonal_shape, orthogonal_coordinates,
                   orthogonal_to_coordinate_map, coordinate_to_orthogonal_map,
                   orthogonal_scale_matrix, orthogonal_reduced_reciprocal_scale_matrix,
                   fractional_momenta)
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
irrep_dimension(sym::TranslationSymmetry, idx::Integer) = 1 # size(first(sym.irreps[idx]), 1)


function get_orbital_permutations(lattice::Lattice,
                                  translation_symmetry::TranslationSymmetry)
    if lattice.hypercube != translation_symmetry.hypercube
        throw(ArgumentError("lattice and translation symmetry not consistent"))
    end
    # n_uc = length(lattice.hypercube.coordinates)
    # n_orb = numorbital(lattice.unitcell)
    permutations = Permutation[]
    for trans_ortho in translation_symmetry.orthogonal_coordinates
        trans_coord = translation_symmetry.orthogonal_to_coordinate_map[trans_ortho]
        push!(permutations, get_orbital_permutation(lattice, trans_coord))
        # p = zeros(Int, n_uc * n_orb)
        # for (orbital_index1, ((orbital_name1, uc_coord1), _)) in enumerate(lattice.supercell.orbitals)
        #     _, uc_coord2 = lattice.hypercube.wrap(uc_coord1 + trans_coord)
        #     orbital_index1 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord1))
        #     orbital_index2 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord2))
        #     p[orbital_index1] = orbital_index2
        # end
        # push!(permutations, Permutation(p))
    end
    return permutations
end


function get_orbital_permutation(lattice::Lattice,
                                 displacement::AbstractVector{<:Integer})
    p = zeros(Int, numorbital(lattice.supercell))
    for (orbital_index1, ((orbital_name1, uc_coord1), _)) in enumerate(lattice.supercell.orbitals)
        _, uc_coord2 = lattice.hypercube.wrap(uc_coord1 + displacement)
        orbital_index1 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord1))
        orbital_index2 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord2))
        p[orbital_index1] = orbital_index2
    end
    return Permutation(p)
end


function generators(lattice::Lattice, tsym::TranslationSymmetry)
    if lattice.hypercube != tsym.hypercube
        throw(ArgumentError("lattice and translation symmetry not consistent"))
    end
    n_uc = length(lattice.hypercube.coordinates)
    n_orb = numorbital(lattice.unitcell)
    dim = dimension(lattice)
    permutations = Permutation[]
    trans_ortho = zeros(Int, dim)
    for d in 1:dimension
        trans_ortho[:] = 0
        trans_ortho[d] = 1
        trans_coord = tsym.orthogonal_to_coordinate_map[trans_ortho]
        p = zeros(Int, n_uc * n_orb)
        for (orbital_index1, ((orbital_name1, uc_coord1), _)) in enumerate(lattice.supercell.orbitals)
            _, uc_coord2 = lattice.hypercube.wrap(uc_coord1 + trans_coord)
            orbital_index1 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord1))
            orbital_index2 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord2))
            p[orbital_index1] = orbital_index2
        end
        push!(permutations, Permutation(p))
    end
    return permutations
end



# function get_irrep_iterator(lattice::Lattice,
#                             tsym::TranslationSymmetry,
#                             tsym_irrep_index::Integer,
#                             tsym_irrep_compo::Integer=1)
#     tsym_permutations = get_orbital_permutations(lattice, tsym)
#     tsym_irrep = irrep(tsym, tsym_irrep_index)
#     tsym_irrep_components = (m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep)
#     return zip(tsym_permutations, tsym_irrep_components)
# end
#
#
# function get_irrep_iterator(lattice::Lattice,
#                             tsym::TranslationSymmetry,
#                             tsym_irrep_index::Integer,
#                             ::Colon)
#     tsym_permutations = get_orbital_permutations(lattice, tsym)
#     tsym_irrep = irrep(tsym, tsym_irrep_index)
#     tsym_irrep_components = view(tsym.character_table, tsym_irrep_index, :)
#     return zip(tsym_permutations, tsym_irrep_components)
# end


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
