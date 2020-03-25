export TranslationSymmetry
export PointSymmetry


struct TranslationSymmetry <:AbstractSymmetry
    hypercube::HypercubicLattice

    group::FiniteAbelianGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{ConjugacyClassType}
    character_table::Matrix{ComplexF64}
    irreducible_representations::Vector{RepresentationType}
    element_names::Vector{String}

    # for quick
    orthogonal_shape::Vector{Int}
    orthogonal_coordinates::Vector{Vector{Int}}

    orthogonal_to_coordinate_map ::Dict{Vector{Int}, Vector{Int}}
    coordinate_to_orthogonal_map ::Dict{Vector{Int}, Vector{Int}}

    function TranslationSymmetry(shape::Matrix{<:Integer}; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(HypercubicLattice(shape))
    end

    function TranslationSymmetry(hypercube::HypercubicLattice; tol::Real=Base.rtoldefault(Float64))
        group = FiniteAbelianGroup(translation_group_multiplication_table(hypercube))
        generators = minimal_generating_set(group)

        orthogonal_shape = [group.period_lengths[g] for g in generators] # in "generator" coordinates
        orthogonal_coordinates = vec([[x...] for x in Iterators.product([0:(d-1) for d in orthogonal_shape]...)])

        orthogonal_to_coordinate_map = Dict{Vector{Int}, Vector{Int}}()
        coordinate_to_orthogonal_map = Dict{Vector{Int}, Vector{Int}}()

        ortho_latvec = hcat(hypercube.coordinates[generators]...)
        for r_ortho in orthogonal_coordinates
            _, i = hypercube.torus_wrap(ortho_latvec * r_ortho)
            r = hypercube.coordinates[i]
            orthogonal_to_coordinate_map[r_ortho] = r
            coordinate_to_orthogonal_map[r] = r_ortho
        end

        element_names = ["$(orthogonal_to_coordinate_map[t])" for t in orthogonal_coordinates]
        conjugacy_classes = [(name=x,) for x in element_names]

        #fractional_momentum(oc::AbstractVector{<:Integer}) = [x/d for (x, d) in zip(oc, orthogonal_shape)]
        momentum(oc::AbstractVector{<:Integer}) = [2π * x / d for (x, d) in zip(oc, orthogonal_shape)]
        character_table = ComplexF64[cis(-dot(momentum(kd), t))
                                     for kd in orthogonal_coordinates,
                                          t in orthogonal_coordinates]

        # clean up values close to zero
        character_table_r = real.(character_table)
        character_table_i = imag.(character_table)
        character_table_r[abs.(character_table_r) .< tol] .= zero(Float64)
        character_table_i[abs.(character_table_i) .< tol] .= zero(Float64)
        character_table[:] = character_table_r[:] + im * character_table_i[:]

        # each element forms a conjugacy class
        irreducible_representations = RepresentationType[]
        for (idx_rep, momentum) in enumerate(orthogonal_coordinates)
            matrices = Matrix{ComplexF64}[]
            for (idx_elem, orthogonal_translation) in enumerate(orthogonal_coordinates)
                push!(matrices, character_table[idx_rep, idx_elem] * ones(ComplexF64, 1, 1))
            end
            push!(irreducible_representations, (name="$momentum", matrices=matrices))
        end

        return new(hypercube, group, generators,
                   conjugacy_classes, character_table, irreducible_representations, element_names,
                   orthogonal_shape, orthogonal_coordinates,
                   orthogonal_to_coordinate_map, coordinate_to_orthogonal_map)
    end
end



#=
struct BravaisPointSymmetry <: AbstractSymmetry
    hypercube::HypercubicLattice

    group::GenericGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{ConjugacyClassType}
    character_table::Matrix{ComplexF64}
    irreducible_representations::Vector{RepresentationType}
    element_names::Vector{String}

    matrix_representations
    schoenflies
    hermann_mauguinn
end
=#
