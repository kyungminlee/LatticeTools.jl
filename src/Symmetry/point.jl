export TranslationSymmetry
export PointSymmetry
export PointSymmetryBravaisRepresentation


struct PointSymmetry
    group::FiniteGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{ConjugacyClassType}
    character_table::Matrix{ComplexF64}
    irreducible_representations::Vector{RepresentationType}
    element_names::Vector{String}

    matrix_representations::Vector{Matrix{Int}}
    schoenflies::String
    hermann_mauguinn::String

    function PointSymmetry(data::AbstractDict)
        schoenflies = data["Schoenflies"]
        hermann_mauguinn = data["HermannMauguinn"]
        element_names = data["ElementNames"]
        multiplication_table = parse_expr(data["MultiplicationTable"])
        conjugacy_classes = [(name=item["Name"],) for item in data["ConjugacyClasses"]]
        character_table = parse_expr(data["CharacterTable"])

        matrix_representations = [parse_expr(x) for x in data["MatrixRepresentations"]]

        irreducible_representations = NamedTuple{(:name, :matrices), Tuple{String, Vector{Matrix{Number}}}}[]
        for item in data["IrreducibleRepresentations"]
            matrices = Matrix{Number}[hcat(parse_expr(elem)) for elem in item["Matrices"]]
            new_item = (name=item["Name"], matrices=matrices)
            push!(irreducible_representations, new_item)
        end

        group = FiniteGroup(multiplication_table)
        generators = minimal_generating_set(group)
        new(group, generators, conjugacy_classes, character_table, irreducible_representations, element_names,
            matrix_representations, schoenflies, hermann_mauguinn)
    end
end



function iscompatible(hypercube::HypercubicLattice, matrix_representation::AbstractMatrix{<:Integer})
    _, elems = hypercube.torus_wrap(matrix_representation * hypercube.scale_matrix)
    all(elems .== 1) # all, since scale_matrix
end

function project(
        matrix_representations::AbstractVector{<:AbstractMatrix{<:Integer}},
        projection::AbstractMatrix{Int})

    isempty(matrix_representations) && throw(ArgumentError("matrix representations cannot be empty"))
    dim = size(matrix_representations[1], 1)
    if any(size(m) != (dim, dim) for m in matrix_representations)
        throw(ArgumentError("matrix representations should be square matrices of the same dimension"))
    end
    size(projection, 2) != dim && throw(ArgumentError("projection does not match matrix_representations dimension"))

    tol = Base.rtoldefault(Float64)
    vals = LinearAlgebra.svdvals(projection)
    if ! all( isapprox(x, 1; atol=tol) || isapprox(x, 0; atol=tol) for x in vals)
        throw(ArgumentError("projection is not projection"))
    end

    out = [projection * m * transpose(projection) for m in matrix_representations]
    !allunique(out) && @warn "projected matrix representations not unique"
    out
end




struct PointSymmetryBravaisRepresentation
    symmetry::PointSymmetry
    hypercube::HypercubicLattice
    matrix_representations::Vector{Matrix{Int}}
    permutations::Vector{Permutation}

    function PointSymmetryBravaisRepresentation(
                    symmetry::PointSymmetry,
                    hypercube::HypercubicLattice,
                    projection::AbstractMatrix{<:Integer})
        size(projection, 1) != dimension(hypercube) && throw(ArgumentError("number of rows of `projection` should match the dimension of hypercube"))
        matrix_representations = project(symmetry.matrix_representations, projection)

        if !all(iscompatible(hypercube, m) for m in matrix_representations)
            throw(ArgumentError("hypercube not compatible with matrix representation. (Origin should remain invariant.)"))
        end

        permutations = Permutation[]
        sizehint!(permutations, length(matrix_representations))
        for (i_elem, mr) in enumerate(matrix_representations)
            p = zeros(Int, length(hypercube.coordinates))
            for (i, r_tilde) in enumerate(hypercube.coordinates)
                _, i_prime = hypercube.torus_wrap(mr * r_tilde)
                p[i_prime] = i  # Important! not the other way around!
            end
            push!(permutations, Permutation(p))
        end
        return new(symmetry, hypercube, matrix_representations, permutations)
    end
end
