export TranslationSymmetry
export dimension
export group, group_order, group_multiplication_table,
       character_table,
       irrep, irreps, irrep_dimension, num_irreps,
       element, elements,
       element_name, element_names,
       generator_elements, generator_indices,
       symmetry_product

export fractional_momentum
export isbragg

export symmetry_name

export translation_symmetry_embedding


struct TranslationSymmetry <: AbstractSymmetry{TranslationOperation{Int}}
    # hypercube::HypercubicLattice
    orthocube::OrthoCube
    
    elements::Vector{TranslationOperation{Int}}
    group::FiniteGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{Vector{Int}}
    character_table::Matrix{ComplexF64}
    irreps::Vector{Vector{Matrix{ComplexF64}}}
    element_names::Vector{String}

    # for quick
    generator_translations::Matrix{Int}
    coordinates::Vector{Vector{Int}}
    orthogonal_shape::Vector{Int}
    orthogonal_coordinates::Vector{Vector{Int}}

    orthogonal_to_coordinate_map::Dict{Vector{Int}, Vector{Int}}
    coordinate_to_orthogonal_map::Dict{Vector{Int}, Vector{Int}}

    orthogonal_shape_matrix::Matrix{Int}
    orthogonal_reduced_reciprocal_shape_matrix::Matrix{Rational{Int}}
    fractional_momenta::Vector{Vector{Rational{Int}}}

    function TranslationSymmetry(shape::Matrix{<:Integer}; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(OrthoCube(shape))
    end

    function TranslationSymmetry(lattice::Lattice; tol::Real=Base.rtoldefault(Float64))
        return TranslationSymmetry(lattice.orthocube)
    end

    function TranslationSymmetry(orthocube::OrthoCube;
                                 tol::Real=Base.rtoldefault(Float64))
        generator_translations = find_generators(orthocube)
        return TranslationSymmetry(orthocube, generator_translations; tol=tol)
    end
   
    function TranslationSymmetry(orthocube::OrthoCube,
                                 generator_translations::AbstractMatrix{<:Integer};
                                 tol::Real=Base.rtoldefault(Float64))
        if ExactLinearAlgebra.determinant(generator_translations) != 1
            throw(ArgumentError("generator translation is not unimodular"))
        end
        
        # group = FiniteGroup(translation_group_multiplication_table(hypercube))
        # ord_group = group_order(group)

        

        # BEGIN Orthogonal
        coordinates = generate_coordinates(orthocube, generator_translations)
        coordinate_indices = Dict(r => i for (i, r) in enumerate(coordinates))

        pbcadd = (x::Vector{Int}, y::Vector{Int}) -> orthocube.wrap(x+y)[2]
        group = FiniteGroup(group_multiplication_table(coordinates, pbcadd))
        @assert isabelian(group)

        ord_group = group_order(group)

        generators = Int[ coordinate_indices[ orthocube.wrap(v)[2] ]
                              for v in eachcol(generator_translations) ]

        orthogonal_shape = Int[group.period_lengths[g] for g in generators] # in "generator" coordinates
        orthogonal_coordinates = vec([[x...] for x in Iterators.product([0:(d-1) for d in orthogonal_shape]...)])

        @assert ord_group == prod(orthogonal_shape) == length(orthogonal_coordinates) == length(coordinates)

        orthogonal_to_coordinate_map = Dict{Vector{Int}, Vector{Int}}()
        coordinate_to_orthogonal_map = Dict{Vector{Int}, Vector{Int}}()

        let ortho_latvec = generator_translations
            for (r, r_ortho) in zip(coordinates, orthogonal_coordinates)
                @assert r == orthocube.wrap(ortho_latvec * r_ortho)[2]
                orthogonal_to_coordinate_map[r_ortho] = r
                coordinate_to_orthogonal_map[r] = r_ortho
            end
        end
        orthogonal_shape_matrix = hcat(
                    [group_order(group, g) * v
                        for (g, v) in zip(generators, eachcol(generator_translations))]...)
        orthogonal_reduced_reciprocal_shape_matrix = ExactLinearAlgebra.inverse(transpose(orthogonal_shape_matrix))
        
        elements = [TranslationOperation(r) for r in coordinates]

        fractional_momenta = let modunit = x -> mod(x, 1)
            [modunit.( orthogonal_reduced_reciprocal_shape_matrix * orthogonal_integer_momentum )
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
        
        let 
            χ = ComplexF64[cis(-2π * dot(k, t)) for k in fractional_momenta, t in coordinates]
            @assert isapprox(character_table, χ; atol=tol)
        end

        # each element forms a conjugacy class
        irreps = Vector{Matrix{ComplexF64}}[]
        for (idx_rep, momentum) in enumerate(orthogonal_coordinates)
            matrices = Matrix{ComplexF64}[]
            for (idx_elem, orthogonal_translation) in enumerate(orthogonal_coordinates)
                push!(matrices, character_table[idx_rep, idx_elem] * ones(ComplexF64, 1, 1))
            end
            push!(irreps, matrices)
        end

        return new(orthocube, elements, group, generators,
                   conjugacy_classes, character_table, irreps, element_names,
                   generator_translations,
                   coordinates,
                   orthogonal_shape, orthogonal_coordinates,
                   orthogonal_to_coordinate_map, coordinate_to_orthogonal_map,
                   orthogonal_shape_matrix, orthogonal_reduced_reciprocal_shape_matrix,
                   fractional_momenta)
    end
end


dimension(sym::TranslationSymmetry) = dimension(sym.orthocube)

group(sym::TranslationSymmetry) = sym.group
group_order(sym::TranslationSymmetry, g...) = group_order(sym.group, g...)
group_multiplication_table(psym::TranslationSymmetry) = group_multiplication_table(psym.group)

import Base.eltype
eltype(sym::TranslationSymmetry) = TranslationOperation{Int}

import Base.valtype
valtype(sym::TranslationSymmetry) = TranslationOperation{Int}


element(sym::TranslationSymmetry, g) = sym.elements[g]
elements(sym::TranslationSymmetry) = sym.elements

element_names(sym::TranslationSymmetry) = sym.element_names
element_name(sym::TranslationSymmetry, g) = sym.element_names[g]

character_table(sym::TranslationSymmetry) = sym.character_table

irreps(sym::TranslationSymmetry) = sym.irreps
irrep(sym::TranslationSymmetry, idx) = sym.irreps[idx]
num_irreps(sym::TranslationSymmetry) = length(sym.irreps)
irrep_dimension(sym::TranslationSymmetry, idx::Integer) = 1 # size(first(sym.irreps[idx]), 1)

generator_indices(sym::TranslationSymmetry) = sym.generators
generator_elements(sym::TranslationSymmetry) = element(sym, sym.generators)


"""
    symmetry_product(tsym::TranslationSymmetry)
"""
function symmetry_product(sym::TranslationSymmetry)
    function product(lhs::TranslationOperation, rhs::TranslationOperation)
        return TranslationOperation(sym.orthocube.wrap(lhs.displacement + rhs.displacement)[2])
    end
    return product
end


import Base.in
in(item::Any, sym::TranslationSymmetry) = false
function in(item::IdentityOperation, sym::TranslationSymmetry)
    dimension(item) == dimension(sym)
end
function in(item::TranslationOperation{<:Integer}, sym::TranslationSymmetry)
    dimension(item) == dimension(sym)
end
function in(item::PointOperation, sym::TranslationSymmetry) 
    dimension(item) == dimension(sym) && istranslation(item)
end
function in(item::SpaceOperation{Tp, <:Integer}, sym::TranslationSymmetry) where {Tp}
    dimension(item) == dimension(sym) && istranslation(item)
end


import Base.iterate
iterate(sym::TranslationSymmetry) = iterate(elements(sym))
iterate(sym::TranslationSymmetry, i) = iterate(elements(sym), i)


import Base.length
length(sym::TranslationSymmetry) = length(elements(sym))


function symmetry_name(sym::TranslationSymmetry)
    n11 = sym.orthocube.shape_matrix[1,1]
    n12 = sym.orthocube.shape_matrix[1,2]
    n21 = sym.orthocube.shape_matrix[2,1]
    n22 = sym.orthocube.shape_matrix[2,2]
    return "TranslationSymmetry[($n11,$n21)x($n12,$n22)]"
end


function fractional_momentum(sym::TranslationSymmetry, g...)
    return sym.fractional_momenta[g...]
end


raw"""
    isbragg(orthogonal_shape, orthogonal_momentum, identity_translation)

Test whether the given orthogonal momentum is compatible with the given identity translation. Since identity translation can only be in the trivial representation
(otherwise the phases cancel), this function tests whether the phase is unity.
The orthogonal momentum is given as an integer vector.

```math
  t(\rho) \vert \psi(k) \rangle
    = exp(i k \cdot \rho ) \vert \psi(k) \rangle
    = \vert \psi(k) \rangle
```

When the orthogonal shape is [n₁, n₂, ...], orthogonal momentum is chosen from
{0, 1, ..., n₁-1} × {0,1, ..., n₂-1} × ....
"""
function isbragg(shape::AbstractVector{<:Integer},
                 integer_momentum::AbstractVector{<:Integer},
                 identity_translation::AbstractVector{<:Integer})
    if length(integer_momentum) != length(shape)
        throw(DimensionMismatch("lengths of momentum and shape do not match"))
    elseif length(shape) != length(identity_translation)
        throw(DimensionMismatch("lengths of shape and translation do not match"))
    elseif any(shape .<= 0)
        throw(ArgumentError("shape needs to be all positive"))
    end
    value = Rational{Int}(0)
    for (i, j, k) in zip(integer_momentum, identity_translation, shape)
        value += i * j // k
    end
    return mod(value, 1) == 0
end


raw"""
    isbragg(orthogonal_shape, orthogonal_momentum, identity_translations)
"""
function isbragg(shape::AbstractVector{<:Integer},
                 integer_momentum::AbstractVector{<:Integer},
                 identity_translations::AbstractVector{<:AbstractVector{<:Integer}})
    return all(isbragg(shape, integer_momentum, t) for t in identity_translations)
end


"""
    isbragg(k, t)

Check for Bragg condition at momentum `k` and translation `t`.
Fractional momentum is normalized to 1, i.e. lies within [0, 1)ⁿ
"""
function isbragg(fractional_momentum::AbstractVector{<:Rational{<:Integer}},
                 translation::AbstractVector{<:Integer})
    if length(fractional_momentum) != length(translation)
        throw(DimensionMismatch("lengths of momentum and translation do not match"))
    end
    value = dot(fractional_momentum, translation)
    return mod(value, 1) == 0
end


raw"""
    isbragg(k, translations)
"""
function isbragg(fractional_momentum::AbstractVector{<:Rational{<:Integer}},
                 translations::AbstractVector{<:AbstractVector{<:Integer}})
    return all(isbragg(fractional_momentum, t) for t in translations)
end


"""
    isbragg(tsym, tsym_irrep_index, translation)
"""
function isbragg(tsym::TranslationSymmetry,
                 tsym_irrep_index::Integer,
                 translation::TranslationOperation{<:Integer})
    fractional_momentum = tsym.fractional_momenta[tsym_irrep_index]
    return isbragg(fractional_momentum, translation.displacement)
end

"""
    isbragg(tsym, tsym_irrep_index, translations)
"""
function isbragg(tsym::TranslationSymmetry,
                 tsym_irrep_index::Integer,
                 translations::AbstractVector{<:TranslationOperation{<:Integer}})
    fractional_momentum = tsym.fractional_momenta[tsym_irrep_index]
    return all(isbragg(fractional_momentum, t.displacement) for t in translations)
end


# """
#     iscompatible(tsym, tsym_irrep_index, identity_translation)

# Test whether the identity translation is compatible with the irreducible representation
# of the translation symmetry. The identity translation is in units of the "generators" of the lattice,
# i.e. orthogonal integer shape.
# """
# function iscompatible(tsym::TranslationSymmetry,
#                       tsym_irrep_index::Integer,
#                       identity_translation::AbstractVector{<:Integer})
#                       #identity_translation::TranslationOperation{<:Integer})
#     orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
#     orthogonal_shape = tsym.orthogonal_shape
#     return isbragg(orthogonal_shape, orthogonal_momentum, identity_translation.displacement)
# end

# function iscompatible(tsym::TranslationSymmetry,
#                       tsym_irrep_index::Integer,
#                       identity_translations::AbstractVector{<:AbstractVector{<:Integer}})
#                     #   identity_translations::AbstractVector{<:TranslationOperation{<:Integer}})
#     orthogonal_shape = tsym.orthogonal_shape
#     orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
#     return all(isbragg(orthogonal_shape, orthogonal_momentum, t.displacement)
#                    for t in identity_translations)
# end


# function iscompatible(# lattice::Lattice,
#                       tsym::TranslationSymmetry,
#                       tsym_irrep_index::Integer,
#                       identity_translation::AbstractVector{<:Integer})
#     # !iscompatible(lattice, tsym) && return false
#     orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
#     orthogonal_shape = tsym.orthogonal_shape
#     return iscompatible(orthogonal_momentum, orthogonal_shape, identity_translation)
# end


# function iscompatible(# lattice::Lattice,
#                       tsym::TranslationSymmetry,
#                       tsym_irrep_index::Integer,
#                       identity_translations::AbstractVector{<:AbstractVector{<:Integer}})
#     orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
#     orthogonal_shape = tsym.orthogonal_shape
#     return all(iscompatible(orthogonal_momentum, orthogonal_shape, t) for t in identity_translations)
# end


# function iscompatible(orthogonal_shape::AbstractVector{<:Integer},
#                       orthogonal_momentum::AbstractVector{<:Integer},
#                       identity_translation::TranslationOperation{<:Integer})
#     value = Rational{Int}(0)
#     for (n, i, j) in zip(orthogonal_shape, orthogonal_momentum, identity_translation.displacement)
#         value += i * j // n
#     end
#     return mod(value, 1) == 0
# end


# """
#     iscompatible(orthogonal_momentum, orthogonal_shape, identity_translations)

# Test whether the given momentum is compatible with *all* the identity translations.
# """
# function iscompatible(orthogonal_shape::AbstractVector{<:Integer},
#                       orthogonal_momentum::AbstractVector{<:Integer},
#                       identity_translations::AbstractVector{<:TranslationOperation{<:Integer}})
#     return all(iscompatible(orthogonal_shape, orthogonal_momentum, t) for t in identity_translations)
# end



# function generators(lattice::Lattice, tsym::TranslationSymmetry)
#     if lattice.hypercube != tsym.hypercube
#         throw(ArgumentError("lattice and translation symmetry not consistent"))
#     end
#     n_uc = length(lattice.hypercube.coordinates)
#     n_orb = numorbital(lattice.unitcell)
#     dim = dimension(lattice)
#     permutations = Permutation[]
#     trans_ortho = zeros(Int, dim)
#     for d in 1:dimension
#         trans_ortho[:] = 0
#         trans_ortho[d] = 1
#         trans_coord = tsym.orthogonal_to_coordinate_map[trans_ortho]
#         p = zeros(Int, n_uc * n_orb)
#         for (orbital_index1, ((orbital_name1, uc_coord1), _)) in enumerate(lattice.supercell.orbitals)
#             _, uc_coord2 = lattice.hypercube.wrap(uc_coord1 + trans_coord)
#             orbital_index1 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord1))
#             orbital_index2 = getorbitalindex(lattice.supercell, (orbital_name1, uc_coord2))
#             p[orbital_index1] = orbital_index2
#         end
#         push!(permutations, Permutation(p))
#     end
#     return permutations
# end




    # function TranslationSymmetry(shape::Matrix{<:Integer}; tol::Real=Base.rtoldefault(Float64))
    #     return TranslationSymmetry(orthogonalize(HypercubicLattice(shape)))
    # end

    # function TranslationSymmetry(lattice::Lattice; tol::Real=Base.rtoldefault(Float64))
    #     return TranslationSymmetry(lattice.hypercube)
    # end

    # function TranslationSymmetry(hypercube::HypercubicLattice; tol::Real=Base.rtoldefault(Float64))
    #     if dimension(hypercube) == 1
    #         generator_translations = ones(Int, (1,1))
    #         return TranslationSymmetry(hypercube, generator_translations; tol=tol)
    #     elseif dimension(hypercube) == 2
    #         generator_translations = decompose_lattice_2d(hypercube)
    #         return TranslationSymmetry(hypercube, generator_translations; tol=tol)
    #     else
    #         error("Currenly only supports 1D and 2D")
    #     end
    # end

    
    # function TranslationSymmetry(hypercube::HypercubicLattice,
    #                              generator_translations::AbstractMatrix{<:Integer};
    #                              tol::Real=Base.rtoldefault(Float64))

    #     if ExactLinearAlgebra.determinant(generator_translations) != 1
    #         throw(ArgumentError("generator translation is not unimodular"))
    #     end

    #     group = FiniteGroup(translation_group_multiplication_table(hypercube))
    #     ord_group = group_order(group)

    #     @assert isabelian(group)
    #     @assert ord_group == length(hypercube.coordinates)

    #     generators = Int[ hypercube.coordinate_indices[ hypercube.wrap(v)[2] ]
    #                          for v in eachcol(generator_translations) ]

    #     # BEGIN Orthogonal
    #     orthogonal_shape = [group.period_lengths[g] for g in generators] # in "generator" coordinates
    #     orthogonal_coordinates = vec([[x...] for x in Iterators.product([0:(d-1) for d in orthogonal_shape]...)])

    #     @assert prod(orthogonal_shape) == ord_group
    #     @assert length(orthogonal_coordinates) == ord_group

    #     orthogonal_to_coordinate_map = Dict{Vector{Int}, Vector{Int}}()
    #     coordinate_to_orthogonal_map = Dict{Vector{Int}, Vector{Int}}()

    #     let ortho_latvec = generator_translations #hcat(hypercube.coordinates[generators]...)
    #         for r_ortho in orthogonal_coordinates
    #             _, r = hypercube.wrap(ortho_latvec * r_ortho)
    #             orthogonal_to_coordinate_map[r_ortho] = r
    #             coordinate_to_orthogonal_map[r] = r_ortho
    #         end
    #     end
    #     orthogonal_shape_matrix = hcat(
    #                 [group_order(group, g) * v
    #                     for (g, v) in zip(generators, eachcol(generator_translations))]...
    #             )
    #     orthogonal_reduced_reciprocal_shape_matrix = ExactLinearAlgebra.inverse(transpose(orthogonal_shape_matrix))
    #     # END Orthogonal
        
    #     elements = [TranslationOperation(orthogonal_to_coordinate_map[r_ortho])
    #                     for r_ortho in orthogonal_coordinates]

    #     @assert hypercube.coordinates == [t.displacement for t in elements] "Remove this when hypercube coordinates convention is fixed"

    #     fractional_momenta = let modunit = x -> mod(x, 1)
    #         [modunit.( orthogonal_reduced_reciprocal_shape_matrix * orthogonal_integer_momentum )
    #              for orthogonal_integer_momentum in orthogonal_coordinates]
    #     end

    #     # each element of an abelian group is a conjugacy class
    #     element_names = ["$(orthogonal_to_coordinate_map[t])" for t in orthogonal_coordinates]
    #     conjugacy_classes = [[i] for (i,x) in enumerate(element_names)]

    #     momentum(oc::AbstractVector{<:Integer}) = [2π * x / d for (x, d) in zip(oc, orthogonal_shape)]
    #     character_table = ComplexF64[cis(-dot(momentum(kd), t))
    #                                  for kd in orthogonal_coordinates,
    #                                       t in orthogonal_coordinates]

    #     character_table = cleanup_number(character_table, tol)

    #     # each element forms a conjugacy class
    #     irreps = Vector{Matrix{ComplexF64}}[]
    #     for (idx_rep, momentum) in enumerate(orthogonal_coordinates)
    #         matrices = Matrix{ComplexF64}[]
    #         for (idx_elem, orthogonal_translation) in enumerate(orthogonal_coordinates)
    #             push!(matrices, character_table[idx_rep, idx_elem] * ones(ComplexF64, 1, 1))
    #         end
    #         push!(irreps, matrices)
    #     end

    #     return new(hypercube, elements, group, generators,
    #                conjugacy_classes, character_table, irreps, element_names,
    #                generator_translations,
    #                orthogonal_shape, orthogonal_coordinates,
    #                orthogonal_to_coordinate_map, coordinate_to_orthogonal_map,
    #                orthogonal_shape_matrix, orthogonal_reduced_reciprocal_shape_matrix,
    #                fractional_momenta)
    # end