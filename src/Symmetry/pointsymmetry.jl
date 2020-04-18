export PointSymmetry

export  group_order,
        element, elements,
        element_name, element_names,
        group_multiplication_table,
        character_table,
        irrep, irreps, num_irreps, irrep_dimension,
        generator_elements, generator_indices

export project

export little_group_elements
export little_group
export little_symmetry, little_symmetry_iso

# export get_orbital_permutation, get_orbital_permutations
export get_irrep_iterator
export read_point_symmetry
export symmetry_name

simplify_name(name::AbstractString) = replace(replace(name, r"<sub>.*?</sub>"=>""), r"<sup>.*?</sup>"=>"")


## PointSymmetry, constructor and related functions

struct PointSymmetry <: AbstractSymmetry
    elements::Vector{PointOperation}
    group::FiniteGroup

    generators::Vector{Int}

    conjugacy_classes::Vector{Vector{Int}}
    character_table::Matrix{ComplexF64}
    irreps::Vector{Vector{Matrix{ComplexF64}}}

    element_names::Vector{String}

    matrix_representations::Vector{Matrix{Int}}
    hermann_mauguinn::String
    schoenflies::String

    function PointSymmetry(
            group::FiniteGroup,
            generators::AbstractVector{<:Integer},
            conjugacy_classes::AbstractVector{<:AbstractVector{<:Integer}},
            character_table::AbstractMatrix{<:Number},
            irreps::AbstractVector{<:AbstractVector{<:AbstractMatrix{<:Number}}},
            element_names::AbstractVector{<:AbstractString},
            matrix_representations::AbstractVector{<:AbstractMatrix{<:Integer}},
            hermann_mauguinn::AbstractString,
            schoenflies::AbstractString)
        tol = Base.rtoldefault(Float64)
        # number counting check
        if length(element_names) != group_order(group)
            throw(ArgumentError("number of elements different from order of group"))
        end
        if length(matrix_representations) != group_order(group)
            throw(ArgumentError("number of matrix representations different from order of group"))
        end
        if !allunique(matrix_representations)
            throw(ArgumentError("matrix representations have to be all unique"))
        end
        # TODO Test matrix representation for isomorphism
        if generate_subgroup(group, generators) != BitSet(1:group_order(group))
            throw(ArgumentError("Generators $(generators) does not generate the group"))
        end
        let prev_set = BitSet()
            for elements in conjugacy_classes
                !isempty(intersect(elements, prev_set)) && throw(ArgumentError("Same element in multiple conjugacy classes"))
                union!(prev_set, elements)
            end
            prev_set != BitSet(1:group_order(group)) && throw(ArgumentError("every element must belong to a conjugacy class"))
        end
        let n = length(conjugacy_classes)
            if size(character_table) != (n, n)
                throw(ArgumentError("size of character table does not match the number of conjugacy classes"))
            end
            if length(irreps) != n
                throw(ArgumentError("number of irreps match the number of conjugacy classes"))
            end
        end
        for rep in irreps
            if length(rep) != group_order(group)
                throw(ArgumentError("wrong number of matrices in irrep $rep"))
            end
            d = size(rep[1], 1)
            if !isapprox(rep[1], Matrix(I, (d,d)); atol=tol)
                throw(ArgumentError("matrix representation of identity should be identity"))
            end
            for m in rep
                if size(m) != (d,d)
                    throw(ArgumentError("matrix representation should all have the same dimension"))
                end
            end
            if !ishomomorphic(group, rep; equal=(x,y) -> isapprox(x,y;atol=tol))
                throw(ArgumentError("matrix representation not homomorphic to group"))
            end
        end
        let dim = size(first(matrix_representations), 1)
            if any(size(m) != (dim, dim) for m in matrix_representations)
                throw(ArgumentError("matrix representations should be square matrices of the same dimension"))
            end
        end

        elements = PointOperation.(matrix_representations)

        return new(elements,
                   group,
                   generators,
                   conjugacy_classes,
                   character_table,
                   irreps,
                   element_names,
                   matrix_representations,
                   hermann_mauguinn,
                   schoenflies)
    end
end


function read_point_symmetry(data::AbstractDict)
    tol = Base.rtoldefault(Float64)
    read_matrix(obj) = collect(transpose(hcat(parse_expr(obj)...)))

    multiplication_table = read_matrix(data["MultiplicationTable"])
    group = FiniteGroup(multiplication_table)
    ord_group = group_order(group)
    generators = data["Generators"]
    conjugacy_classes = [item for item in data["ConjugacyClasses"]]
    character_table = cleanup_number(read_matrix(data["CharacterTable"]), tol)
    irreps = [ Matrix{ComplexF64}[cleanup_number(read_matrix(elem), tol) for elem in item]
                   for item in data["IrreducibleRepresentations"] ]
    element_names = data["ElementNames"]
    matrix_representations = Matrix{Int}[read_matrix(x) for x in data["MatrixRepresentations"]]
    hermann_mauguinn = data["HermannMauguinn"]
    schoenflies = data["Schoenflies"]
    PointSymmetry(group, generators,
                  conjugacy_classes, character_table, irreps,
                  element_names, matrix_representations, hermann_mauguinn, schoenflies)
end


function project(psym::PointSymmetry,
                 projection::AbstractMatrix{<:Integer};
                 tol::Real=Base.rtoldefault(Float64))
    dim = size(psym.matrix_representations[1], 1)
    size(projection, 2) != dim && throw(ArgumentError("projection does not match matrix_representations dimension"))

    vals = LinearAlgebra.svdvals(projection)

    new_matrix_representations = [projection * m * transpose(projection) for m in psym.matrix_representations]

    return PointSymmetry(psym.group,
                         psym.generators,
                         psym.conjugacy_classes,
                         psym.character_table,
                         psym.irreps,
                         psym.element_names,
                         new_matrix_representations,
                         psym.hermann_mauguinn,
                         psym.schoenflies)
end


## Basic properties

element(sym::PointSymmetry, g) = sym.elements[g]
elements(sym::PointSymmetry) = sym.elements

element_name(sym::PointSymmetry, g) = sym.element_names[g]
element_names(sym::PointSymmetry) = sym.element_names

group(psym::PointSymmetry) = psym.group
group_order(psym::PointSymmetry) = group_order(psym.group)
group_multiplication_table(psym::PointSymmetry) = group_multiplication_table(psym.group)

character_table(sym::PointSymmetry) = sym.character_table
irreps(sym::PointSymmetry) = sym.irreps
irrep(sym::PointSymmetry, idx::Integer) = sym.irreps[idx]
num_irreps(sym::PointSymmetry) = length(sym.irreps)
irrep_dimension(sym::PointSymmetry, idx::Integer) = size(first(irrep(sym, idx)), 2)

generator_indices(sym::PointSymmetry) = sym.generators
generator_elements(sym::PointSymmetry) = element(sym, sym.generators)

symmetry_name(sym::PointSymmetry) = "PointSymmetry[$(sym.hermann_mauguinn)]"



# """
#     get_orbital_permutation(lattice, matrix_representation, orbital_map)
#
# Get a list of `Permutation` which represents the element of the point group
# specified by the `matrix_representation` and `orbital_map`, which respectively
# contain information about how the Bravais lattice transforms, and how the
# basis sites transforms.
# """
# function get_orbital_permutation(
#             lattice::Lattice,
#             matrix_representation::AbstractMatrix{<:Integer},
#             orbital_map::AbstractVector{<:Tuple{<:Integer, <:AbstractVector{<:Integer}}})
#     p = zeros(Int, numorbital(lattice.supercell))
#     for (i, (j, dR)) in enumerate(orbital_map)
#         namei = getorbitalname(lattice.unitcell, i)
#         namej = getorbitalname(lattice.unitcell, j)
#         for Ri in lattice.hypercube.coordinates
#             _, Rj = lattice.hypercube.wrap(matrix_representation * Ri + dR)
#             i_super = lattice.supercell.orbitalindices[(namei, Ri)]
#             j_super = lattice.supercell.orbitalindices[(namej, Rj)]
#             p[i_super] = j_super
#         end
#     end
#     return Permutation(p)
# end


# function get_orbital_permutations(lattice::Lattice, point_symmetry::PointSymmetry)
#     # orbitalmap contains how orbitals of the "unitcell" transform under
#     # every point group operations. Uses heuristics
#     orbitalmap = findorbitalmap(lattice.unitcell, point_symmetry)

#     permutations = Permutation[]
#     sizehint!(permutations, group_order(point_symmetry))

#     for (i_elem, (matrep, orbmap)) in enumerate(zip(point_symmetry.matrix_representations, orbitalmap))
#         p = zeros(Int, numorbital(lattice.supercell))
#         for (i, (j, dR)) in enumerate(orbmap)
#             namei = getorbitalname(lattice.unitcell, i)
#             namej = getorbitalname(lattice.unitcell, j)
#             for Ri in lattice.hypercube.coordinates
#                 _, Rj = lattice.hypercube.wrap(matrep * Ri + dR)
#                 i_super = lattice.supercell.orbitalindices[(namei, Ri)]
#                 j_super = lattice.supercell.orbitalindices[(namej, Rj)]
#                 p[i_super] = j_super
#             end
#         end
#         push!(permutations, Permutation(p))
#     end
#     return permutations
# end


# function get_irrep_iterator(lattice::Lattice,
#                             tsym::TranslationSymmetry,
#                             psym::PointSymmetry,
#                             tsym_irrep_index::Integer,
#                             psym_irrep_index::Integer,
#                             tsym_irrep_compo::Integer=1,
#                             psym_irrep_compo::Integer=1,
#                             tol::Real=Base.rtoldefault(Float64)
#                             )
#
#     if tsym_irrep_compo != 1 || psym_irrep_compo != 1
#         @warn "Currently only supports Gamma point, trivial point irrep"
#     end
#     #@assert iscompatible(tsym, psym)
#
#     tsym_permutations = get_orbital_permutations(lattice, tsym)
#     tsym_irrep = irrep(tsym, tsym_irrep_index)
#     tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep]
#
#     psym_permutations = get_orbital_permutations(lattice, psym)
#
#     psym_irrep = irrep(psym, psym_irrep_index)
#     psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep]
#
#     return (
#         (psym_perm * tsym_perm ,  psym_phase * tsym_phase)
#         for (tsym_perm, tsym_phase) in zip(tsym_permutations,
#                                            tsym_irrep_components)
#         for (psym_perm, psym_phase) in zip(psym_permutations,
#                                            psym_irrep_components)
#         #if abs(tsym_phase)>tol && abs(psym_phase) > tol
#     )
# end


# function parse_seitz(name::AbstractString)
#     m = match(r"([12346m])(<sup>([-+])</sup>)?(<sub>(.*)?</sub>)?", name)
#     main_str, chiral_str, axis_str = m[1], m[3], m[5]
#     axis = Int[]
#     while !isempty(axis_str)
#         m2 = match(r"^(-?\d)(.*)$", axis_str)
#         push!(axis, parse(Int, m2[1]))
#         axis_str = m2[2]
#     end
#     if !isempty(axis) && first(axis[axis!=0]) < 0
#         axis[:] = -axis
#     end
#     return (main_str, isnothing(chiral_str) ? "" : chiral_str, axis)
# end

# export get_hermann_mauguinn
# function get_hermann_mauguinn(group::FiniteGroup,
#                               element_names::AbstractVector{<:AbstractString})
#     # 1. find principal
#     is_rotation = [startswith(name, r"-?[2346]") for name in element_names]
#     is_reflection = [startswith(name, r"m") for name in element_names]
#     elements = parse_seitz.(element_names)
#
#     if any(is_rotation)
#         _, principal_rotation_index = minimum((-pl, i) for (i, pl) in enumerate(group.period_lengths) if is_rotation[i])
#
#
#     else
#         if any(is_reflection)
#             conjugacy_elements = [x[1] for x in elements[ [first(cc) for cc in group.conjugacy_classes[2:end]] ]]
#             return join(conjugacy_elements)
#         else
#             return "1"
#         end
#     end
# end
