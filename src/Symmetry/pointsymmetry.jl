export PointSymmetry

export group_order,
       group_multiplication_table,
       element_names,
       element_name,
       character_table,
       irreps,
       irrep,
       num_irreps,
       irrep_dimension

export iscompatible
export findorbitalmap
export project

export little_group_elements
export little_group
export little_symmetry

export get_orbital_permutations
export get_irrep_iterator
export read_point_symmetry


struct PointSymmetry <: AbstractSymmetry
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
        return new(group,
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
    # multiplication_table = transpose(hcat(parse_expr(data["MultiplicationTable"])...))
    # character_table = cleanup_number(transpose(hcat(parse_expr(data["CharacterTable"])...)), tol)
    # irreps = [ Matrix{ComplexF64}[transpose(hcat(parse_expr(elem)...)) for elem in item]
    #             for item in data["IrreducibleRepresentations"] ]
    # matrix_representations = Matrix{Int}[transpose(hcat(parse_expr(x)...)) for x in data["MatrixRepresentations"]]
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

group_order(psym::PointSymmetry) = group_order(psym.group)
group_multiplication_table(psym::PointSymmetry) = group_multiplication_table(psym.group)

element_name(sym::PointSymmetry, g) = sym.element_names[g]
element_names(sym::PointSymmetry) = sym.element_names

character_table(sym::PointSymmetry) = sym.character_table
irreps(sym::PointSymmetry) = sym.irreps
irrep(sym::PointSymmetry, idx::Integer) = sym.irreps[idx]
num_irreps(sym::PointSymmetry) = length(sym.irreps)
irrep_dimension(sym::PointSymmetry, idx::Integer) = size(first(irrep(sym, idx)), 2)


function iscompatible(hypercube::HypercubicLattice, matrix_representation::AbstractMatrix{<:Integer})
    _, elems = hypercube.wrap(matrix_representation * hypercube.scale_matrix)
    all(iszero(elems)) # all, since scale_matrix
end


function iscompatible(hypercube::HypercubicLattice, psym::PointSymmetry)
    return all(iscompatible(hypercube, m) for m in psym.matrix_representations)
end


iscompatible(tsym::TranslationSymmetry, psym::PointSymmetry) = iscompatible(tsym.hypercube, psym)


function findorbitalmap(unitcell::UnitCell, psym_matrep::AbstractMatrix{<:Integer})
    norb = numorbital(unitcell)
    map = Tuple{Int, Vector{Int}}[]
    for (orbname, orbfc) in unitcell.orbitals
        j, Rj = findorbitalindex(unitcell, psym_matrep * orbfc)
        j <= 0 && return nothing
        push!(map, (j, Rj))
    end
    return map
end


function findorbitalmap(unitcell::UnitCell, psym::PointSymmetry)
    return [findorbitalmap(unitcell, m) for m in psym.matrix_representations]
end


function project(psym::PointSymmetry,
                 projection::AbstractMatrix{<:Integer};
                 tol::Real=Base.rtoldefault(Float64))
    dim = size(psym.matrix_representations[1], 1)
    size(projection, 2) != dim && throw(ArgumentError("projection does not match matrix_representations dimension"))

    vals = LinearAlgebra.svdvals(projection)
    if ! all( isapprox(x, 1; atol=tol) || isapprox(x, 0; atol=tol) for x in vals)
        throw(ArgumentError("projection is not projection"))
    end

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


function get_orbital_permutations(lattice::Lattice, point_symmetry::PointSymmetry)
    # orbitalmap contains how orbitals of the "unitcell" transform under
    # every point group operations. Uses heuristics
    orbitalmap = findorbitalmap(lattice.unitcell, point_symmetry)

    permutations = Permutation[]
    sizehint!(permutations, group_order(point_symmetry))

    for (i_elem, (matrep, orbmap)) in enumerate(zip(point_symmetry.matrix_representations, orbitalmap))
        p = zeros(Int, numorbital(lattice.supercell))
        for (i, (j, dR)) in enumerate(orbmap)
            namei = getorbitalname(lattice.unitcell, i)
            namej = getorbitalname(lattice.unitcell, j)
            for Ri in lattice.hypercube.coordinates
                _, Rj = lattice.hypercube.wrap(matrep * Ri + dR)
                i_super = lattice.supercell.orbitalindices[(namei, Ri)]
                j_super = lattice.supercell.orbitalindices[(namej, Rj)]
                p[i_super] = j_super
            end
        end
        push!(permutations, Permutation(p))
    end
    return permutations
end


function little_group_elements(tsym::TranslationSymmetry,
                               tsym_irrep_index::Integer,
                               psym::PointSymmetry)

    k_o2 = tsym.orthogonal_coordinates[tsym_irrep_index]
    lg = Int[]
    for (i_elem, matrep) in enumerate(psym.matrix_representations)
        k_o1 = tsym.orthogonal_to_coordinate_map[k_o2]
        k_r1 = tsym.hypercube.wrap(matrep * k_o1)[2]
        k_r2 = tsym.coordinate_to_orthogonal_map[k_r1]
        k_r2 == k_o2 && push!(lg, i_elem)
    end
    return lg
end


function little_group(tsym::TranslationSymmetry, tsym_irrep_index::Integer, psym::PointSymmetry)
    return little_group(tsym, psym, little_group_elements(tsym, tsym_irrep_index, psym))
end


function iscompatible(tsym::TranslationSymmetry,
                      tsym_irrep_index::Integer,
                      psym::PointSymmetry)
    ! iscompatible(tsym, psym) && return false
    return little_group_elements(tsym, tsym_irrep_index, psym) == 1:group_order(psym)
end


"""
Generate a little group with given elements.
The elements of the little group, which may be sparse, are compressed into consecutive integers.
"""
function little_group(tsym::TranslationSymmetry,
                      psym::PointSymmetry,
                      elements::AbstractVector{<:Integer})
    element_lookup = Dict(x=>i for (i, x) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[ group_product(psym.group, elements[i], elements[j]) ]
    end
    return FiniteGroup(mtab)
end


function little_symmetry(tsym::TranslationSymmetry, tsym_irrep::Integer, psym::PointSymmetry)
    tsym_irrep == 1 && return psym
    (lg_raw, lg_matrep_raw, lg_element_names_raw) = let
        lg_elements = little_group_elements(tsym, tsym_irrep, psym)
        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        lg_element_names_raw = psym.element_names[lg_elements]
        (lg_raw, lg_matrep_raw, lg_element_names_raw)
    end

    group_index = PointSymmetryDatabase.find(lg_element_names_raw)
    psym2 = PointSymmetryDatabase.get(group_index)
    ϕ = group_isomorphism(lg_raw, psym2.group)
    isnothing(ϕ) && error("Group not isomorphic")

    lg_matrep = lg_matrep_raw[ϕ]
    lg_element_names = lg_element_names_raw[ϕ]

    PointSymmetry(psym2.group,
                  psym2.generators,
                  psym2.conjugacy_classes,
                  psym2.character_table,
                  psym2.irreps,
                  lg_element_names,
                  lg_matrep,
                  psym2.hermann_mauguinn,
                  psym2.schoenflies)
end


# function little_symmetry_iso(tsym::TranslationSymmetry, tsym_irrep::Integer, psym::PointSymmetry)
#     tsym_irrep == 1 && return psym
#     (lg_irrep, lg_matrep, lg_element_names) = let
#         lg_elements = little_group_elements(tsym, tsym_irrep, psym)
#
#         lg_raw = little_group(tsym, psym, lg_elements)
#         lg_matrep_raw = psym.matrix_representations[lg_elements]
#         lg_element_names_raw = psym.element_names[lg_elements]
#
#         (lg_irrep, ϕ) = IrrepDatabase.find(lg_raw)
#
#         lg_matrep = lg_matrep_raw[ϕ]
#         lg_element_names = lg_element_names_raw[ϕ]
#         (lg_irrep, lg_matrep, lg_element_names)
#     end
#
#     generators = minimal_generating_set(lg_irrep.group)
#     hermann_mauguinn = join(lg_element_names[generators])
#     schoenflies = "unknown"
#
#     simple_element_names = sort(simplify_name(lg_element_names))
#     for i in 1:32
#         psym = PointSymmetryDatabase.get(i)
#         if ( sort(simplify_name(psym.element_names)) == simple_element_names )
#              #&& !isnothing(group_isomorphism(lg_irrep.group, psym.group)) )
#             hermann_mauguinn = psym.hermann_mauguinn
#             break
#         end
#     end
#
#     PointSymmetry(lg_irrep.group,
#                   generators,
#                   lg_irrep.conjugacy_classes,
#                   lg_irrep.character_table,
#                   lg_irrep.irreps,
#                   lg_element_names,
#                   lg_matrep,
#                   hermann_mauguinn,
#                   schoenflies)
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

simplify_name(name::AbstractString) = replace(replace(name, r"<sub>.*?</sub>"=>""), r"<sup>.*?</sup>"=>"")
# simplify_name(names::AbstractVector{<:AbstractString}) = simplify_name.(names)

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
