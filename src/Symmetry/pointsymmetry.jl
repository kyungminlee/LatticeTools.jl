export PointSymmetry

#export PointSymmetryBravaisRepresentation
#export SymmorphicSpaceSymmetry
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

struct PointSymmetry
    group::FiniteGroup

    generators::Vector{Int}
    conjugacy_classes::Vector{ConjugacyClassType}
    character_table::Matrix{ComplexF64}
    irreps::Vector{IrrepType}
    element_names::Vector{String}

    matrix_representations::Vector{Matrix{Int}}
    hermann_mauguinn::String

    function PointSymmetry(
            group::FiniteGroup,
            generators::AbstractVector{<:Integer},
            conjugacy_classes::AbstractVector{ConjugacyClassType},
            character_table::AbstractMatrix{<:Number},
            irreps::AbstractVector{IrrepType},
            element_names::AbstractVector{<:AbstractString},
            matrix_representations::AbstractVector{<:AbstractMatrix{<:Integer}},
            hermann_mauguinn::AbstractString)

        dim = size(first(matrix_representations), 1)
        if any(size(m) != (dim, dim) for m in matrix_representations)
            throw(ArgumentError("matrix representations should be square matrices of the same dimension"))
        end
        return new(group,
                   generators,
                   conjugacy_classes,
                   character_table,
                   irreps,
                   element_names,
                   matrix_representations,
                   hermann_mauguinn)
    end
end

function read_point_symmetry(data::AbstractDict)
    tol = Base.rtoldefault(Float64)
    multiplication_table = transpose(hcat(parse_expr(data["MultiplicationTable"])...))

    group = FiniteGroup(multiplication_table)
    ord_group = group_order(group)

    generators = data["Generators"]
    if generate_subgroup(group, generators) != BitSet(1:group_order(group))
        throw(ArgumentError("Generators $(generators) does not generate the group"))
    end

    conjugacy_classes = [(name=item["Name"], elements=item["Elements"]) for item in data["ConjugacyClasses"]]
    let prev_set = BitSet()
        for (name, elements) in conjugacy_classes
            !isempty(intersect(elements, prev_set)) && throw(ArgumentError("Same element in multiple conjugacy classes"))
            union!(prev_set, elements)
        end
        prev_set != BitSet(1:ord_group) && throw(ArgumentError("every element must belong to a conjugacy class"))
    end

    character_table = cleanup_number(transpose(hcat(parse_expr(data["CharacterTable"])...)), tol)
    let nc = length(conjugacy_classes)
        size(character_table) != (nc, nc) && throw(ArgumentError("character table has wrong size"))
    end

    irreps = IrrepType[]
    for item in data["IrreducibleRepresentations"]
        matrices = Matrix{Number}[transpose(hcat(parse_expr(elem)...)) for elem in item["Matrices"]]
        new_item = (name=item["Name"], matrices=cleanup_number(matrices, tol))
        push!(irreps, new_item)
    end

    element_names = data["ElementNames"]
    matrix_representations = [transpose(hcat(parse_expr(x)...)) for x in data["MatrixRepresentations"]]
    hermann_mauguinn = data["HermannMauguinn"]

    PointSymmetry(group, generators,
                  conjugacy_classes, character_table, irreps,
                  element_names, matrix_representations, hermann_mauguinn)
end

group_order(psym::PointSymmetry) = group_order(psym.group)
group_multiplication_table(psym::PointSymmetry) = group_multiplication_table(psym.group)

element_name(sym::PointSymmetry, g) = sym.element_names[g]
element_names(sym::PointSymmetry) = sym.element_names

character_table(sym::PointSymmetry) = sym.character_table
irreps(sym::PointSymmetry) = sym.irreps
irrep(sym::PointSymmetry, idx::Integer) = sym.irreps[idx]
num_irreps(sym::PointSymmetry) = length(sym.irreps)
irrep_dimension(sym::PointSymmetry, idx::Integer) = size(sym.irreps[idx].matrices[1], 2)


function iscompatible(hypercube::HypercubicLattice, matrix_representation::AbstractMatrix{<:Integer})
    #_, elems = hypercube.torus_wrap(matrix_representation * hypercube.scale_matrix)
    #all(elems .== 1) # all, since scale_matrix
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
        @assert j > 0
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
    !allunique(new_matrix_representations) && @warn "projected matrix representations not unique"

    return PointSymmetry(psym.group,
                         psym.generators,
                         psym.conjugacy_classes,
                         psym.character_table,
                         psym.irreps,
                         psym.element_names,
                         new_matrix_representations,
                         psym.hermann_mauguinn,
                         )
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
                #_, idx_Rj = supercell.hypercube.torus_wrap(matrep * Ri + dR)
                #Rj = supercell.hypercube.coordinates[idx_Rj]
                _, Rj = lattice.hypercube.wrap(matrep * Ri + dR)
                #Rj = supercell.hypercube.coordinates[idx_Rj]

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
                               momentum_index::Integer,
                               psym::PointSymmetry)

    k_o2 = tsym.orthogonal_coordinates[momentum_index]
    lg = Int[]
    for (i_elem, matrep) in enumerate(psym.matrix_representations)
        k_o1 = tsym.orthogonal_to_coordinate_map[k_o2]
        k_r1 = tsym.hypercube.wrap(matrep * k_o1)[2]
        k_r2 = tsym.coordinate_to_orthogonal_map[k_r1]
        k_r2 == k_o2 && push!(lg, i_elem)
    end
    return lg
end


function little_group(tsym::TranslationSymmetry, irrep_index::Integer, psym::PointSymmetry)
    return little_group(tsym, psym, little_group_elements(tsym, irrep_index, psym))
end


function iscompatible(tsym::TranslationSymmetry,
                      tsym_irrep::Integer,
                      psym::PointSymmetry)
    ! iscompatible(tsym, psym) && return false

    return little_group_elements(tsym, tsym_irrep, psym) == 1:group_order(psym)

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
    (lg_raw, lg_matrep_raw) = let
        lg_elements = little_group_elements(tsym, tsym_irrep, psym)
        lg_raw = little_group(tsym, psym, lg_elements)
        lg_matrep_raw = psym.matrix_representations[lg_elements]
        (lg_raw, lg_matrep_raw)
    end
    little_symmetry_candidates = Tuple{PointSymmetry, Vector{Int}}[]
    for i in 1:32
        ps = PointSymmetryDatabase.get(i)
        ϕ = group_isomorphism(lg_raw, ps.group)
        if !isnothing(ϕ)
            push!(little_symmetry_candidates, (ps, ϕ))
            break
        end
    end

    #if length(little_symmetry_candidates) > 1
    #    @warn ("more than one matching point symmetry: " *
    #           join([ps.hermann_mauguinn for (ps,ϕ) in little_symmetry_candidates], ", "))
    #end

    (psym2, ϕ) = first(little_symmetry_candidates)

    lg_matrep = lg_matrep_raw[ϕ]

    # element names may be wrong.
    PointSymmetry(psym2.group,
                  psym2.generators,
                  psym2.conjugacy_classes,
                  psym2.character_table,
                  psym2.irreps,
                  psym2.element_names,
                  lg_matrep,
                  psym2.hermann_mauguinn)
end



function get_irrep_iterator(lattice::Lattice,
    tsym::TranslationSymmetry,
    psym::PointSymmetry,
    tsym_irrep_index::Integer,
    psym_irrep_index::Integer,
    tsym_irrep_compo::Integer=1,
    psym_irrep_compo::Integer=1,
    tol::Real=Base.rtoldefault(Float64)
    )

    if tsym_irrep_compo != 1 || psym_irrep_compo != 1
        @warn "Currently only supports Gamma point, trivial point irrep"
    end
    #@assert iscompatible(tsym, psym)

    tsym_permutations = get_orbital_permutations(lattice)
    tsym_irrep = irrep(tsym, tsym_irrep_index)
    tsym_irrep_components = [m[tsym_irrep_compo, tsym_irrep_compo] for m in tsym_irrep.matrices]

    psym_permutations = get_orbital_permutations(lattice, psym)
    #psym_little_group_generators = little_group(tsym, tsym_irrep_index, psym)

    psym_irrep = irrep(psym, psym_irrep_index)
    psym_irrep_components = [m[psym_irrep_compo, psym_irrep_compo] for m in psym_irrep.matrices]

    return [
        (psym_perm * tsym_perm ,  psym_phase * tsym_phase)
        for (tsym_perm, tsym_phase) in zip(tsym_permutations,
                                           tsym_irrep_components)
        for (psym_perm, psym_phase) in zip(psym_permutations,
                                           psym_irrep_components)
        if abs(tsym_phase)>tol && abs(psym_phase) > tol
    ]
end
