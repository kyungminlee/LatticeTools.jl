export IrrepDatabase
module IrrepDatabase

using YAML
using LinearAlgebra

import ..TightBindingLattice: FiniteGroup, cleanup_number, parse_expr, group_isomorphism, group_multiplication_table, simplify_name


struct IrrepData
    group::FiniteGroup
    conjugacy_classes::Vector{Vector{Int}}
    character_table::Matrix
    irreps::Vector{Vector{Matrix{ComplexF64}}}
end

IRREP_DATABASE = IrrepData[]

function __init__()
    tol = Base.rtoldefault(Float64)
    global IRREP_DATABASE
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "Irreps"))
    database_raw = YAML.load_file(joinpath(data_directory, "irrep_database.yaml"))

    out = IrrepData[]
    for data in database_raw
        group = FiniteGroup(transpose(hcat(data["MultiplicationTable"]...)))
        conjugacy_classes = data["ConjugacyClasses"]

        @assert group.conjugacy_classes == conjugacy_classes

        character_table = cleanup_number(transpose(hcat(parse_expr(data["CharacterTable"])...)), tol)
        let nc = length(conjugacy_classes)
            size(character_table) != (nc, nc) && throw(ArgumentError("character table has wrong size"))
        end

        irreps = Vector{Matrix{ComplexF64}}[]
        for item in data["IrreducibleRepresentations"]
            matrices = Matrix{ComplexF64}[cleanup_number(transpose(hcat(parse_expr(elem)...)), tol) for elem in item]
            push!(irreps, matrices)
        end
        push!(out, IrrepData(group, conjugacy_classes, character_table, irreps))
    end
    IRREP_DATABASE = out
end


function find(group::FiniteGroup)
    global IRREP_DATABASE
    for irrep_item in IRREP_DATABASE
        group_found = irrep_item.group
        emap = group_isomorphism(group_found, group)
        if !isnothing(emap)
            # which_cc1 = zeros(Int, group_order(group))
            # for (icc, cc) in group.conjugacy_classes
            #     which_cc1[cc] .= icc
            # end
            # which_cc2 = zeros(Int, group_order(group))
            # for (icc, cc) in group_found.conjugacy_classes
            #     which_cc1[cc] .= icc
            # end
            #
            # cmap = zeros(Int, length(group.conjugacy_classes))
            # for (i, j) in enumerate(emap)
            #     cmap[which_cc1[i]] = which_cc2[j]
            #     @assert something
            # end
            # character_table = irrep_item.character_table[:, cmap]
            return (irrep=irrep_item, element_map=emap)
        end
    end
    return nothing
end


function find_irrep(matrep::AbstractVector{<:AbstractMatrix{<:Integer}})
    group = FiniteGroup(group_multiplication_table(matrep))
    (irrep, ϕ) = IrrepDatabase.find(group)
    matrep_new = matrep[ϕ]
    return (irrep, matrep_new, ϕ)
end


end
