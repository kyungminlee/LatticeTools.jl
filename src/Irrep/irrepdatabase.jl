export IrrepDatabase
module IrrepDatabase

using YAML
using LinearAlgebra

using ..TightBindingLattice: FiniteGroup, IrrepData
using ..TightBindingLattice: cleanup_number, parse_expr, group_isomorphism, group_multiplication_table, simplify_name

IRREP_DATABASE = IrrepData[]

function __init__()
    tol = Base.rtoldefault(Float64)

    global IRREP_DATABASE
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "Irreps"))
    database_raw = YAML.load_file(joinpath(data_directory, "irrep_database.yaml"))

    for data in database_raw
        group = FiniteGroup(transpose(hcat(data["MultiplicationTable"]...)))
        conjugacy_classes = data["ConjugacyClasses"]

        if group.conjugacy_classes != conjugacy_classes
            error("Conjugacy class not canonical")
        end

        character_table = cleanup_number(transpose(hcat(parse_expr(data["CharacterTable"])...)), tol)
        let nc = length(conjugacy_classes)
            size(character_table) != (nc, nc) && error("character table has wrong size")
        end

        irreps = Vector{Matrix{ComplexF64}}[]
        for item in data["IrreducibleRepresentations"]
            matrices = Matrix{ComplexF64}[cleanup_number(transpose(hcat(parse_expr(elem)...)), tol) for elem in item]
            push!(irreps, matrices)
        end
        push!(IRREP_DATABASE, IrrepData(group, conjugacy_classes, character_table, irreps))
    end
end


function find(group::FiniteGroup)
    global IRREP_DATABASE
    for irrep_item in IRREP_DATABASE
        group_found = irrep_item.group
        emap = group_isomorphism(group_found, group)
        if !isnothing(emap)
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
