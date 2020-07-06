module IrrepDatabase

import Pkg
import Serialization
using LinearAlgebra

using YAML

using ..TightBindingLattice: FiniteGroup, IrrepData
using ..TightBindingLattice: cleanup_number, parse_expr, group_isomorphism, group_multiplication_table, simplify_name


IRREP_DATABASE = IrrepData[]


function load_yaml()
    tol = Base.rtoldefault(Float64)
    global IRREP_DATABASE
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "Irreps"))
    database_raw = YAML.load_file(joinpath(data_directory, "irrep_database.yaml"))

    for data in database_raw
        group = FiniteGroup(transpose(hcat(data["MultiplicationTable"]...)))
        conjugacy_classes = data["ConjugacyClasses"]

        if group.conjugacy_classes != conjugacy_classes
            error("Conjugacy class not canonical") # COV_EXCL_LINE
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


function load()
    mutable_artifacts_toml = joinpath(@__DIR__, "..", "..", "MutableArtifacts.toml")
    try
        cache_hash = Pkg.Artifacts.artifact_hash("IrrepDatabase", mutable_artifacts_toml)
        if !isnothing(cache_hash) && Pkg.Artifacts.artifact_exists(cache_hash)
            cache_filepath = joinpath(Pkg.Artifacts.artifact_path(cache_hash), "IrrepDatabase.cache")
            global IRREP_DATABASE = Serialization.deserialize(cache_filepath)
            return
        end
    catch e
        @warn "Failed to load IrrepDatabase."
        @warn "$e"
        @warn "Loading YAML"
    end
    load_yaml()
    cache_hash = Pkg.Artifacts.create_artifact() do cache_directory
        global IRREP_DATABASE
        cache_filepath = joinpath(cache_directory, "IrrepDatabase.cache")
        Serialization.serialize(cache_filepath, IRREP_DATABASE)
    end
    Pkg.Artifacts.bind_artifact!(mutable_artifacts_toml, "IrrepDatabase", cache_hash; force=true)
end


function __init__()
    load()
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
