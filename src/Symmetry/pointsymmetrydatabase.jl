
module PointSymmetryDatabase

using YAML
using JLD2
using FileIO

import ..TightBindingLattice: PointSymmetry, read_point_symmetry, simplify_name

NUM_POINT_SYMMETRIES = 32
POINT_SYMMETRY_DATABASE = Vector{PointSymmetry}(undef, NUM_POINT_SYMMETRIES)
POINT_SYMMETRY_LOOKUP = Dict{Vector{String}, Int}()

function __init__()
    for i in 1:NUM_POINT_SYMMETRIES
        psym = get(i)
        POINT_SYMMETRY_LOOKUP[sort(simplify_name.(psym.element_names))] = i
    end
end

function load_group(group_index::Integer)
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "PointGroup3D"))
    file_path = joinpath(data_directory, "PointGroup3D-$group_index.yaml")
    data_yaml = YAML.load_file(file_path)
    point_symmetry = read_point_symmetry(data_yaml)
    return point_symmetry
end

function get(group_index::Integer)
    (group_index < 1 || group_index > NUM_POINT_SYMMETRIES) && throw(ArgumentError("Point group 3D #$group_index not found"))
    if !isassigned(POINT_SYMMETRY_DATABASE, group_index)
        POINT_SYMMETRY_DATABASE[group_index] = load_group(group_index)
    end
    return POINT_SYMMETRY_DATABASE[group_index]
end

function find(group_name::AbstractString)
    for i in 1:32
        psym = get(i)
        if group_name == psym.hermann_mauguinn
            return psym
        end
    end
    return nothing
end

function find(element_names::AbstractVector{<:AbstractString})
    simple_names = sort(simplify_name.(element_names))
    return POINT_SYMMETRY_LOOKUP[simple_names]
end

end # module PointSymmetryDatabase
