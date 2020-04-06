
module PointSymmetryDatabase

using YAML
using JLD2
using FileIO

import ..TightBindingLattice: PointSymmetry, read_point_symmetry, simplify_name

NUM_POINT_SYMMETRIES = 32
POINT_SYMMETRY_DATABASE = Vector{PointSymmetry}(undef, NUM_POINT_SYMMETRIES)
POINT_SYMMETRY_LOOKUP = Dict{Vector{String}, Int}()

function __init__()
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "PointGroup3D"))
    if isfile(joinpath(data_directory, "PointGroup3D.jld2"))
        cache_database = load(joinpath(data_directory, "PointGroup3D.jld2"), "POINT_SYMMETRY_DATABASE")
        for gn in 1:NUM_POINT_SYMMETRIES
            POINT_SYMMETRY_DATABASE[gn] = cache_database[gn]
        end
    end
	# println("$(@__MODULE__)::__init")
    for i in 1:NUM_POINT_SYMMETRIES
        psym = get(i)
        POINT_SYMMETRY_LOOKUP[sort(simplify_name.(psym.element_names))] = i
    end
end

function load_group(groupnum)
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "PointGroup3D"))
    file_path = joinpath(data_directory, "PointGroup3D-$groupnum.yaml")
    data_yaml = YAML.load_file(file_path)
    point_symmetry = read_point_symmetry(data_yaml)
    return point_symmetry
end

function get(groupnum::Integer)
    (groupnum < 1 || groupnum > NUM_POINT_SYMMETRIES) && throw(ArgumentError("group #$groupnum not found"))
    if !isassigned(POINT_SYMMETRY_DATABASE, groupnum)
        POINT_SYMMETRY_DATABASE[groupnum] = load_group(groupnum)
    end
    return POINT_SYMMETRY_DATABASE[groupnum]
end

function find(element_names::AbstractVector{<:AbstractString})
    simple_names = sort(simplify_name.(element_names))
    return POINT_SYMMETRY_LOOKUP[simple_names]
end

end # module PointSymmetryDatabase
