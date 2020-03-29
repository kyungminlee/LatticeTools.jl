module PointSymmetryDatabase
using YAML

import ..TightBindingLattice: PointSymmetry

NUM_POINT_SYMMETRIES = 32
POINT_SYMMETRY_DATABASE = Vector{PointSymmetry}(undef, NUM_POINT_SYMMETRIES)

function load(groupnum)
    data_directory = abspath(joinpath(@__DIR__, "..", "..", "data", "PointGroup3D"))
    file_path = joinpath(data_directory, "PointGroup3D-$groupnum.yaml")
    data_yaml = YAML.load_file(file_path)
    point_symmetry = PointSymmetry(data_yaml)
    return point_symmetry
end

#=
function get(name::String)
    idx = find(x -> x.hermann_mauguinn == name, POINT_SYMMETRY_DATABASE)
    isnothing(idx) && throw(ArgumentError("group with name \"$name\" not found"))
    return get(idx)
end
=#

function get(groupnum::Integer)
    (groupnum < 1 || groupnum > NUM_POINT_SYMMETRIES) && throw(ArgumentError("group #$groupnum not found"))
    if !isassigned(POINT_SYMMETRY_DATABASE, groupnum)
        POINT_SYMMETRY_DATABASE[groupnum] = load(groupnum)
    end
    return POINT_SYMMETRY_DATABASE[groupnum]
end

end # module PointSymmetryDatabase
