using TightBindingLattice

#data_yaml = YAML.load_file(joinpath(data_directory, "PointGroup3D-13.yaml"))
point_symmetry = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])

for (ename, mat) in zip(element_names(point_symmetry), irrep(point_symmetry, 5).matrices)
    println(ename)
    display(Int.(real.(mat)))
    println()
    println()
end
