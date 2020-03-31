using TightBindingLattice

all_elements = Dict()
for i in 1:32
    println("Reading group #$i")
    psym = PointSymmetryDatabase.get(i)
    for (elname, matrep) in zip(psym.element_names, psym.matrix_representations)
        if !haskey(all_elements, elname)
            all_elements[elname] = []
        end
        push!(all_elements[elname], (matrep, psym.hermann_mauguinn))
    end
end

println("All elements:")
for (k, v) in all_elements
    println("- Element name: $k")
    println("  Matrix representations:")
    for v2 in v
        println("  - $v2")
    end
end
