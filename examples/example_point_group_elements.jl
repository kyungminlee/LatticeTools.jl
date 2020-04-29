
# Search elements by their Seitz name, and show their matrix representations

using TightBindingLattice

all_elements = Dict()
print("# Reading groups")
for i in 1:32
    print(" $i")
    psym = PointSymmetryDatabase.get(i)
    for (elname, matrep) in zip(psym.element_names, psym.matrix_representations)
        if !haskey(all_elements, elname)
            all_elements[elname] = []
        end
        push!(all_elements[elname], (matrep, psym.hermann_mauguinn))
    end
end
println()

println("elements:")
for (k, v) in all_elements
    println("- element_name: $k")
    println("  matrix_representations:")
    for v2 in v
        m = [collect(x) for x in eachrow(v2[1])]
        g = v2[2]
        println("  - { matrix: $(m), group: \"$g\" }")
    end
end
