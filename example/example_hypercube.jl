# # Test Hypercubic Lattice
#
using TightBindingLattice
using Plots
#-----------------------------------------------------------------------------
# Draw function
function draw_hypercube(hypercube ::HypercubicLattice)
    coordinates = hcat(hypercube.coordinates...)
    xlim = (minimum(coordinates[1,:]) - 3.5, maximum(coordinates[1,:]) + 3.5)
    ylim = (minimum(coordinates[2,:]) - 3.5, maximum(coordinates[2,:]) + 3.5)

    boundary = let
        r0 = [0,0]
        r1 = hypercube.scale_matrix[:,1]
        r2 = hypercube.scale_matrix[:,2]
        hcat(r0, r1, r1 .+ r2, r2, r0)
    end
    p = plot(boundary[1,:], boundary[2,:], label="")

    for Ri in Iterators.product(-1:1, -1:1)
        alpha = (Ri == (0,0)) ? 1.0 : 0.1
        Rr = hypercube.scale_matrix * [Ri...]
        scatter!(coordinates[1,:] .+ Rr[1], coordinates[2,:] .+ Rr[2], aspect_ratio=1, markersize=12, xlim=xlim, ylim=ylim, markeralpha=alpha, label="")
        annotation = []
        for i in 1:size(coordinates, 2)
            x = coordinates[1,i] + Rr[1]
            y = coordinates[2,i] + Rr[2]
            if xlim[1] < x < xlim[2] && ylim[1] < y < ylim[2]
                push!(annotation, (x, y, text("$i", 8, :black, :center)))
            end
        end
        if !isempty(annotation)
            annotation = [annotation...]
            annotate!(annotation)
        end
    end
    p
end
#-----------------------------------------------------------------------------
# ### (4,-4) x (4,4)
#-----------------------------------------------------------------------------
size_matrix = [ 4 4;
               -4 4]
hypercube = HypercubicLattice(size_matrix)
draw_hypercube(hypercube)
#-----------------------------------------------------------------------------
# ## All elements
for (i, c) in enumerate(hypercube.coordinates)
    println("$i : $c")
end
#-----------------------------------------------------------------------------
# ## Generators
translation_elements = [translation_element(hypercube, t) for t in hypercube.coordinates]
unit_translations = get_generators(hypercube)
println("Generators")
for it in unit_translations
    println("t(", it, ") = ", hypercube.coordinates[it], " (order = ", translation_elements[it].order, ")")
end
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ### (4,-1) x (1,4)
#-----------------------------------------------------------------------------
size_matrix = [ 4 2;
               -2 4]
hypercube = HypercubicLattice(size_matrix)
draw_hypercube(hypercube)
#-----------------------------------------------------------------------------
# ## All elements
for (i, c) in enumerate(hypercube.coordinates)
    println("$i : $c")
end
#-----------------------------------------------------------------------------
# ## Generators
translation_elements = [translation_element(hypercube, t) for t in hypercube.coordinates]
unit_translations = get_generators(hypercube)
for it in unit_translations
    println("t(", it, ") = ", hypercube.coordinates[it], " (order = ", translation_elements[it].order, ")")
end
#-----------------------------------------------------------------------------
