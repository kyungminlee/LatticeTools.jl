# # Orthocube Examples (Bravais Lattice)

using TightBindingLattice
using Plots

function draw_orthocube(orthocube::OrthoCube, coordinates::AbstractMatrix{<:Integer})
    xlim = (minimum(coordinates[1,:]) - 3.5, maximum(coordinates[1,:]) + 3.5)
    ylim = (minimum(coordinates[2,:]) - 3.5, maximum(coordinates[2,:]) + 3.5)
    boundary = let
        r0 = [0,0]
        r1 = orthocube.shape_matrix[:,1]
        r2 = orthocube.shape_matrix[:,2]
        hcat(r0, r1, r1 .+ r2, r2, r0)
    end
    fig = plot(boundary[1,:], boundary[2,:], label="", size=(400, 400))
    for Ri in Iterators.product(-1:1, -1:1)
        alpha = (Ri == (0,0)) ? 1.0 : 0.1
        Rr = orthocube.shape_matrix * [Ri...]
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
    fig
end

# ## (4,-4) x (4,4)
size_matrix = [ 4 4; -4 4]
orthocube = OrthoCube(size_matrix)
generator_translations = find_generators(orthocube)
coordinates = generate_coordinates(orthocube, generator_translations)
coordmat = hcat(coordinates...)
println("All elements")
for (i, c) in enumerate(coordinates)
    println("$i : $c")
end
println("Generator translations")
for (it, t) in enumerate(eachcol(generator_translations))
    println("t($it) = $t")
end
draw_orthocube(orthocube, coordmat)

# ## (2,-2) x (2,4)
size_matrix = [ 2 2; -2 4]
orthocube = OrthoCube(size_matrix)
generator_translations = find_generators(orthocube)
coordinates = generate_coordinates(orthocube, generator_translations)
coordmat = hcat(coordinates...)
draw_orthocube(orthocube, coordmat)
println("All elements")
for (i, c) in enumerate(coordinates)
    println("$i : $c")
end
println("Generator translations")
for (it, t) in enumerate(eachcol(generator_translations))
    println("t($it) = $t")
end
draw_orthocube(orthocube, coordmat)
