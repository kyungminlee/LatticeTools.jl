# # Honeycomb lattice

# ## Preamble
using TightBindingLattice
using LinearAlgebra
using Plots

simplifyname(s::AbstractString) = (
                    s |> (x-> replace(x, r"<sup>(.+?)</sup>" => s"\1"))
                      |> (x-> replace(x, r"<sub>(.+?)</sub>" => s"[\1]"))
                )
mkpath("example_honeycomb_symmetry")
extent = [-2.5, 2.5, -2.5, 2.5]
within(r) = (extent[1] <= r[1] <= extent[2] && extent[3] <= r[2] <= extent[4]);

# ## Set up lattice and symmetry
scale_matrix = [2 2; -2 4]
@assert det(scale_matrix) ≈ 12

unitcell = make_unitcell([1 -0.5; 0 0.5*sqrt(3.0)])
addorbital!(unitcell, "A", carte2fract(unitcell, [0.5, 0.5/sqrt(3)]))
addorbital!(unitcell, "B", carte2fract(unitcell, [0.5,-0.5/sqrt(3)]))

lattice = make_lattice(unitcell, scale_matrix)

tsymbed = translation_symmetry_embedding(lattice)
tsym = symmetry(tsymbed)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
psymbed = embed(lattice, psym)

print("Point group: ", psym.hermann_mauguin)



# ## Permutations
println("Translation permutations")
println("------------------------")

for t in tsymbed
    println(t)
end
println()

println("Point permutations")
println("------------------")
for p in psymbed
    println(p)
end
println()


# ## Plots for translation symmetry embeddings

for (i_elem, perm) in enumerate(elements(tsymbed))
    elname = element_name(tsym, i_elem)
    fig = plot(title=elname, aspect=1, size=(200, 250), grid=false, titlefont=Plots.font("sans-serif", pointsize=8))
    orbcoords = []
    orbnames = []

    for iorb in eachindex(lattice.supercell.orbitals)
        orbfc = getorbitalcoord(lattice.supercell, perm(iorb))
        orbcc = fract2carte(lattice.supercell, orbfc)
        push!(orbnames, "$iorb")
        push!(orbcoords, orbcc)
    end
    orbcoords = hcat(orbcoords...)
    L = lattice.supercell.latticevectors
    for R1 in -2:2, R2 in -2:2
        R = [R1, R2]
        LR = L * R
        idx_filt = [i for (i, r) in enumerate(eachcol(orbcoords)) if within(r .+ LR)]
        scatter!(orbcoords[1,idx_filt] .+ LR[1],
                 orbcoords[2,idx_filt] .+ LR[2],
                 color="blue",
                 markerstrokecolor="blue",
                 series_annotations=[Plots.text(x, 6, :left, :bottom) for x in orbnames[idx_filt]], label=nothing)
    end
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(fig, "example_honeycomb_symmetry/translation_symmetry-$i_elem.svg")
end

# ![](example_honeycomb_symmetry/translation_symmetry-1.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-2.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-3.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-4.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-5.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-6.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-7.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-8.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-9.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-10.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-11.svg)
# ![](example_honeycomb_symmetry/translation_symmetry-12.svg)


# ## Plots for point symmetry embeddings

for (i_elem, perm) in enumerate(elements(psymbed))
    elname = element_name(psym, i_elem)
    fig = plot(title=simplifyname(elname), aspect=1, size=(200, 250), grid=false, titlefont=Plots.font("sans-serif", pointsize=8))
    orbcoords = []
    orbnames = []

    for iorb in eachindex(lattice.supercell.orbitals)
        orbfc = getorbitalcoord(lattice.supercell, perm(iorb))
        orbcc = fract2carte(lattice.supercell, orbfc)
        push!(orbnames, "$iorb")
        push!(orbcoords, orbcc)
    end
    orbcoords = hcat(orbcoords...)
    L = lattice.supercell.latticevectors
    for R1 in -2:2, R2 in -2:2
        R = [R1, R2]
        LR = L * R
        idx_filt = [i for (i, r) in enumerate(eachcol(orbcoords)) if within(r .+ LR)]
        scatter!(orbcoords[1,idx_filt] .+ LR[1],
                 orbcoords[2,idx_filt] .+ LR[2],
                 color="blue",
                 markerstrokecolor="blue",
                 series_annotations=[Plots.text(x, 6, :left, :bottom) for x in orbnames[idx_filt]], label=nothing)
    end
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(fig, "example_honeycomb_symmetry/point_symmetry-$i_elem.svg")
end

# ![](example_honeycomb_symmetry/point_symmetry-1.svg)
# ![](example_honeycomb_symmetry/point_symmetry-2.svg)
# ![](example_honeycomb_symmetry/point_symmetry-3.svg)
# ![](example_honeycomb_symmetry/point_symmetry-4.svg)
# ![](example_honeycomb_symmetry/point_symmetry-5.svg)
# ![](example_honeycomb_symmetry/point_symmetry-6.svg)
# ![](example_honeycomb_symmetry/point_symmetry-7.svg)
# ![](example_honeycomb_symmetry/point_symmetry-8.svg)
# ![](example_honeycomb_symmetry/point_symmetry-9.svg)
# ![](example_honeycomb_symmetry/point_symmetry-10.svg)
# ![](example_honeycomb_symmetry/point_symmetry-11.svg)
# ![](example_honeycomb_symmetry/point_symmetry-12.svg)