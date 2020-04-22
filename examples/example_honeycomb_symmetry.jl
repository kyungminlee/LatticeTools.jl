using TightBindingLattice

using YAML
using LinearAlgebra
using Plots


scale_matrix = [4 -2; 2 2]
@assert det(scale_matrix) ≈ 12

unitcell = make_unitcell([1 -0.5; 0 0.5*sqrt(3.0)])
addorbital!(unitcell, "A", carte2fract(unitcell, [0.5, 0.5/sqrt(3)]))
addorbital!(unitcell, "B", carte2fract(unitcell, [0.5,-0.5/sqrt(3)]))

lattice = make_lattice(unitcell, scale_matrix)

tsymbed = translation_symmetry_embedding(lattice)
tsym = symmetry(tsymbed)
psym = project(PointSymmetryDatabase.get(19), [1 0 0; 0 1 0])
psymbed = embed(lattice, psym)

println("# Translation permutations")
for t in tsymbed
    println(t)
end
println()

println("# Point permutations")
for p in psymbed
    println(p)
end
println()


# ## Plots

extent = [-2.5, 2.5, -2.5, 2.5]
within(r) = (extent[1] <= r[1] <= extent[2] && extent[3] <= r[2] <= extent[4])

println("# Plotting translation symmetry embeddings")

for (i_elem, perm) in enumerate(elements(tsymbed))
    elname = element_name(tsym, i_elem)
    fig = plot(title=elname, aspect=1, size=(500, 500), grid=false)
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
                 series_annotations=[Plots.text(x, 12, :left, :bottom) for x in orbnames[idx_filt]], label=nothing)
    end
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(fig, "translation_symmetry-$i_elem.html")
end

println("# Plotting point symmetry embeddings")

for (i_elem, perm) in enumerate(elements(psymbed))
    elname = element_name(psym, i_elem)
    fig = plot(title=elname, aspect=1, size=(500, 500), grid=false)
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
                 series_annotations=[Plots.text(x, 12, :left, :bottom) for x in orbnames[idx_filt]], label=nothing)
    end
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(fig, "point_symmetry-$i_elem.html")
end