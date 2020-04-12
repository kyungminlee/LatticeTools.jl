using LinearAlgebra
using Printf
using Plots

using TightBindingLattice

include("Kagome.jl")

kagome = make_kagome_lattice([4 -2; 2 2])
#kagome = make_kagome_lattice([2 -1; 1 1])

#=
kagome.unitcell
kagome.lattice
kagome.translation_symmetry
kagome.point_symmetry
kagome.nearest_neighbor_bonds
kagome.next_nearest_neighbor_bonds
kagome.nearest_neighbor_triangles
=#

println("Point group: ", kagome.point_symmetry.hermann_mauguinn)


# ## Orbital map
orbital_map = findorbitalmap(kagome.lattice.unitcell, kagome.point_symmetry)

println("=== Orbital map ===")
for (n, map) in zip(kagome.point_symmetry.element_names, orbital_map)
    @printf("%32s:", n)
    for (i_elem, (j_elem, R)) in enumerate(map)
        print("\t$i_elem ↦ $j_elem, $R")
    end
    println()
end
println()


permutations_direct = Permutation[]
for (n, mat, map) in zip(kagome.point_symmetry.element_names,
                         kagome.point_symmetry.matrix_representations,
                         orbital_map)
    push!(permutations_direct, get_orbital_permutation(kagome.lattice, mat, map))
end
permutations = get_orbital_permutations(kagome.lattice, kagome.point_symmetry)

@show permutations_direct == permutations




# ## Angle

@printf("%32s", "name")
for (iorb, (orbname, orbfc)) in enumerate(kagome.unitcell.orbitals)
    @printf("\t%9s", orbname)
end
println()
println("---------------------------------------------------------------------------------")
for (n, m) in zip(kagome.point_symmetry.element_names, orbital_map)
    @printf("%32s", n)
    for (i_elem, (j_elem, R)) in enumerate(m)
        ri = fract2carte(kagome.lattice.unitcell, getorbitalcoord(kagome.lattice.unitcell, i_elem))
        rj = fract2carte(kagome.lattice.unitcell, getorbitalcoord(kagome.lattice.unitcell, j_elem) + R)
        # Rotation angle
        θi, θj = atan(ri[2], ri[1]), atan(rj[2], rj[1])
        Δθij = round(Int, mod((θj-θi) * 180 / π, 360))
        @printf("\t%3.0f", Δθij)

        # Mirror direction
        r = ri + rj
        ρ = ri - rj
        if !isapprox(norm(r), 0)
            ϕ = round(Int, mod(atan(r[2], r[1]) * 180 / π + 90, 180))
        elseif !isapprox(norm(ρ), 0)
            ϕ = round(Int, mod(atan(ρ[2], ρ[1]) * 180 / π, 180))
        else
            ϕ = NaN
        end
        @printf(" | %3.0f", ϕ)
    end
    println()
end


# ## Point Symmetry

permutations = get_orbital_permutations(kagome.lattice, kagome.point_symmetry)
extent = [-2, 2, -2, 2]
within = (r) -> extent[1] <= r[1] <= extent[2] && extent[3] <= r[2] <= extent[4]

for (i_elem, perm) in enumerate(permutations)
    elname = element_name(kagome.point_symmetry, i_elem)
    fig = plot(title=elname, aspect=1, size=(300, 300), grid=false)
    orbcoords = []
    orbnames = []

    for iorb in eachindex(kagome.lattice.supercell.orbitals)
        orbfc = getorbitalcoord(kagome.lattice.supercell, perm.map[iorb])
        orbcc = fract2carte(kagome.lattice.supercell, orbfc)
        push!(orbnames, "$iorb")
        push!(orbcoords, orbcc)
    end
    orbcoords = hcat(orbcoords...)
    L = kagome.lattice.supercell.latticevectors
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
    xlims!(-2, 2)
    ylims!(-2, 2)
    # savefig(fig, "point_symmetry-$i_elem.png")
end




# ## Translation Symmetry

permutations = get_orbital_permutations(kagome.lattice, kagome.translation_symmetry)
extent = [-2, 2, -2, 2]
within = (r) -> extent[1] <= r[1] <= extent[2] && extent[3] <= r[2] <= extent[4]

for (i_elem, perm) in enumerate(permutations)
    elname = element_name(kagome.translation_symmetry, i_elem)
    fig = plot(title=elname, aspect=1, size=(300, 300), grid=false)
    orbcoords = []
    orbnames = []

    for iorb in eachindex(kagome.lattice.supercell.orbitals)
        orbfc = getorbitalcoord(kagome.lattice.supercell, perm.map[iorb])
        orbcc = fract2carte(kagome.lattice.supercell, orbfc)
        push!(orbnames, "$iorb")
        push!(orbcoords, orbcc)
    end
    orbcoords = hcat(orbcoords...)
    L = kagome.lattice.supercell.latticevectors
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
    xlims!(-2, 2)
    ylims!(-2, 2)
    #savefig(fig, "translation_symmetry-$i_elem.png")
end
