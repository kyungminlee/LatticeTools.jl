# # Little Symmetry

using LinearAlgebra
using Plots
using TightBindingLattice

simplifyname(s::AbstractString) = (
                    s |> (x-> replace(x, r"<sup>(.+?)</sup>" => s"\1"))
                      |> (x-> replace(x, r"<sub>(.+?)</sub>" => s"[\1]"))
                )
function makewithin(extent::AbstractVector{<:Real})
    a, b, c, d = extent
    (x::Real, y::Real) -> ((a <= x <= b) && (c <= y <= d))
end
mkpath("example_little_symmetry_kspace_square")


# ## Honeycome lattice

# ### Define Unitcell
latticevectors = [1.0 0.0; 0.0 1.0];
unitcell = make_unitcell(latticevectors; OrbitalType=String);
addsite!(unitcell, "A", carte2fract(unitcell, [0.0, 0.0]));

# # ### Plot lattice
# let bravais_lattice = [], sites = []
#     extent = [-2, 2, -2, 2]
#     within = makewithin(extent)
#     for i1 in -5:5, i2 in -5:5
#         R = latticevectors * [i1, i2]
#         within(R...) && push!(bravais_lattice, R)
#         for (orb_name, orb_fc) in unitcell.sites
#             orb_cc = fract2carte(unitcell, orb_fc)
#             r = R + orb_cc
#             within(r...) && push!(sites, r)
#         end
#     end
#     bravais_lattice = hcat(bravais_lattice...)
#     sites = hcat(sites...)
    
#     img = plot(aspect_ratio=1)
#     scatter!(bravais_lattice[1,:], bravais_lattice[2,:], color="black", markershape=:circle, markersize=3, label="Bravais")
#     scatter!(sites[1,:], sites[2,:], color="red", markerstrokewidth=0, markersize=5, label="A")
#     xlims!(extent[1], extent[2])
#     ylims!(extent[3], extent[4])
#     savefig(img, "example_little_symmetry_kspace_square/realspace.svg")
# end

# ![](example_little_symmetry_kspace_square/realspace.svg)

# ## 4 × 4 supercell

# ### Make Superlattice
shape = [4 0; 0 4]
lattice = make_lattice(unitcell, shape)

# ### Plot lattice
let bravais_lattice = [], sites = []
    extent = [-2, 6, -2, 6]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.sites
            orb_cc = fract2carte(lattice.supercell, orb_fc)
            r = R + orb_cc
            within(r...) && push!(sites, r)
        end
    end
    bravais_lattice = hcat(bravais_lattice...)
    sites = hcat(sites...)
    
    img = plot(aspect_ratio=1)
    scatter!(bravais_lattice[1,:], bravais_lattice[2,:],
             markershape=:circle, markersize=7,
             markercolor=RGBA(1,1,1,0), markerstrokecolor=RGBA(1,0,0,1),
             markerstrokewidth=2,
             label="Bravais (supercell)",
            )
    scatter!(sites[1,:], sites[2,:], color="black", markerstrokewidth=0, markersize=4, label="Sites")
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(img, "example_little_symmetry_kspace_square/realspace-4x4.svg")
end

# ![](example_little_symmetry_kspace_square/realspace-4x4.svg)


# ### Symmetries
tsym = TranslationSymmetry(lattice)
# psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
psym = PointSymmetryDatabase.find2d("4mm")
println("point symmetry (Hermann Mauguin): $(psym.hermann_mauguin)")
println("translation symmetry compatible with point symmetry? ", iscompatible(tsym, psym))

println("irreducible representations:")
for tsym_irrep_index in 1:num_irreps(tsym)
    kf = tsym.fractional_momenta[tsym_irrep_index]
    kc = lattice.unitcell.reducedreciprocallatticevectors * kf
    psym_little = little_symmetry(tsym, tsym_irrep_index, psym)
    print("- irrep index: $tsym_irrep_index\n")
    print("  fractional momentum: [")
    for (i, x) in enumerate(kf)
        i != 1 && print(", ")
        print(x)
    end
    print("]\n")
    print("  momentum: [")
    for (i, x) in enumerate(kc)
        i != 1 && print(", ")
        print(x)
    end
    print("]\n")
    print("  is point symmetry compatible with momentum: $(iscompatible(tsym, tsym_irrep_index, psym))\n")
    print("  little point symmetry: $(psym_little.hermann_mauguin)\n")
end


# ### Plot momentum space
let
    extent = [-1.3, 1.3, -1.3, 1.3]
    within = makewithin(extent)

    reciprocallatticepoints = []
    kpoints = []
    littlegroupnames = String[]
    generatornames = String[]
    kpointnames = String[]

    for tsym_irrep_index in 1:num_irreps(tsym)
        psym_little = little_symmetry(tsym, tsym_irrep_index, psym)
        kf = tsym.fractional_momenta[tsym_irrep_index]
        kc = lattice.unitcell.reducedreciprocallatticevectors * kf
        for i1 in -2:2, i2 in -2:2
            G = lattice.unitcell.reducedreciprocallatticevectors * [i1, i2]
            within(G...) && push!(reciprocallatticepoints, G)
            k = G + kc
            if within(k...)
                push!(kpoints, k)
                push!(littlegroupnames, psym_little.hermann_mauguin)
                push!(generatornames, join(simplifyname.(psym_little.element_names[psym_little.generators]), "\n"))
                push!(kpointnames, "$tsym_irrep_index")
            end
        end
    end
    reciprocallatticepoints = hcat(reciprocallatticepoints...)
    kpoints = hcat(kpoints...)

    img = plot(size=(800, 800), xlims=(extent[1]-0.1, extent[2]+0.1), ylims=(extent[3]-0.1, extent[4]+0.1), aspect_ratio=1)
    scatter!(reciprocallatticepoints[1,:], reciprocallatticepoints[2,:],
             markershape=:circle, markersize=6, markercolor=RGBA(1,1,1,0), markerstrokecolor=RGBA(1,0,0,1), label="reciprocal lattice")
    scatter!(kpoints[1,:], kpoints[2,:],
             markershape=:circle, markersize=4, markercolor=RGBA(0,0,1,0.5), markerstrokecolor=RGBA(0,0,1,0.5),
             series_annotations=[Plots.text("k[$x]", 8, :left, :bottom) for x in kpointnames], label="momentum points")
    scatter!(kpoints[1,:], kpoints[2,:],
             markershape=:circle, markersize=0, color=RGBA(0,0,0,0),
             series_annotations=[Plots.text("$x ", 8, :right, :top) for x in littlegroupnames], label="")
    scatter!(kpoints[1,:], kpoints[2,:],
             markershape=:circle, markersize=0, color=RGBA(0,0,0,0),
             series_annotations=[Plots.text("$x", 6, :left, :top) for x in generatornames], label="")
    savefig(img, "example_little_symmetry_kspace_square/momentumspace-4x4.svg")
end

# ![](example_little_symmetry_kspace_square/momentumspace-4x4.svg)
