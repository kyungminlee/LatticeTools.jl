# # Honeycomb Lattice in k-space

using LinearAlgebra
using Plots
using TightBindingLattice

simplifyname(s::AbstractString) = (
    s
    |> (x-> replace(x, r"<sup>(.+?)</sup>" => s"\1"))
    |> (x-> replace(x, r"<sub>(.+?)</sub>" => s"[\1]"))
)

function makewithin(extent::AbstractVector{<:Real})
    a, b, c, d = extent
    (x::Real, y::Real) -> ((a <= x <= b) && (c <= y <= d))
end
mkpath("example_little_symmetry_kspace_honeycomb")


# ## Honeycome Lattice

# ### Define UnitCell
latticevectors = [1 -0.5; 0 sqrt(3)*0.5];
unitcell = makeunitcell(latticevectors; SiteType=String);
addsite!(unitcell, "A", carte2fract(unitcell, [0.5, 0.5 / sqrt(3)]));
addsite!(unitcell, "B", carte2fract(unitcell, [0.5, -0.5 / sqrt(3)]));

# ### Plot Lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-2, 2, -2, 2]
    within = makewithin(extent)
    for i1 in -5:5, i2 in -5:5
        R = latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in unitcell.sites
            orb_cc = fract2carte(unitcell, orb_fc)
            r = R + orb_cc
            within(r...) && push!(orb_name == "A" ? a_sites : b_sites, r)
        end
    end
    bravais_lattice = hcat(bravais_lattice...)
    a_sites = hcat(a_sites...)
    b_sites = hcat(b_sites...)

    img = plot(aspect_ratio=1)
    scatter!(bravais_lattice[1,:], bravais_lattice[2,:], color="black", markershape=:circle, markersize=3, label="Bravais")
    scatter!(a_sites[1,:], a_sites[2,:], color="red", markerstrokewidth=0, markersize=5, label="A")
    scatter!(b_sites[1,:], b_sites[2,:], color="blue", markerstrokewidth=0, markersize=5, label="B")
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(img, "example_little_symmetry_kspace_honeycomb/realspace.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/realspace.svg)

# ## √3 × √3 Supercell

# ### Make Superlattice
shape = [2 -1; 1 1]
lattice = makelattice(unitcell, shape)

# ### Plot Lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-4, 4, -4, 4]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.sites
            orb_cc = fract2carte(lattice.supercell, orb_fc)
            r = R + orb_cc
            within(r...) && push!(orb_name[1] == "A" ? a_sites : b_sites, r)
        end
    end
    bravais_lattice = hcat(bravais_lattice...)
    a_sites = hcat(a_sites...)
    b_sites = hcat(b_sites...)

    img = plot(aspect_ratio=1)
    scatter!(bravais_lattice[1,:], bravais_lattice[2,:], color="black", markershape=:circle, markersize=3, label="Bravais (supercell)")
    scatter!(a_sites[1,:], a_sites[2,:], color="red", markerstrokewidth=0, markersize=5, label="A")
    scatter!(b_sites[1,:], b_sites[2,:], color="blue", markerstrokewidth=0, markersize=5, label="B")
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(img, "example_little_symmetry_kspace_honeycomb/realspace-root3xroot3.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/realspace-root3xroot3.svg)


# ### Symmetries
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
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


# ### Plot Momentum Space
let
    extent = [-1.1, 1.1, -1.4, 1.4]
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

    img = plot(size=(400, 500), xlims=(extent[1]-0.1, extent[2]+0.1), ylims=(extent[3]-0.1, extent[4]+0.1), aspect_ratio=1)
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
    savefig(img, "example_little_symmetry_kspace_honeycomb/momentumspace-root3xroot3.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/momentumspace-root3xroot3.svg)

# ## 2√3 × 2√3 Supercell

# ### Make Superlattice
shape = [2 2; -2 4]
lattice = makelattice(unitcell, shape)

# ### Plot Lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-8, 8, -8, 8]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.sites
            orb_cc = fract2carte(lattice.supercell, orb_fc)
            r = R + orb_cc
            within(r...) && push!(orb_name[1] == "A" ? a_sites : b_sites, r)
        end
    end
    bravais_lattice = hcat(bravais_lattice...)
    a_sites = hcat(a_sites...)
    b_sites = hcat(b_sites...)

    img = plot(aspect_ratio=1)
    scatter!(bravais_lattice[1,:], bravais_lattice[2,:], color="black", markershape=:circle, markersize=5, label="Bravais (supercell)")
    scatter!(a_sites[1,:], a_sites[2,:], color="red", markerstrokewidth=0, markersize=3, label="A")
    scatter!(b_sites[1,:], b_sites[2,:], color="blue", markerstrokewidth=0, markersize=3, label="B")
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(img, "example_little_symmetry_kspace_honeycomb/realspace-2root3x2root3.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/realspace-2root3x2root3.svg)


# ### Symmetries
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
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


# ### Plot Momentum Space
let
    extent = [-1.2, 1.2, -1.6, 1.6]
    extent = [-1.1, 1.1, -1.4, 1.4]
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

    img = plot(size=(600, 750), xlims=(extent[1]-0.1, extent[2]+0.1), ylims=(extent[3]-0.1, extent[4]+0.1), aspect_ratio=1)
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
    savefig(img, "example_little_symmetry_kspace_honeycomb/momentumspace-2root3x2root3.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/momentumspace-2root3x2root3.svg)

# ## 6 × 6 Supercell

# ### Make Superlattice
shape = [6 0; 0 6]
lattice = makelattice(unitcell, shape)

# ### Plot Lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-8, 8, -8, 8]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.sites
            orb_cc = fract2carte(lattice.supercell, orb_fc)
            r = R + orb_cc
            within(r...) && push!(orb_name[1] == "A" ? a_sites : b_sites, r)
        end
    end
    bravais_lattice = hcat(bravais_lattice...)
    a_sites = hcat(a_sites...)
    b_sites = hcat(b_sites...)

    img = plot(aspect_ratio=1)
    scatter!(bravais_lattice[1,:], bravais_lattice[2,:], color="black", markershape=:circle, markersize=5, label="Bravais (supercell)")
    scatter!(a_sites[1,:], a_sites[2,:], color="red", markerstrokewidth=0, markersize=3, label="A")
    scatter!(b_sites[1,:], b_sites[2,:], color="blue", markerstrokewidth=0, markersize=3, label="B")
    xlims!(extent[1], extent[2])
    ylims!(extent[3], extent[4])
    savefig(img, "example_little_symmetry_kspace_honeycomb/realspace-6x6.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/realspace-6x6.svg)


# ### Symmetries
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
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


# ### Plot Momentum Space
let
    extent = [-1.2, 1.2, -1.6, 1.6]
    extent = [-1.1, 1.1, -1.4, 1.4]
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

    img = plot(size=(800, 1000), xlims=(extent[1]-0.1, extent[2]+0.1), ylims=(extent[3]-0.1, extent[4]+0.1), aspect_ratio=1)
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
    savefig(img, "example_little_symmetry_kspace_honeycomb/momentumspace-6x6.svg")
end

# ![](example_little_symmetry_kspace_honeycomb/momentumspace-6x6.svg)
