# # Little Symmetry

using LinearAlgebra
using Plots
pyplot()
import DisplayAs

using TightBindingLattice

latexify(s::AbstractString) = "\$" * (
                    s |> (x-> replace(x, r"<sup>(.+?)</sup>" => s"^{\1}"))
                      |> (x-> replace(x, r"<sub>(.+?)</sub>" => s"_{\1}"))
                      |> (x-> replace(x, r"-(\d)" => s"\\overline{\1}"))
                )* "\$"

function makewithin(extent::AbstractVector{<:Real})
    a, b, c, d = extent
    (x::Real, y::Real) -> ((a <= x <= b) && (c <= y <= d))
end

# ## Honeycome lattice

# ### Define Unitcell
latticevectors = [1 -0.5; 0 sqrt(3)*0.5];
unitcell = make_unitcell(latticevectors; OrbitalType=String);
addorbital!(unitcell, "A", carte2fract(unitcell, [0.5, 0.5 / sqrt(3)]));
addorbital!(unitcell, "B", carte2fract(unitcell, [0.5, -0.5 / sqrt(3)]));

# ### Plot lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-2, 2, -2, 2]
    within = makewithin(extent)
    for i1 in -5:5, i2 in -5:5
        R = latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in unitcell.orbitals
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
    img
end

# ## √3 × √3 supercell

# ### Make Superlattice
shape = [2 -1; 1 1]
lattice = make_lattice(unitcell, shape)

# ### Plot lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-4, 4, -4, 4]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.orbitals
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
    img
end


# ### Symmetries
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
println("point symmetry (Hermann Mauguinn): $(psym.hermann_mauguinn)")
println("translation symmetry compatible with point symmetry? ", iscompatible(tsym, psym))

println("irreducible represnetations:")
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
    print("  little point symmetry: $(psym_little.hermann_mauguinn)\n")
end


# ### Plot momentum space
let
    extent = [-1.2, 1.2, -1.6, 1.6]
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
                push!(littlegroupnames, psym_little.hermann_mauguinn)
                push!(generatornames, join(latexify.(psym_little.element_names[psym_little.generators]), "\n"))
                push!(kpointnames, "$tsym_irrep_index")
            end
        end
    end
    reciprocallatticepoints = hcat(reciprocallatticepoints...)
    kpoints = hcat(kpoints...)

    img = plot(aspect_ratio=1, size=(300, 400))
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=4, color="blue",markerstrokecolor="blue", 
             series_annotations=[Plots.text("$x ", 8, :right, :top) for x in littlegroupnames], label="momentum points")
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=0,
             series_annotations=[Plots.text("$x", 8, :left, :top) for x in generatornames], label="")
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=0,
             series_annotations=[Plots.text("\$\\mathbf{k}_{$x}\$", 12, :left, :bottom) for x in kpointnames], label="")
    scatter!(reciprocallatticepoints[1,:], reciprocallatticepoints[2,:], markershape=:circle, markersize=6, color="red", markerstrokecolor="red", label="reciprocal lattice")
    xlims!(extent[1] - 0.1, extent[2] + 0.1)
    ylims!(extent[3] - 0.1, extent[4] + 0.1)
    savefig(img, "momentumspace-root3xroot3.pdf")
end



# ## 2√3 × 2√3 supercell

# ### Make Superlattice
shape = [4 -2; 2 2]
lattice = make_lattice(unitcell, shape)

# ### Plot lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-8, 8, -8, 8]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.orbitals
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
end


# ### Symmetries
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
println("point symmetry (Hermann Mauguinn): $(psym.hermann_mauguinn)")
println("translation symmetry compatible with point symmetry? ", iscompatible(tsym, psym))

println("irreducible represnetations:")
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
    print("  little point symmetry: $(psym_little.hermann_mauguinn)\n")
end


# ### Plot momentum space
let
    extent = [-1.2, 1.2, -1.6, 1.6]
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
                push!(littlegroupnames, psym_little.hermann_mauguinn)
                push!(generatornames, join(latexify.(psym_little.element_names[psym_little.generators]), "\n"))
                push!(kpointnames, "$tsym_irrep_index")
            end
        end
    end
    reciprocallatticepoints = hcat(reciprocallatticepoints...)
    kpoints = hcat(kpoints...)

    img = plot(aspect_ratio=1, size=(600, 800))
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=4, color="blue",markerstrokecolor="blue", 
             series_annotations=[Plots.text("$x ", 8, :right, :top) for x in littlegroupnames], label="momentum points")
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=0,
             series_annotations=[Plots.text("$x", 8, :left, :top) for x in generatornames], label="")
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=0,
             series_annotations=[Plots.text("\$\\mathbf{k}_{$x}\$", 12, :left, :bottom) for x in kpointnames], label="")
    scatter!(reciprocallatticepoints[1,:], reciprocallatticepoints[2,:], markershape=:circle, markersize=6, color="red", markerstrokecolor="red", label="reciprocal lattice")
    xlims!(extent[1]-0.1, extent[2]+0.1)
    ylims!(extent[3]-0.1, extent[4]+0.1)
    savefig(img, "momentumspace-2root3x2root3.pdf")
end


# ## 6 × 6 supercell

# ### Make Superlattice
shape = [6 0; 0 6]
lattice = make_lattice(unitcell, shape)

# ### Plot lattice
let bravais_lattice = [], a_sites = [], b_sites = []
    extent = [-8, 8, -8, 8]
    within = makewithin(extent)

    for i1 in -5:5, i2 in -5:5
        R = lattice.supercell.latticevectors * [i1, i2]
        within(R...) && push!(bravais_lattice, R)
        for (orb_name, orb_fc) in lattice.supercell.orbitals
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
    img
end


# ### Symmetries
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
println("point symmetry (Hermann Mauguinn): $(psym.hermann_mauguinn)")
println("translation symmetry compatible with point symmetry? ", iscompatible(tsym, psym))

println("irreducible represnetations:")
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
    print("  little point symmetry: $(psym_little.hermann_mauguinn)\n")
end


# ### Plot momentum space
let
    extent = [-1.2, 1.2, -1.6, 1.6]
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
                push!(littlegroupnames, psym_little.hermann_mauguinn)
                push!(generatornames, join(latexify.(psym_little.element_names[psym_little.generators]), "\n"))
                push!(kpointnames, "$tsym_irrep_index")
            end
        end
    end
    reciprocallatticepoints = hcat(reciprocallatticepoints...)
    kpoints = hcat(kpoints...)

    img = plot(aspect_ratio=1, size=(900, 1200))
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=4, color="blue",markerstrokecolor="blue", 
             series_annotations=[Plots.text("$x ", 8, :right, :top) for x in littlegroupnames], label="momentum points")
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=0,
             series_annotations=[Plots.text("$x", 8, :left, :top) for x in generatornames], label="")
    scatter!(kpoints[1,:], kpoints[2,:], markershape=:circle, markersize=0,
             series_annotations=[Plots.text("\$\\mathbf{k}_{$x}\$", 12, :left, :bottom) for x in kpointnames], label="")
    scatter!(reciprocallatticepoints[1,:], reciprocallatticepoints[2,:], markershape=:circle, markersize=6, color="red", markerstrokecolor="red", label="reciprocal lattice")
    xlims!(extent[1] - 0.1, extent[2] + 0.1)
    ylims!(extent[3] - 0.1, extent[4] + 0.1)
    savefig(img, "momentumspace-6x6.pdf")
end

