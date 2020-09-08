using Luxor
using JSON
using LatticeTools

include("Kagome.jl")

lattice_data = JSON.parsefile("lattice-3.json")
sites = [(sitetype, [R...]) for (sitetype, R) in lattice_data["sites"]]
kagome = Kagome.make_kagome_obc(sites)
site_coordinates = []
site_types = String[]
for ((sitetype, sitesubscript), sitefc) in kagome.unitcell.sites
    sitecc = fract2carte(kagome.unitcell, sitefc)
    push!(site_coordinates, sitecc)
    push!(site_types, sitetype)
end
site_coordinates = hcat(site_coordinates...)
site_coordinates[2,:] = -site_coordinates[2,:]

Drawing(800, 300, "logo.png")

origin()
translate(-450, 130)

let ratio = 100, radius = 13, radius_large = 18, offset = (0,0),
    sitecolor = Dict("A" => "#389826", "B" => "#9558B2", "C" => "#CB3C33")

    sethue("gray")
    setline(2)
    for (i, j) in kagome.nearest_neighbor_bonds
        ri = site_coordinates[:,i]
        rj = site_coordinates[:,j]
        line(
            Point(ri[1]*ratio + offset[1], ri[2]*ratio + offset[2]),
            Point(rj[1]*ratio + offset[1], rj[2]*ratio + offset[2]),
            :stroke
        )
    end
    for (sitetype, r) in zip(site_types, eachcol(site_coordinates))
        sethue("white")
        circle(Point(r[1]*ratio + offset[1], r[2]*ratio + offset[2]), radius_large, :fill)
        sethue(sitecolor[sitetype])
        circle(Point(r[1]*ratio + offset[1], r[2]*ratio + offset[2]), radius, :fill)
    end
end

translate(360, -160)
sethue("black")
fontsize(144)
fontface("Tamil MN")
texttrack("Lattice", Point(0, 0), -50, 144)
translate(0, 140)
texttrack("Tools.jl", Point(0, 0), -50, 144)

setcolor(sethue("#4063D8")..., 1.0)
#ellipse(286.3, -239.6, 25, 22.6, :fill)
ellipse(286, -238, 32, 32, :fill)

finish()
#preview()
