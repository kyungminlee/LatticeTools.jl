using Luxor
using JSON
using LatticeTools

include("Kagome.jl")

lattice_data = JSON.parsefile("lattice.json")
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

Drawing(840, 300, "logo.png")

origin()
translate(-460, 110)

let ratio = 88, radius = 10, radius_large = 14, offset = (0,0)
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

    #sethue("white")
    #for sc in values(site_coordinates), r in eachcol(sc)
    #    circle(Point(r[1]*ratio + offset[1], r[2]*ratio + offset[2]), radius_large, :fill)
    #end
    sitecolor = Dict("A" => "brown3", "B" => "forestgreen", "C" => "mediumorchid3")
    for (sitetype, r) in zip(site_types, eachcol(site_coordinates))
        sethue("white")
        circle(Point(r[1]*ratio + offset[1], r[2]*ratio + offset[2]), radius_large, :fill)
        sethue(sitecolor[sitetype])
        circle(Point(r[1]*ratio + offset[1], r[2]*ratio + offset[2]), radius, :fill)
    end
end

translate(400, -136)
sethue("black")
fontsize(144)
fontface("Tamil MN")
texttrack("Lattice", Point(0, 0), -50, 144)
translate(0, 140)
texttrack("Tools.jl", Point(0, 0), -50, 144)

setcolor(sethue("royalblue")..., 1.0)
ellipse(286.3, -239.6, 25, 22.6, :fill)

finish()
#preview()
