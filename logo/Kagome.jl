module Kagome

using LatticeTools


NN_BOND_TYPES = [
    ([ 0, 0], "A", [ 0, 0], "B", 1),
    ([ 0, 0], "B", [ 0, 0], "C", 1),
    ([ 0, 0], "C", [ 0, 0], "A", 1),

    ([ 1, 1], "A", [ 1, 0], "B",-1),
    ([ 1, 0], "B", [ 0, 0], "C",-1),
    ([ 0, 0], "C", [ 1, 1], "A",-1),
]


NNN_BOND_TYPES = [
    ([ 0, 0], "A", [ 1, 0], "B", 1), # ◁
    ([ 1, 0], "B", [ 0,-1], "C", 1),
    ([ 0,-1], "C", [ 0, 0], "A", 1),
    ([ 0, 0], "C", [ 1, 0], "A",-1),
    ([ 1, 0], "A", [ 0,-1], "B",-1), # ▷
    ([ 0,-1], "B", [ 0, 0], "C",-1),
]


export make_kagome_lattice
export make_kagome_obc


function make_kagome_lattice(size_matrix ::AbstractMatrix{<:Integer}; compute_symmetry::Bool=false)
    latticevectors = [1 -0.5; 0 0.5*sqrt(3.0)];
    unitcell = makeunitcell(latticevectors, SiteType=String)
    addsite!(unitcell, "A", carte2fract(unitcell, [0.5, 0.0]))
    addsite!(unitcell, "B", carte2fract(unitcell, [0.25, 0.25*sqrt(3.0)]))
    addsite!(unitcell, "C", carte2fract(unitcell, [0.5+0.25, 0.25*sqrt(3.0)]))

    lattice = make_lattice(unitcell, size_matrix)
    hypercube = lattice.hypercube
    supercell = lattice.supercell
    tsym = FiniteTranslationSymmetry(lattice)
    psym = little_symmetry(tsym, project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0]))
    tsymbed = embed(lattice, tsym)
    psymbed = embed(lattice, psym)
    ssymbed = SymmorphicSymmetry(tsymbed, psymbed)

    nnbonds = []
    nnnbonds = []

    for r in lattice.bravais_coordinates
        for (rowvec, rowsite, colvec, colsite, bondsign) in NN_BOND_TYPES
            R_row, r_row = hypercube.wrap(r .+ rowvec)
            R_col, r_col = hypercube.wrap(r .+ colvec)
            rowsite_super = (rowsite, r_row)
            colsite_super = (colsite, r_col)
            irow = get(supercell.siteindices, rowsite_super, -1)
            icol = get(supercell.siteindices, colsite_super, -1)
            push!(nnbonds, ((irow, icol), R_col-R_row, bondsign))
        end
        for (rowvec, rowsite, colvec, colsite, bondsign) in NNN_BOND_TYPES
            R_row, r_row = hypercube.wrap(r .+ rowvec)
            R_col, r_col = hypercube.wrap(r .+ colvec)
            rowsite_super = (rowsite, r_row)
            colsite_super = (colsite, r_col)
            irow = get(supercell.siteindices, rowsite_super, -1)
            icol = get(supercell.siteindices, colsite_super, -1)
            push!(nnnbonds, ((irow, icol), R_col-R_row, bondsign))
        end
    end

    nn_triangles = []
    for r in lattice.bravais_coordinates
      triangle = []
      for (rowvec, roworb, colvec, colorb, bondsign) in NN_BOND_TYPES[1:3]
        R_row, r_row = hypercube.wrap(r .+ rowvec)
        R_col, r_col = hypercube.wrap(r .+ colvec)
        roworb_super = (roworb, r_row)
        colorb_super = (colorb, r_col)
        irow = get(supercell.siteindices, roworb_super, -1)
        icol = get(supercell.siteindices, colorb_super, -1)
        push!(triangle, ((irow, icol), R_col-R_row))
      end
      push!(nn_triangles, (triangle, 1))

      triangle = []
      for (rowvec, roworb, colvec, colorb, bondsign) in NN_BOND_TYPES[4:6]
        R_row, r_row = hypercube.wrap(r .+ rowvec)
        R_col, r_col = hypercube.wrap(r .+ colvec)
        roworb_super = (roworb, r_row)
        colorb_super = (colorb, r_col)
        irow = get(supercell.siteindices, roworb_super, -1)
        icol = get(supercell.siteindices, colorb_super, -1)
        push!(triangle, ((irow, icol), R_col-R_row))
      end
      push!(nn_triangles, (triangle, -1))
    end


    nnn_triangles = []
    for r in lattice.bravais_coordinates
      triangle = []
      for (rowvec, roworb, colvec, colorb, bondsign) in NNN_BOND_TYPES[1:3]
        R_row, r_row = hypercube.wrap(r .+ rowvec)
        R_col, r_col = hypercube.wrap(r .+ colvec)
        roworb_super = (roworb, r_row)
        colorb_super = (colorb, r_col)
        irow = get(supercell.siteindices, roworb_super, -1)
        icol = get(supercell.siteindices, colorb_super, -1)
        push!(triangle, ((irow, icol), R_col-R_row))
      end
      push!(nnn_triangles, (triangle, 1))

      triangle = []
      for (rowvec, roworb, colvec, colorb, bondsign) in NNN_BOND_TYPES[4:6]
        R_row, r_row = hypercube.wrap(r .+ rowvec)
        R_col, r_col = hypercube.wrap(r .+ colvec)
        roworb_super = (roworb, r_row)
        colorb_super = (colorb, r_col)
        irow = get(supercell.siteindices, roworb_super, -1)
        icol = get(supercell.siteindices, colorb_super, -1)
        push!(triangle, ((irow, icol), R_col-R_row))
      end
      push!(nnn_triangles, (triangle, -1))
    end


    return (unitcell=unitcell,
            lattice=lattice,
            space_symmetry_embedding=ssymbed,
            nearest_neighbor_bonds=nnbonds,
            next_nearest_neighbor_bonds=nnnbonds,
            nearest_neighbor_triangles=nn_triangles,
            next_nearest_neighbor_triangles=nnn_triangles,)
end




function make_kagome_obc(sites::AbstractVector{<:Tuple})
    site_set = Set(sites)

    min_R1 = minimum(R[1] for (sitetype, R) in sites)
    max_R1 = maximum(R[1] for (sitetype, R) in sites)
    min_R2 = minimum(R[2] for (sitetype, R) in sites)
    max_R2 = maximum(R[2] for (sitetype, R) in sites)

    n1 = max_R1 - min_R1 + 1
    n2 = max_R2 - min_R2 + 1

    big_kagome = make_kagome_lattice([n1+1 0; 0 n2+1])

    site_folding = Dict{Int, Int}()

    new_kagome_supercell = makeunitcell(
        big_kagome.lattice.supercell.latticevectors,
        SiteType=Tuple{String, Vector{Int}}
    )

    for (siteindex, ((sitetype, siteuc), sitecoord)) in enumerate(big_kagome.lattice.supercell.sites)
        if (sitetype, siteuc + [min_R1, min_R2]) in site_set
            addsite!(new_kagome_supercell, (sitetype, siteuc + [min_R1, min_R2]), sitecoord)
            site_folding[siteindex] = numsite(new_kagome_supercell)
        end
    end


    nn_bonds = Tuple{Int, Int}[]
    for ((i, j), R, _) in big_kagome.nearest_neighbor_bonds
        if !haskey(site_folding, i) || !haskey(site_folding, j) || !iszero(R)
            continue
        end
        ip = site_folding[i]
        jp = site_folding[j]
        if ip < jp
            push!(nn_bonds, (ip, jp))
        else
            push!(nn_bonds, (jp, ip))
        end
    end

    #@show big_kagome.nearest_neighbor_triangles

    nn_triangles = Tuple{Int, Int, Int}[]
    for (tri, sgn) in big_kagome.nearest_neighbor_triangles
        triangle_vertices = Set{Int}()
        complete_triangle = true
        for ((i, j), R) in tri
            if !iszero(R)
                complete_triangle = false
                break
            end
            ip = get(site_folding, i, -1)
            jp = get(site_folding, j, -1)

            if ip <= 0 || jp <= 0
                complete_triangle = false
                break
            end

            push!(triangle_vertices, ip)
            push!(triangle_vertices, jp)
        end

        !complete_triangle && continue

        push!(nn_triangles, tuple(sort(collect(triangle_vertices))...))
    end

    return (unitcell=new_kagome_supercell,
            nearest_neighbor_bonds=nn_bonds,
            nearest_neighbor_triangles=nn_triangles)

end


end # module Kagome
