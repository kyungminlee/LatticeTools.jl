using LatticeTools

function make_kagome_lattice(size_matrix ::AbstractMatrix{<:Integer}; compute_symmetry::Bool=false)
    latticevectors = [1 -0.5; 0 0.5*sqrt(3.0)];
    unitcell = make_unitcell(latticevectors, SiteType=String)
    addsite!(unitcell, "A", carte2fract(unitcell, [0.5, 0.0]))
    addsite!(unitcell, "B", carte2fract(unitcell, [0.25, 0.25*sqrt(3.0)]))
    addsite!(unitcell, "C", carte2fract(unitcell, [0.5+0.25, 0.25*sqrt(3.0)]))

    nnbondtypes = [
        ([ 0, 0], "A", [ 0, 0], "B", 1),
        ([ 0, 0], "A", [ 0, 0], "C", 1),
        ([ 0, 0], "B", [ 0, 0], "C", 1),
        ([ 1, 1], "A", [ 1, 0], "B",-1),
        ([ 1, 0], "B", [ 0, 0], "C",-1),
        ([ 0, 0], "C", [ 1, 1], "A",-1),
    ]

    nnnbondtypes = [
        ([ 0, 0], "A", [ 1, 0], "B", 1), # ◁
        ([ 1, 0], "B", [ 0,-1], "C", 1),
        ([ 0,-1], "C", [ 0, 0], "A", 1),
        ([ 0, 0], "C", [ 1, 0], "A",-1),
        ([ 1, 0], "A", [ 0,-1], "B",-1), # ▷
        ([ 0,-1], "B", [ 0, 0], "C",-1),
    ]

    lattice = make_lattice(unitcell, size_matrix)
    orthocube = lattice.orthocube
    supercell = lattice.supercell
    tsym = TranslationSymmetry(lattice)
    psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])

    nnbonds = []
    nnnbonds = []

    for r in lattice.bravais_coordinates
        for (rowvec, roworb, colvec, colorb, bondsign) in nnbondtypes
            R_row, r_row = orthocube.wrap(r .+ rowvec)
            R_col, r_col = orthocube.wrap(r .+ colvec)
            roworb_super = (roworb, r_row)
            colorb_super = (colorb, r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(nnbonds, ((irow, icol), R_col-R_row, bondsign))
        end
        for (rowvec, roworb, colvec, colorb, bondsign) in nnnbondtypes
            R_row, r_row = orthocube.wrap(r .+ rowvec)
            R_col, r_col = orthocube.wrap(r .+ colvec)
            roworb_super = (roworb, r_row)
            colorb_super = (colorb, r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(nnnbonds, ((irow, icol), R_col-R_row, bondsign))
        end
    end

    nn_triangles = []
    for r in lattice.bravais_coordinates
      triangle = []
      for (rowvec, roworb, colvec, colorb, bondsign) in nnbondtypes[1:3]
        R_row, r_row = orthocube.wrap(r .+ rowvec)
        R_col, r_col = orthocube.wrap(r .+ colvec)
        roworb_super = (roworb, r_row)
        colorb_super = (colorb, r_col)
        irow = get(supercell.siteindices, roworb_super, -1)
        icol = get(supercell.siteindices, colorb_super, -1)
        push!(triangle, ((irow, icol), R_col-R_row))
      end
      push!(nn_triangles, (triangle, 1))

      triangle = []
      for (rowvec, roworb, colvec, colorb, bondsign) in nnbondtypes[4:6]
        R_row, r_row = orthocube.wrap(r .+ rowvec)
        R_col, r_col = orthocube.wrap(r .+ colvec)
        roworb_super = (roworb, r_row)
        colorb_super = (colorb, r_col)
        irow = get(supercell.siteindices, roworb_super, -1)
        icol = get(supercell.siteindices, colorb_super, -1)
        push!(triangle, ((irow, icol), R_col-R_row))
      end
      push!(nn_triangles, (triangle, -1))
    end

    return (unitcell=unitcell,
            lattice=lattice,
            translation_symmetry=tsym,
            point_symmetry=psym,
            nearest_neighbor_bonds=nnbonds,
            next_nearest_neighbor_bonds=nnnbonds,
            nearest_neighbor_triangles=nn_triangles)
end
