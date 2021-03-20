using Test
using LatticeTools


@testset "Lattice1D" begin
    unitcell = make_unitcell(2.0; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, 3)
    scale_matrix = 3 * ones(Int, (1,1))
    @test lattice.unitcell == unitcell
    #@test lattice.hypercube == HypercubicLattice(scale_matrix)
    @test lattice.hypercube == Hypercube(scale_matrix)
    @test numsite(lattice.supercell) == 3

    @test_throws DimensionMismatch make_lattice(unitcell, [3 0])
    @test_throws DimensionMismatch make_lattice(unitcell, [3 0; 0 3])
end


@testset "Lattice" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
    addsite!(unitcell, "d", FractCoord([0,0], [0.0, 0.0]))
    addsite!(unitcell, "px", FractCoord([0,0], [0.5, 0.0]))
    addsite!(unitcell, "py", FractCoord([0,0], [0.0, 0.5]))

    @test_throws DimensionMismatch make_lattice(unitcell, 3)
    @test_throws DimensionMismatch make_lattice(unitcell, [1 0 0; 0 1 0; 0 0 1])
    @test_throws DimensionMismatch make_lattice(unitcell, [1 0 0; 0 1 0; 0 0 1], [1 0; 0 1])

    @test_throws DimensionMismatch make_lattice(unitcell, [3 0; 0 3], ones(Int, 1, 1))
    @test_throws DimensionMismatch make_lattice(unitcell, [3 0; 0 3], [1 0 0; 0 1 0; 0 0 1])
    @test_throws ArgumentError make_lattice(unitcell, [4 0; 0 3], [1 1; 0 1]) # generators not orthogonal
    @test_throws ArgumentError make_lattice(unitcell, [4 0; 0 4], [1 0; 0 2]) # generator not unimodular

    lattice = make_lattice(unitcell, [2 0; 0 2])

    @test lattice.unitcell == unitcell
    @test lattice.hypercube.shape_matrix == [2 0; 0 2]
    @test numsite(lattice.supercell) == 3 * 2 * 2
    @test lattice.supercell.latticevectors == [2.0 0.0; 0.0 2.0]

    @test lattice == make_lattice(unitcell, [2 0; 0 2])
    @test lattice != make_lattice(unitcell, [3 0; 0 2])

    let lattice2 = make_lattice(unitcell, [2 0; 0 2], [1 0; 0 1])
        @test lattice == lattice2
    end
    let lattice3 = make_lattice(unitcell, [2 0; 0 2], [0 1; 1 0])
        @test lattice != lattice3
        @test lattice.unitcell == lattice3.unitcell
        @test lattice.hypercube == lattice3.hypercube
        @test lattice.bravais_coordinates != lattice3.bravais_coordinates
    end
    let unitcell2 = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
        addsite!(unitcell2, "d", FractCoord([0,0], [0.0, 0.0]))
        @test lattice != make_lattice(unitcell2, [2 0; 0 2])
    end
end


@testset "LatticeWithOrbital" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String, OrbitalType=String)
    addsite!(unitcell, "Cu", FractCoord([0,0], [0.0, 0.0]))
    addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
    addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

    addorbital!(unitcell, "d", "Cu")
    addorbital!(unitcell, "px", "Ox")
    addorbital!(unitcell, "py", 3)

    @test_throws DimensionMismatch make_lattice(unitcell, 3)
    @test_throws DimensionMismatch make_lattice(unitcell, [1 0 0; 0 1 0; 0 0 1])
    @test_throws DimensionMismatch make_lattice(unitcell, [1 0 0; 0 1 0; 0 0 1], [1 0; 0 1])

    @test_throws DimensionMismatch make_lattice(unitcell, [3 0; 0 3], ones(Int, 1, 1))
    @test_throws DimensionMismatch make_lattice(unitcell, [3 0; 0 3], [1 0 0; 0 1 0; 0 0 1])
    @test_throws ArgumentError make_lattice(unitcell, [4 0; 0 3], [1 1; 0 1]) # generators not orthogonal
    @test_throws ArgumentError make_lattice(unitcell, [4 0; 0 4], [1 0; 0 2]) # generator not unimodular

    lattice = make_lattice(unitcell, [2 0; 0 2])

    @test lattice.unitcell == unitcell
    @test lattice.hypercube.shape_matrix == [2 0; 0 2]
    @test numsite(lattice.supercell) == 3 * 2 * 2
    @test lattice.supercell.latticevectors == [2.0 0.0; 0.0 2.0]

    @test lattice == make_lattice(unitcell, [2 0; 0 2])
    @test lattice != make_lattice(unitcell, [3 0; 0 2])

    let lattice2 = make_lattice(unitcell, [2 0; 0 2], [1 0; 0 1])
        @test lattice == lattice2
    end
    let lattice3 = make_lattice(unitcell, [2 0; 0 2], [0 1; 1 0])
        @test lattice != lattice3
        @test lattice.unitcell == lattice3.unitcell
        @test lattice.hypercube == lattice3.hypercube
        @test lattice.bravais_coordinates != lattice3.bravais_coordinates
    end
    let unitcell2 = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String, OrbitalType=String)
        addsite!(unitcell2, "d", FractCoord([0,0], [0.0, 0.0]))
        @test lattice != make_lattice(unitcell2, [2 0; 0 2])
    end

    @test numsite(lattice.supercell) == 12
    for (siteindex, (sitename, sitecoord)) in enumerate(lattice.supercell.sites)
        @test hassite(lattice.unitcell, sitename[1])
        @test isa(sitename[2], Vector{Int})
        @test getsite(lattice.supercell, siteindex) == (sitename, sitecoord)
        @test getsitename(lattice.supercell, siteindex) == sitename
        @test getsitecoord(lattice.supercell, siteindex) == sitecoord
    end

    @test numorbital(lattice.supercell) == 12
    for (orbitalindex, (orbitalname, siteindex)) in enumerate(lattice.supercell.orbitals)
        @test hasorbital(lattice.unitcell, orbitalname[1])
        @test isa(orbitalname[2], Vector{Int})
        small_orbitalindex = getorbitalindex(lattice.unitcell, orbitalname[1])
        small_siteindex = getorbitalsiteindex(lattice.unitcell, small_orbitalindex)
        small_sitename = getsitename(lattice.unitcell, small_siteindex)
        @test small_sitename == getsitename(lattice.supercell, siteindex)[1]
        @test getsitename(lattice.supercell, siteindex) == (small_sitename, orbitalname[2]) 
    end
end