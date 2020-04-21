using Test
using TightBindingLattice


@testset "Lattice1D" begin
    unitcell = make_unitcell(2.0; OrbitalType=String)
    addorbital!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, 3)
    scale_matrix = 3 * ones(Int, (1,1))
    @test lattice.unitcell == unitcell
    #@test lattice.hypercube == HypercubicLattice(scale_matrix)
    @test lattice.orthocube == OrthoCube(scale_matrix)
    @test numorbital(lattice.supercell) == 3

    @test_throws DimensionMismatch make_lattice(unitcell, [3 0])
    @test_throws DimensionMismatch make_lattice(unitcell, [3 0; 0 3])
end


@testset "Lattice" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
    addorbital!(unitcell, "d", FractCoord([0,0], [0.0, 0.0]))
    addorbital!(unitcell, "px", FractCoord([0,0], [0.5, 0.0]))
    addorbital!(unitcell, "py", FractCoord([0,0], [0.0, 0.5]))

    @test_throws DimensionMismatch make_lattice(unitcell, 3)
    @test_throws DimensionMismatch make_lattice(unitcell, [1 0 0; 0 1 0; 0 0 1])
    @test_throws ArgumentError make_lattice(unitcell, [4 0; 0 3], [1 1; 0 1]) # generators not orthogonal
    @test_throws ArgumentError make_lattice(unitcell, [4 0; 0 4], [1 0; 0 2]) # generator not unimodular

    lattice = make_lattice(unitcell, [2 0; 0 2])

    @test lattice.unitcell == unitcell
    @test lattice.orthocube.shape_matrix == [2 0; 0 2]
    @test numorbital(lattice.supercell) == 3 * 2 * 2
    @test lattice.supercell.latticevectors == [2.0 0.0; 0.0 2.0]

    @test lattice == make_lattice(unitcell, [2 0; 0 2])
    @test lattice != make_lattice(unitcell, [3 0; 0 2])

    let lattice2 = make_lattice(unitcell, [2 0; 0 2], [1 0; 0 1])
        @test lattice == lattice2
    end
    let lattice3 = make_lattice(unitcell, [2 0; 0 2], [0 1; 1 0])
        @test lattice != lattice3
        @test lattice.unitcell == lattice3.unitcell
        @test lattice.orthocube == lattice3.orthocube
        @test lattice.bravais_coordinates != lattice3.bravais_coordinates
    end
    let unitcell2 = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addorbital!(unitcell2, "d", FractCoord([0,0], [0.0, 0.0]))
        @test lattice != make_lattice(unitcell2, [2 0; 0 2])
    end
end
