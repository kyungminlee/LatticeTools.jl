using Test
using TightBindingLattice

@testset "Lattice" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
    addorbital!(unitcell, "d", FractCoord([0,0], [0.0, 0.0]))
    addorbital!(unitcell, "px", FractCoord([0,0], [0.5, 0.0]))
    addorbital!(unitcell, "py", FractCoord([0,0], [0.0, 0.5]))
    lattice = make_lattice(unitcell, [2 0; 0 2])

    @test lattice.unitcell == unitcell
    @test lattice.hypercube.scale_matrix == [2 0; 0 2]
    @test numorbital(lattice.supercell) == 3 * 2 * 2
    @test lattice.supercell.latticevectors == [2.0 0.0; 0.0 2.0]
end