using Test
using TightBindingLattice


@testset "orbitalmap" begin
    unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
    addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
    addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))
    lattice = make_lattice(unitcell, [4 0; 0 4])
    
    tsym = TranslationSymmetry(lattice)
    psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
    
    @test findorbitalmap(unitcell, TranslationOperation([1,0])) == [(1,[0,0]), (2,[0,0])]
    @test findorbitalmap(unitcell, TranslationOperation([2,0])) == [(1,[0,0]), (2,[0,0])]
    
    @test findorbitalmap(unitcell, tsym) == collect(
        Iterators.repeated([(1,[0,0]), (2,[0,0])], 16)
    )

    rotC4 = element(psym, 3)
    @test findorbitalmap(unitcell, rotC4) == [ (2, [0,0]), (1, [-1,0]) ]
end
