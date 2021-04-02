using Test
using LatticeTools


@testset "sitemap" begin
    unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
    addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
    addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))
    lattice = make_lattice(unitcell, [4 0; 0 4])

    tsym = FiniteTranslationSymmetry(lattice)
    psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])

    @test findsitemap(unitcell, TranslationOperation([1,0])) == [(1,[0,0]), (2,[0,0])]
    @test findsitemap(unitcell, TranslationOperation([2,0])) == [(1,[0,0]), (2,[0,0])]

    @test findsitemap(unitcell, tsym) == collect(
        Iterators.repeated([(1,[0,0]), (2,[0,0])], 16)
    )

    rotC4 = element(psym, 3)
    @test findsitemap(unitcell, rotC4) == [ (2, [0,0]), (1, [-1,0]) ]
end
