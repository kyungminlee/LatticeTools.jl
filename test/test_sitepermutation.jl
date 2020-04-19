using Test
using TightBindingLattice


@testset "sitepermutation" begin
    op1 = SitePermutation([3,1,2]) # 1=>3, 2=>1, 3=>2
    op1p = SitePermutation([3,1,2])
    op2 = SitePermutation([2,1,3]) # 1=>2, 2=>1, 3=>3
    @test op2 * op1 != SitePermutation([1,2,3])
    @test op2 * op1 == SitePermutation([3,2,1])

    @test op1 !== op1p
    @test op1 == op1p
    @test hash(op1) == hash(op1p)
    @test isequal(op1, op1p)

    @test op1 != op2
    @test hash(op1) != hash(op2)
    @test !isequal(op1, op2)

    @test inv(op1) == SitePermutation([2,3,1])

    @test op1^-2 == SitePermutation([3,1,2])
    @test op1^-1 == SitePermutation([2,3,1])
    @test op1^0 == SitePermutation([1,2,3])
    @test op1^1 == SitePermutation([3,1,2])
    @test op1^2 == SitePermutation([2,3,1])
    @test op1^3 == SitePermutation([1,2,3])

    let
        n = -2
        @test op1^n == SitePermutation([3,1,2])
        n = -1
        @test op1^n == SitePermutation([2,3,1])
        n = 0
        @test op1^n == SitePermutation([1,2,3])
        n = 1
        @test op1^n == SitePermutation([3,1,2])
        n = 2
        @test op1^n == SitePermutation([2,3,1])
        n = 3
        @test op1^n == SitePermutation([1,2,3])
    end

    @testset "kagome" begin
        include("Kagome.jl")
        kagome = make_kagome_lattice([3 0; 0 3])

        # TranslationOperation
        t1 = TranslationOperation([1, 0])
        t1e = embed(kagome.lattice, t1)
        t2 = TranslationOperation([0, 1])
        t2e = embed(kagome.lattice, t2)

        # PointOperation
        psym = project(PointSymmetryDatabase.get(25), [1 0 0; 0 1 0])
        p6 = psym.elements[6]
        p6e = embed(kagome.lattice, p6)
        @test_throws ArgumentError embed(kagome.lattice, PointOperation([1 2; -1 1]))

        # SpaceOperation
        embed(kagome.lattice, t1 * p6)
        embed(kagome.lattice, p6 * t1)
    end

end