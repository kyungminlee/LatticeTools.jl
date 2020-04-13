using Test
using TightBindingLattice

@testset "symmetryoperations" begin

    latticevectors = [1.0 -0.5; 0.0 sqrt(3.0)*0.5]

    @testset "translation" begin
        t1 = TranslationOperation([1,0])
        t2 = TranslationOperation([2,4])
        t3 = t1*t2
        @test t3.displacement == [3,4]
        @test t3 == TranslationOperation([3, 4])
        @test inverse(t3) == TranslationOperation([-3, -4])
        @test apply_symmetry(t3, [5,6]) == [8,10]
        # TODO: Exceptions
    end

    @testset "point" begin
        p1 = PointOperation([0 -1; 1 -1])
        p2 = PointOperation([0 1; 1 0])
        
        @test p1*p2 == PointOperation([0 -1; 1 -1] * [0 1; 1 0])
        @test inverse(p1) * p1 == PointOperation([1 0; 0 1])
        @test p1 * inverse(p1) == PointOperation([1 0; 0 1])

        @test canonize(p1*inverse(p1)) == IdentityOperation()
        @test apply_symmetry(p2, [5, 6]) == [6, 5]
        # TODO: Exceptions
    end

    @testset "product" begin
        t = TranslationOperation([2, 4])
        p = PointOperation([0 -1; 1 -1])

        tp = t * p
        tpc = canonize(tp)

        pt = p * t
        ptc = canonize(pt)

        @test tpc.factors[1] == ptc.factors[1]
        @test isa(tpc.factors[1], PointOperation) && isa(tpc.factors[2], TranslationOperation)
        @test isa(ptc.factors[1], PointOperation) && isa(ptc.factors[2], TranslationOperation)
        @test apply_symmetry(tp, [5,0]) == apply_symmetry(tpc, [5,0])
        @test apply_symmetry(pt, [5,0]) == apply_symmetry(ptc, [5,0])

        @show tpc
        @show inverse(tpc)
        @test iscanonical(tpc)
        @test iscanonical(ptc)
    end
end