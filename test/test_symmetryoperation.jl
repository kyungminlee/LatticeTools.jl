using Test
using TightBindingLattice

@testset "symmetryoperations" begin

    latticevectors = [1.0 -0.5; 0.0 sqrt(3.0)*0.5]
    iden = IdentityOperation()

    @testset "identity" begin
        @test iden * iden == iden
        @test inv(iden) == iden
        @test apply_operation(iden, [3,2,1]) == [3,2,1]
        @test iden([3,2,1]) == [3,2,1]
        @test combinable(iden, iden)
        @test canonize(iden) == iden
        @test iscanonical(iden)
        @test scalartype(iden) == Bool
        @test dimension(iden) == 0
    end

    @testset "translation" begin
        t0 = TranslationOperation([0,0])
        t1 = TranslationOperation([1, -1])
        t1p = TranslationOperation{Int}([1, -1])
        t2 = TranslationOperation([2,4])

        
        @test hash(t1) == hash(t1p)
        @test t1 == t1p
        @test t1 !== t1p
        @test isequal(t1, t1p)

        @test dimension(t2) == 2
        @test t2^2 == TranslationOperation([4,8])

        @test t1 * iden == t1
        @test iden * t1 == t1

        t3 = t1*t2
        @test t3.displacement == [3, 3]
        @test t3 == TranslationOperation([3, 3])
        @test inv(t3) == TranslationOperation([-3, -3])
        @test apply_operation(t3, [5,6]) == [8, 9]

        @test combinable(iden, t1)
        @test combinable(t1, iden)
        @test combinable(t1, t1)

        @test !iscanonical(t0)
        @test iscanonical(t1)
        @test canonize(t0) == IdentityOperation()
        @test canonize(t1) == t1

        @test scalartype(t1) == Int
        @test dimension(t1) == 2

        t4 = TranslationOperation([1,2,3,4])
        @test_throws DimensionMismatch t1 * t4
        @test_throws DimensionMismatch apply_operation(t1, [1,2,3])
        @test_throws DimensionMismatch t1([1,2,3])
    end

    @testset "point" begin
        p0 = PointOperation([1 0; 0 1])
        p1 = PointOperation([0 -1; 1 -1])
        p1p = PointOperation{Int}([0 -1; 1 -1])
        p2 = PointOperation([0 1; 1 0])

        @test p1 == p1p
        @test p1 !== p1p
        @test isequal(p1, p1p)
        @test hash(p1) == hash(p1p)

        @test dimension(p1) == 2

        @test p1*p2 == PointOperation([0 -1; 1 -1] * [0 1; 1 0])
        @test inv(p1) * p1 == PointOperation([1 0; 0 1])
        @test p1 * inv(p1) == PointOperation([1 0; 0 1])
        
        @test p1^3 == p1 * p1 * p1
        @test p1^(-1) == inv(p1)
        @test p1^(-2) == inv(p1)*inv(p1)
        @test p2^2 == PointOperation([1 0; 0 1])

        n = 3
        @test p1^n == p1 * p1 * p1
        n = -1
        @test p1^n == inv(p1)
        n = -2
        @test p1^n == inv(p1)*inv(p1)
        n = 2
        @test p2^n == PointOperation([1 0; 0 1])

        @test combinable(iden, p1)
        @test combinable(p1, iden)
        @test combinable(p1, p1)

        t1 = TranslationOperation([1,0])
        @test !combinable(p1, t1)
        @test !combinable(t1, p1)
                
        @test canonize(p1*inv(p1)) == IdentityOperation()
        @test apply_operation(p2, [5, 6]) == [6, 5]
        @test p2([5, 6]) == [6, 5]
        # TODO: Exceptions

        p3 = PointOperation([0 1 0; 1 0 0; 0 0 1])
        @test_throws DimensionMismatch p1*p3
        @test_throws DimensionMismatch apply_operation(p1, [1,2,3])

        @test !iscanonical(p0)
        @test iscanonical(p1)
        @test canonize(p0) == IdentityOperation()
        @test canonize(p1) == p1

        @test scalartype(p0) == Int
        @test dimension(p0) == 2
    end

    @testset "product" begin
        r0 = ProductOperation()

        t = TranslationOperation([2, 4])
        p = PointOperation([0 -1; 1 -1])

        rt = ProductOperation(t)
        rp = ProductOperation{Int}(p)

        @test r0 * t == rt
        @test t * r0 == rt

        tp = t * p
        pt = p * t
        
        tpc = canonize(tp)
        ptc = canonize(pt)

        @test tp^3 == t * p * t * p * t * p
        @test tpc.factors[1] == ptc.factors[1]
        @test isa(tpc.factors[1], PointOperation) && isa(tpc.factors[2], TranslationOperation)
        @test isa(ptc.factors[1], PointOperation) && isa(ptc.factors[2], TranslationOperation)
        @test apply_operation(tp, [5,0]) == apply_operation(tpc, [5,0])
        @test apply_operation(pt, [5,0]) == apply_operation(ptc, [5,0])
        @test tp([5,0]) == tpc([5,0])
        @test pt([5,0]) == ptc([5,0])

        @test tp^0 == ProductOperation()
        @test tp^1 == tp
        @test tp^2 == t * p * t * p
        @test tp^(-2) == inv(p) * inv(t) * inv(p) * inv(t)

        n = 0
        @test tp^n == ProductOperation()
        n = 1
        @test tp^n == tp
        n = 2
        @test tp^n == t * p * t * p
        n = -2
        @test tp^n == inv(p) * inv(t) * inv(p) * inv(t)


        @test canonize(tp^3 * inv(tp^3)) == IdentityOperation()
        @test iscanonical(pt)
        @test !iscanonical(tp)
        @test iscanonical(tpc)
        @test iscanonical(ptc)

        pp = ProductOperation(p, p)
        @test pp != p * p
        @test canonize(pp) == p*p

        @test scalartype(pp) == Int
        @test dimension(pp) == 2
    end
end