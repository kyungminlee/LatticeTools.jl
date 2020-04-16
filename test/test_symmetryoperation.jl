using Test
using TightBindingLattice

@testset "symmetryoperation" begin
    iden = IdentityOperation(Int, 2)

    @testset "IdentityOperation" begin
        @test iden == IdentityOperation{Int}(2)
       
        @test iden != IdentityOperation(Int, 3)
        @test iden * iden == iden
        n=0;  @test iden^n == iden
        n=4;  @test iden^n == iden
        n=-3; @test iden^n == iden
        @test inv(iden) == iden
        @test isidentity(iden)
        @test domaintype(iden) == Int
        @test dimension(iden) == 2
        @test hash(iden) == hash(IdentityOperation(Int, 2))
        @test hash(iden) != hash(IdentityOperation(Int, 3))
        
        @test apply_operation(iden, [3,2]) == [3,2]
        @test iden([3,2]) == [3,2]

        @test_throws DimensionMismatch iden * IdentityOperation(Int, 3)
        @test_throws DimensionMismatch iden([3,2,1])
        @test_throws DimensionMismatch apply_operation(iden, [3,2,1])
    end

    @testset "TranslationOperation" begin
        t0  = TranslationOperation([0,0])
        t1  = TranslationOperation([1, -1])
        t1p = TranslationOperation{Int}([1, -1])
        t2  = TranslationOperation([2,4])

        @testset "properties" begin
            @test domaintype(t1) == Int
            @test dimension(t1) == 2
            @test dimension(t2) == 2
            @test isidentity(t0)
            @test !isidentity(t1)
        end

        @testset "equality" begin
            @test t1 == t1p
            @test t1 !== t1p
            @test isequal(t1, t1p)
            @test hash(t1) == hash(t1p)
            @test hash(t1) != hash(t2)
        end

        @testset "equality-heterotype" begin
            @test t0 == IdentityOperation(Int, 2)
            @test IdentityOperation(Int, 2) == t0
        end
        
        @testset "operator" begin
            @test t1 * iden == t1
            @test iden * t1 == t1

            @test_throws DimensionMismatch t1 * IdentityOperation(Int, 5)
            @test_throws DimensionMismatch IdentityOperation(Int, 5) * t1

            @test t2^2 == TranslationOperation([4,8])

            t3 = t1*t2
            @test t3.displacement == [3, 3]
            @test t3 == TranslationOperation([3, 3])
            @test inv(t3) == TranslationOperation([-3, -3])

            t4 = TranslationOperation([1,2,3,4])
            @test_throws DimensionMismatch t1 * t4
            @test_throws DimensionMismatch apply_operation(t1, [1,2,3])
            @test_throws DimensionMismatch t1([1,2,3])
        end

        @testset "apply" begin
            @test apply_operation(t1, [5,6]) == [6, 5]
            @test t1([5,6]) == [6, 5]
        end

        @testset "promotion" begin
            @test promote_type(typeof(t1), typeof(iden)) == TranslationOperation{Int}
            arr = TranslationOperation{Int}[]
            push!(arr, IdentityOperation(Int, 2))
            @test isa(arr[1], TranslationOperation{Int})
            @test arr[1].displacement == [0,0]
        end
    end

    @testset "PointOperation" begin
        p0 = PointOperation([1 0; 0 1])
        p1 = PointOperation([0 -1; 1 -1])
        p1p = PointOperation{Int}([0 -1; 1 -1])
        p2 = PointOperation([0 1; 1 0])

        @testset "properties" begin
            @test domaintype(p0) == Int
            @test dimension(p0) == 2
            @test isidentity(p0)
            @test !isidentity(p1)
            @test dimension(p1) == 2
        end

        @testset "equality" begin
            @test p1 == p1p
            @test p1 !== p1p
            @test isequal(p1, p1p)
            @test hash(p1) == hash(p1p)
            @test hash(p1) != hash(p2)
        end

        @testset "equality-heterotype" begin
            @test p0 == iden
            @test iden == p0
            @test p1 != iden
            @test iden != p1
            @test p1 != TranslationOperation([1,0])
            @test p0 == TranslationOperation([0,0])
            @test TranslationOperation([1,0]) != p1
            @test TranslationOperation([0,0]) == p0
        end

        @testset "operator" begin
            @test p1*p2 == PointOperation([0 -1; 1 -1] * [0 1; 1 0])
            @test inv(p1) * p1 == PointOperation([1 0; 0 1])
            @test p1 * inv(p1) == PointOperation([1 0; 0 1])
            
            @test p1^3 == p1 * p1 * p1
            @test p1^(-1) == inv(p1)
            @test p1^(-2) == inv(p1)*inv(p1)
            @test p2^2 == PointOperation([1 0; 0 1])

            n =  3;  @test p1^n == p1 * p1 * p1
            n = -1;  @test p1^n == inv(p1)
            n = -2;  @test p1^n == inv(p1)*inv(p1)
            n =  2;  @test p2^n == PointOperation([1 0; 0 1])
        end

        @testset "apply" begin
            @test apply_operation(p2, [5, 6]) == [6, 5]
            @test p2([5, 6]) == [6, 5]

            p3 = PointOperation([0 1 0; 1 0 0; 0 0 1])
            @test_throws DimensionMismatch p1*p3
            @test_throws DimensionMismatch apply_operation(p1, [1,2,3])
        end

        @testset "promotion" begin
            @test promote_type(typeof(p1), typeof(iden)) == PointOperation{Int}
            arr = PointOperation{Int}[]
            push!(arr, IdentityOperation(Int, 2))
            @test isa(arr[1], PointOperation{Int})
            @test arr[1].matrix == [1 0; 0 1]
        end
    end

    @testset "SpaceOperation" begin
        @test_throws DimensionMismatch SpaceOperation([1 0 0; 0 1 0], [1,2,3])
        @test_throws DimensionMismatch SpaceOperation([1 0; 0 1], [1,2,3])

        let sop = SpaceOperation{Int}(2)
            @test sop.matrix == [1 0; 0 1]
            @test sop.displacement == [0, 0]
        end

        c4p = PointOperation([ 0 -1;  1  0])
        c4m = PointOperation([ 0  1; -1  0])
        m10 = PointOperation([-1  0;  0  1])
        t10 = TranslationOperation([ 1,  0])

        @testset "properties" begin
            s0 = SpaceOperation(Int, 2)
            s1 = SpaceOperation(c4p, t10)
            @test domaintype(s0) == Int
            @test dimension(s0) == 2
            @test isidentity(s0)
            @test !isidentity(s1)
            @test dimension(s1) == 2
        end

        @testset "equality" begin
            sop = SpaceOperation(c4p)
            @test sop.matrix == [0 -1; 1 0] && sop.displacement == [0,0]
            
            sop = SpaceOperation(t10)
            @test sop.matrix == [1 0; 0 1] && sop.displacement == [1,0]
            
            sop = SpaceOperation(c4p, t10)
            @test sop.matrix == [0 -1; 1 0] && sop.displacement == [1, 0]
            @test sop == SpaceOperation([0 -1; 1 0], [1, 0])
            @test sop != SpaceOperation([0 -1; 1 0], [0, 1])
            @test sop != SpaceOperation([0  1; 1 0], [0, 1])

            @test hash(sop) == hash(SpaceOperation([0 -1; 1 0], [1,0]))
            @test hash(sop) != hash(SpaceOperation([0  1; 1 0], [1,0]))
            @test hash(sop) != hash(SpaceOperation([0 -1; 1 0], [0,0]))
        end

        @testset "product" begin
            c4p_t10 = c4p*t10
            t10_c4p = t10*c4p

            @test c4p_t10 == SpaceOperation([0 -1; 1 0], [1,  0])
            @test t10_c4p == SpaceOperation([0 -1; 1 0], [0, -1])
            
            @test c4p * c4p_t10 == SpaceOperation([-1 0; 0 -1], [ 1,  0])
            @test c4p * t10_c4p == SpaceOperation([-1 0; 0 -1], [ 0, -1])
            @test c4p_t10 * c4p == SpaceOperation([-1 0; 0 -1], [ 0, -1])
            @test t10_c4p * c4p == SpaceOperation([-1 0; 0 -1], [-1,  0])

            # associativity
            for A in [t10, m10, c4p, c4m],
                B in [t10, m10, c4p, c4m],
                C in [t10, m10, c4p, c4m]
                @test (A * B) * C == A * (B * C)
            end
            # inverse
            for A in [t10, m10, c4p, c4m], B in [t10, m10, c4p, c4m]
                @test inv(A * B) == inv(B) * inv(A)
            end
        end

        @testset "comparison with other types" begin
            @test SpaceOperation(c4p) == c4p
            @test c4p == SpaceOperation(c4p)
            @test isidentity(SpaceOperation(c4p * c4m))
            @test SpaceOperation(c4p * c4m) == IdentityOperation(Int, 2)
            @test IdentityOperation(Int, 2) == SpaceOperation(c4p * c4m)

            @test m10 * m10 * t10 == t10
            @test c4p * t10 * c4m == TranslationOperation([0, 1])
            @test TranslationOperation([0, 1]) == c4p * t10 * c4m
            @test t10 * c4p * c4p * t10 == PointOperation([-1 0; 0 -1])
            @test PointOperation([-1 0; 0 -1]) == t10 * c4p * c4p * t10
        end

        @testset "promotion" begin
            @test promote_type(IdentityOperation{Int}, SpaceOperation{Int}) == SpaceOperation{Int}
            @test promote_type(TranslationOperation{Int}, SpaceOperation{Int}) == SpaceOperation{Int}
            @test promote_type(PointOperation{Int}, SpaceOperation{Int}) == SpaceOperation{Int}

            arr = SpaceOperation{Int}[]
            push!(arr, PointOperation([0 1; 1 0]))
            push!(arr, TranslationOperation([1, 0]))
            push!(arr, IdentityOperation(Int, 2))

            @test isa(arr[1], SpaceOperation{Int})
            @test arr[1].matrix == [0 1; 1 0] && arr[1].displacement == [0,0]
            @test isa(arr[2], SpaceOperation{Int})
            @test arr[2].matrix == [1 0; 0 1] && arr[2].displacement == [1,0]
            @test isa(arr[3], SpaceOperation{Int})
            @test arr[3].matrix == [1 0; 0 1] && arr[3].displacement == [0,0]
        end

        @testset "apply" begin
            t = TranslationOperation([1, 0])
            p = PointOperation([0 1; 1 0])
            tp = t * p
            pt = p * t

            @test tp^3 == t * p * t * p * t * p
            @test apply_operation(tp, [5,0]) == [1, 5]
            @test apply_operation(pt, [5,0]) == [0, 6]
            @test tp([5,0]) == [1, 5]
            @test pt([5,0]) == [0, 6]
            # @test tp([5,0]) == tpc([5,0])
            # @test pt([5,0]) == ptc([5,0])
            @test tp^0 == IdentityOperation(Int, 2)
            @test tp^1 == tp
            @test tp^2 == t * p * t * p
            @test tp^(-2) == inv(p) * inv(t) * inv(p) * inv(t)

            n = 0;   @test tp^n == IdentityOperation(Int, 2)
            n = 1;   @test tp^n == tp
            n = 2;   @test tp^n == t * p * t * p
            n = -2;  @test tp^n == inv(p) * inv(t) * inv(p) * inv(t)
        end
    end

    # @testset "product" begin
    #     r0 = ProductOperation()

    #     t = TranslationOperation([2, 4])
    #     p = PointOperation([0 -1; 1 -1])

    #     rt = ProductOperation(t)
    #     rp = ProductOperation{Int}(p)

    #     @test r0 * t == rt
    #     @test t * r0 == rt

    #     tp = t * p
    #     pt = p * t
        
    #     tpc = canonize(tp)
    #     ptc = canonize(pt)

    #     @test tp^3 == t * p * t * p * t * p
    #     @test tpc.factors[1] == ptc.factors[1]
    #     @test isa(tpc.factors[1], PointOperation) && isa(tpc.factors[2], TranslationOperation)
    #     @test isa(ptc.factors[1], PointOperation) && isa(ptc.factors[2], TranslationOperation)
    #     @test apply_operation(tp, [5,0]) == apply_operation(tpc, [5,0])
    #     @test apply_operation(pt, [5,0]) == apply_operation(ptc, [5,0])
    #     @test tp([5,0]) == tpc([5,0])
    #     @test pt([5,0]) == ptc([5,0])

    #     @test tp^0 == ProductOperation()
    #     @test tp^1 == tp
    #     @test tp^2 == t * p * t * p
    #     @test tp^(-2) == inv(p) * inv(t) * inv(p) * inv(t)

    #     n = 0
    #     @test tp^n == ProductOperation()
    #     n = 1
    #     @test tp^n == tp
    #     n = 2
    #     @test tp^n == t * p * t * p
    #     n = -2
    #     @test tp^n == inv(p) * inv(t) * inv(p) * inv(t)


    #     @test canonize(tp^3 * inv(tp^3)) == IdentityOperation()
    #     @test iscanonical(pt)
    #     @test !iscanonical(tp)
    #     @test iscanonical(tpc)
    #     @test iscanonical(ptc)

    #     pp = ProductOperation(p, p)
    #     @test pp != p * p
    #     @test canonize(pp) == p*p

    #     @test domaintype(pp) == Int
    #     @test dimension(pp) == 2
    # end
end