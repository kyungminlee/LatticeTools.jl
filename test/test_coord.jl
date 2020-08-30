using Test
using LatticeTools


@testset "FractCoord" begin
    @testset "Explicit Construction" begin
        fractcoord = FractCoord([-2, 1], [0.8, 0.5])
        @test fractcoord.whole == [-2, 1]
        @test fractcoord.fraction == [0.8, 0.5]
    end

    @testset "Construction by Raw" begin
        fractcoord = FractCoord([-1.2, 1.5])
        @test fractcoord.whole == [-2, 1]
        @test isapprox(fractcoord.fraction, [0.8, 0.5])
    end

    @testset "Zero Construction" begin
        fractcoord = FractCoord(2)
        @test fractcoord.whole == [0,0]
        @test fractcoord.fraction == [0.0, 0.0]
    end

    @testset "methods" begin
        fractcoord = FractCoord([-2, 1], [0.8, 0.5])
        @test dimension(fractcoord) == 2
    end

    @testset "failures" begin
        @test_throws ArgumentError FractCoord(-1)
        @test_throws ArgumentError FractCoord(0)
        @test_throws ArgumentError FractCoord([-1, 1], [0.0])
        @test_throws ArgumentError FractCoord([-1, 1], [1.0, 2.0])
        @test_throws ArgumentError FractCoord([-1, 1], [0.0, 0.0]; tol=-1)
    end

    @testset "isapprox" begin
        fractcoord1 = FractCoord([-2, 1], [0.8, 0.5])
        fractcoord2 = FractCoord([-2, 1], [0.8, 0.5 + 1E-10])
        fractcoord3 = FractCoord([-1, 1], [0.8, 0.5])
        fractcoord4 = FractCoord([-1, 1], [0.8, 0.5+1E-2])

        @test isapprox(fractcoord1, fractcoord2)
        @test !isapprox(fractcoord1, fractcoord3)
        @test !isapprox(fractcoord1, fractcoord4)
    end

    @testset "round" begin
        fractcoord1 = FractCoord([1, 0], [0.0, 0.0])
        fractcoord2 = FractCoord([0, 0], [1.0-1E-12, 0.0]; tol=1E-6)
        @test isapprox(fractcoord1, fractcoord2)
    end

    @testset "Operators" begin
        lvs = [0.5 0.0; 0.0 1.0]
        fc1 = FractCoord([-2, 1], [0.8, 0.5])
        fc2 = FractCoord([1, 2], [0.1, 0.8])
        w = [3, 1]
        @test (+fc1) == fc1
        @test isapprox(+fc1, fc1)
        @test isapprox(-fc1, FractCoord([1, -2], [ 0.2, 0.5]))
        @test isapprox(fc1 + fc2, FractCoord([-1, 4], [0.9, 0.3]))
        @test isapprox(fc1 - fc2, FractCoord([-3, -2], [0.7, 0.7]))
        @test isapprox(fc1 + w, FractCoord([1, 2], [0.8, 0.5]))
        @test isapprox(fc1 - w, FractCoord([-5, 0], [0.8, 0.5]))
    end

    @testset "show" begin
        fc1 = FractCoord([-2, 1], [0.8, 0.5])
        s = string(fc1)
        @test string(fc1) == "FractCoord([-2, 1] + [0.8, 0.5])"
    end

    @testset "fract2carte exceptions" begin
        fc = FractCoord([-2, 1], [0.8, 0.5])

        @test_throws ArgumentError fract2carte([0.5 0.0;], fc)
        @test_throws ArgumentError fract2carte([0.5 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], fc)
    end

    @testset "fract2carte/carte2fract" begin
        latticevectors = [0.5 0.0; 0.0 1.0]
        rawfractcoord = [-1.2, 1.5]

        correctfractcoord = FractCoord([-2, 1], [0.8, 0.5])
        correctcartecoord = [-0.6, 1.5]

        fractcoord = FractCoord(rawfractcoord)
        cartecoord = fract2carte(latticevectors, fractcoord)
        newfractcoord = carte2fract(latticevectors, cartecoord)
        @test isapprox(fractcoord, correctfractcoord)
        @test isapprox(cartecoord, correctcartecoord)
        @test isapprox(newfractcoord, correctfractcoord)

        @testset "tolerance" begin
            fc1 = carte2fract(latticevectors, [0.5 - 1E-9, 1.0])
            fc2 = carte2fract(latticevectors, [0.5 - 1E-9, 1.0]; tol=0.0)
            @test fc1.whole == [1, 1]
            @test isapprox(fc1.fraction, [0.0, 0.0])
            @test fc2.whole == [0, 1]
            @test isapprox(fc2.fraction, [1.0, 0.0])
        end
    end

    @testset "multiplication" begin
        fc = FractCoord([-2, 1], [0.8, 0.5])
        @test isapprox(FractCoord([1, -2], [0.5, 0.8]), [0 1; 1 0] * fc)
    end

end
