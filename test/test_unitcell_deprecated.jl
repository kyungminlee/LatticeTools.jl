using Test
using TightBindingLattice


@testset "UnitCell" begin
    latticevectors = [0.5 0.0; 0.0 1.0]

    @testset "Constructors" begin
        uc = make_unitcell(latticevectors)
        @test isapprox(uc.latticevectors, latticevectors)
        @test isapprox(uc.reducedreciprocallatticevectors, [2.0 0.0; 0.0 1.0])
        @test isapprox(uc.reciprocallatticevectors, [4*pi 0.0; 0.0 2*pi])

        uc2 = make_unitcell([[0.5,0.0], [0.0,1.0]])

        @test_throws ArgumentError make_unitcell(latticevectors; OrbitalType=Int)
        @test_throws ArgumentError make_unitcell([1.0 0.0;]; OrbitalType=String)
        @test_throws ArgumentError make_unitcell([1.0 0.0; 1.0 0.0]; OrbitalType=String)
        @test_throws ArgumentError make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String, tol=-1E-8)

        uc1d = make_unitcell(2.0)
        uc1d.latticevectors == 2.0 * ones(Float64, 1, 1)
    end

    @testset "Equality" begin
        uc1 = make_unitcell([2.0 0.0; 1.0 2.0]; OrbitalType=String)
        uc2 = make_unitcell([2.0 0.0; 1.0 2.0]; OrbitalType=String)
        uc3 = make_unitcell([2.0 0.0; 1.0 2.0]; OrbitalType=Tuple{String, Int})
        uc4 = make_unitcell([2.0 0.0; 0.0 2.0]; OrbitalType=String)
        @test uc1 == uc2
        @test uc1 != uc3
        @test uc1 != uc4
    end

    @testset "AngledLatticeVectors" begin
        uc1 = make_unitcell([2.0 0.0; 1.0 2.0])
        uc2 = make_unitcell([[2.0, 1.0], [0.0, 2.0]])
        @test uc1 == uc2
        @test isapprox(uc1.reciprocallatticevectors / π, [1.0 -0.5; 0.0 1.0])
        @test isapprox(uc2.reciprocallatticevectors / π, [1.0 -0.5; 0.0 1.0])
    end

    @testset "Methods" begin
        uc = make_unitcell(latticevectors)
        fc1 = FractCoord([0, 0], [0.5, 0.0])
        fc2 = FractCoord([0, 0], [0.0, 0.5])
        index1 = addorbital!(uc, "Ox", fc1)
        index2 = addorbital!(uc, "Oy", fc2)
        @test index1 == 1
        @test index2 == 2

        @test dimension(uc) == 2
        @test numorbital(uc) == 2
        @test hasorbital(uc, "Ox")
        @test hasorbital(uc, "Oy")
        @test !hasorbital(uc, "Oz")
        @test getorbitalindex(uc, "Ox") == 1
        @test getorbitalindex(uc, "Oy") == 2
        @test getorbital(uc, "Ox") == ("Ox", fc1)
        @test getorbital(uc, "Oy") == ("Oy", fc2)
        @test getorbitalcoord(uc, "Ox") == fc1
        @test getorbitalcoord(uc, "Oy") == fc2
        @test getorbitalindexcoord(uc, "Ox") == (1, fc1)
        @test getorbitalindexcoord(uc, "Oy") == (2, fc2)

        @test getorbital(uc, 1) == ("Ox", fc1)
        @test getorbital(uc, 2) == ("Oy", fc2)
        @test getorbitalname(uc, 1) == "Ox"
        @test getorbitalname(uc, 2) == "Oy"
        @test getorbitalcoord(uc, 1) == fc1
        @test getorbitalcoord(uc, 2) == fc2
    end

    @testset "Methods Exceptions" begin
        uc = make_unitcell(latticevectors; OrbitalType=String)
        addorbital!(uc, "Ox", FractCoord([0, 0], [0.0, 0.0]))
        @test_throws ArgumentError addorbital!(uc, "Oy", FractCoord([0], [0.5]))
        @test_throws ArgumentError addorbital!(uc, "Ox", FractCoord([0, 0], [0.0, 0.0]))
    end

    @testset "Type" begin
        uc = make_unitcell(latticevectors; OrbitalType=String)
        fc = Dict(
            "Ox" => FractCoord([0, 0], [0.5, 0.0]),
            "Oy" => FractCoord([0, 0], [0.0, 0.5]),
        )
        for orbital in ["Ox", "Oy"]
            addorbital!(uc, orbital, fc[orbital])
        end
        @test getorbitalname(uc, 1) == "Ox"
        @test getorbitalname(uc, 2) == "Oy"
    end

    @testset "fract2carte/carte2fract" begin
        latticevectors = [0.5 0.0; 0.0 1.0]
        uc = make_unitcell(latticevectors)

        rawfractcoord = [-1.2, 1.5]

        correctfractcoord = FractCoord([-2, 1], [0.8, 0.5])
        correctcartecoord = [-0.6, 1.5]

        fractcoord = FractCoord(rawfractcoord)
        cartecoord = fract2carte(uc, fractcoord)
        newfractcoord = carte2fract(uc, cartecoord)
        @test isapprox(fractcoord, correctfractcoord)
        @test isapprox(cartecoord, correctcartecoord)
        @test isapprox(newfractcoord, correctfractcoord)

        @testset "tolerance" begin
            latticevectors = [1.0 0.0; 0.0 1.0]
            uc = make_unitcell(latticevectors)

            fc1 = carte2fract(uc, [1.0 - 1E-9, 1.0])
            fc2 = carte2fract(uc, [1.0 - 1E-9, 1.0]; tol=0.0)
            @test fc1.whole == [1, 1]
            @test isapprox(fc1.fraction, [0.0, 0.0])
            @test fc2.whole == [0, 1]
            @test isapprox(fc2.fraction, [1.0, 0.0])
        end
    end

    @testset "fract2carte/carte2fract exception" begin
        latticevectors = [0.5 0.0; 0.0 1.0]
        uc1 = make_unitcell(latticevectors)

        carte2fract(uc1, [1.0, 2.0])
        @test_throws ArgumentError carte2fract(uc1, [1.0])
        @test_throws ArgumentError carte2fract(uc1, [1.0, 2.0, 3.0])

        fract2carte(uc1, FractCoord([1,2], [0.1, 0.2]))
        @test_throws ArgumentError fract2carte(uc1, FractCoord([1], [0.1]))
        @test_throws ArgumentError fract2carte(uc1, FractCoord([1,2,3], [0.1, 0.2, 0.3]))
    end

    @testset "momentumgrid" begin
        uc = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)

        @test_throws ArgumentError momentumgrid(uc, [2,])
        @test_throws ArgumentError momentumgrid(uc, [2,3,4])
        @test_throws ArgumentError momentumgrid(uc, [-2,3])

        kg = momentumgrid(uc, [2,3])
        for i in 1:2, j in 1:3
            @test isapprox(kg[i,j], [2π * (i-1) / 2, 2π * (j-1) / 3])
        end
    end

    @testset "whichunitcell" begin
        uc = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addorbital!(uc, "A", FractCoord([0,0], [0.1, 0.1]))
        addorbital!(uc, "B", FractCoord([0,0], [0.2, 0.2]))
        @test whichunitcell(uc, "A", [1.1, 2.1]) == [1, 2]
        @test whichunitcell(uc, "A", FractCoord([1,2], [0.1, 0.1])) == [1, 2]
        @test_throws ArgumentError whichunitcell(uc, "B", [1.1, 2.1])
        @test_throws ArgumentError whichunitcell(uc, "B", FractCoord([1,2], [0.1, 0.1]))
    end

    @testset "findorbitalindex" begin
        uc = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
        addorbital!(uc, "A", FractCoord([0, 0], [0.0, 0.0]))
        @test findorbitalindex(uc, FractCoord([0,0], [0.0, 0.0])) == (1, [0,0])
        @test findorbitalindex(uc, FractCoord([1,0], [0.0, 0.0])) == (1, [1,0])

        @test findorbitalindex(uc, FractCoord([0,0], [0.2, 0.2])) == (-1, Int[])
    end

end
