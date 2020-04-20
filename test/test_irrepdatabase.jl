using Test
using TightBindingLattice


@testset "irrepdatabase" begin
    # D<sub>10</sub> group
    @testset "find" begin
        group1 = FiniteGroup([
            1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20;
            2  3  4  5  6  7  8  9  10  1  20  11  12  13  14  15  16  17  18  19;
            3  4  5  6  7  8  9  10  1  2  19  20  11  12  13  14  15  16  17  18;
            4  5  6  7  8  9  10  1  2  3  18  19  20  11  12  13  14  15  16  17;
            5  6  7  8  9  10  1  2  3  4  17  18  19  20  11  12  13  14  15  16;
            6  7  8  9  10  1  2  3  4  5  16  17  18  19  20  11  12  13  14  15;
            7  8  9  10  1  2  3  4  5  6  15  16  17  18  19  20  11  12  13  14;
            8  9  10  1  2  3  4  5  6  7  14  15  16  17  18  19  20  11  12  13;
            9  10  1  2  3  4  5  6  7  8  13  14  15  16  17  18  19  20  11  12;
            10  1  2  3  4  5  6  7  8  9  12  13  14  15  16  17  18  19  20  11;
            11  12  13  14  15  16  17  18  19  20  1  2  3  4  5  6  7  8  9  10;
            12  13  14  15  16  17  18  19  20  11  10  1  2  3  4  5  6  7  8  9;
            13  14  15  16  17  18  19  20  11  12  9  10  1  2  3  4  5  6  7  8;
            14  15  16  17  18  19  20  11  12  13  8  9  10  1  2  3  4  5  6  7;
            15  16  17  18  19  20  11  12  13  14  7  8  9  10  1  2  3  4  5  6;
            16  17  18  19  20  11  12  13  14  15  6  7  8  9  10  1  2  3  4  5;
            17  18  19  20  11  12  13  14  15  16  5  6  7  8  9  10  1  2  3  4;
            18  19  20  11  12  13  14  15  16  17  4  5  6  7  8  9  10  1  2  3;
            19  20  11  12  13  14  15  16  17  18  3  4  5  6  7  8  9  10  1  2;
            20  11  12  13  14  15  16  17  18  19  2  3  4  5  6  7  8  9  10  1
        ])
        @test isnothing(IrrepDatabase.find(group1))
        group2 = FiniteGroup([1 2; 2 1])
        (irrep_data, elmap) = IrrepDatabase.find(group2)
        @test irrep_data.group == group2
        @test irrep_data.irreps == [[ones(ComplexF64, 1, 1), ones(ComplexF64, 1, 1)],
                                    [ones(ComplexF64, 1, 1), -ones(ComplexF64, 1, 1)]]
        @test irrep_data.conjugacy_classes == [[1], [2]]
        @test irrep_data.character_table == [1 1; 1 -1]
        @test elmap == [1,2]
    end

    @testset "find_irrep" begin
        (irrep, matrep_new, ϕ) = IrrepDatabase.find_irrep([ [1 0; 0 1], [1 0; 0 -1]])
        @test irrep.group == FiniteGroup([1 2; 2 1])
    end
end
