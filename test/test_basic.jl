using Test

using TightBindingLattice

@testset "Basic" begin
    using TightBindingLattice: parse_expr
    @test 1 == parse_expr(1)
    @test 1 == parse_expr("1")
    @test 1.5 == parse_expr("1.5")
    @test 1.5im == parse_expr("1.5im")
    @test 1.5im == parse_expr("1.5i")
    @test [1,2,3] == parse_expr("[1,2,3]")
    @test [1,2,3im] == parse_expr([1, "2", "3i"])


    using TightBindingLattice: cleanup_number

    @test cleanup_number(42, 1E-8) == 42
    @test cleanup_number(1.5 + 1E-12, 1E-8) == 1.5
    @test cleanup_number([1.0 + 1E-12, 2.0 - 1E-12, 3.0], 1E-8) == [1.0, 2.0, 3.0]
    @test cleanup_number(0.1234567, 1E-8) == 0.1234567
end
