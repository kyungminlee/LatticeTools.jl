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
end
