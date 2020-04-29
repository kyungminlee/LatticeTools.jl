using Test
using TightBindingLattice

using TightBindingLattice: parse_expr


@testset "parse_expr" begin
    @test 1 == parse_expr(1)
    @test 1 == parse_expr("1")
    @test 1.5 == parse_expr("1.5")
    @test 1.5im == parse_expr("1.5im")
    @test 1.5im == parse_expr("1.5i")
    @test [1,2,3] == parse_expr("[1,2,3]")
    @test [1,2,3im] == parse_expr([1, "2", "3i"])

    @test parse_expr("2") == 2
    @test parse_expr("i") == im
    @test parse_expr("cis(2)") == cis(2)
    @test parse_expr("exp(2i)") == exp(2im)
    @test parse_expr("[10, π, exp(1)]") == [10, π, exp(1)]
    @test parse_expr("[10 π]") == [10 π]
    @test parse_expr("[1 2; 3 4]") == [1 2; 3 4]
    @test parse_expr("(1, 2.0)") == (1, 2.0)

    @test parse_expr("2+3") == 2+3
    @test parse_expr("2-3") == 2-3
    @test parse_expr("2*3") == 2*3
    @test parse_expr("2/3") == 2/3
    @test parse_expr("2//3") == 2//3
    @test parse_expr("2\\3") == 2\3
    @test parse_expr("2^3") == 2^3

    for unary_ftn in [cos, sin, tan, exp, cis, cospi, sinpi, sqrt, log, abs, abs2, sign, conj, conj, real, imag, angle]
        @test parse_expr("$unary_ftn(3.0 + 5.0i)") == unary_ftn(3.0 + 5.0im)
    end
    @test_throws ErrorException parse_expr("j")
    @test_throws ErrorException parse_expr("true ? π : -1")
end

@testset "cleanup" begin
    using TightBindingLattice: cleanup_number
    @test cleanup_number(42, 1E-8) == 42
    @test cleanup_number(1.5 + 1E-12, 1E-8) == 1.5
    @test cleanup_number([1.0 + 1E-12, 2.0 - 1E-12, 3.0], 1E-8) == [1.0, 2.0, 3.0]
    @test cleanup_number(0.1234567, 1E-8) == 0.1234567
end
