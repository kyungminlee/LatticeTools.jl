using Test
using TightBindingLattice

using TightBindingLattice: parse_expr
@testset "parse_expr" begin
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
    @test parse_expr("2\\3") == 2\3
    @test parse_expr("2^3") == 2^3

    for unary_ftn in [cos, sin, tan, exp, cis, cospi, sinpi, sqrt, log, abs, abs2, sign, conj, conj, real, imag, angle]
        @test parse_expr("$unary_ftn(3.0 + 5.0i)") == unary_ftn(3.0 + 5.0im)
    end
    @test_throws ErrorException parse_expr("j")
    @test_throws ErrorException parse_expr("true ? π : -1")
end
