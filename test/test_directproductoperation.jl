using Test
using LatticeTools

struct Phase{T<:Number} <: AbstractSymmetryOperation
    value::T
    function Phase(value::T) where {T<:Number}
        if !isapprox(abs(value), one(value))
            throw(ArgumentError("value $value not of unit length"))
        end
        return new{T}(value)
    end
end

LatticeTools.isidentity(arg::Phase) = isone(arg)
Base.:(*)(lhs::Phase, rhs::Phase) = Phase(lhs.value * rhs.value)
Base.:(^)(lhs::Phase, rhs::Real) = Phase(lhs.value^rhs)
Base.inv(lhs::Phase) = Phase(conj(lhs.value))
Base.isapprox(lhs::Phase, rhs::Phase; kwargs...) = isapprox(lhs.value, rhs.value; kwargs...)

@testset "DirectProductOperation" begin
    α = Phase(-1)
    θ = Phase(cis(2π/3))
    t = TranslationOperation([1, 0])
    let
        θt = directproduct(θ, t)
        @test θt != t ×ˢ θ
        @test θt^3 == directproduct(θ^3, t^3)
        @test inv(θt) == directproduct(inv(θ), inv(t))
    end
    let
        θt = θ ×ˢ t
        @test θt^3 == directproduct(θ^3, t^3)
        @test inv(θt) == directproduct(inv(θ), inv(t))
    end
    let 
        αθt = (α ×ˢ θ) ×ˢ t
        @test αθt == α ×ˢ (θ ×ˢ t)
        @test αθt == directproduct(directproduct(α, θ), t)
        @test αθt == directproduct(α, directproduct(θ, t))
        @test αθt == directproduct(α, directproduct(θ, t))
    end
    @test directproduct([α, α^2], [t, t^2, t^3]) == [
        (α ×ˢ t)   (α ×ˢ t^2)   (α ×ˢ t^3) ;
        (α^2 ×ˢ t) (α^2 ×ˢ t^2) (α^2 ×ˢ t^3)
    ]
    @test [α, α^2] ×ˢ [t, t^2, t^3] == [
        (α ×ˢ t)   (α ×ˢ t^2)   (α ×ˢ t^3) ;
        (α^2 ×ˢ t) (α^2 ×ˢ t^2) (α^2 ×ˢ t^3)
    ]
end