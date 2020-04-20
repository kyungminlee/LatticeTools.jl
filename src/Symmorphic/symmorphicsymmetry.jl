#  extension of N {\displaystyle N} N by H {\displaystyle H} H).
#  If any of these statements holds (and hence all of them hold, by their equivalence), we say G is the semidirect product of N and H, written
#      G = N ⋊ H {\displaystyle G=N\rtimes H} {\displaystyle G=N\rtimes H} or G = H ⋉ N , {\displaystyle G=H\ltimes N,} {\displaystyle G=H\ltimes N,}
#  or that G splits over N; one also says that G is a semidirect product of H acti

export SymmorphicSymmetry
export SymmorphicIrrepComponent
export ⋊, ⋉
export symmetry_product
export group, group_order, group_multiplication_table,
       element, elements, element_name, element_names

export generator_indices, generator_elements


"""
    SymmorphicSymmetry{S1, S2, E}
"""
struct SymmorphicSymmetry{S1<:SymmetryOrEmbedding,
                          S2<:SymmetryOrEmbedding,
                          E<:AbstractSymmetryOperation}<:AbstractSymmetry{E}
    normal::S1   # e.g.) Translation Symmetry
    rest::S2     # e.g.) Point Symmetry
    elements::Array{E}
    element_names::Array{String}

    function SymmorphicSymmetry(normal::S1, rest::S2) where {S1<:SymmetryOrEmbedding, S2<:SymmetryOrEmbedding}
        if !iscompatible(normal, rest)
            throw(ArgumentError("symmetries $normal and $rest are not compatible"))
        end
        E = promote_type(eltype(normal), eltype(rest))
        if !isconcretetype(E)
            throw(ArgumentError("element type $E is not a concrete type"))
        end
        for g in generator_elements(rest), n in generator_elements(normal)
            n2 = g * n * inv(g)
            if n2 ∉ normal
                throw(ArgumentError("first argument must be normal"))
            end
        end
        elems = E[x*y for y in elements(normal), x in elements(rest)]
        elem_names = String["{ $y | $x }"  for y in element_names(normal), x in element_names(rest)]
        @assert allunique(elems)

        return new{S1, S2, E}(normal, rest, elems, elem_names)
    end
end

function ⋊(normal::AbstractSymmetry, rest::AbstractSymmetry)
    return SymmorphicSymmetry(normal, rest)
end
function ⋉(rest::AbstractSymmetry, normal::AbstractSymmetry)
    return SymmorphicSymmetry(normal, rest)
end

import Base.==
function (==)(lhs::SymmorphicSymmetry{S1, S2, E}, rhs::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E}
    return lhs.normal == rhs.normal && lhs.rest == rhs.rest
end

"""
    symmetry_product(sym::SymmorphicSymmetry{TranslationSymmetry,PointSymmetry,S}) where {S<:SpaceOperation{<:Integer}}
"""
function symmetry_product(sym::SymmorphicSymmetry{TranslationSymmetry,PointSymmetry,S}
            ) where {S<:SpaceOperation{<:Integer}}
    function product(lhs::SpaceOperation{<:Integer}, rhs::SpaceOperation{<:Integer})
        foo = lhs * rhs
        return S(foo.matrix, sym.normal.orthocube.wrap(foo.displacement)[2])
    end
    return product
end

group(sym::SymmorphicSymmetry) = FiniteGroup(group_multiplication_table(vec(sym.elements), symmetry_product(sym)))
group_order(sym::SymmorphicSymmetry) = group_order(sym.normal) * group_order(sym.rest)
group_multiplication_table(sym::SymmorphicSymmetry) = group_multiplication_table(vec(sym.elements), symmetry_product(sym))

import Base.eltype
eltype(sym::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E} = E
import Base.valtype
valtype(sym::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E} = E

element(sym::SymmorphicSymmetry, g...) = sym.elements[g...]
elements(sym::SymmorphicSymmetry) = sym.elements

element_name(sym::SymmorphicSymmetry, g...) = sym.element_names[g...]
element_names(sym::SymmorphicSymmetry) = sym.element_names

# character_table(sym::SymmorphicSymmetry) = sym.character_table

# irreps(sym::SymmorphicSymmetry) = sym.irreps
# irrep(sym::SymmorphicSymmetry, idx) = sym.irreps[idx]
# num_irreps(sym::SymmorphicSymmetry) = length(sym.irreps)
# irrep_dimension(sym::SymmorphicSymmetry, idx::Integer) = 1 # size(first(sym.irreps[idx]), 1)


"""
    generator_indices(sym::SymmorphicSymmetry)
"""
function generator_indices(sym::SymmorphicSymmetry)
    gn, gr = sym.normal.generators, sym.rest.generators
    return vcat([CartesianIndex(x..., 1) for x in gn], [CartesianIndex(1, x...) for x in gr])
end


"""
    generator_elements(sym::SymmorphicSymmetry)
"""
function generator_elements(sym::SymmorphicSymmetry) 
    g = generator_indices(sym)
    return element(sym, g)
end
