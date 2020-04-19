using TightBindingLattice
import TightBindingLattice.symmetry_product

#  extension of N {\displaystyle N} N by H {\displaystyle H} H).
#  If any of these statements holds (and hence all of them hold, by their equivalence), we say G is the semidirect product of N and H, written
#      G = N ⋊ H {\displaystyle G=N\rtimes H} {\displaystyle G=N\rtimes H} or G = H ⋉ N , {\displaystyle G=H\ltimes N,} {\displaystyle G=H\ltimes N,}
#  or that G splits over N; one also says that G is a semidirect product of H acti

struct SymmorphicSymmetry{S1<:AbstractSymmetry, S2<:AbstractSymmetry, E<:AbstractSymmetryOperation}
    normal::S1   # e.g.) Translation Symmetry
    rest::S2     # e.g.) Point Symmetry
    elements::Array{E}

    function SymmorphicSymmetry(normal::S1, rest::S2) where {S1<:AbstractSymmetry, S2<:AbstractSymmetry}
        E = promote_type(elementtype(normal), elementtype(rest))
        for g in rest, n in normal
            n2 = g * n * inv(g)
            if n2 ∉ normal
                throw(ArgumentError("first argument must be normal"))
            end
        end
        elems = E[x*y for y in elements(normal), x in elements(rest)]
        @assert allunique(elems)
        return new{S1, S2, E}(normal, rest, elems)
    end
end

#import Base.⋊

export ⋊
function ⋊(normal::AbstractSymmetry, rest::AbstractSymmetry)
    return SymmorphicSymmetry(normal, rest)
end


function symmetry_product(sym::SymmorphicSymmetry{TranslationSymmetry, PointSymmetry, SpaceOperation{Int}})
    function product(lhs::SpaceOperation, rhs::SpaceOperation)
        foo = lhs * rhs
        return SpaceOperation(foo.matrix, sym.normal.orthocube.wrap(foo.displacement)[2])
    end
    return product
end

group(sym::SymmorphicSymmetry) = FiniteGroup(group_multiplication_table(vec(sym.elements), symmetry_product(sym)))
group_order(sym::SymmorphicSymmetry) = group_order(normal) * group_order(rest)
group_multiplication_table(sym::SymmorphicSymmetry) = group_multiplication_table(vec(sym.elements), symmetry_product(sym))

element(sym::SymmorphicSymmetry, g) = sym.elements[g]
elements(sym::SymmorphicSymmetry) = sym.elements

# element_names(sym::SymmorphicSymmetry) = sym.element_names
# element_name(sym::SymmorphicSymmetry, g) = sym.element_names[g]

# character_table(sym::SymmorphicSymmetry) = sym.character_table

# irreps(sym::SymmorphicSymmetry) = sym.irreps
# irrep(sym::SymmorphicSymmetry, idx) = sym.irreps[idx]
# num_irreps(sym::SymmorphicSymmetry) = length(sym.irreps)
# irrep_dimension(sym::SymmorphicSymmetry, idx::Integer) = 1 # size(first(sym.irreps[idx]), 1)

# generator_indices(sym::SymmorphicSymmetry) = sym.generators
# generator_elements(sym::SymmorphicSymmetry) = element(sym, sym.generators)

tsym = TranslationSymmetry([3 0; 0 3])
psym = project(PointSymmetryDatabase.find("4mm"), [1 0 0; 0 1 0])

ssym = tsym ⋊ psym

@show ssym
