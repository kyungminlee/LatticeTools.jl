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

export fractional_momentum
export generator_indices, generator_elements


"""
    SymmorphicSymmetry{S1, S2, E}

Symmorphic symmetry of symmetry `S1` ⋊ `S2`, whose elements are of type `E`.
For symmorphic space group, the normal symmetry `S1` is translation symmetry, and `S2` is
point symmetry. `S1` and `S2` can either be `AbstractSymmetry` or `SymmetryEmbedding`.

# Fields
* `normal::S1`: normal subgroup
* `rest::S2`
* `elements::Array{E}`
* `element_names::Array{String}`
"""
struct SymmorphicSymmetry{
    S1<:SymmetryOrEmbedding,
    S2<:SymmetryOrEmbedding,
    E<:AbstractSpaceSymmetryOperation
}<:AbstractSymmetry{E}

    normal::S1   # e.g.) Translation Symmetry
    rest::S2     # e.g.) Point Symmetry
    elements::Array{E}
    element_names::Array{String}

    @doc """
        SymmorphicSymmetry(normal, rest)

    Construct a symmorphic symmetry `normal ⋊ rest`.
    """
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


"""
    ⋊(normal::AbstractSymmetry, rest::AbstractSymmetry)

Return symmorphic symmetry `normal ⋊ rest`.
"""
function ⋊(normal::AbstractSymmetry, rest::AbstractSymmetry)
    return SymmorphicSymmetry(normal, rest)
end


"""
    ⋉(rest::AbstractSymmetry, normal::AbstractSymmetry)

Return symmorphic symmetry `normal ⋊ rest`.
"""
function ⋉(rest::AbstractSymmetry, normal::AbstractSymmetry)
    return SymmorphicSymmetry(normal, rest)
end

function Base.:(==)(lhs::SymmorphicSymmetry{S1, S2, E}, rhs::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E}
    return lhs.normal == rhs.normal && lhs.rest == rhs.rest
end


"""
    symmetry_product(sym)

# Arguments
* `sym::SymmorphicSymmetry{TranslationSymmetry,PointSymmetry,S}) where {S<:SpaceOperation{<:Integer}}`
"""
function symmetry_product(
    sym::SymmorphicSymmetry{TranslationSymmetry,PointSymmetry,S},
) where {S<:SpaceOperation{<:Integer}}
    function product(lhs::SpaceOperation{<:Integer}, rhs::SpaceOperation{<:Integer})
        foo = lhs * rhs
        return S(foo.matrix, sym.normal.orthocube.wrap(foo.displacement)[2])
    end
    return product
end


"""
    group(sym::SymmorphicSymmetry)

Return `FiniteGroup` of the elements.
"""
function group(sym::SymmorphicSymmetry)
    return FiniteGroup(group_multiplication_table(vec(sym.elements), symmetry_product(sym)))
end


"""
    group_order(sym::SymmorphicSymmetry)

Group order of `sym`, which is the product of the subgroup and the normal subgroup.
"""
function group_order(sym::SymmorphicSymmetry)
    return group_order(sym.normal) * group_order(sym.rest)
end

"""
    group_multiplication_table(sym::SymmorphicSymmetry)

Return group multiplication table of `sym`.
"""
function group_multiplication_table(sym::SymmorphicSymmetry)
    return group_multiplication_table(vec(sym.elements), symmetry_product(sym))
end

Base.eltype(sym::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E} = E
Base.valtype(sym::SymmorphicSymmetry{S1, S2, E}) where {S1, S2, E} = E

"""
    element(sym::SymmorphicSymmetry, g...)

Return the `i`th element of the symmorphic symmetry.
"""
element(sym::SymmorphicSymmetry, g...) = sym.elements[g...]

"""
    elements(sym::SymmorphicSymmetry)

Return the elements of the symmorphic symmetry.
"""
elements(sym::SymmorphicSymmetry) = sym.elements


"""
    element_names(sym::SymmorphicSymmetry)

Return the element names of the symmorphic symmetry.
"""
element_names(sym::SymmorphicSymmetry) = sym.element_names

"""
    element_name(sym::SymmorphicSymmetry, g...)

Return the name of the `g`th element of the symmorphic symmetry.
"""
element_name(sym::SymmorphicSymmetry, g...) = sym.element_names[g...]


# character_table(sym::SymmorphicSymmetry) = sym.character_table
# irreps(sym::SymmorphicSymmetry) = sym.irreps
# irrep(sym::SymmorphicSymmetry, idx) = sym.irreps[idx]
# num_irreps(sym::SymmorphicSymmetry) = length(sym.irreps)
# irrep_dimension(sym::SymmorphicSymmetry, idx::Integer) = 1 # size(first(sym.irreps[idx]), 1)


"""
    fractional_momentum(sym, g)

Return the `g`th fractional momentum of the normal (translation) symmetry.

# Arguments
* `sym::SymmorphicSymmetry{<:TranslationSymmetry, S2, E}`
"""
function fractional_momentum(
    sym::SymmorphicSymmetry{<:TranslationSymmetry, S2, E},
    args...
) where {S2, E}
    return fractional_momentum(sym.normal, args...)
end


"""
    generator_indices(sym::SymmorphicSymmetry)

Return the indices of the generators, which is a union of the generators of normal and rest.
"""
function generator_indices(sym::SymmorphicSymmetry)
    gn, gr = sym.normal.generators, sym.rest.generators
    return vcat(
        [CartesianIndex(x..., 1) for x in gn],
        [CartesianIndex(1, x...) for x in gr],
    )
end


"""
    generator_elements(sym::SymmorphicSymmetry)

Return the generating elements
"""
function generator_elements(sym::SymmorphicSymmetry)
    g = generator_indices(sym)
    return element(sym, g)
end
