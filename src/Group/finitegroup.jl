export FiniteGroup

export group_order
export period_length
export group_product
export group_inverse

export generate_subgroup
export issubgroup
export isabelian
export minimal_generating_set
export group_multiplication_table

export element, elements
export element_name, element_names

export group_isomorphism
export ishomomorphic


"""
    FiniteGroup

Finite group, with elements {1, 2, 3,..., n}. The identity element is always 1.
Can be constructed using `FiniteGroup(multiplication_table)`

# Fields
* `multiplication_table::Matrix{Int}`: multiplication table
* `period_lengths::Vector{Int}`: period length (order) of every element
* `inverses::Vector{Int}`: inverse of every element
* `conjugacy_classes::Vector{Vector{Int}}`: conjugacy classes

# Examples
```jldoctest
julia> using LatticeTools

julia> FiniteGroup([1 2; 2 1])
FiniteGroup([1 2; 2 1], [1, 2], [1, 2], [[1], [2]])
```
"""
struct FiniteGroup <: AbstractGroup
    multiplication_table::Matrix{Int}

    period_lengths::Vector{Int}
    inverses::Vector{Int}
    conjugacy_classes::Vector{Vector{Int}}

    @doc """
        FiniteGroup(multiplication_table)
    """
    function FiniteGroup(mtab::AbstractMatrix{<:Integer})
        n_elem = size(mtab, 1)
        if size(mtab, 2) != n_elem
            throw(ArgumentError("Multiplication table should be a square matrix"))
        end
        # check identity and closure
        if (mtab[:,1] != 1:n_elem) || (mtab[1,:] != 1:n_elem)
            throw(ArgumentError("Element 1 should be identity"))
        end
        elements = BitSet(1:n_elem)
        # check closure and inverse (sudoku)
        for i in 2:n_elem
            if (BitSet(mtab[:,i]) != elements) || (BitSet(mtab[i,:]) != elements)
                throw(ArgumentError("Multiplication not a group"))
            end
        end
        # check associativity
        for i in 1:n_elem, j in 1:n_elem, k in 1:n_elem
            if mtab[mtab[i, j], k] != mtab[i, mtab[j, k]]
                throw(ArgumentError("Multiplication not associative"))
            end
        end
        # compute cycles
        period_lengths = zeros(Int, n_elem)
        for idx in 1:n_elem
            jdx = idx
            for i in 1:n_elem
                if jdx == 1
                    period_lengths[idx] = i
                    break
                end
                jdx = mtab[jdx, idx]
            end
        end
        # compute inverses
        inverses = zeros(Int, n_elem)
        for i in 1:n_elem
            for j in i:n_elem
                if mtab[i,j] == 1
                    inverses[i] = j
                    inverses[j] = i
                    break
                end
            end
        end
        # compute conjugacy classes
        conjugacy_classes = Vector{Int}[]
        let
            adjacency = [Int[] for i in 1:n_elem]
            for i in 1:n_elem, j in 1:n_elem
                k = mtab[ mtab[j,i], inverses[j] ]
                push!(adjacency[i], k)
                push!(adjacency[k], i)
            end
            for i in eachindex(adjacency)
                sort!(adjacency[i])
                unique!(adjacency[i])
            end
            visited = falses(n_elem)
            for i in 1:n_elem
                visited[i] && continue
                visited[adjacency[i]] .= true
                push!(conjugacy_classes, adjacency[i])
            end
        end

        @assert all(
            period_lengths[first(cc)] == period_lengths[i]
            for cc in conjugacy_classes for i in cc
        )

        new(mtab, period_lengths, inverses, conjugacy_classes)
    end
end


function Base.:(==)(lhs::FiniteGroup, rhs::FiniteGroup)
    return (lhs.multiplication_table == rhs.multiplication_table)
end

"""
    element(group, idx)

Return the element of index `idx`. For `FiniteGroup`, this is somewhat meaningless
since the `idx`th element is `idx`. The sole purpose of this function is the bounds checking.
"""
element(group::FiniteGroup, idx) = Base.OneTo(group_order(group))[idx]

"""
    elements(group)

Return the elements of the group.
"""
elements(group::FiniteGroup) = Base.OneTo(group_order(group))


"""
    element_name(group, idx)

Return the name of element at index `idx`, which is just the string of `idx`.
"""
element_name(group::FiniteGroup, idx) = string.(element(group, idx))


"""
    element_names(group)

Return the names of element.
"""
element_names(group::FiniteGroup) = string.(elements(group))


"""
    group_order(group)

Order of group (i.e. number of elements)
"""
group_order(group::FiniteGroup) = size(group.multiplication_table, 1)


"""
    group_order(group, g)

Order of group element (i.e. period length)
"""
group_order(group::FiniteGroup, g) = group.period_lengths[g]


"""
    period_length(group, g)

Order of group element (i.e. period length)
"""
period_length(group::FiniteGroup, g) = group.period_lengths[g]


"""
    group_multiplication_table(group)

Return multiplcation table of the group.
"""
group_multiplication_table(group::FiniteGroup) = group.multiplication_table


"""
    isabelian(group)

Check if the group is abelian.
"""
function isabelian(group::FiniteGroup)
    return group.multiplication_table == transpose(group.multiplication_table)
end


"""
    group_product(group)

Return a function which computes the group product.
"""
function group_product(group::FiniteGroup) # a bit like currying
    function product(lhs::Integer, rhs::Integer)
        return group.multiplication_table[lhs, rhs]
    end
    function product(lhs::Integer, rhs::AbstractSet{<:Integer})
        return BitSet([group.multiplication_table[lhs, x] for x in rhs])
    end
    function product(lhs::AbstractSet{<:Integer}, rhs::Integer)
        return BitSet([group.multiplication_table[x, rhs] for x in lhs])
    end
    function product(lhs::AbstractSet{<:Integer}, rhs::AbstractSet{<:Integer})
        return BitSet([group.multiplication_table[x, y] for x in lhs for y in rhs])
    end
    return product
end


"""
    group_product(group, lhs, rhs)

Return the result of group multiplication of `lhs` and `rhs`.
If both `lhs` and `rhs` are integers, return an integer.
If either of them is a set (`AbstractSet`) of integers, then return a `BitSet`.
"""
function group_product(group::FiniteGroup, lhs::Integer, rhs::Integer)
    return group.multiplication_table[lhs, rhs]
end


function group_product(group::FiniteGroup, lhs::AbstractSet{<:Integer}, rhs::Integer)
    return BitSet([group_product(group, x, rhs) for x in lhs])
end


function group_product(group::FiniteGroup, lhs::Integer, rhs::AbstractSet{<:Integer})
    return BitSet([group_product(group, lhs, x) for x in rhs])
end


function group_product(
    group::FiniteGroup,
    lhs::AbstractSet{<:Integer},
    rhs::AbstractSet{<:Integer},
)
    return BitSet([group_product(group, x, y) for x in lhs for y in rhs])
end


"""
    group_inverse(group)

Get a function which gives inverse.
"""
function group_inverse(group::FiniteGroup)
    inverse(idx::Integer) = group.inverses[idx]
    inverse(idx::AbstractVector{<:Integer}) = group.inverses[idx]
    return inverse
end


"""
    group_inverse(group, g)

Get inverse of element/elements `g`.
"""
group_inverse(group::FiniteGroup, g::Integer) = group.inverses[g]
group_inverse(group::FiniteGroup, g::AbstractVector{<:Integer}) = group.inverses[g]


"""
    conjugacy_class(group::FiniteGroup, i::Integer)

Conjugacy class of the element `i`.
"""
function conjugacy_class(group::FiniteGroup, i::Integer)
    return findfirst(
        c -> let j = searchsortedfirst(c, i)
            j <= length(c) && c[j] == i
        end,
        group.conjugacy_classes
    )
end


"""
    generate_subgroup(group::FiniteGroup, idx::Integer)

subgroup generated by `generators`. ⟨ {`idx`} ⟩
"""
function generate_subgroup(group::FiniteGroup, idx::Integer)
    out = BitSet()
    sizehint!(out, group_order(group, idx))
    jdx = 1
    for i in 1:group_order(group, idx)
        push!(out, jdx)
        jdx = group_product(group, jdx, idx)
    end
    @assert jdx == 1
    return out
end


"""
    generate_subgroup(group::FiniteGroup, generators)

subgroup generated by `generators`. ⟨ S ⟩
"""
function generate_subgroup(
    group::FiniteGroup,
    generators::G
) where {G<:Union{<:AbstractSet{<:Integer}, <:AbstractVector{<:Integer}}}
    change = true
    subgroup = BitSet(generators)
    push!(subgroup, 1)
    while change
        change = false
        for g1 in generators, g2 in subgroup
            g3 = group_product(group, g1, g2)
            if !(g3 in subgroup)
                change = true
                push!(subgroup, g3)
            end
        end
    end
    return subgroup
end


"""
    issubgroup(group, subset)

Check whether the given subset is a subgroup of `group`.
"""
function issubgroup(group::FiniteGroup, subset::AbstractSet{<:Integer})
    return all(group_product(group, x, y) in subset for x in subset for y in subset)
end


"""
    minimal_generating_set(group)

Get minimally generating set of the finite group.
"""
function minimal_generating_set(group::FiniteGroup)
    ord_group::Int = group_order(group)
    element_queue::Vector{Tuple{Int, Int}} = collect(enumerate(group.period_lengths))
    sort!(element_queue, by=item->(-item[2], item[1]))

    queue_begin = 1
    span = BitSet([1])
    generators = Int[]
    while queue_begin <= ord_group && length(span) < ord_group
        next_index = queue_begin
        next_elem = element_queue[queue_begin][1]
        next_span = generate_subgroup(group, group_product(group, span, next_elem))
        for i in (queue_begin+1):ord_group
            (g, _) = element_queue[i]
            new_span = generate_subgroup(group, group_product(group, span, g))
            if length(new_span) > length(next_span)
                next_index = i
                next_elem = g
                next_span = new_span
            end
        end
        queue_begin = next_index + 1
        span = next_span
        push!(generators, next_elem)
    end
    return generators
end


"""
    group_isomorphism(group1, group2)

Find the isomorphism ϕ: G₁ → G₂. Return nothing if not isomorphic.
"""
function group_isomorphism(group1::FiniteGroup, group2::FiniteGroup)
    group_order(group1) != group_order(group2) && return nothing
    sort(group1.period_lengths) != sort(group2.period_lengths) && return nothing
    let cl1 = sort(length.(group1.conjugacy_classes)),
        cl2 = sort(length.(group1.conjugacy_classes))
        cl1 != cl2 && return nothing
    end
    ord_group = group_order(group1)

    pl1_list = group1.period_lengths
    pl2_list = group2.period_lengths
    cci1_list = zeros(Int, ord_group)

    for (cc_index, cc) in enumerate(group1.conjugacy_classes), i in cc
        cci1_list[i] = cc_index
    end

    element_map = zeros(Int, ord_group)
    class_map = zeros(Int, length(group1.conjugacy_classes))

    function suggest(i::Integer)
        # @assert element_map[i] == 0
        cci1 = cci1_list[i]
        cc1 = group1.conjugacy_classes[cci1]
        if class_map[cci1] != 0
            return (j for j in group2.conjugacy_classes[class_map[cci1]]
                    if !(j in element_map) && (pl2_list[j] == pl1_list[i]))
        else
            return (j for (icc2, cc2) in enumerate(group2.conjugacy_classes)
                          if length(cc2) == length(cc1)
                      for j in cc2
                          if !(j in element_map) && (pl2_list[j] == pl1_list[i]))
        end
    end

    function dfs(i::Integer)
        i > ord_group && return true
        element_map[i] != 0 && return dfs(i+1)
        for j in suggest(i)
            cset = class_map[cci1_list[i]] == 0
            element_map[i] = j
            class_map[cci1_list[i]] = conjugacy_class(group2, j)
            if dfs(i+1)
                return true
            else
                element_map[i] = 0
                if cset
                    class_map[cci1_list[i]] = 0
                end
            end
        end
        return false
    end

    if dfs(1)
        return element_map
    else
        return nothing
    end
end


"""
    group_multiplication_table(elements, product=(*))

Generate a multiplication table from elements with product.
"""
function group_multiplication_table(
    elements::AbstractVector{ElementType},
    product::Function=Base.:(*)
) where {ElementType}
    element_lookup = Dict(k=>i for (i, k) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[ product(elements[i], elements[j]) ]
    end
    return mtab
end


"""
    ishomomorphic(group, representation; product=(*), equal=(==))

Check whether `representation` is homomorphic to `group` under `product` and `equal`,
order preserved.
"""
function ishomomorphic(
    group::FiniteGroup,
    representation::AbstractVector;
    product::Function=Base.:(*),
    equal::Function=Base.:(==)
)
    ord_group = group_order(group)
    if length(representation) != ord_group
        return false
    end
    for i in 1:ord_group, j in 1:ord_group
        if !equal( product(representation[i], representation[j]),
                   representation[ group_product(group, i, j)] )
            return false
        end
    end
    return true
end


# function group_isomorphism_naive(group1::FiniteGroup, group2::FiniteGroup)
#     group_order(group1) != group_order(group2) && return nothing
#     sort(group1.period_lengths) != sort(group2.period_lengths) && return nothing
#
#     ord_group = group_order(group1)
#     element_groups1 = Dict{Tuple{Int, Int}, Vector{Int}}() # group by period lengths and conjugacy class size
#     element_groups2 = Dict{Tuple{Int, Int}, Vector{Int}}() # group by period lengths
#
#     for i in 1:group_order(group1)
#         pl = group1.period_lengths[i]
#         cc = length(group1.conjugacy_classes[conjugacy_class(group1, i)])
#         if !haskey(element_groups1, (pl, cc))
#             element_groups1[(pl, cc)] = Int[]
#         end
#         push!(element_groups1[(pl, cc)], i)
#     end
#
#     for i in 1:group_order(group2)
#         pl = group2.period_lengths[i]
#         cc = length(group2.conjugacy_classes[conjugacy_class(group2, i)])
#         if !haskey(element_groups2, (pl, cc))
#             element_groups2[(pl, cc)] = Int[]
#         end
#         push!(element_groups2[(pl, cc)], i)
#     end
#
#     #q = sort([(pl, i) for (pl, els) in element_groups1 for i in els], rev=true)
#     plccs = sort(collect(keys(element_groups1)))
#     mapping = zeros(Int, ord_group)
#     for perm_set in Iterators.product([
#             permutations(1:length(element_groups1[plcc]), length(element_groups1[plcc]))
#             for plcc in plccs
#         ]...)
#         #@assert length(pls) == length(perm_set)
#         for (ipl, (plcc, perm)) in enumerate(zip(plccs, perm_set))
#             elg1, elg2 = element_groups1[plcc], element_groups2[plcc]
#             for j in 1:length(perm)
#                 mapping[ element_groups1[plcc][j] ] = element_groups2[plcc][ perm[j] ]
#             end
#         end
#         #mtab1p = zeros(Int, (ord_group, ord_group))
#         isiso = true
#         for i in 1:ord_group, j in 1:ord_group
#             #mtab1p[mapping[i], mapping[j]] = mapping[group1.multiplication_table[i, j]]
#             if group2.multiplication_table[mapping[i], mapping[j]] != mapping[group1.multiplication_table[i, j]]
#                 isiso = false
#                 break
#             end
#         end
#         isiso && return mapping
#         #!isiso && continue
#         #mtab1p == group2.multiplication_table && return mapping
#     end
#     return nothing
# end
