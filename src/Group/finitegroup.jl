export FiniteGroup

export group_order
export group_product
export generate_subgroup
export issubgroup
export isabelian
export minimal_generating_set
export group_multiplication_table

export group_isomorphism

using Combinatorics

struct FiniteGroup <: AbstractGroup
    multiplication_table::Matrix{Int}
    period_lengths::Vector{Int}

    function FiniteGroup(mtab::AbstractMatrix{<:Integer})
        n_elem = size(mtab, 1)
        size(mtab, 2) != n_elem && throw(ArgumentError("Multiplication table should be a square matrix"))

        elements = BitSet(1:n_elem)
        mtab[:,1] != 1:n_elem && throw(ArgumentError("Element 1 should be identity"))
        mtab[1,:] != 1:n_elem && throw(ArgumentError("Element 1 should be identity"))
        for i in 2:n_elem
            BitSet(mtab[:,i]) != elements && throw(ArgumentError("Multiplication not a group"))
            BitSet(mtab[i,:]) != elements && throw(ArgumentError("Multiplication not a group"))
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
        new(mtab, period_lengths)
    end
end


"""
    group_order(group)

Order of group (i.e. number of elements)
"""
group_order(group::FiniteGroup) = size(group.multiplication_table, 1)


"""
    group_order(group, g)

Order of group element (i.e. period length)
"""
group_order(group::FiniteGroup, g::Integer) = group.period_lengths[g]


isabelian(group::FiniteGroup) = group.multiplication_table == transpose(group.multiplication_table)


group_multiplication_table(group::FiniteGroup) = group.multiplication_table


"""
    group_product(group, lhs, rhs)
"""
group_product(group::FiniteGroup, lhs::Integer, rhs::Integer) = group.multiplication_table[lhs, rhs]


function group_product(group::FiniteGroup,
                       lhs::AbstractSet{<:Integer}, rhs::Integer)
    return BitSet([group_product(group, x, rhs) for x in lhs])
end


function group_product(group::FiniteGroup,
                       lhs::Integer, rhs::AbstractSet{<:Integer})
    return BitSet([group_product(group, lhs, x) for x in rhs])
end


function group_product(group::FiniteGroup,
                       lhs::AbstractSet{<:Integer},
                       rhs::AbstractSet{<:Integer})
    return BitSet([group_product(group, x, y) for x in lhs for y in rhs])
end


"""
    generate_subgroup(group, idx)

subgroup generated by `generators`. ⟨ {g} ⟩
"""
function generate_subgroup(group::FiniteGroup, idx::Integer)
    out = BitSet([1])
    jdx = idx
    for i in 1:group_order(group)
        jdx == 1 && return out
        push!(out, jdx)
        jdx = group_product(group, jdx, idx)
    end
    assert(false)
    return BitSet()
end


"""
    generate_subgroup(group, generators)

subgroup generated by `generators`. ⟨ S ⟩
"""
function generate_subgroup(group::FiniteGroup,
                           generators::G) where {G<:Union{<:AbstractSet{<:Integer}, <:AbstractVector{<:Integer}}}
    change = true
    subgroup = BitSet(generators)
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
    issubgroup(mtab, subset)
"""
function issubgroup(group::FiniteGroup,
                    subset::AbstractSet{<:Integer})
    return all(group_product(group, x, y) in subset for x in subset for y in subset)
end


"""
    minimal_generating_set
"""
function minimal_generating_set(group::FiniteGroup)
    ord_group ::Int = group_order(group)
    element_queue ::Vector{Tuple{Int, Int}} = collect(enumerate(group.period_lengths))
    sort!(element_queue, by=item->(-item[2], item[1]))

    function factorize(generators::Vector{Int}, span::BitSet, queue_begin::Int)::Bool
        ord_span = length(span)
        ord_span == ord_group && return true
        for i in queue_begin:ord_group
            (g, pl) = element_queue[i]
            if ord_group % (ord_span * pl) == 0
                new_span = generate_subgroup(group, group_product(group, span, g))
                if length(new_span) == ord_span * pl
                    push!(generators, g)
                    factorize(generators, new_span, i+1) && return true
                    pop!(generators)
                end
            end
        end
        return false
    end

    generators = Int[]
    sizehint!(generators, ord_group)
    factorize(generators, BitSet([1]), 1)
    return generators
end







#=
function orthogonalize(group::FiniteGroup)
    !isabelian(group) && throw(ArgumentError("orthogonalize only works for abelian groups"))

    generators = minimal_generating_set(group)
    generator_orders = group.period_lengths[generators]

    elements = Int[1]

    for multiindex in Iterator.product([0:(x-1) for x in generator_orders]...)
        multiindex_vec = [multiindex...]
    end
end
=#




function group_isomorphism(group1::FiniteGroup, group2::FiniteGroup)
    group_order(group1) != group_order(group2) && return nothing
    sort(group1.period_lengths) != sort(group2.period_lengths) && return nothing

    ord_group = group_order(group1)
    element_groups1 = Dict{Int, Vector{Int}}() # group by period lengths
    element_groups2 = Dict{Int, Vector{Int}}() # group by period lengths

    for i in 1:group_order(group1)
        pl = group1.period_lengths[i]
        if !haskey(element_groups1, pl)
            element_groups1[pl] = Int[]
        end
        push!(element_groups1[pl], i)
    end

    for i in 1:group_order(group2)
        pl = group2.period_lengths[i]
        if !haskey(element_groups2, pl)
            element_groups2[pl] = Int[]
        end
        push!(element_groups2[pl], i)
    end

    #q = sort([(pl, i) for (pl, els) in element_groups1 for i in els], rev=true)
    pls = sort(collect(keys(element_groups1)))
    for perm_set in Iterators.product([
            permutations(1:length(element_groups1[pl]), length(element_groups1[pl]))
            for pl in pls
        ]...)
        mapping = zeros(Int, ord_group)
        @assert length(pls) == length(perm_set)
        for (ipl, (pl, perm)) in enumerate(zip(pls, perm_set))

            elg1 = element_groups1[pl]
            elg2 = element_groups2[pl]

            for j in 1:length(perm)
                mapping[ element_groups1[pl][j] ] = element_groups2[pl][perm[j]]
            end
        end

        mtab1p = zeros(Int, (ord_group, ord_group))

        for i in 1:ord_group, j in 1:ord_group
            mtab1p[mapping[i], mapping[j]] = mapping[group1.multiplication_table[i, j]]
        end

        if mtab1p == group2.multiplication_table
            return mapping
        end
    end
    return nothing
end


function group_multiplication_table(elements::AbstractVector{ElementType}) where {ElementType}
    element_lookup = Dict(k=>i for (i, k) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[ elements[i] * elements[j] ]
    end
    return mtab
end
