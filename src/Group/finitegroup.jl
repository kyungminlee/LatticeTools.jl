export FiniteGroup

export group_order
export group_product
export generate_subgroup
export issubgroup
export isabelian
export minimal_generating_set
export group_multiplication_table

export group_isomorphism
export ishomomorphic

using Combinatorics

struct FiniteGroup <: AbstractGroup
    multiplication_table::Matrix{Int}

    period_lengths::Vector{Int}
    inverses::Vector{Int}
    conjugacy_classes::Vector{Vector{Int}}

    function FiniteGroup(mtab::AbstractMatrix{<:Integer})
        n_elem = size(mtab, 1)
        size(mtab, 2) != n_elem && throw(ArgumentError("Multiplication table should be a square matrix"))

        elements = BitSet(1:n_elem)
        # check identity and closure
        mtab[:,1] != 1:n_elem && throw(ArgumentError("Element 1 should be identity"))
        mtab[1,:] != 1:n_elem && throw(ArgumentError("Element 1 should be identity"))
        # check closure and inverse
        for i in 2:n_elem
            BitSet(mtab[:,i]) != elements && throw(ArgumentError("Multiplication not a group"))
            BitSet(mtab[i,:]) != elements && throw(ArgumentError("Multiplication not a group"))
        end
        for i in 1:n_elem, j in 1:n_elem, k in 1:n_elem
            mtab[mtab[i, j], k] != mtab[i, mtab[j, k]] && throw(Argument("Multiplication not associative"))
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
            # adjacency = [BitSet() for i in 1:n_elem]
            # for i in 1:n_elem, j in 1:n_elem
            #     k = mtab[ mtab[j,i], inverses[j] ]
            #     push!(adjacency[i], k)
            #     push!(adjacency[k], i)
            # end
            # visited = falses(n_elem)
            # for i in 1:n_elem
            #     visited[i] && continue
            #     for j in adjacency[i]
            #         visited[j] = true
            #     end
            #     push!(conjugacy_classes, sort(collect(adjacency[i])))
            # end
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

        @assert all(period_lengths[first(cc)] == period_lengths[i] for cc in conjugacy_classes for i in cc)

        new(mtab, period_lengths, inverses, conjugacy_classes)
    end
end


import Base.==

function ==(lhs::FiniteGroup, rhs::FiniteGroup)
    return lhs.multiplication_table == rhs.multiplication_table
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


group_inverse(group::FiniteGroup, idx::Integer) = group.inverses[idx]


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
    conjugacy_class(group, i)

Conjugacy class of the element i.
"""
function conjugacy_class(group::FiniteGroup, i::Integer)
    return findfirst(c -> let j = searchsortedfirst(c, i)
        j <= length(c) && c[j] == i
    end, group.conjugacy_classes)
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
    @assert false "algorithm error with generate_subgroup"
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


function group_isomorphism_naive(group1::FiniteGroup, group2::FiniteGroup)
    group_order(group1) != group_order(group2) && return nothing
    sort(group1.period_lengths) != sort(group2.period_lengths) && return nothing

    ord_group = group_order(group1)
    element_groups1 = Dict{Tuple{Int, Int}, Vector{Int}}() # group by period lengths and conjugacy class size
    element_groups2 = Dict{Tuple{Int, Int}, Vector{Int}}() # group by period lengths

    for i in 1:group_order(group1)
        pl = group1.period_lengths[i]
        cc = length(group1.conjugacy_classes[conjugacy_class(group1, i)])
        if !haskey(element_groups1, (pl, cc))
            element_groups1[(pl, cc)] = Int[]
        end
        push!(element_groups1[(pl, cc)], i)
    end

    for i in 1:group_order(group2)
        pl = group2.period_lengths[i]
        cc = length(group2.conjugacy_classes[conjugacy_class(group2, i)])
        if !haskey(element_groups2, (pl, cc))
            element_groups2[(pl, cc)] = Int[]
        end
        push!(element_groups2[(pl, cc)], i)
    end

    #q = sort([(pl, i) for (pl, els) in element_groups1 for i in els], rev=true)
    plccs = sort(collect(keys(element_groups1)))
    mapping = zeros(Int, ord_group)
    for perm_set in Iterators.product([
            permutations(1:length(element_groups1[plcc]), length(element_groups1[plcc]))
            for plcc in plccs
        ]...)
        #@assert length(pls) == length(perm_set)
        for (ipl, (plcc, perm)) in enumerate(zip(plccs, perm_set))
            elg1, elg2 = element_groups1[plcc], element_groups2[plcc]
            for j in 1:length(perm)
                mapping[ element_groups1[plcc][j] ] = element_groups2[plcc][ perm[j] ]
            end
        end
        #mtab1p = zeros(Int, (ord_group, ord_group))
        isiso = true
        for i in 1:ord_group, j in 1:ord_group
            #mtab1p[mapping[i], mapping[j]] = mapping[group1.multiplication_table[i, j]]
            if group2.multiplication_table[mapping[i], mapping[j]] != mapping[group1.multiplication_table[i, j]]
                isiso = false
                break
            end
        end
        isiso && return mapping
        #!isiso && continue
        #mtab1p == group2.multiplication_table && return mapping
    end
    return nothing
end

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
    # cci2_list = zeros(Int, ord_group)

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
            return [j for j in group2.conjugacy_classes[class_map[cci1]]
                    if !(j in element_map) && (pl2_list[j] == pl1_list[i])]
        else
            return [j for (icc2, cc2) in enumerate(group2.conjugacy_classes)
                          if length(cc2) == length(cc1)
                      for j in cc2
                          if !(j in element_map) && (pl2_list[j] == pl1_list[i])]
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



function group_multiplication_table(elements::AbstractVector{ElementType}) where {ElementType}
    element_lookup = Dict(k=>i for (i, k) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[ elements[i] * elements[j] ]
    end
    return mtab
end


function ishomomorphic(group::FiniteGroup, representation::AbstractVector; tol=Base.rtoldefault(Float64))
    ord_group = group_order(group)
    if length(representation) != ord_group
        return false
    end

    for i in 1:ord_group, j in 1:ord_group
        if !isapprox( representation[i] * representation[j],
                      representation[ group_product(group, i, j)] )
            return false
        end
    end
    return true
end
