
export HypercubicLattice
export dimension
export isequiv
export translation_group_multiplication_table
export orthogonalize

import LinearAlgebra

struct HypercubicLattice
    scale_matrix ::Matrix{Int}

    inverse_scale_matrix ::Matrix{Rational{Int}}
    coordinates ::Vector{Vector{Int}}
    coordinate_indices ::Dict{Vector{Int}, Int}
    wrap ::Function

    function HypercubicLattice(scale_matrix::AbstractMatrix{<:Integer})
        n, m = size(scale_matrix)
        n != m && throw(DimensionMismatch("scale_matrix is not square: dimensions are ($n, $m)"))
        d = abs( ExactLinearAlgebra.determinant(scale_matrix) )
        d == 0 && throw(ArgumentError("scale matrix null"))

        inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)

        function wrap(r::AbstractArray{<:Integer}, mode::RoundingMode=RoundDown)
            rnd = (x) -> round(Int, x, mode)
            R = rnd.(inverse_scale_matrix * r)
            r2 = r - scale_matrix * R
            return R, r2
        end

        max_range = sum((abs.(scale_matrix[:, i]) for i in 1:n))
        coords = Vector{Int}[]
        sizehint!(coords, abs(d))
        for g in Iterators.product([0:2*x+1 for x in max_range]...)
            r1 = collect(g)
            _, r2 = wrap(r1)
            if !(r2 in coords)
                push!(coords, r2)
            end
        end

        @assert length(coords) == d
        @assert all(wrap(r) == (zeros(n), r) for (i, r) in enumerate(coords))
        coord_indices = Dict{Vector{Int}, Int}(r => i for (i, r) in enumerate(coords))

        return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, wrap)
    end


    function HypercubicLattice(scale_matrix ::AbstractMatrix{<:Integer},
                               coords::AbstractVector{<:AbstractVector{<:Integer}})
        n, m = size(scale_matrix)
        n != m && throw(DimensionMismatch("scale_matrix is not square: dimensions are ($n, $m)"))
        d = abs( ExactLinearAlgebra.determinant(scale_matrix) )
        d == 0 && throw(ArgumentError("scale matrix null"))

        !allunique(coords) && throw(ArgumentError("coordinates not unique"))
        length(coords) != d && throw(ArgumentError("too few coordinates"))

        inverse_scale_matrix = ExactLinearAlgebra.inverse(scale_matrix)
        coord_indices = Dict{Vector{Int}, Int}(r => i for (i, r) in enumerate(coords))

        function wrap(r::AbstractArray{<:Integer}, mode::RoundingMode=RoundDown)
            rnd = (x) -> round(Int, x, mode)
            R = rnd.(inverse_scale_matrix * r)
            r2 = r - scale_matrix * R
            return R, r2
        end

        if !all(wrap(r) == (zeros(n), r) for (i, r) in enumerate(coords))
            throw(ArgumentError("coordinates are not that of the hypercubic lattice with size $(scale_matrix)"))
        end

        return new(scale_matrix, inverse_scale_matrix, coords, coord_indices, wrap)
    end
end


dimension(hypercube::HypercubicLattice) = size(hypercube.scale_matrix, 1)


import Base.(==)

function (==)(lhs::HypercubicLattice, rhs::HypercubicLattice)
    (lhs.scale_matrix == rhs.scale_matrix) && (lhs.coordinates == rhs.coordinates)
end

function isequiv(lhs::HypercubicLattice, rhs::HypercubicLattice)
    det_lhs = TightBindingLattice.ExactLinearAlgebra.determinant(lhs.scale_matrix)
    det_rhs = TightBindingLattice.ExactLinearAlgebra.determinant(rhs.scale_matrix)
    det_lhs != det_rhs && return false

    inv_lhs = TightBindingLattice.ExactLinearAlgebra.inverse(lhs.scale_matrix)
    inv_rhs = TightBindingLattice.ExactLinearAlgebra.inverse(rhs.scale_matrix)

    return all(isinteger.(inv_lhs * rhs.scale_matrix)) && all(isinteger.(inv_rhs * lhs.scale_matrix))
end


function translation_group_multiplication_table(hypercube::HypercubicLattice)
    n = length(hypercube.coordinates)
    mtab = zeros(Int, (n, n))
    wrap = hypercube.wrap
    for (i, ri) in enumerate(hypercube.coordinates), (j, rj) in enumerate(hypercube.coordinates)
        _, rk = wrap(ri.+rj)
        k = hypercube.coordinate_indices[rk]
        mtab[i,j] = k
    end
    return mtab
end


# """
#         orthogonalize

# Find orthogonal generators of the hypercube, and reorder the coordinates accordingly
# """
# function orthogonalize(hypercube::HypercubicLattice)
#     first_group = FiniteGroup(translation_group_multiplication_table(hypercube))
#     ord_group = group_order(first_group)

#     @assert isabelian(first_group)
#     @assert ord_group == length(hypercube.coordinates)

#     first_generators = minimal_generating_set(first_group)
#     orthogonal_shape = [first_group.period_lengths[g] for g in first_generators]
#     orthogonal_coordinates = vec([[x...] for x in Iterators.product([0:(d-1) for d in orthogonal_shape]...)])
#     @assert prod(orthogonal_shape) == ord_group
#     @assert length(orthogonal_coordinates) == ord_group

#     coordinates = Vector{Int}[]
#     let ortho_latvec = hcat(hypercube.coordinates[first_generators]...)
#         for r_ortho in orthogonal_coordinates
#             _, r = hypercube.wrap(ortho_latvec * r_ortho)
#             push!(coordinates, r)
#         end
#     end

#     return HypercubicLattice(hypercube.scale_matrix, coordinates)
# end



function decompose_lattice_2d(hypercube::HypercubicLattice)
    dimension(hypercube) != 2 && throw(ArgumentError("shape matrix should be 2x2"))
    group = FiniteGroup(translation_group_multiplication_table(hypercube))
    ord_group = group_order(group)
    allowed_pairs = Vector{Int}[[1,0]]
    sizehint!(allowed_pairs, ord_group^2*3÷4)
    for y in 1:ord_group, x in 0:ord_group
        if gcd(x, y) == 1
            push!(allowed_pairs, [x, y])
        end
    end
    for r1p in allowed_pairs
        if r1p[1] == 0
            continue
        end
        r1 = [r1p[1], -r1p[2]]
        j1 = hypercube.coordinate_indices[ hypercube.wrap(r1)[2] ]
        n1 = group_order(group, j1)
        for r2 in allowed_pairs
            j2 = hypercube.coordinate_indices[ hypercube.wrap(r2)[2] ]
            # condition 1/2: unimodular
            if ExactLinearAlgebra.determinant(hcat(r1, r2)) != 1
                continue
            end
            n2 = group_order(group, j2)
            nprod = length(generate_subgroup(group, [j1, j2]))
            # condition 2/2: orthogonal generators
            n1 * n2 == nprod && return hcat(r1, r2)
        end
    end
    error("Failed to find decomposition")
end


"""
        orthogonalize

Find orthogonal generators of the hypercube, and reorder the coordinates accordingly
"""
function orthogonalize(hypercube::HypercubicLattice)
    if dimension(hypercube) == 1
        return hypercube
    elseif dimension(hypercube) == 2
        generator_translations = decompose_lattice_2d(hypercube)

        R1 = Vector{Int}[]
        let ρ = generator_translations[:,1]
            r = [0,0]
            for i in 1:length(hypercube.coordinates)
                push!(R1, r)
                r = hypercube.wrap(r + ρ)[2]
                iszero(r) && break
            end
        end
        R2 = Vector{Int}[]
        let ρ = generator_translations[:,2]
            r = [0,0]
            for i in 1:length(hypercube.coordinates)
                push!(R2, r)
                r = hypercube.wrap(r + ρ)[2]
                iszero(r) && break
            end
        end

        coordinates = Vector{Int}[]
        sizehint!(coordinates, length(R1)*length(R2))
        for r2 in R2, r1 in R1
            r = hypercube.wrap(r1 + r2)[2]
            push!(coordinates, r)
        end
        return HypercubicLattice(hypercube.scale_matrix, coordinates)
    else
        error("currently only 1d and 2d are supported")
    end

end