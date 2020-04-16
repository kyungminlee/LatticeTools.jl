using TightBindingLattice

# Little group

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])
#idx_C4 = 3

lattice = make_lattice(unitcell, [4 0; 0 4])
perms = get_orbital_permutations(lattice, psym)

tsym = TranslationSymmetry(lattice)

lge = little_group_element_indices(tsym, 2, psym)
lg = little_group(tsym, 2, psym)
@show lg

lg_matrep = psym.matrix_representations[lge]
@show lg_matrep

little_symmetry_candidates = Tuple{PointSymmetry, Vector{Int}}[]
for i in 1:32
    ps = PointSymmetryDatabase.get(i)
    ϕ = group_isomorphism(lg, ps.group)
    if !isnothing(ϕ)
        push!(little_symmetry_candidates, (ps, ϕ))
    end
end
(psym2, ϕ) = first(little_symmetry_candidates)

lg_matrep2 = lg_matrep[ϕ]

# import TightBindingLattice.group_multiplication_table
# function group_multiplication_table(elements::AbstractVector{ElementType}) where {ElementType}
#     element_lookup = Dict(k=>i for (i, k) in enumerate(elements))
#     ord_group = length(elements)
#     mtab = zeros(Int, (ord_group, ord_group))
#     for i in 1:ord_group, j in 1:ord_group
#         mtab[i,j] = element_lookup[ elements[i] * elements[j] ]
#     end
#     return mtab
# end

@show group_multiplication_table(psym2)
@show group_multiplication_table(lg_matrep)
@show group_multiplication_table(lg_matrep2)


for idx in 1:num_irreps(tsym)
    psym_little = little_symmetry(tsym, idx, psym)
    @show idx
    @show tsym.hypercube.coordinates[idx]
    @show psym_little.hermann_mauguinn, group_order(psym_little)
    @show iscompatible(tsym, idx, psym)
    @show iscompatible(tsym, idx, psym_little)
end
