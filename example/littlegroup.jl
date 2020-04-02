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

lge = little_group_elements(tsym, 2, psym)
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

#=
group::FiniteGroup
generators::Vector{Int}
conjugacy_classes::Vector{ConjugacyClassType}
character_table::Matrix{ComplexF64}
irreps::Vector{IrrepType}
element_names::Vector{String}
matrix_representations::Vector{Matrix{Int}}
hermann_mauguinn::String
=#

# lg_matrep[ϕ]
# for each i: lg_matrep2[i] = lg_matrep[ϕ[i]]
lg_matrep2 = lg_matrep[ϕ]

import TightBindingLattice.group_multiplication_table
function group_multiplication_table(elements::AbstractVector{ElementType}) where {ElementType}
    element_lookup = Dict(k=>i for (i, k) in enumerate(elements))
    ord_group = length(elements)
    mtab = zeros(Int, (ord_group, ord_group))
    for i in 1:ord_group, j in 1:ord_group
        mtab[i,j] = element_lookup[ elements[i] * elements[j] ]
    end
    return mtab
end

@show group_multiplication_table(psym2)
@show group_multiplication_table(lg_matrep)
@show group_multiplication_table(lg_matrep2)
