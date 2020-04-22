# # Little Group Example

using TightBindingLattice
using Plots


# ## Lattice and symmetry setup

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
addorbital!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]))
addorbital!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]))

lattice = make_lattice(unitcell, [4 0; 0 4])
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])


# ## Little Group

lge = little_group_elements(tsym, 2, psym)
lg = little_group(tsym, 2, psym)
println("Little group: $lg")

lg_matrep = psym.matrix_representations[lge]
println("Matrix representations: $lg_matrep")


# ## Finding isomorphic point groups

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
println("Matrix representations (iso): $lg_matrep2")


# ## Multiplication Tables

@show group_multiplication_table(psym2)
@show group_multiplication_table(lg_matrep)
@show group_multiplication_table(lg_matrep2)

# ## Irreps and Little Groups

println("Irreps and Little Groups")
for tsic in get_irrep_components(tsym)
    idx = tsic.irrep_index
    kf = tsym.fractional_momenta[idx]
    k = lattice.unitcell.reducedreciprocallatticevectors * kf
    psym_little = little_symmetry(tsym, idx, psym)
    println("- irrep_index: $(idx)")
    println("  momentum: $(k)")
    println("  little_point_group: { name: \"$(psym_little.hermann_mauguinn)\", order: $(group_order(psym_little)) }")
    println("  is_psym_compatible: $(iscompatible(tsym, idx, psym))")
    println("  is_psym_little_compatible: $(iscompatible(tsym, idx, psym_little))")
end
