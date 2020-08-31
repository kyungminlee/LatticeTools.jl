# # Little Group in 2D

# ## Preamble
using LatticeTools
using Formatting
using Plots

function display_matrix(io::IO, matrix::AbstractMatrix; prefix::AbstractString="")
    width = ceil(Int, maximum(length("$item") for item in matrix)/4)*4
    for row in eachrow(matrix)
        for (icol, col) in enumerate(row)
            if icol == 1
                print(io, prefix)
                printfmt(io, "{:>$(width)s}", "$col")
            else
                printfmt(io, " {:>$(width)s}", "$col")
            end
        end
        println(io)
    end
end

# ## Set up Lattice and Symmetry
unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String);
addsite!(unitcell, "Ox", FractCoord([0,0], [0.5, 0.0]));
addsite!(unitcell, "Oy", FractCoord([0,0], [0.0, 0.5]));

lattice = makelattice(unitcell, [4 0; 0 4]);
tsym = TranslationSymmetry(lattice);
psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0]);


# ## Little Group
lge = little_group_elements(tsym, 2, psym)
lg = little_group(tsym, 2, psym)
println("Little Group")
println("------------")
display_matrix(stdout, group_multiplication_table(lg))


lg_matrep = psym.matrix_representations[lge]
println("Matrix Representations: $lg_matrep")


# ## Finding Point Groups Isomorphic to the Little Group
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
println("Matrix Representations (Isomorphic): $lg_matrep2")


# ## Multiplication Tables
println("Parent Point Group")
println("------------------")
display_matrix(stdout, group_multiplication_table(psym2))
println("Little Group")
println("------------")
display_matrix(stdout, group_multiplication_table(lg_matrep))
println("Isomorphic Little Group")
println("-----------------------")
display_matrix(stdout, group_multiplication_table(lg_matrep2))

# ## Irreps and Little Groups

println("Irreps and Little Groups")
println("------------------------")
for tsic in get_irrep_components(tsym)
    idx = tsic.irrep_index
    kf = tsym.fractional_momenta[idx]
    k = lattice.unitcell.reducedreciprocallatticevectors * kf
    psym_little = little_symmetry(tsym, idx, psym)
    println("- irrep_index: $(idx)")
    println("  momentum: $(k)")
    println("  little_point_group: { name: \"$(psym_little.hermann_mauguin)\", order: $(group_order(psym_little)) }")
    println("  is_psym_compatible: $(iscompatible(tsym, idx, psym))")
    println("  is_psym_little_compatible: $(iscompatible(tsym, idx, psym_little))")
end
