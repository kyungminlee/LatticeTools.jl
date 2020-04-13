using LinearAlgebra
using Printf
using Plots

using TightBindingLattice

include("Kagome.jl")

kagome = make_kagome_lattice([4 -2; 2 2])
#kagome = make_kagome_lattice([2 -1; 1 1])

tsym_embed = embed(kagome.lattice, kagome.translation_symmetry)
psym_embed = embed(kagome.lattice, kagome.point_symmetry)
ssym_embed = embed(kagome.lattice, kagome.translation_symmetry, kagome.point_symmetry)

println()
println()
