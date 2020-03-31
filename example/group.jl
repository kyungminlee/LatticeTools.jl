using TightBindingLattice
using Combinatorics

println("# Group Isomorphism")
println()

group1 = FiniteGroup([1 2 3 4;
                      2 1 4 3;
                      3 4 2 1;
                      4 3 1 2])

group2 = let mtab1 = group1.multiplication_table,
             mtab2 = zeros(Int, (4, 4)),
             mapping = [1, 3, 2, 4]
             for i in 1:4, j in 1:4
                 mtab2[mapping[i], mapping[j]] = mapping[mtab1[i,j]]
             end
             FiniteGroup(mtab2)
         end

println("## Group G₁")
display(group_multiplication_table(group1))
println()
println()

println("## Group G₂")
display(group_multiplication_table(group2))
println()
println()

println("## Group isomorphism  ϕ: G₁ → G₂")
ϕ = group_isomorphism(group1, group2)
for g in 1:group_order(group1)
    println("  $g ↦ $(ϕ[g])")
end
println()

mtab2 = zeros(Int, (group_order(group1), group_order(group1)))
for g in 1:group_order(group1), h in 1:group_order(group1)
    mtab2[ϕ[g], ϕ[h]] = ϕ[ group_product(group1, g, h) ]
end

println("# Multiplication table of ϕ(G₁)")
println("  ϕ(g)⋅ϕ(h) = ϕ(g⋅h)")
println("  ϕ(G₁) = G₂ should hold.")

display(mtab2)
println()
