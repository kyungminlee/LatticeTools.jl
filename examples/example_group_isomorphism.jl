# # Group Isomorphism

# ## Set up Basic Functions

using LatticeTools
using Formatting

function display_matrix(io::IO, matrix::AbstractMatrix; prefix::AbstractString="")
    width = ceil(Int, maximum(length("$item") for item in matrix)/4)*4
    for row in eachrow(matrix)
        for (icol, col) in enumerate(row)
            print(io, icol == 1 ? prefix : " ")
            printfmt(io, "{:>$(width)s}", "$col")
        end
        println(io)
    end
end

# ## Generate Two Isomorphic Groups

group1 = FiniteGroup([
    1 2 3 4;
    2 1 4 3;
    3 4 2 1;
    4 3 1 2;
]);
group2 = let mtab1 = group1.multiplication_table,
             mtab2 = zeros(Int, (4, 4)),
             mapping = [1, 3, 2, 4]
             for i in 1:4, j in 1:4
                 mtab2[mapping[i], mapping[j]] = mapping[mtab1[i,j]]
             end
             FiniteGroup(mtab2)
         end;


# ## Group Multiplication Tables

println("Multiplication Table of G₁")
println("--------------------------")
display_matrix(stdout, group_multiplication_table(group1))
println()

println("Multiplication Table of G₂")
println("--------------------------")
display_matrix(stdout, group_multiplication_table(group2))
println()


# ## Group Isomorphism

println("Group Isomorphism  ϕ: G₁ → G₂")
println("-----------------------------")
ϕ = group_isomorphism(group1, group2)
for g in 1:group_order(group1)
    println("  ϕ($g) = $(ϕ[g])")
end
println()

mtab2 = zeros(Int, (group_order(group1), group_order(group1)))
for g in 1:group_order(group1), h in 1:group_order(group1)
    mtab2[ϕ[g], ϕ[h]] = ϕ[ group_product(group1, g, h) ]
end

println("Multiplication Table of ϕ(G₁)")
println("-----------------------------")
display_matrix(stdout, mtab2)
println()
println("  ϕ(g)⋅ϕ(h) = ϕ(g⋅h)")
println("  ϕ(G₁) ≡ G₂ should hold.")
