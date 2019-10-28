export Permutation
export generate_group

"""
    Permutation(perms ::AbstractVector{Int}; max_order=2048)

Create a permutation of integers from 1 to n.
`perms` should be a permutation of `1:n`.

# Arguments
- `perms`: an integer vector containing a permutation of integers from 1 to n
- `max_order`: maximum order
"""
struct Permutation <: AbstractSymmetryOperation
  map ::Vector{Int}
  order ::Int
  function Permutation(perms ::AbstractVector{Int}; max_order=2048)
    n = length(perms)
    map = Vector{Int}(perms)
    let # check for duplicates
      duplicates = Set{Int}()
      for j in perms
        (1 <= j <= n) || throw(ArgumentError("argument not a proper permutation (target != universe)"))
        j in duplicates && throw(ArgumentError("argument not a proper permutation (contains duplicates)"))
        push!(duplicates, j)
      end
    end

    order = 1
    let # compute cycle length
      current = Int[map[x] for x in 1:n]
      while !issorted(current) && order <= max_order
        current = Int[map[x] for x in current]
        order += 1
      end
      if order > max_order
        throw(OverflowError("cycle length exceeds maximum value (max = $max_order)"))
      end
    end
    return new(map, order)
  end
end


import Base.*
"""
    *(p1 ::Permutation, p2 ::Permutation)

Multiply the two permutation.
Return `[p2.map[x] for x in p1.map]`.

# Examples
```jldoctest
julia> using TightBindingLattice

julia> Permutation([1,3,2]) * Permutation([2,1,3])
Permutation([2, 3, 1], 3)

julia> Permutation([2,1,3]) * Permutation([1,3,2])
Permutation([3, 1, 2], 3)
```
"""
function *(p1 ::Permutation, p2 ::Permutation)
  if length(p1.map) != length(p2.map)
    throw(ArgumentError("permutations of different universes"))
  end
  return Permutation(Int[p2.map[x] for x in p1.map])
end

import Base.*
function *(lhs ::Set{Permutation}, rhs::Permutation)
    return generate_group(lhs..., rhs)
end

function *(lhs ::Permutation, rhs::Set{Permutation})
    return generate_group(lhs, rhs...)
end

function *(lhs ::Set{Permutation}, rhs::Set{Permutation})
    return generate_group(lhs..., rhs...)
end

import Base.^
"""
    ^(perm ::Permutation, pow ::Integer)

Exponentiate the permutation.

# Examples
```jldoctest
julia> using TightBindingLattice

julia> Permutation([2,3,4,1])^2
Permutation([3, 4, 1, 2], 2)
```
"""
function ^(perm ::Permutation, pow ::Integer)
  p = mod(pow, perm.order)
  out = collect(1:length(perm.map))
  for i in 1:p
    out = collect(perm.map[x] for x in out)
  end
  return Permutation(out)
end

export inverse
function inverse(perm ::Permutation)
  out = zeros(Int, length(perm.map))
  for (i, x) in enumerate(perm.map)
    out[x] = i
  end
  return Permutation(out)
end

import Base.==
function ==(p1 ::Permutation, p2::Permutation)
  return p1.map == p2.map
end

import Base.isless
function isless(p1 ::Permutation, p2::Permutation)
  return isless(p1.order, p2.order) || (isequal(p1.order, p2.order) && isless(p1.map, p2.map))
end

import Base.isequal
function isequal(p1 ::Permutation, p2::Permutation)
  return isequal(p1.map, p2.map)
end

import Base.hash
hash(p ::Permutation) = hash(p.map)


function generate_group(generators ::Permutation...)
  change = true
  group = Set{Permutation}([generators...])
  while change
    change = false
    for g1 in generators, g2 in group
      g3 = g1 * g2
      if !(g3 in group)
        change = true
        push!(group, g3)
      end
    end
  end
  return group
end

#
#   shape = [g.order for g in generators]
#   translations = vcat( collect( Iterators.product([0:g.order-1 for g in generators]...) )...)
#   translations = [ [x...] for x in translations]
#   elements = [prod(gen^d for (gen, d) in zip(generators, dist)) for (ig, dist) in enumerate(translations)]
#   return Set(elements)
# end



function groupmod(numer::Permutation, denominators ::AbstractVector{Permutation})
    min_perm = numer
    for denom in denominators
        g = numer * denom
        while g != numer
            if g < min_perm
                min_perm = g
            end
            g = g * denom
        end
    end
    return min_perm
end
