export Permutation

"""
    Permutation(perms ::AbstractVector{Int}; max_cycle=2048)

Create a permutation of integers from 1 to n.
`perms` should be a permutation of `1:n`.

# Arguments
- `perms`: an integer vector containing a permutation of integers from 1 to n
- `max_cycle`: maximum length of cycle
"""
struct Permutation <: AbstractSymmetryOperation
  map ::Vector{Int}
  cycle_length ::Int
  function Permutation(perms ::AbstractVector{Int}; max_cycle=2048)
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

    cycle_length = 1
    let # compute cycle length
      current = Int[map[x] for x in 1:n]
      while !issorted(current) && cycle_length <= max_cycle
        current = Int[map[x] for x in current]
        cycle_length += 1
      end
      if cycle_length > max_cycle
        throw(OverflowError("cycle length exceeds maximum value (max = $max_cycle)"))
      end
    end
    return new(map, cycle_length)
  end
end


import Base.*
"""
    *(p1 ::Permutation, p2 ::Permutation)

Multiply the two permutation.
Return `[p2.map[x] for x in p1.map]`.

# Examples
```jldoctest
julia> using ExactDiagonalization

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

import Base.^
"""
    ^(perm ::Permutation, pow ::Integer)

Exponentiate the permutation.

# Examples
```jldoctest
julia> using ExactDiagonalization

julia> Permutation([2,3,4,1])^2
Permutation([3, 4, 1, 2], 2)
```
"""
function ^(perm ::Permutation, pow ::Integer)
  p = mod(pow, perm.cycle_length)
  out = collect(1:length(perm.map))
  for i in 1:p
    out = collect(perm.map[x] for x in out)
  end
  return Permutation(out)
end

import Base.==
function ==(p1 ::Permutation, p2::Permutation)
  return p1.map == p2.map
end

import Base.hash
hash(p ::Permutation) = hash(p.map)
