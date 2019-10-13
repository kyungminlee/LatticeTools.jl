
export Permutation

struct Permutation <: AbstractSymmetryOperation
  map ::Vector{Int}
  cycle_length ::Int
  function Permutation(perms ::AbstractVector{Int}; max_cycle=2048)
    n = length(perms)
    for (i, j) in enumerate(perms)
      if ! (1 <= i <= n)
        throw(ArgumentError("argument is not a proper permutation (domain != universe)"))
      elseif !(1 <= j <= n)
        throw(ArgumentError("argument is not a proper permutation (target != universe)"))
      end
    end
    map = Vector{Int}(perms)
    
    cycle_length = 1
    current = Int[map[x] for x in 1:n]
    while !issorted(current) && cycle_length < max_cycle
      current = Int[map[x] for x in current]
      cycle_length += 1
    end
    if cycle_length == max_cycle
      throw(OverflowError("cycle length exceeds maximum value (max = $max_cycle)"))
    end
    return new(map, cycle_length)
  end
end


import Base.*
function *(p1 ::Permutation, p2 ::Permutation)
  if length(p1.map) != length(p2.map)
    throw(ArgumentError("permutations of different universes"))
  end
  p3 = Int[p2.map[ p1.map[x] ] for x in 1:length(p1.map)]
  return Permutation(p3)
end

import Base.^
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

