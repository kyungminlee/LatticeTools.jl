function display(unitcell ::UnitCell{O}) where O
  lvs = unitcell.latticevectors
  arrow(0,0, lvs[:,1]...)
  arrow(0,0, lvs[:,2]...)

  orbital_ccs = transpose(hcat([fract2carte(unitcell, orbcoord) for (orbname, orbcoord) in unitcell.orbitals]...))
  plot(orbital_ccs[:,1], orbital_ccs[:,2], "o")
end

