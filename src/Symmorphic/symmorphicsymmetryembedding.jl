

embed(lattice::Lattice, ssym::SymmorphicSymmetry) = SymmetryEmbedding(lattice, ssym)



function get_irrep_components(sym::SymmorphicSymmetry)
    (SymmorphicIrrepComponent(normal_sic, rest_sic)
        for normal_sic in get_irrep_components(sym.normal)
        for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest)))
end

function get_irrep_components(sym::SymmetryEmbedding{<:SymmorphicSymmetry})
    (SymmorphicIrrepComponent(normal_sic, rest_sic)
        for normal_sic in get_irrep_components(sym.normal)
        for rest_sic in get_irrep_components(little_symmetry(normal_sic, sym.rest)))
end



