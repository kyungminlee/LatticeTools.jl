# Irreducible Representations

Given a symmetry of a system, its Hamiltonian can be broken into sectors corresponding to the irreducible representations, or *irreps*, of the symmetry group.
(See [Using Symmetry in Exact Diagonalization](http://kyungminlee.org/symmetry-for-numerics).)


## IrrepComponent

A symmetry group in general has multiple irreps, and some of those irreps may be of dimension higher than one.
The data structure that represents the component of an irrep of a symmetry is [`IrrepComponent`](@ref), which has three fields:
`symmetry`, the underlying symmetry, `irrep_index`, the irrep it refers to, and `irrep_component`, the component in a multidimensional irrep.


## Irreps of Translation Symmetry

Translation symmetry forms an Abelian group.
All irreps of a translation symmetry are, therefore, one-dimensional, corresponding to *Fourier modes*, or *momentum sectors*.
Given a Bravais lattice, together with its generating translations, irreps can be computed simply as
```math
\Phi_{\mathbf{k}} = \exp \left( - 2 \pi i \sum_{i} k_{i} x_{i} / L_{i} \right)
```
where $k_{i} \in \{0, 1, \ldots, L_{i}-1 \}$ labels the irrep,
$x_{i}$ labels the translation operations in units of the generators,
and $L_{i}$ is the order (i.e. period length) of the *i*th generator.

Also since all irreps are one-dimensional, `IrrepComponent` of translation symmetry must, therefore, always have `irrep_component=1`.


## Irreps of Point Symmetry

Point symmetry, on the other hand, is not always Abelian.
Unlike with translation symmetry which is abelian, **LatticeTools.jl** does not compute the irreps of the point group.
Instead, it keeps a database of the point symmetries in two and three dimensions, and their irreps in [`PointSymmetryDatabase`](@ref).
(* This may later be replaced by [`IrrepDatabase`](@ref).)
The representation of the point operation, and their irreps are taken from the [Bilbao Crystallographic Server](https://www.cryst.ehu.es).


[Group representation](https://en.wikipedia.org/wiki/Group_representation)
[Irreducible representation](https://en.wikipedia.org/wiki/Irreducible_representation)
