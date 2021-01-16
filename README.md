<img src="https://kyungminlee.org/LatticeTools.jl/stable/assets/logo.png" width="400px">

# LatticeTools.jl

| **Documentation** | **Build Status** | **Code Coverage** |
|:-----------------:|:----------------:|:-----------------:|
| [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] | [![Build][githubaction-img]][githubaction-url] | [![Code Coverage][codecov-img]][codecov-url] |

[LatticeTools.jl](https://github.com/kyungminlee/LatticeTools.jl) is a Julia package that provides functionalities to define lattices and perform symmetry analyses useful for studying interacting quantum many-body Hamiltonians.


## Installation

[LatticeTools.jl](https://github.com/kyungminlee/LatticeTools.jl) is currently not included in Julia's default package registry. To install, add the package registry [KyungminLeeRegistry](https://github.com/kyungminlee/KyungminLeeRegistry.jl) and then install the package:

```julia-repl
(@v1.5) pkg> registry add https://github.com/kyungminlee/KyungminLeeRegistry.git
(@v1.5) pkg> add LatticeTools
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://kyungminlee.org/LatticeTools.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://kyungminlee.org/LatticeTools.jl/dev

[githubaction-img]: https://github.com/kyungminlee/LatticeTools.jl/workflows/Build/badge.svg
[githubaction-url]: https://github.com/kyungminlee/LatticeTools.jl/actions?query=workflow%3ABuild

[codecov-img]: https://codecov.io/gh/kyungminlee/LatticeTools.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/LatticeTools.jl
