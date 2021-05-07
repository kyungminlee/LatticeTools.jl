<img src="https://kyungminlee.org/LatticeTools.jl/stable/assets/logo.svg" width="400px">

# LatticeTools.jl

| **Documentation** | **Build Status** | **Code Coverage** |
|:-----------------:|:----------------:|:-----------------:|
| [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] | [![Build][githubaction-img]][githubaction-url] | [![Code Coverage][codecov-main-img]][codecov-url] [![Code Coverage][codecov-dev-img]][codecov-url] |

[LatticeTools.jl](https://github.com/kyungminlee/LatticeTools.jl) is a Julia package that provides functionalities to define lattices and perform symmetry analyses useful for studying interacting quantum many-body Hamiltonians.


## Installation

[LatticeTools.jl](https://github.com/kyungminlee/LatticeTools.jl) is currently not included in Julia's default package registry. To install, add the package registry [KyungminLeeRegistry](https://github.com/kyungminlee/KyungminLeeRegistry.jl) and then install the package:

```julia-repl
pkg> registry add https://github.com/kyungminlee/KyungminLeeRegistry.git
pkg> add LatticeTools
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://kyungminlee.org/LatticeTools.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://kyungminlee.org/LatticeTools.jl/dev

[githubaction-img]: https://github.com/kyungminlee/LatticeTools.jl/workflows/Build/badge.svg
[githubaction-url]: https://github.com/kyungminlee/LatticeTools.jl/actions?query=workflow%3ABuild

[codecov-main-img]: https://codecov.io/gh/kyungminlee/LatticeTools.jl/branch/main/graph/badge.svg?token=GSYT8B0CVI
[codecov-dev-img]: https://codecov.io/gh/kyungminlee/LatticeTools.jl/branch/dev/graph/badge.svg?token=GSYT8B0CVI
[codecov-url]: https://codecov.io/gh/kyungminlee/LatticeTools.jl


