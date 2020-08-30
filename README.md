<img src="https://kyungminlee.org/LatticeTools.jl/stable/assets/logo.png" width="400px">

# LatticeTools.jl

| **Documentation** | **Build Status** | **Code Coverage** |
|:-----------------:|:----------------:|:-----------------:|
| [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] | [![Build Status][travis-img]][travis-url] [![Build Status][appveyor-img]][appveyor-url] | [![Code Coverage][codecov-img]][codecov-url] [![Code Coverage][coveralls-img]][coveralls-url] |

[LatticeTools.jl](https://github.com/kyungminlee/LatticeTools.jl) is a Julia package which provides functionalities to define lattices and perform symmetry analysis, with focus on interacting many-body quantum Hamiltonians.


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

[travis-img]: https://travis-ci.org/kyungminlee/LatticeTools.jl.svg?branch=master
[travis-url]: https://travis-ci.org/kyungminlee/LatticeTools.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/1yrosfyjvn4u61nw?svg=true
[appveyor-url]: https://ci.appveyor.com/project/kyungminlee/LatticeTools-jl

[codecov-img]: https://codecov.io/gh/kyungminlee/LatticeTools.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/LatticeTools.jl

[coveralls-img]: https://coveralls.io/repos/github/kyungminlee/LatticeTools.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/kyungminlee/LatticeTools.jl?branch=master
