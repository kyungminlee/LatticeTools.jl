# TightBindingLattice.jl

| **Documentation** | **Build Status** | **Code Coverage** |
|:-----------------:|:----------------:|:-----------------:|
| [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] | [![Build Status][travis-img]][travis-url] [![Build Status][appveyor-img]][appveyor-url] | [![Code Coverage][codecov-img]][codecov-url] [![Code Coverage][coveralls-img]][coveralls-url] |

[TightBindingLattice.jl](https://github.com/kyungminlee/TightBindingLattice.jl) is a Julia package which provides functionalities to define lattices and perform symmetry analysis, with focus on interacting many-body quantum Hamiltonians.


## Installation

[TightBindingLattice.jl](https://github.com/kyungminlee/TightBindingLattice.jl) is currently not included in Julia's default package registry. To install, add the package registry [KyungminLeeRegistry](https://github.com/kyungminlee/KyungminLeeRegistry.jl) and then install the package:

```julia-repl
(@v1.5) pkg> registry add https://github.com/kyungminlee/KyungminLeeRegistry.git
(@v1.5) pkg> add TightBindingLattice
```


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://kyungminlee.org/TightBindingLattice.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: http://kyungminlee.org/TightBindingLattice.jl/dev

[travis-img]: https://travis-ci.org/kyungminlee/TightBindingLattice.jl.svg?branch=master
[travis-url]: https://travis-ci.org/kyungminlee/TightBindingLattice.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/1yrosfyjvn4u61nw?svg=true
[appveyor-url]: https://ci.appveyor.com/project/kyungminlee/tightbindinglattice-jl

[codecov-img]: https://codecov.io/gh/kyungminlee/TightBindingLattice.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/TightBindingLattice.jl

[coveralls-img]: https://coveralls.io/repos/github/kyungminlee/TightBindingLattice.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/kyungminlee/TightBindingLattice.jl?branch=master
