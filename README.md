# FEMSparse.jl

[![][travis-img]][travis-url]
[![][pkg-0.7-img]][pkg-0.7-url]
[![][pkg-1.0-img]][pkg-1.0-url]
[![][coveralls-img]][coveralls-url]
[![][issues-img]][issues-url]

FEMSparse package contains sparse matrix operations spesifically designed for
finite element simulations. In particular, we aim to provide support for
sparse matrices which are fast to fill with dense local element matrices.
In literature, this is called to *finite element assembly procedure*, where
element local degrees of freedom are connected to the global degrees of freedom
of model. Typically this procedure looks something similar to above:

```julia
K = zeros(N, N)
Ke = [1.0 -1.0; -1.0 1.0]
dofs1 = [4, 5]
dofs2 = [4, 5]
K[dofs1, dofs2] += Ke
```

[gitter-url]: https://gitter.im/JuliaFEM/JuliaFEM.jl

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://juliafem.github.io/FEMSparse.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliafem.github.io/FEMSparse.jl/stable

[travis-img]: https://travis-ci.org/JuliaFEM/FEMSparse.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaFEM/FEMSparse.jl

[coveralls-img]: https://coveralls.io/repos/github/JuliaFEM/FEMSparse.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaFEM/FEMSparse.jl?branch=master

[issues-img]: https://img.shields.io/github/issues/JuliaFEM/FEMSparse.jl.svg
[issues-url]: https://github.com/JuliaFEM/FEMSparse.jl/issues

[pkg-0.7-img]: http://pkg.julialang.org/badges/FEMSparse_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=FEMSparse&ver=0.7
[pkg-1.0-img]: http://pkg.julialang.org/badges/FEMSparse_1.0.svg
[pkg-1.0-url]: http://pkg.julialang.org/?pkg=FEMSparse&ver=1.0
