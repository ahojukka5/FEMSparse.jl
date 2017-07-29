# FEMSparse.jl

Sparse matrix operations spesifically designed for FE calculations.

[![Build Status](https://travis-ci.org/JuliaFEM/FEMSparse.jl.svg?branch=master)](https://travis-ci.org/JuliaFEM/FEMSparse.jl)[![Coverage Status](https://coveralls.io/repos/github/JuliaFEM/FEMSparse.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaFEM/FEMSparse.jl?branch=master)[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliafem.github.io/FEMSparse.jl/stable)[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliafem.github.io/FEMSparse.jl/latest)

Design goals:
- fast add of local matrices / vector to global matrices / vectors
- fast removal of rows / columns
- possible to assemble using threads
