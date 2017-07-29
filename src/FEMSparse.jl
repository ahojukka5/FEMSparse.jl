# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMSparse.jl/blob/master/LICENSE

module FEMSparse

include("sparsematrixcoo.jl")
include("sparsevectorcoo.jl")
include("sparsevectordok.jl")

SparseMatrixFEM{Tv,Ti<:Integer} = SparseMatrixCOO{Tv,Ti}
SparseVectorFEM{Tv,Ti<:Integer} = SparseVectorCOO{Tv,Ti}
#include("assembly.jl")

export SparseMatrixFEM,
       add!,
       get_nonzero_rows,
       get_nonzero_columns,
       remove_row!,
       remove_column!

end
