# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMSparse.jl/blob/master/LICENSE

"""
Container holding all global matrices.
"""
type Assembly{Tv,Ti<:Integer}
    M :: SparseMatrixFEM{Tv, Ti} 
    K :: SparseMatrixFEM{T, Ti}
    Kg :: SparseMatrixFEM{T, Ti}
    f_int :: SparseVectorFEM{T, Ti}
    f_ext :: SparseVectorFEM{T, Ti}
    C1 :: SparseMatrixFEM{T, Ti}
    C2 :: SparseMatrixFEM{T, Ti}
    D :: SparseMatrixFEM{T, Ti}
    g :: SparseVectorFEM{T, Ti}
end

function Assembly(T=Float64)
    M = SparseMatrixFEM(Matrix{T}())
    K = SparseMatrixFEM(Matrix{T}())
    Kg = SparseMatrixFEM(Matrix{T}())
    f_int = SparseVectorFEM(Vector{T}())
    f_ext = SparseVectorFEM(Vector{T}())
    C1 = SparseMatrixFEM(Matrix{T}())
    C2 = SparseMatrixFEM(Matrix{T}())
    D = SparseMatrixFEM(Matrix{T}())
    g = SparseVectorFEM(Vector{T}())
    assembly = Assembly(M, K, Kg, f_int, f_ext, C1, C2, D, g)
    return assembly
end

function empty!(assembly::Assembly)
    empty!(assembly.M)
    empty!(assembly.K)
    empty!(assembly.Kg)
    empty!(assembly.f_int)
    empty!(assembly.f_ext)
    empty!(assembly.C1)
    empty!(assembly.C2)
    empty!(assembly.D)
    empty!(assembly.g)
end

function isempty(assembly::Assembly)
    T = isempty(assembly.M)
    T &= isempty(assembly.K)
    T &= isempty(assembly.Kg)
    T &= isempty(assembly.f_int)
    T &= isempty(assembly.f_ext)
    T &= isempty(assembly.C1)
    T &= isempty(assembly.C2)
    T &= isempty(assembly.D)
    T &= isempty(assembly.g)
    return T
end
