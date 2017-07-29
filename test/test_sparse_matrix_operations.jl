# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMSparse.jl/blob/master/LICENSE

using Base.Test
using FEMSparse

@testset "add to sparse matrix" begin
    A1 = SparseMatrixFEM()
    A2 = zeros(5, 5)
    k = [1.0 -1.0; -1.0 1.0]
    for i=1:4
        # TODO: support unit range
        add!(A1, collect(i:i+1), collect(i:i+1), k)
        A2[i:i+1,i:i+1] += k
    end
    @test isapprox(full(A1), A2)
end
