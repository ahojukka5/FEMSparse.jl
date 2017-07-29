# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMSparse.jl/blob/master/LICENSE

using Base.Test
using TimerOutputs
const to = TimerOutput()

test_files = String[]
push!(test_files, "test_sparse_matrix_operations.jl")

@testset "FEMSparse.jl" begin
    for fn in test_files
        timeit(to, fn) do
            include(fn)
        end
    end
end
println()
println("Test statistics:")
println(to)
