using FEMSparse, SparseArrays, Test

K = sparse(Float64[1 0 1 1;
                   0 1 0 1;
                   1 0 1 0;
                   1 1 0 1];)

fill!(K, 0.0)

dofs1 = [1, 3]
dofs2 = [2, 4]
dofs3 = [1, 4]
Ke1 = ones(2, 2)
Ke2 = ones(2, 2)
Ke3 = ones(2, 2)
assembler = FEMSparse.start_assemble(K)
for (dofs, Ke) in zip([dofs1, dofs2, dofs3], [Ke1, Ke2, Ke3])
    FEMSparse.assemble_local_matrix!(assembler, dofs, Ke)
end

@test Matrix(K) ≈  [2.0  0.0  1.0  1.0;
                    0.0  1.0  0.0  1.0;
                    1.0  0.0  1.0  0.0;
                    1.0  1.0  0.0  2.0]

