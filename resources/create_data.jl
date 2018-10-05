# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMSparse.jl/blob/master/LICENSE

using JuliaFEM, HDF5

function create_data()
    mesh = abaqus_read_mesh("EIFFEL_TOWER_TET10_220271.inp")
    tower = Problem(Elasticity, "tower", 3)
    tower_elements = create_elements(mesh, "TOWER")
    update!(tower_elements, "youngs modulus", 210.0E3)
    update!(tower_elements, "poissons ratio", 0.3)
    update!(tower_elements, "density", 7.85E-9)
    update!(tower_elements, "displacement load 3", -9810.0)
    add_elements!(tower, tower_elements)
    initialize!(tower, 0.0)
    assemble!(tower, 0.0)
    # contains 99342 x Tet10 elements
    # 30x30 matrices => length of I, J, V = 89407800

    dofs = zeros(Int, 30, length(tower_elements))
    for (i,element) in enumerate(tower_elements)
        for (j,dof) in enumerate(get_gdofs(tower, element))
            dofs[j,i] = dof
        end
    end

    K = copy(reshape(tower.assembly.K.V, 900, length(tower_elements)))

    return K, dofs
end

K, dofs = create_data()

# save = false
# if save
#     h5open("data.h5", "w") do fid
#         fid["dofs", "blosc", 9] = dofs
#         fid["stiffness", "blosc", 9] = K
#     end
# end

using SparseArrays, TimerOutputs

function create_stiffness_matrix_1(stiffness, dofs, N)
    K = SparseMatrixCOO()
    for j in 1:N
        @info("create sparse matrix #$j")
        timeit("create sparse stiffness matrix $j") do
            ndofs, nelements = size(dofs)
            @timeit "assemble" for i=1:nelements
                gdofs = @view dofs[:,i]
                Ke = @view stiffness[:,i]
                add!(K, gdofs, gdofs, Ke)
            end
            @timeit "convert to CSC" Ks = sparse(K)
            empty!(K)
        end
    end
end

function FEMBase.add!(A::SparseMatrixCSC, dofs1, dofs2, K)
    rows = rowvals(A)
    vals = nonzeros(A)
    k = 1
    for j in dofs2
        for tmp in nzrange(A, j)
            if in(rows[tmp], dofs1)
                vals[tmp] += K[k]
                k += 1
            end
        end
    end
end

function create_stiffness_matrix_2(stiffness, dofs, N)
    K1 = SparseMatrixCOO()
    local K2
    for j in 1:N
        @info("create sparse matrix #$j")
        timeit("create sparse stiffness matrix $j") do
            ndofs, nelements = size(dofs)
            if j == 1
                @timeit "assemble" for i=1:nelements
                    gdofs = @view dofs[:,i]
                    Ke = @view stiffness[:,i]
                    add!(K1, gdofs, gdofs, Ke)
                end
                @timeit "convert to CSC" K2 = sparse(K1)
                empty!(K1)
            else
                @timeit "assemble" for i=1:nelements
                    gdofs = @view dofs[:,i]
                    Ke = @view stiffness[:,i]
                    add!(K2, gdofs, gdofs, Ke)
                end
            end
        end
    end
end

@info("Warm-up")

create_stiffness_matrix_1(K, dofs, 2)
create_stiffness_matrix_2(K, dofs, 2)
using Profile
Profile.clear_malloc_data()

@info("Starting performance test")

reset_timer!()
create_stiffness_matrix_1(K, dofs, 10)
println()
println("ALGORITHM 1")
println()
print_timer()
println()
println()

reset_timer!()
create_stiffness_matrix_2(K, dofs, 10)
println()
println("ALGORITHM 2")
println()
print_timer()
println()
println()
