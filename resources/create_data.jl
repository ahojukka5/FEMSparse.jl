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

save = false
if save
    h5open("data.h5", "w") do fid
        fid["dofs", "blosc", 9] = dofs
        fid["stiffness", "blosc", 9] = K
    end
end

using SparseArrays, TimerOutputs

function create_stiffness_matrix(stiffness, dofs, N)
    @timeit "initialize sparse matrix" K = SparseMatrixCOO()
    for _ in 1:N
        @timeit "create sparse stiffness matrix" begin
            empty!(K)
            ndofs, nelements = size(dofs)
            Ke = zeros(30,30)
            gdofs = zeros(Int, 30)
            @timeit "assemble" for i=1:nelements
                gdofs[:] .= dofs[:,i]
                Ke[:] .= stiffness[:,i]
                add!(K, gdofs, gdofs, Ke)
            end
            @timeit "convert to CSC" Ks = sparse(K)
        end
    end
end

reset_timer!()
create_stiffness_matrix(K, dofs, 10)
print_timer()
