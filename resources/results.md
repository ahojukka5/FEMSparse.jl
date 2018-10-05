# Constructing 10 CSC Sparse Matrices by assembling to SparseMatrixCOO and convert

## Initial version

### Implementation

```julia
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
```

### Results

```
─────────────────────────────────────────────────────────────────────────────────────────
                                                 Time                   Allocations
                                         ──────────────────────   ───────────────────────
            Tot / % measured:                 57.4s / 96.6%           29.5GiB / 100%

Section                          ncalls     time   %tot     avg     alloc   %tot      avg
─────────────────────────────────────────────────────────────────────────────────────────
create sparse stiffness matrix       10    55.4s   100%   5.54s   29.4GiB  100%   2.94GiB
  assemble                           10    33.0s  59.5%   3.30s   9.45GiB  32.1%  0.94GiB
  convert to CSC                     10    22.4s  40.5%   2.24s   20.0GiB  67.9%  2.00GiB
initialize sparse matrix              1   3.30μs  0.00%  3.30μs         -  0.00%        -
─────────────────────────────────────────────────────────────────────────────────────────
```


## Update 1

- Assemble COO first time, then just assemble directly to CSC form.
- Use `@view`s in implementations to reduce memory allocations.

### Implementations

```julia
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
```

### Results

First algorithm is the same than before, but now using `@view`. Memory allocations 29.5G -> 22.4G.

```
 ────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                   Allocations
                                             ──────────────────────   ───────────────────────
              Tot / % measured:                   49.3s / 96.8%           22.4GiB / 99.5%

 Section                             ncalls     time   %tot     avg     alloc   %tot      avg
 ────────────────────────────────────────────────────────────────────────────────────────────
 create sparse stiffness matrix 1         1    7.56s  15.8%   7.56s   4.17GiB  18.7%  4.17GiB
   assemble                               1    5.37s  11.3%   5.37s   2.17GiB  9.76%  2.17GiB
   convert to CSC                         1    2.18s  4.57%   2.18s   2.00GiB  8.99%  2.00GiB
 create sparse stiffness matrix 10        1    4.48s  9.39%   4.48s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.31s  4.83%   2.31s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.18s  4.56%   2.18s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 8         1    4.48s  9.39%   4.48s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.29s  4.81%   2.29s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.19s  4.58%   2.19s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 2         1    4.48s  9.38%   4.48s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.32s  4.86%   2.32s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.16s  4.53%   2.16s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 7         1    4.48s  9.38%   4.48s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.29s  4.81%   2.29s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.18s  4.57%   2.18s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 9         1    4.46s  9.34%   4.46s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.30s  4.82%   2.30s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.16s  4.52%   2.16s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 3         1    4.46s  9.34%   4.46s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.29s  4.79%   2.29s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.17s  4.54%   2.17s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 6         1    4.46s  9.33%   4.46s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.29s  4.79%   2.29s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.17s  4.55%   2.17s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 5         1    4.46s  9.33%   4.46s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.29s  4.80%   2.29s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.17s  4.54%   2.17s   9.10MiB  0.04%  9.10MiB
 create sparse stiffness matrix 4         1    4.44s  9.30%   4.44s   2.01GiB  9.03%  2.01GiB
   convert to CSC                         1    2.28s  4.77%   2.28s   2.00GiB  8.99%  2.00GiB
   assemble                               1    2.16s  4.53%   2.16s   9.10MiB  0.04%  9.10MiB
 ────────────────────────────────────────────────────────────────────────────────────────────
```

Second algorithm does not allocate after first assembly, but the total time is bigger
probably due to the very simple implementation. Better implementation of `add!` and finding
bandwidth-minimizing permutation is very likely to going make algorithm faster.

```
 ────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                   Allocations
                                             ──────────────────────   ───────────────────────
              Tot / % measured:                   59.8s / 100%            4.25GiB / 100%

 Section                             ncalls     time   %tot     avg     alloc   %tot      avg
 ────────────────────────────────────────────────────────────────────────────────────────────
 create sparse stiffness matrix 1         1    7.25s  12.2%   7.25s   4.17GiB  98.1%  4.17GiB
   assemble                               1    5.07s  8.51%   5.07s   2.17GiB  51.1%  2.17GiB
   convert to CSC                         1    2.18s  3.65%   2.18s   2.00GiB  47.0%  2.00GiB
 create sparse stiffness matrix 9         1    5.85s  9.83%   5.85s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.85s  9.83%   5.85s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 7         1    5.85s  9.81%   5.85s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.85s  9.81%   5.85s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 2         1    5.83s  9.79%   5.83s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.83s  9.79%   5.83s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 6         1    5.82s  9.78%   5.82s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.82s  9.78%   5.82s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 10        1    5.81s  9.75%   5.81s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.81s  9.75%   5.81s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 4         1    5.79s  9.73%   5.79s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.79s  9.73%   5.79s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 8         1    5.79s  9.72%   5.79s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.79s  9.72%   5.79s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 5         1    5.79s  9.72%   5.79s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.79s  9.72%   5.79s   9.10MiB  0.21%  9.10MiB
 create sparse stiffness matrix 3         1    5.78s  9.70%   5.78s   9.10MiB  0.21%  9.10MiB
   assemble                               1    5.78s  9.70%   5.78s   9.10MiB  0.21%  9.10MiB
 ────────────────────────────────────────────────────────────────────────────────────────────
```