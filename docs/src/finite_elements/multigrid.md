# Multigrid Solver

## Intergrid Operations

```@docs
prolongate!
prolongate_E!
prolongate_vector_Nxe
restrict!
restrict_V!
restrict_E!
```

## MGM-CN Solver 

```@docs
cn_mgm
cn_mgm_cycle!
cn_coarse_grid_correction!
cn_coarse_grid_correction_solve!
cn_smooth_jacobi!
cn_fill_jacobi_edges!
cn_fill_jacobi_vertices!
cn_weighted_degree!
cn_matrix_free_jacobi!
cn_mat_mul!
cn_mat_mul_VV!
cn_mat_mul_VE!
cn_mat_mul_EV_EE!
cn_mat_mul_e!
```

## MG Nested Iteration Solver (Elliptic)

```@docs
solve_mgm
MGM_matrix_free_jacobi!
mat_mul_M!
mat_mul_M_VE!
mat_mul_M_EE_EV!
mat_mul_M_EE!
mat_mul_H!
mat_mul_H_VE!
mat_mul_H_EE_EV!
mat_mul_H_EE!
smooth_jacobi!
matrix_free_jacobi!
fill_jacobi_vertices!
fill_jacobi_edges!
```