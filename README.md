# FEM Solver for Saint-Venant Warping Problem

## Overview

This project presents a finite element implementation of the Saint-Venant torsion problem for a rectangular cross-section. 

The objective is to compute the warping function using linear triangular finite elements (CST — Constant Strain Triangle). 

The solver is written in Fortran 90 and uses:

- Linear triangular finite elements (P1/CST)
- Variational formulation of the Laplace equation
- LU factorization from LAPACK
- Gmsh mesh input/output
- Python post-processing for visualization

---

# Physical Problem

For a beam subjected to torsion, the Saint-Venant formulation introduces a scalar warping function:

$$\omega(y,z)$$

defined on the cross-section $\Omega$.

The governing problem is:

$$
\begin{cases}
\Delta \omega = 0 & \text{in } \Omega \\
\dfrac{\partial \omega}{\partial n} = z n_y - y n_z & \text{on } \partial\Omega
\end{cases}
$$

This corresponds to a pure Neumann problem.

---

# Numerical Method

The domain is discretized using:

- 3-node triangular elements (CST)
- Linear interpolation of the warping field
- Constant gradient per element

For each element:

$$K^e = S^e B^T B$$

$$q^e = S^e B^T \begin{pmatrix} -z_c \\ y_c \end{pmatrix}$$

The global system is assembled as:

$$K \omega = f$$

A single nodal constraint is added to remove the singularity associated with the pure Neumann problem.

---

# Project Structure

```text
.
├── mesh.f90           # Mesh reader and utilities
├── TP_elas.f90        # Main FEM solver
├── MAIL.msh           # Gmsh mesh file
├── resu.msh           # FEM results exported to Gmsh
└── rapport.pdf        # Final report
```
---

# Compilation

The code requires:

- gfortran
- LAPACK
- BLAS/OpenBLAS

Compilation:


```bash
gfortran -c mesh.f90
gfortran -c TP_elas.f90
gfortran mesh.f90 TP_elas.f90 -llapack -lblas -o solveur

```
Or:
```bash
gfortran *.f90 -llapack -lopenblas -o solveur
```

On HPC clusters using modules:

```bash
module load OpenBLAS
gfortran -c mesh.f90
gfortran -c TP_elas.f90
gfortran *.f90 -llapack -lopenblas -o solveur
```

---

# Execution

Run the solver with:

```bash
./solveur
```

The program generates:

```text
resu.msh
```

containing:

- mesh
- nodal values of the warping function
- FEM solution exported in Gmsh format

---

# Main Results

The numerical solution reproduces the expected behavior:

- antisymmetry of the warping function
- saddle-shaped 3D profile
- dominant shear stress along the long direction
- physically consistent Saint-Venant torsion behavior

---

# Author

**Oumkalthoum M'HAMDI**

Finite Element Methods Project  
Academic Year 2024–2025
