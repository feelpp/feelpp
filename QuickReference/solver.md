Algebraic solutions
====================


#  Definitions

##  Matrices

A \b matrix is a linear transformation between finite dimensional vector spaces.

Assembling a matrix
\b Assembling a matrix means defining it's action in terms of entries
(usually stored in a sparse format).

Symmetric matrix
$$A = A^T$$


Definite (resp. semi-definite) positive matrix
All eigenvalue are $$>0$$ (resp $$\geq 0$$) or $$x^TAx >0\, \forall x$$ (resp. $$x^TAx
\geq 0\, \forall x$$)

Definite (resp. semi-negative matrix)
All eigenvalue are $$<0$$ (resp. $$\leq 0$$) or $$x^TAx <0 \forall x$$ (resp $$x^TAx
\leq 0\, \forall x)$$

Indefinite matrix
There exists positive and negative eigenvalue (Stokes, Helmholtz) or there
exists $$x,y$$ such that $$x^TAx > 0 > y^T A y$$

## Preconditioners

### Preconditioning

Preconditioning improves the conditioning of the Krylov operator.

Left preconditioning
$$
  \begin{gather*}
  (P^{-1} A) x = P^{-1} b \\
  \{ P^{-1} b, (P^{-1}A) P^{-1} b, (P^{-1}A)^2 P^{-1} b, \dotsc \}
  \end{gather*}
$$

Right preconditioning
$$
  \begin{gather*}
  (A P^{-1}) P x = b \\
  \{ b, (P^{-1}A)b, (P^{-1}A)^2b, \dotsc \}
  \end{gather*}
$$
Note that the product $$P^{-1}A$$ or $$A P^{-1}$$ is \b never formed.

### Preconditioner

Definition
A \b preconditioner $$\mathcal{P}$$ is a method for constructing a
matrix (just a linear function, not assembled!)  $$P^{-1} = \mathcal{P}(A,A_p)$$
using a matrix $$A$$ and extra information $$A_p$$, such that the spectrum
of $$P^{-1}A$$ (or $$A P^{-1}$$) is well-behaved.


  - $$P^{-1}$$ is dense, $$P$$ is often not available and is not needed
  - $$A$$ is rarely used by $$\mathcal{P}$$, but $$A_p = A$$ is common
  - $$A_p$$ is often a sparse matrix, the \e preconditioning  \e matrix
  - Matrix-based: Jacobi, Gauss-Seidel, SOR, ILU(k), LU
  - Parallel: Block-Jacobi, Schwarz, Multigrid, FETI-DP, BDDC
  - Indefinite: Schur-complement, Domain Decomposition, Multigrid

  Various preconditioning strategies are presented in
  - \subpage PreconditionerPage

# Principles 

  Feel++ abstracts the PETSc library and provides a subset (sufficient in most
  cases) to the PETSc features. It interfaces with the following PETSc
  libraries: `Mat` , `Vec` , `KSP` , `PC` , `SNES.` 
  - `Vec`  Vector handling library
  - `Mat`  Matrix handling library
  - `KSP`  Krylov SubSpace library implements various iterative solvers
  - `PC`  Preconditioner library implements various  preconditioning strategies
  - `SNES`  Nonlinear solver library implements various  nonlinear solve strategies

  All linear algebra are encapsulated within backend which interfaces PETSc (`petsc)` ,
  Eigen sparse (`eigen)`  or dense (`eigen_dense)` 


  The \b Default \b backend is `petsc.` 

# Examples

## Laplacian

  We start with the quickstart Laplacian example, recall that we wish to, given
  a domain $$\Omega$$, find $$u$$ such that

  $$ -\nabla \cdot (k \nabla u) = f \mbox{ in } \Omega \subset \mathbb{R}^{2}, u = g \mbox{ on } \partial \Omega $$

  Monitoring KSP solvers

  ```cpp
  feelpp_qs_laplacian --ksp-monitor=true
  ```

  Viewing KSP solvers


  ```cpp{.sh}
shell> mpirun -np 2 feelpp_qs_laplacian --ksp-monitor=1  --ksp-view=1
  0 KSP Residual norm 8.953261456448e-01
  1 KSP Residual norm 7.204431786960e-16
KSP Object: 2 MPI processes
  type: gmres
    GMRES: restart=30, using Classical (unmodified) Gram-Schmidt
     Orthogonalization with no iterative refinement
    GMRES: happy breakdown tolerance 1e-30
  maximum iterations=1000
  tolerances:  relative=1e-13, absolute=1e-50, divergence=100000
  left preconditioning
  using nonzero initial guess
  using PRECONDITIONED norm type for convergence test
PC Object: 2 MPI processes
  type: shell
    Shell:
  linear system matrix = precond matrix:
  Matrix Object:   2 MPI processes
    type: mpiaij
    rows=525, cols=525
    total: nonzeros=5727, allocated nonzeros=5727
    total number of mallocs used during MatSetValues calls =0
      not using I-node (on process 0) routines

  ```

### Solvers and preconditioners

## Stokes

  We now turn to the quickstart Stokes example, recall that we wish to, given a
  domain $$\Omega$$, find $$(\mathbf{u},p) $$ such that

  $$
  \begin{eqnarray*}
  -\Delta \mathbf{u} + \nabla p &= \mathbf{ f},\\
  \nabla \cdot \mathbf{u} &=    0 \mbox{ in } \Omega\\
  \mathbf{u} &= \mathbf{g} \mbox{ on } \partial \Omega
  \end{eqnarray*}     
  $$

  This problem is indefinite, \ref SolverDefinitionsMatrices.

  - \b Strategies : Uzawa, penalty(techniques from optimisation), augmented
  lagrangian approach (Glowinski,Le Tallec)

  Note that The Inf-sup condition must be satisfied. In particular for a multigrid strategy, the smoother needs to preserve it

### General approach for saddle point problems

  - \b Solvers: MINRES, GMRES
  - \b Preconditioners : look first at the saddle point matrix $$M$$ and its
    block factorization $$M = LDL^T$$:
      $$
        \label{eq:2}
        M =
        \begin{pmatrix}
          A & B\\
          B^T & 0
        \end{pmatrix}
        =
        \begin{pmatrix}
          I & 0\\
          B^T C & I
        \end{pmatrix}
        \begin{pmatrix}
          A & 0\\
          0 & - B^T A^{-1} B
        \end{pmatrix}
        \begin{pmatrix}
          I & A^{-1} B\\
          0 & I
        \end{pmatrix}
        $$
        - Elman, Silvester and Wathen propose 3 preconditioners:
    $$
      \label{eq:3}
      P_1 =
      \begin{pmatrix}
        \tilde{A}^{-1} & B\\
        B^T & 0
      \end{pmatrix}, \quad
      P_2 =
      \begin{pmatrix}
        \tilde{A}^{-1} & 0\\
        0 & \tilde{S}
      \end{pmatrix},\quad
      P_3 =
      \begin{pmatrix}
        \tilde{A}^{-1} & B\\
        0 & \tilde{S}
      \end{pmatrix}
      $$
    where $$\tilde{S} \approx S^{-1} = B^T A^{-1} B$$ and  $$\tilde{A}^{-1}
    \approx A^{-1}$$

# Options

