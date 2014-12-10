/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page Solver Algebraic solutions


  \tableofcontents

  \section SolverDefinitions Definitions

  \subsection SolverDefinitionsMatrices Matrices
  \par Matrix
  A \b matrix is a linear transformation between finite dimensional vector spaces.

  \par Assembling a matrix
  \b Assembling a matrix means defining it's action in terms of entries
  (usually stored in a sparse format).

  \par Symmetric matrix
  \f$A = A^T\f$


  \par Definite (resp. semi-definite) positive matrix
  All eigenvalue are \f$>0\f$ (resp \f$\geq 0\f$) or \f$x^TAx >0\, \forall x\f$ (resp. \f$x^TAx
  \geq 0\, \forall x\f$)

  \par Definite (resp. semi-negative matrix)
  All eigenvalue are \f$<0\f$ (resp. \f$\leq 0\f$) or \f$x^TAx <0 \forall x\f$ (resp \f$x^TAx
  \leq 0\, \forall x)\f$

  \par Indefinite matrix
  There exists positive and negative eigenvalue (Stokes, Helmholtz) or there
  exists \f$x,y\f$ such that \f$x^TAx > 0 > y^T A y\f$

  \subsection SolverDefinitionsPreconditioners Preconditioners

  \subsubsection SolverDefinitionsPreconditioning Preconditioning

  Preconditioning improves the conditioning of the Krylov operator.

  \par Left preconditioning
  \f[
  \begin{gather*}
  (P^{-1} A) x = P^{-1} b \\
  \{ P^{-1} b, (P^{-1}A) P^{-1} b, (P^{-1}A)^2 P^{-1} b, \dotsc \}
  \end{gather*}
  \f]

  \par Right preconditioning
  \f[
  \begin{gather*}
  (A P^{-1}) P x = b \\
  \{ b, (P^{-1}A)b, (P^{-1}A)^2b, \dotsc \}
  \end{gather*}
  \f]
  \note The product \f$P^{-1}A\f$ or \f$A P^{-1}\f$ is \b never formed.

  \subsubsection SolverDefinitionsPreconditioner Preconditioner

  \par Definition
  A \b preconditioner \f$\mathcal{P}\f$ is a method for constructing a
  matrix (just a linear function, not assembled!)  \f$P^{-1} = \mathcal{P}(A,A_p)\f$
  using a matrix \f$A\f$ and extra information \f$A_p\f$, such that the spectrum
  of \f$P^{-1}A\f$ (or \f$A P^{-1}\f$) is well-behaved.


  - \f$P^{-1}\f$ is dense, \f$P\f$ is often not available and is not needed
  - \f$A\f$ is rarely used by \f$\mathcal{P}\f$, but \f$A_p = A\f$ is common
  - \f$A_p\f$ is often a sparse matrix, the \e preconditioning  \e matrix
  - Matrix-based: Jacobi, Gauss-Seidel, SOR, ILU(k), LU
  - Parallel: Block-Jacobi, Schwarz, Multigrid, FETI-DP, BDDC
  - Indefinite: Schur-complement, Domain Decomposition, Multigrid

  Various preconditioning strategies are presented in
  - \subpage PreconditionerPage

  \section SolverPrinciples Principles

  \feel abstracts the PETSc library and provides a subset (sufficient in most
  cases) to the PETSc features. It interfaces with the following PETSc
  libraries: \c Mat, \c Vec, \c KSP, \c PC, \c SNES.
  - \c Vec: Vector handling library
  - \c Mat: Matrix handling library
  - \c KSP: Krylov SubSpace library implements various iterative solvers
  - \c PC: Preconditioner library implements various  preconditioning strategies
  - \c SNES: Nonlinear solver library implements various  nonlinear solve strategies

  All linear algebra are encapsulated within backend which interfaces PETSc (\c petsc),
  Eigen sparse (\c eigen) or dense (\c eigen_dense)


  The \b Default \b backend is \c petsc.

  \section SolverExamples Examples

  \subsection SolverExamplesLaplacian Laplacian

  We start with the quickstart Laplacian example, recall that we wish to, given
  a domain \f$\Omega\f$, find \f$u\f$ such that

  \f[ -\nabla \cdot (k \nabla u) = f \mbox{ in } \Omega \subset \mathbb{R}^{2}, u = g \mbox{ on } \partial \Omega \f]

  \par Monitoring KSP solvers

  \code
  feelpp_qs_laplacian --ksp-monitor=true
  \endcode

  \par Viewing KSP solvers


  \code{.sh}
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

  \endcode

  \subsubsection SolverExamplesLaplacianSolvers Solvers and preconditioners

  \subsection SolverExamplesStokes Stokes

  We now turn to the quickstart Stokes example, recall that we wish to, given a
  domain \f$\Omega\f$, find \f$(\mathbf{u},p) \f$ such that

  \f{eqnarray*}{
  -\Delta \mathbf{u} + \nabla p &= \mathbf{ f},\\
	\nabla \cdot \mathbf{u} &=    0 \mbox{ in } \Omega\\
	 \mathbf{u} &= \mathbf{g} \mbox{ on } \partial \Omega
  \f}

  This problem is indefinite, \ref SolverDefinitionsMatrices.

  - \b Strategies : Uzawa, penalty(techniques from optimisation), augmented
  lagrangian approach (Glowinski,Le Tallec)

  \note The Inf-sup condition must be satisfied. In particular for a multigrid strategy, the smoother needs to preserve it

  \subsubsection SolverExamplesStokesApproach General approach for saddle point problems

  - \b Solvers: MINRES, GMRES
  - \b Preconditioners : look first at the saddle point matrix \f$M\f$ and its
    block factorization \f$M = LDL^T\f$:
      \f[
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
        \f]
        - Elman, Silvester and Wathen propose 3 preconditioners:
    \f[
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
      \f]
    where \f$\tilde{S} \approx S^{-1} = B^T A^{-1} B\f$ and  \f$\tilde{A}^{-1}
    \approx A^{-1}\f$

  \section SolverOptions Options

*/
}
