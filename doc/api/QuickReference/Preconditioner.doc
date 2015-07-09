/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*!
   \page PreconditionerPage  Preconditioner strategies

   \tableofcontents

   <hr>

   \section PreconditionerRelaxation Relaxation

   Split into lower, diagonal, upper parts: \f$ A = L + D + U \f$

   \subsection PreconditionerRelaxationJacobi Jacobi

   Cheapest preconditioner: \f$P^{-1} = D^{-1}\f$

   \subsection PreconditionerRelaxationSOR    Successive over-relaxation (SOR)

   \f[
   \begin{gather*}
     \left(L + \frac 1 \omega D\right) x_{n+1} = \left[\left(\frac 1\omega-1\right)D - U\right] x_n + \omega b \\
     P^{-1} = \text{$k$ iterations starting with $x_0=0$}
   \end{gather*}
   \f]

   \li Implemented as a sweep
   \li \f$\omega = 1\f$ corresponds to Gauss-Seidel
   \li Very effective at removing high-frequency components of residual

   <hr>
   \section PreconditionerFactorization Factorization
   Two phases

   - symbolic factorization: find where fill occurs, only uses sparsity pattern
   - numeric factorization: compute factors

   \subsection PreconditionerFactorizationLU LU decomposition

   - preconditioner
   - Expensive, for \f$m\times m\f$ sparse matrix with bandwidth \f$b\f$, traditionally requires \f$\mathcal{O}(mb^2)\f$ time and \f$\mathcal{O}(mb)\f$ space.
     - Bandwidth scales as \f$m^{\frac{d-1}{d}}\f$ in \f$d\f$-dimensions
     - Optimal in 2D: \f$\mathcal{O}(m \cdot \log m)\f$ space, \f$\mathcal{O}(m^{3/2})\f$ time
     - Optimal in 3D: \f$\mathcal{O}(m^{4/3})\f$ space, \f$\mathcal{O}(m^2)\f$ time
   - Symbolic factorization is problematic in parallel

   \subsection PreconditionerFactorizationILU Incomplete LU

   - Allow a limited number of levels of fill: ILU(\f$k\f$)
   - Only allow fill for entries that exceed threshold: ILUT
   - Usually poor scaling in parallel
   - No guarantees

   \section PreconditionerDD 1-level Domain decomposition

   Domain size \f$L\f$, subdomain size \f$H\f$, element size \f$h\f$

   \par Overlapping/Schwarz
    - Solve Dirichlet problems on overlapping subdomains
    - No overlap: \f$\textit{its} \in \mathcal{O}\big( \frac{L}{\sqrt{Hh}} \big)\f$
    - Overlap \f$\delta\f$: \f$\textit{its} \in \big( \frac L {\sqrt{H\delta}} \big)\f$

   \par Neumann-Neumann

    - Solve Neumann problems on non-overlapping subdomains
    - \f$\textit{its} \in \mathcal{O}\big( \frac{L}{H}(1+\log\frac H h) \big)\f$
    - Tricky null space issues (floating subdomains)
    - Need subdomain matrices, not globally assembled matrix.

   \note Multilevel variants knock off the leading \f$\frac L H\f$
   \note Both overlapping and nonoverlapping with this bound

  \par BDDC and FETI-DP
     - Neumann problems on subdomains with coarse grid correction
     - \f$\textit{its} \in \mathcal{O}\big(1 + \log\frac H h \big)\f$


   \section PreconditionerMultigrid Multigrid

   \subsection PreconditionerMultigridIntro Introduction

   \par Hierarchy: Interpolation and restriction operators}
   \f[ \Pi^\uparrow : X_{\text{coarse}} \to X_{\text{fine}} \qquad
   \Pi^\downarrow :  X_{\text{fine}} \to X_{\text{coarse}} \f]
   - Geometric: define problem on multiple levels, use grid to compute hierarchy
   - Algebraic: define problem only on finest level, use matrix structure to build hierarchy

   \par Galerkin approximation
   Assemble this matrix: \f[A_{\text{coarse}} = \Pi^\downarrow A_{\text{fine}} \Pi^\uparrow\f]

   \par Application of multigrid preconditioner (\f$ V \f$-cycle)
    - Apply pre-smoother on fine level (any preconditioner)
    - Restrict residual to coarse level with \f$\Pi^\downarrow\f$
    - Solve on coarse level \f$A_{\text{coarse}} x = r\f$
    - Interpolate result back to fine level with \f$\Pi^\uparrow\f$
    - Apply post-smoother on fine level (any preconditioner)


 \subsection PreconditionerMultigridProps Multigrid convergence properties
  - Textbook: \f$P^{-1}A\f$ is spectrally equivalent to identity
    - Constant number of iterations to converge up to discretization error
  - Most theory applies to SPD systems
    - variable coefficients (\eg discontinuous): low energy interpolants
    - mesh- and/or physics-induced anisotropy: semi-coarsening/line smoothers
    - complex geometry: difficult to have meaningful coarse levels
  - Deeper algorithmic difficulties
    - nonsymmetric (\eg advection, shallow water, Euler)
    - indefinite (\eg incompressible flow, Helmholtz)
  - Performance considerations
    - Aggressive coarsening is critical in parallel
    - Most theory uses SOR smoothers, ILU often more robust
    - Coarsest level usually solved semi-redundantly with direct solver
  - Multilevel Schwarz is essentially the same with different language
    - assume strong smoothers, emphasize aggressive coarsening
*/
}
