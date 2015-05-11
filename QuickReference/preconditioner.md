Preconditioner strategies
=========================

# Relaxation

   Split into lower, diagonal, upper parts: $$ A = L + D + U $$

## Jacobi

   Cheapest preconditioner: $$P^{-1} = D^{-1}$$

##    Successive over-relaxation (SOR)

   $$
   \begin{gather*}
     \left(L + \frac 1 \omega D\right) x_{n+1} = \left[\left(\frac 1\omega-1\right)D - U\right] x_n + \omega b \\
     P^{-1} = \text{$k$ iterations starting with $x_0=0$}
   \end{gather*}
   $$

   * Implemented as a sweep
   * $$\omega = 1$$ corresponds to Gauss-Seidel
   * Very effective at removing high-frequency components of residual

   <hr>
# Factorization
   Two phases

   - symbolic factorization: find where fill occurs, only uses sparsity pattern
   - numeric factorization: compute factors

## LU decomposition

   - preconditioner
   - Expensive, for $$m\times m$$ sparse matrix with bandwidth $$b$$, traditionally requires $$\mathcal{O}(mb^2)$$ time and $$\mathcal{O}(mb)$$ space.
     - Bandwidth scales as $$m^{\frac{d-1}{d}}$$ in $$d$$-dimensions
     - Optimal in 2D: $$\mathcal{O}(m \cdot \log m)$$ space, $$\mathcal{O}(m^{3/2})$$ time
     - Optimal in 3D: $$\mathcal{O}(m^{4/3})$$ space, $$\mathcal{O}(m^2)$$ time
   - Symbolic factorization is problematic in parallel

## Incomplete LU

   - Allow a limited number of levels of fill: ILU($$k$$)
   - Only allow fill for entries that exceed threshold: ILUT
   - Usually poor scaling in parallel
   - No guarantees

# 1-level Domain decomposition

   Domain size $$L$$, subdomain size $$H$$, element size $$h$$

 * Overlapping/Schwarz
    - Solve Dirichlet problems on overlapping subdomains
    - No overlap: $$\textit{its} \in \mathcal{O}\big( \frac{L}{\sqrt{Hh}} \big)$$
    - Overlap $$\delta$$: $$\textit{its} \in \big( \frac L {\sqrt{H\delta}} \big)$$

 * Neumann-Neumann

    - Solve Neumann problems on non-overlapping subdomains
    - $$\textit{its} \in \mathcal{O}\big( \frac{L}{H}(1+\log\frac H h) \big)$$
    - Tricky null space issues (floating subdomains)
    - Need subdomain matrices, not globally assembled matrix.

   \note Multilevel variants knock off the leading $$\frac L H$$
   \note Both overlapping and nonoverlapping with this bound

 * BDDC and FETI-DP
     - Neumann problems on subdomains with coarse grid correction
     - $$\textit{its} \in \mathcal{O}\big(1 + \log\frac H h \big)$$


# Multigrid

## Introduction

   Hierarchy: Interpolation and restriction operators}
   $$ \Pi^\uparrow : X_{\text{coarse}} \to X_{\text{fine}} \qquad
   \Pi^\downarrow :  X_{\text{fine}} \to X_{\text{coarse}} $$
   - Geometric: define problem on multiple levels, use grid to compute hierarchy
   - Algebraic: define problem only on finest level, use matrix structure to build hierarchy

   Galerkin approximation
   Assemble this matrix: $$A_{\text{coarse}} = \Pi^\downarrow A_{\text{fine}} \Pi^\uparrow$$

   Application of multigrid preconditioner ($$ V $$-cycle)
    - Apply pre-smoother on fine level (any preconditioner)
    - Restrict residual to coarse level with $$\Pi^\downarrow$$
    - Solve on coarse level $$A_{\text{coarse}} x = r$$
    - Interpolate result back to fine level with $$\Pi^\uparrow$$
    - Apply post-smoother on fine level (any preconditioner)


## Multigrid convergence properties
  - Textbook: $$P^{-1}A$$ is spectrally equivalent to identity
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
