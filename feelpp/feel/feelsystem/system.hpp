/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_SYSTEM_HPP
#define FEELPP_FEELSYSTEM_SYSTEM_HPP 1

/**
 * @file system.hpp
 * @brief Umbrella header for the Feel++ Problem/System API.
 *
 * The `Feel::systems` namespace separates mathematical problem descriptors from
 * assembled or evaluable algebraic systems.
 *
 * The API provides:
 *
 * - tags such as `linear_t`, `elliptic_t`, `parabolic_t`, and
 *   `generalized_eigen_t`;
 * - concepts such as `ProblemLike`, `LinearProblemLike`,
 *   `NonlinearProblemLike`, and `GeneralizedEigenvalueProblemLike`;
 * - `Problem<Space, Tags...>` descriptors with named term roles;
 * - first-class constraint terms through `.constraint()`, `.bc()`, and
 *   `.dirichlet(...)`;
 * - natural boundary terms through `.neumann(...)`, `.robin(...)`, and
 *   `.fourier(...)`;
 * - stationary `LinearSystem` assembly and linear solves;
 * - residual/Jacobian nonlinear systems backed by `backend()->nlSolve()`;
 * - Picard fixed-point systems with per-iteration linear assembly callbacks;
 * - SLEPc-backed `EigenSystem` assembly and eigenvalue solves;
 * - BDF/Newmark `TimeStepSystem` assembly through `problem.evolve()`;
 * - block/mixed product-space descriptors backed by `blockform2` and
 *   `blockform1`;
 * - generic block roles in `Feel::systems::roles`;
 * - HDG role-space descriptors through `HdgSpaces` and
 *   `Feel::hdg::poissonSpaces`;
 * - HDG-aware `HdgProblem`/`HdgSystem` facades with role-indexed block access
 *   and typed static-condensation policy defaults;
 * - `assemble()` and `solve()` customization-point style facades.
 *
 * The fluent syntax is intentionally layered:
 *
 * @code
 * using namespace Feel;
 * using namespace Feel::vf;
 * auto uh = Xh->element( "u" );
 * auto u = trial( Xh );
 * auto v = test( Xh );
 *
 * auto p =
 *     linearProblem( _space = Xh, _name = "laplacian" )
 *         .lhs( integrate( _range = elements( mesh ),
 *                          _expr = inner( grad( u ), grad( v ) ) ) )
 *         .rhs( integrate( _range = elements( mesh ),
 *                          _expr = f * v ) )
 *         .dirichlet( _range = boundaryfaces( mesh ),
 *                     _element = uh,
 *                     _expr = g );
 *
 * auto robin =
 *     linearProblem( _space = Xh, _name = "heat-exchange" )
 *         .lhs( integrate( _range = elements( mesh ),
 *                          _expr = inner( grad( u ), grad( v ) ) ) )
 *         .fourier( _range = boundaryfaces( mesh ),
 *                   _trial = u,
 *                   _test = v,
 *                   _alpha = h,
 *                   _expr = uExt );
 *
 * auto sys = p.assemble();
 * auto result = sys.solve( _solution = uh );
 *
 * auto oneShot = p.solve( _solution = uh );
 * @endcode
 *
 * `Problem` stores linear and eigen form contributions as type-erased assembly
 * actions. Each action retains its variational expression and applies it to
 * `form2` or `form1` when the problem is assembled or reassembled.
 *
 * `Problem::assemble()` and `Problem::solve()` are convenience facades over the
 * same `systems::assemble` and `systems::solve` customization points. They do
 * not make the problem own assembled state.
 *
 * Backend ownership is explicit when needed. By default, `problem.assemble()`
 * creates or reuses a backend named from `problem.name()`, so solver CLI options
 * are naturally scoped by the problem name. Pass `_name = "other-prefix"` to
 * select a distinct derived backend, or pass `_backend = backend(...)` when an
 * existing backend and its preconditioner/cache state must be reused exactly.
 * Reassembly of a `LinearSystem` preserves its backend and algebraic object
 * handles. Applications that want CLI help, option validation, and config-file
 * parsing for a problem-derived backend prefix should register
 * `backend_options(problemName)` in `Environment` before parsing.
 *
 * Mixed/product-space systems use a block-native descriptor:
 *
 * @code
 * auto Xh = product( Vh, Qh );
 * auto U = Xh.element();
 *
 * auto p =
 *     mixedProblem( _space = Xh, _name = "stokes" )
 *         .block( 0_c, 0_c, integrate( _range = elements( mesh ),
 *                                      _expr = inner( grad( u ), grad( v ) ) ) )
 *         .block( 0_c, 1_c, integrate( _range = elements( mesh ),
 *                                      _expr = -p * div( v ) ) )
 *         .row( 0_c, integrate( _range = elements( mesh ),
 *                               _expr = inner( f, v ) ) );
 *
 * auto r = p.solve( _solution = U );
 * @endcode
 *
 * HDG and static-condensation workflows can keep their local/trace algebra
 * visible while using role-aware block selectors and a typed solve policy:
 *
 * @code
 * namespace role = Feel::systems::roles;
 * auto H = hdg::poissonSpaces<1>( mesh );
 * auto U = H.element();
 *
 * auto hdg =
 *     hdgProblem( _spaces = H, _name = "hdg-laplacian" )
 *         .condensation( hdg::condensation()
 *                            .enabled( boption( "sc.condense" ) )
 *                            .condenser( condenser_poisson() ) )
 *         .lhs( [&]( auto const& X, auto& a )
 *         {
 *             a( role::flux, role::flux ) += integrate( ... );
 *             a( role::trace, role::potential ) += integrate( ... );
 *         } )
 *         .rhs( [&]( auto const& X, auto& l )
 *         {
 *             l( role::potential ) += integrate( ... );
 *         } );
 *
 * auto result = hdg.solve( _solution = U );
 * @endcode
 *
 * `LinearSystem::matrixPtr()` and `rhsPtr()` expose the assembled base forms.
 * `activeMatrixPtr()`, `activeVectorPtr()`, and `activeSystem()` expose the
 * algebraic pair that should be solved. When constraints are present this pair
 * is produced by the same deferred constraint machinery as `BilinearFormBase`.
 *
 * \example feelsystem/qs_sys_linear.cpp
 * \example feelsystem/qs_sys_dirichlet.cpp
 * \example feelsystem/qs_sys_nonlinear.cpp
 * \example feelsystem/qs_sys_eigen.cpp
 * \example feelsystem/qs_sys_evolution.cpp
 * \example feelsystem/qs_sys_mixed.cpp
 * \example feelsystem/qs_sys_model.cpp
 */

#include <feel/feelsystem/assembly.hpp>
#include <feel/feelsystem/auxiliary.hpp>
#include <feel/feelsystem/block.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/eigen.hpp>
#include <feel/feelsystem/evolution.hpp>
#include <feel/feelsystem/execution.hpp>
#include <feel/feelsystem/hdg.hpp>
#include <feel/feelsystem/hdgspaces.hpp>
#include <feel/feelsystem/linear.hpp>
#include <feel/feelsystem/nonlinear.hpp>
#include <feel/feelsystem/picard.hpp>
#include <feel/feelsystem/problem.hpp>
#include <feel/feelsystem/result.hpp>
#include <feel/feelsystem/roles.hpp>
#include <feel/feelsystem/solve.hpp>
#include <feel/feelsystem/strategy.hpp>
#include <feel/feelsystem/tags.hpp>
#include <feel/feelsystem/terms.hpp>

// Compatibility only. New code should use Feel::systems.
#include <feel/feelsystem/systemlegacy.hpp>

#endif // FEELPP_FEELSYSTEM_SYSTEM_HPP
