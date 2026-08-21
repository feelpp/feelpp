/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#ifndef FEELPP_FEELOPT_SOLVEROPTIMIZATIONPETSC_HPP
#define FEELPP_FEELOPT_SOLVEROPTIMIZATIONPETSC_HPP 1

#include <exception>

#include <feel/feelopt/solveroptimization.hpp>

#if defined( FEELPP_HAS_PETSC_TAO )

#include <petsctao.h>

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>

namespace Feel
{

/** PETSc TAO implementation of the Feel++ optimization solver contract. */
template<typename T, typename SizeT = uint32_type>
class SolverOptimizationPetsc : public SolverOptimization<T, SizeT>
{
public:
    using super = SolverOptimization<T, SizeT>;
    using self_type = SolverOptimizationPetsc<T, SizeT>;
    using typename super::real_type;
    using typename super::monitor_record_type;
    using typename super::objective_gradient_data_type;
    using typename super::result_type;
    using typename super::size_type;
    using typename super::hessian_data_type;
    using typename super::sparse_matrix_ptrtype;
    using typename super::sparse_matrix_type;
    using typename super::vector_type;
    using typename super::vector_ptrtype;
    using datamap_type = typename vector_type::datamap_type;
    using datamap_ptrtype = typename vector_type::datamap_ptrtype;

    /**
     * Create and own one TAO object on @p worldComm.
     *
     * @param prefix PETSc options prefix; a trailing underscore is added when needed
     * @param worldComm communicator used to create the TAO object
     * @param vm Feel++ program-options values
     */
    explicit SolverOptimizationPetsc(
        std::string const& prefix = {},
        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
        po::variables_map const& vm = Environment::vm() );

    /** Destroy the owned TAO object exactly once. */
    ~SolverOptimizationPetsc() override;

    /** TAO ownership cannot be copied. */
    SolverOptimizationPetsc( SolverOptimizationPetsc const& ) = delete;

    /** Moving is disabled because PETSc retains the callback-context address. */
    SolverOptimizationPetsc( SolverOptimizationPetsc&& ) = delete;

    /** TAO ownership cannot be copy-assigned. */
    SolverOptimizationPetsc& operator=( SolverOptimizationPetsc const& ) = delete;

    /** Move assignment is disabled because PETSc retains the callback-context address. */
    SolverOptimizationPetsc& operator=( SolverOptimizationPetsc&& ) = delete;

    using super::solve;

    /** Apply programmatic defaults followed by prefixed PETSc runtime options. */
    void setFromOptions() override;

    /** Destroy the owned TAO object; repeated calls are safe. */
    void clear() noexcept override;

    /** Solve in place using a borrowed PETSc vector owned by @p variable. */
    [[nodiscard]] result_type solve( vector_type& variable ) override;

    /**
     * Return the borrowed native TAO handle.
     *
     * The caller must not destroy the returned object. A null handle indicates
     * that clear() has already released the solver.
     *
     * @return non-owning native TAO handle
     */
    [[nodiscard]] Tao nativeTao() const noexcept { return M_tao; }

    /**
     * Return the borrowed KSP used by the active TAO method, when available.
     *
     * The caller must not destroy the returned object. Methods without an
     * inner linear solver may return null.
     *
     * @return non-owning native KSP handle or null
     */
    [[nodiscard]] KSP nativeKsp() const;

private:
    /** State passed through PETSc's C callback interface for one solve. */
    struct CallbackContext
    {
        self_type* solver = nullptr; ///< Non-owning solver receiving the callback.
        worldcomm_ptr_t worldComm; ///< Feel++ communicator retained from the algebraic variable.
        datamap_ptrtype map; ///< Distribution map used by borrowed vector views.
        vector_type* variable = nullptr; ///< Non-owning backend variable retained during solve().
        sparse_matrix_ptrtype hessian; ///< Backend-created exact Hessian storage.
        sparse_matrix_ptrtype preconditioner; ///< Backend-created Hessian preconditioner storage.
        std::exception_ptr exception; ///< First exception captured across the C boundary.
    };

    /** PETSc C callback exposing a borrowed state Vec to the objective callable. */
    static PetscErrorCode objectiveBridge(
        Tao tao, Vec state, PetscReal* objective, void* context ) noexcept;

    /** PETSc C callback exposing borrowed state and gradient Vec objects. */
    static PetscErrorCode gradientBridge(
        Tao tao, Vec state, Vec gradient, void* context ) noexcept;

    /** PETSc C callback that exposes borrowed Vec objects as temporary Feel++ views. */
    static PetscErrorCode objectiveGradientBridge(
        Tao tao, Vec state, PetscReal* objective, Vec gradient, void* context ) noexcept;

    /** PETSc C callback exposing borrowed state and Hessian Mat objects. */
    static PetscErrorCode hessianBridge(
        Tao tao, Vec state, Mat hessian, Mat preconditioner, void* context ) noexcept;

    /** PETSc C monitor forwarding structured iteration data to C++ callables. */
    static PetscErrorCode monitorBridge( Tao tao, void* context ) noexcept;

    /**
     * Create a non-owning Feel++ view of a PETSc vector for the active WorldComm.
     *
     * @param vector borrowed native PETSc vector
     * @param context retained algebraic callback state
     * @return sequential or MPI-aware Feel++ vector view
     */
    static vector_ptrtype makeVectorView( Vec vector, CallbackContext const& context );

    /**
     * Create a non-owning Feel++ view of a PETSc matrix for the active WorldComm.
     *
     * @param matrix borrowed native PETSc matrix
     * @param rowMap Feel++ row distribution map
     * @param columnMap Feel++ column distribution map
     * @param context retained algebraic callback state
     * @return sequential or MPI-aware Feel++ sparse-matrix view
     */
    static sparse_matrix_ptrtype makeMatrixView(
        Mat matrix, datamap_ptrtype const& rowMap,
        datamap_ptrtype const& columnMap, CallbackContext const& context );

    /** Normalize a Feel++ prefix to PETSc's trailing-underscore convention. */
    static std::string normalizePrefix( std::string prefix );

    /** Map a native TAO reason to the backend-independent status enumeration. */
    static OptimizationStatus normalizeReason( TaoConvergedReason reason );

    /** Convert a native TAO convergence reason to a stable diagnostic string. */
    static std::string reasonString( TaoConvergedReason reason );

    /** Apply all explicitly configured solver defaults to the native TAO object. */
    void applyProgrammaticConfiguration();

    /** Install the selected callback mode and clear stale native callback slots. */
    void configureCallbacks();

    /** Reinstall the C++ bridge before PETSc adds option-driven monitors. */
    void configureMonitors();

    /** Install or clear the optional exact Hessian callback and matrices. */
    void configureHessian( vector_type const& variable );

    /** Install explicit or materialized variable bounds on the native TAO object. */
    void configureVariableBounds( vector_type const& variable );

    /** Validate that a bound vector has exactly the solution distribution. */
    void validateBoundLayout( vector_type const& variable, vector_type const& bound,
                              char const* role ) const;

    /** Materialize a constant or missing one-sided bound using a Feel++ vector clone. */
    Vec materializeBound( vector_type const& variable, vector_ptrtype& storage,
                          real_type value, char const* role );

    /** @return true when two Feel++ vectors have identical distributed layouts */
    static bool haveSameLayout( vector_type const& left, vector_type const& right );

    /** @return true when two data maps have identical distributed layouts */
    static bool haveSameLayout( datamap_type const& left, datamap_type const& right );

    /** Validate that a borrowed PETSc object communicator matches the solver communicator. */
    void checkCommunicator( PetscObject object, char const* role ) const;

    /** Extract a borrowed Vec without allocating or copying. */
    Vec nativeVector( vector_type const& vector, char const* role ) const;

    /** Extract a borrowed Mat without allocating or copying. */
    Mat nativeMatrix( sparse_matrix_type const& matrix, char const* role ) const;

    /** Throw a descriptive C++ exception for a nonzero PETSc error code. */
    static void checkPetsc( PetscErrorCode error, char const* operation );

    Tao M_tao = nullptr; ///< Owned native TAO handle.
    std::string M_petscPrefix; ///< Normalized PETSc runtime-options prefix.
    CallbackContext M_context; ///< Instance-local callback and exception state.
    vector_ptrtype M_materializedLowerBound; ///< Cached constant or implicit lower bound.
    vector_ptrtype M_materializedUpperBound; ///< Cached constant or implicit upper bound.
    bool M_hessianConfigured = false; ///< True after an exact Hessian was installed in TAO.
};

extern template class SolverOptimizationPetsc<double>;

} // namespace Feel

#endif // FEELPP_HAS_PETSC_TAO

#endif // FEELPP_FEELOPT_SOLVEROPTIMIZATIONPETSC_HPP
