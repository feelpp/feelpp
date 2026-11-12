/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_HDG_HPP
#define FEELPP_FEELSYSTEM_HDG_HPP 1

#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

#include <feel/feelalg/backend.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelsystem/block.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/hdgpolicy.hpp>
#include <feel/feelsystem/hdgspaces.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/roles.hpp>

namespace Feel::systems
{

namespace detail
{
/**
 * @brief Role-aware view over a block bilinear form.
 *
 * The view maps semantic roles to the underlying product-space indices and
 * otherwise forwards to the raw block form. It does not own algebraic state.
 */
template <typename SpacesType, typename FormType>
class HdgBilinearFormView
{
  public:
    using spaces_type = std::remove_cvref_t<SpacesType>;
    using form_type = std::remove_cvref_t<FormType>;

    /**
     * @brief Bind a role layout and a live block form.
     */
    HdgBilinearFormView( spaces_type const& spaces, form_type& form )
        : M_spaces( std::addressof( spaces ) ),
          M_form( std::addressof( form ) )
    {
    }

    /**
     * @brief Raw block form for expert operations not wrapped by this view.
     */
    form_type& form() noexcept { return *M_form; }

    /**
     * @brief Raw block form for expert operations not wrapped by this view.
     */
    form_type const& form() const noexcept { return *M_form; }

    /**
     * @brief Implicit raw-form access for existing helper callbacks.
     */
    operator form_type&() noexcept { return *M_form; }

    /**
     * @brief Access a block by semantic row and column roles.
     */
    template <BlockRole RowRole, BlockRole ColumnRole>
    decltype( auto ) operator()( RowRole row, ColumnRole column )
    {
        return ( *M_form )( M_spaces->index( row ), M_spaces->index( column ) );
    }

    /**
     * @brief Access a block by semantic row role and existing column selector.
     */
    template <BlockRole RowRole, typename Column>
        requires( !BlockRole<Column> )
    decltype( auto ) operator()( RowRole row, Column&& column )
    {
        return ( *M_form )( M_spaces->index( row ), std::forward<Column>( column ) );
    }

    /**
     * @brief Access a block by existing row selector and semantic column role.
     */
    template <typename Row, BlockRole ColumnRole>
        requires( !BlockRole<Row> )
    decltype( auto ) operator()( Row&& row, ColumnRole column )
    {
        return ( *M_form )( std::forward<Row>( row ), M_spaces->index( column ) );
    }

    /**
     * @brief Forward existing non-role block selectors unchanged.
     */
    template <typename Row, typename Column>
        requires( !BlockRole<Row> && !BlockRole<Column> )
    decltype( auto ) operator()( Row&& row, Column&& column )
    {
        return ( *M_form )( std::forward<Row>( row ), std::forward<Column>( column ) );
    }

    /**
     * @brief Access a row constraint view by semantic role.
     */
    template <BlockRole RowRole>
    decltype( auto ) row( RowRole row )
    {
        return M_form->row( M_spaces->index( row ) );
    }

    /**
     * @brief Forward an existing row selector unchanged.
     */
    template <typename Row>
        requires( !BlockRole<Row> )
    decltype( auto ) row( Row&& row )
    {
        return M_form->row( std::forward<Row>( row ) );
    }

    /**
     * @brief Forward whole-form additions to the raw block form.
     */
    template <typename Expr>
    decltype( auto ) operator+=( Expr&& expr )
    {
        return ( *M_form += std::forward<Expr>( expr ) );
    }

  private:
    /** @brief Role/product-space descriptor used to map roles to block indices. */
    spaces_type const* M_spaces = nullptr;
    /** @brief Non-owning pointer to the live block bilinear form. */
    form_type* M_form = nullptr;
};

/**
 * @brief Role-aware view over a block linear form.
 */
template <typename SpacesType, typename FormType>
class HdgLinearFormView
{
  public:
    using spaces_type = std::remove_cvref_t<SpacesType>;
    using form_type = std::remove_cvref_t<FormType>;

    /**
     * @brief Bind a role layout and a live block linear form.
     */
    HdgLinearFormView( spaces_type const& spaces, form_type& form )
        : M_spaces( std::addressof( spaces ) ),
          M_form( std::addressof( form ) )
    {
    }

    /**
     * @brief Raw block form for expert operations not wrapped by this view.
     */
    form_type& form() noexcept { return *M_form; }

    /**
     * @brief Raw block form for expert operations not wrapped by this view.
     */
    form_type const& form() const noexcept { return *M_form; }

    /**
     * @brief Implicit raw-form access for existing helper callbacks.
     */
    operator form_type&() noexcept { return *M_form; }

    /**
     * @brief Access a vector block by semantic role.
     */
    template <BlockRole RowRole>
    decltype( auto ) operator()( RowRole row )
    {
        return ( *M_form )( M_spaces->index( row ) );
    }

    /**
     * @brief Forward existing non-role selectors unchanged.
     */
    template <typename Row>
        requires( !BlockRole<Row> )
    decltype( auto ) operator()( Row&& row )
    {
        return ( *M_form )( std::forward<Row>( row ) );
    }

    /**
     * @brief Zero a vector block selected by semantic role.
     */
    template <BlockRole RowRole>
    void zeroBlock( RowRole row )
    {
        M_form->zeroBlock( int( M_spaces->index( row ) ) );
    }

    /**
     * @brief Forward whole-form additions to the raw block form.
     */
    template <typename Expr>
    decltype( auto ) operator+=( Expr&& expr )
    {
        return ( *M_form += std::forward<Expr>( expr ) );
    }

  private:
    /** @brief Role/product-space descriptor used to map roles to block indices. */
    spaces_type const* M_spaces = nullptr;
    /** @brief Non-owning pointer to the live block linear form. */
    form_type* M_form = nullptr;
};

/**
 * @brief Invoke an HDG matrix assembly action using the best supported signature.
 */
template <typename Action, typename SpacesType, typename ViewType, typename RawFormType>
void invokeHdgMatrixAction( Action& action, SpacesType const& spaces, ViewType& view, RawFormType& rawForm )
{
    if constexpr ( requires { action( spaces, view ); } )
        action( spaces, view );
    else if constexpr ( requires { action( view ); } )
        action( view );
    else if constexpr ( !ExpressionLike<Action> && requires { action( rawForm ); } )
        action( rawForm );
    else
        rawForm += action;
}

/**
 * @brief Invoke an HDG vector assembly action using the best supported signature.
 */
template <typename Action, typename SpacesType, typename ViewType, typename RawFormType>
void invokeHdgVectorAction( Action& action, SpacesType const& spaces, ViewType& view, RawFormType& rawForm )
{
    if constexpr ( requires { action( spaces, view ); } )
        action( spaces, view );
    else if constexpr ( requires { action( view ); } )
        action( view );
    else if constexpr ( !ExpressionLike<Action> && requires { action( rawForm ); } )
        action( rawForm );
    else
        rawForm += action;
}

/**
 * @brief Invoke an HDG constraint action using role-aware or raw block forms.
 */
template <typename Constraint, typename SpacesType,
          typename MatrixViewType, typename VectorViewType,
          typename RawMatrixFormType, typename RawVectorFormType>
void invokeHdgConstraintAction( Constraint& constraint, SpacesType const& spaces,
                                MatrixViewType& matrixView, VectorViewType& vectorView,
                                RawMatrixFormType& rawMatrix, RawVectorFormType& rawVector )
{
    if constexpr ( requires { constraint( spaces, matrixView, vectorView ); } )
        constraint( spaces, matrixView, vectorView );
    else if constexpr ( requires { constraint( matrixView, vectorView ); } )
        constraint( matrixView, vectorView );
    else
        constraint( rawMatrix, rawVector );
}
} // namespace detail

template <typename SpacesType, typename CondensationPolicyType = HdgCondensationPolicy<>>
class HdgSystem;

/**
 * @brief HDG-aware linear block problem descriptor.
 *
 * `HdgProblem` owns semantic HDG spaces, a typed condensation policy, and an
 * internal `BlockProblem` that stores replayable assembly callbacks. It makes
 * block roles and solve defaults explicit, but the user still writes the HDG
 * weak form block by block.
 *
 * @tparam SpacesType `HdgSpaces`-like role/product-space descriptor.
 * @tparam CondensationPolicyType Static-condensation policy type.
 */
template <typename SpacesType, typename CondensationPolicyType = HdgCondensationPolicy<>>
class HdgProblem
{
  public:
    using spaces_type = std::remove_cvref_t<SpacesType>;
    using condensation_policy_type = std::remove_cvref_t<CondensationPolicyType>;
    using product_space_type = typename spaces_type::product_space_type;
    using space_type = product_space_type;
    using value_type = typename spaces_type::value_type;
    using element_type = typename spaces_type::element_type;
    using backend_type = Backend<value_type>;
    using backend_ptrtype = std::shared_ptr<backend_type>;
    using block_problem_type = BlockProblem<product_space_type,
                                            linear_t,
                                            stationary_t,
                                            block_t,
                                            hdg_t>;
    using block_bilinear_form_type = typename block_problem_type::block_bilinear_form_type;
    using block_linear_form_type = typename block_problem_type::block_linear_form_type;
    using term_container_type = typename block_problem_type::term_container_type;
    using system_type = HdgSystem<spaces_type, condensation_policy_type>;

    template <typename Tag>
    static constexpr bool has_tag = block_problem_type::template has_tag<Tag>;

    /**
     * @brief Construct an HDG problem from spaces and a problem/backend name.
     */
    explicit HdgProblem( spaces_type spaces, std::string name = {} )
        : M_spaces( std::move( spaces ) ),
          M_name( std::move( name ) ),
          M_blockProblem( M_spaces.product(), M_name )
    {
    }

    /**
     * @brief Construct an HDG problem while preserving an existing block body.
     */
    HdgProblem( spaces_type spaces,
                std::string name,
                block_problem_type blockProblem,
                condensation_policy_type policy )
        : M_spaces( std::move( spaces ) ),
          M_name( std::move( name ) ),
          M_blockProblem( std::move( blockProblem ) ),
          M_condensation( std::move( policy ) )
    {
        M_blockProblem.strategy( M_condensation.strategy() );
    }

    /**
     * @brief Semantic HDG spaces with role metadata.
     */
    spaces_type const& spaces() const noexcept { return M_spaces; }

    /**
     * @brief Mutable semantic HDG spaces.
     */
    spaces_type& spaces() noexcept { return M_spaces; }

    /**
     * @brief Product space used by the underlying block forms.
     */
    product_space_type const& space() const noexcept { return M_blockProblem.space(); }

    /**
     * @brief Mutable product space used by the underlying block forms.
     */
    product_space_type& space() noexcept { return M_blockProblem.space(); }

    /**
     * @brief Problem/backend option prefix.
     */
    std::string const& name() const noexcept { return M_name; }

    /**
     * @brief Internal block problem used for assembly storage.
     */
    block_problem_type const& blockProblem() const noexcept { return M_blockProblem; }

    /**
     * @brief Mutable internal block problem for expert workflows.
     */
    block_problem_type& blockProblem() noexcept { return M_blockProblem; }

    /**
     * @brief Static-condensation policy.
     */
    condensation_policy_type const& condensationPolicy() const noexcept { return M_condensation; }

    /**
     * @brief Default block assembly strategy.
     */
    solve::strategy strategy() const noexcept { return M_blockProblem.strategy(); }

    /**
     * @brief Set the default block assembly strategy.
     */
    HdgProblem& strategy( solve::strategy value ) noexcept
    {
        M_condensation.strategy( value );
        M_blockProblem.strategy( value );
        return *this;
    }

    /**
     * @brief Enable or disable static condensation and update block strategy.
     */
    HdgProblem& condense( bool value ) noexcept
    {
        M_condensation.condense( value );
        M_blockProblem.strategy( M_condensation.strategy() );
        return *this;
    }

    /**
     * @brief Set the final condensed backend name.
     */
    HdgProblem& backendName( std::string name )
    {
        M_condensation.backendName( std::move( name ) );
        return *this;
    }

    /**
     * @brief Return a copy using a different condenser tag.
     */
    template <CondenserTag NewCondenserType>
    [[nodiscard]] auto condenser( NewCondenserType condenser ) const
    {
        return this->condensation( M_condensation.condenser( std::move( condenser ) ) );
    }

    /**
     * @brief Return a copy using a full condensation policy.
     */
    template <typename PolicyType>
    [[nodiscard]] auto condensation( PolicyType&& policy ) const
    {
        using new_policy_type = std::remove_cvref_t<PolicyType>;
        return HdgProblem<spaces_type, new_policy_type>(
            M_spaces, M_name, M_blockProblem, std::forward<PolicyType>( policy ) );
    }

    /**
     * @brief Final backend name used by the condensed trace/global solve.
     */
    std::string condensedBackendName() const
    {
        return M_condensation.backendName( M_name );
    }

    /**
     * @brief Default backend name used by HDG postprocessing solves.
     */
    std::string postBackendName() const
    {
        return M_condensation.postBackendName( M_name );
    }

    /**
     * @brief Assemble this descriptor into an `HdgSystem`.
     */
    template <typename... Args>
    auto assemble( Args&&... args ) const
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        auto strategy = namedArgs.get_else( _strategy, this->strategy() );
        auto system = system_type( *this, std::move( backend ), std::move( backendName ), strategy );
        system.assemble();
        return system;
    }

    /**
     * @brief Assemble and solve this HDG descriptor.
     */
    template <typename... Args>
    auto solve( Args&&... args ) const
    {
        auto system = this->assemble( args... );
        return system.solve( std::forward<Args>( args )... );
    }

    /**
     * @brief Add a whole-block matrix assembly action.
     *
     * Supported callback signatures are `(spaces, a)`, `(a)`, and the raw
     * block-form `(rawA)`. Non-callback form expressions are forwarded to the
     * raw block form.
     */
    template <typename Action>
    HdgProblem& lhs( Action&& action )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        M_blockProblem.lhs(
            [spaces = M_spaces, storedAction]( block_bilinear_form_type& rawA )
            {
                detail::HdgBilinearFormView<spaces_type, block_bilinear_form_type> a( spaces, rawA );
                detail::invokeHdgMatrixAction( *storedAction, spaces, a, rawA );
            } );
        return *this;
    }

    /**
     * @brief Add a whole-block vector assembly action.
     */
    template <typename Action>
    HdgProblem& rhs( Action&& action )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        M_blockProblem.rhs(
            [spaces = M_spaces, storedAction]( block_linear_form_type& rawL )
            {
                detail::HdgLinearFormView<spaces_type, block_linear_form_type> l( spaces, rawL );
                detail::invokeHdgVectorAction( *storedAction, spaces, l, rawL );
            } );
        return *this;
    }

    /**
     * @brief Add a constraint action over live role-aware block forms.
     */
    template <typename Constraint>
    HdgProblem& constraint( Constraint&& constraint )
    {
        auto storedConstraint = std::make_shared<std::decay_t<Constraint>>( std::forward<Constraint>( constraint ) );
        M_blockProblem.constraint(
            [spaces = M_spaces, storedConstraint]( block_bilinear_form_type& rawA, block_linear_form_type& rawL )
            {
                detail::HdgBilinearFormView<spaces_type, block_bilinear_form_type> a( spaces, rawA );
                detail::HdgLinearFormView<spaces_type, block_linear_form_type> l( spaces, rawL );
                detail::invokeHdgConstraintAction( *storedConstraint, spaces, a, l, rawA, rawL );
            } );
        return *this;
    }

    /**
     * @brief Alias for `.constraint(...)`.
     */
    template <typename Constraint>
    HdgProblem& bc( Constraint&& constraint )
    {
        return this->constraint( std::forward<Constraint>( constraint ) );
    }

    /**
     * @brief Add a contribution to one matrix block selected by roles.
     */
    template <BlockRole RowRole, BlockRole ColumnRole, typename Term>
    HdgProblem& block( RowRole row, ColumnRole column, Term&& term )
    {
        M_blockProblem.block( M_spaces.index( row ), M_spaces.index( column ), std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add a contribution to one matrix block selected by existing indices.
     */
    template <typename Row, typename Column, typename Term>
        requires( !BlockRole<Row> || !BlockRole<Column> )
    HdgProblem& block( Row&& row, Column&& column, Term&& term )
    {
        M_blockProblem.block( std::forward<Row>( row ), std::forward<Column>( column ), std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add a contribution to one rhs block selected by role.
     */
    template <BlockRole RowRole, typename Term>
    HdgProblem& row( RowRole row, Term&& term )
    {
        M_blockProblem.row( M_spaces.index( row ), std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add a contribution to one rhs block selected by existing index.
     */
    template <typename Row, typename Term>
        requires( !BlockRole<Row> )
    HdgProblem& row( Row&& row, Term&& term )
    {
        M_blockProblem.row( std::forward<Row>( row ), std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Matrix term metadata.
     */
    auto lhsTerms() const noexcept { return M_blockProblem.lhsTerms(); }

    /**
     * @brief Rhs term metadata.
     */
    auto rhsTerms() const noexcept { return M_blockProblem.rhsTerms(); }

    /**
     * @brief Constraint term metadata.
     */
    auto constraintTerms() const noexcept { return M_blockProblem.constraintTerms(); }

    /**
     * @brief All stored term descriptors.
     */
    term_container_type const& terms() const noexcept { return M_blockProblem.terms(); }

    /**
     * @brief Structured summaries of all stored HDG terms without assembly.
     */
    std::vector<TermSummary> termSummaries() const { return M_blockProblem.termSummaries(); }

    /**
     * @brief Conservative metadata consistency diagnostics.
     */
    std::vector<TermConsistencyIssue> checkTerms() const { return M_blockProblem.checkTerms(); }

    /**
     * @brief Write a side-effect-free diagnostic description of this HDG problem.
     */
    std::ostream& describe( std::ostream& os ) const { return systems::describeTerms( os, M_name, this->terms() ); }

    /**
     * @brief Return a side-effect-free diagnostic description of this HDG problem.
     */
    std::string describe() const { return systems::describeTerms( M_name, this->terms() ); }

    /**
     * @brief Side-effect-free assembly plan derived from HDG term metadata.
     */
    AssemblyPlan assemblyPlan() const
    {
        AssemblyPlanOptions options;
        options.hdgContext = true;
        return systems::buildAssemblyPlan( this->terms(), options );
    }

    /**
     * @brief Write a side-effect-free HDG assembly-plan description.
     */
    std::ostream& describeAssemblyPlan( std::ostream& os ) const
    {
        return systems::describeAssemblyPlan( os, M_name, this->assemblyPlan() );
    }

    /**
     * @brief Return a side-effect-free HDG assembly-plan description.
     */
    std::string describeAssemblyPlan() const
    {
        return systems::describeAssemblyPlan( M_name, this->assemblyPlan() );
    }

    /**
     * @brief Scheduler-ready task graph derived from HDG term metadata.
     */
    AssemblyTaskGraph assemblyTaskGraph( AssemblyExecutionOptions const& options = {} ) const
    {
        return systems::buildAssemblyTaskGraph( this->assemblyPlan(), options );
    }

    /**
     * @brief Write a side-effect-free HDG task graph description.
     */
    std::ostream& describeAssemblyTaskGraph( std::ostream& os,
                                             AssemblyExecutionOptions const& options = {} ) const
    {
        return systems::describeAssemblyTaskGraph( os, M_name, this->assemblyTaskGraph( options ) );
    }

    /**
     * @brief Return a side-effect-free HDG task graph description.
     */
    std::string describeAssemblyTaskGraph( AssemblyExecutionOptions const& options = {} ) const
    {
        return systems::describeAssemblyTaskGraph( M_name, this->assemblyTaskGraph( options ) );
    }

    /**
     * @brief Number of all stored terms.
     */
    std::size_t nTerms() const noexcept { return M_blockProblem.nTerms(); }

    /**
     * @brief Replay all stored actions into live block forms.
     */
    void assembleBlock( block_bilinear_form_type& a, block_linear_form_type& l ) const
    {
        M_blockProblem.assembleBlock( a, l );
    }

  private:
    /** @brief Semantic HDG spaces with role metadata. */
    spaces_type M_spaces;
    /** @brief Problem/backend option prefix. */
    std::string M_name;
    /** @brief Underlying block problem that stores replayable actions. */
    block_problem_type M_blockProblem;
    /** @brief Static-condensation and backend policy. */
    condensation_policy_type M_condensation;
};

/**
 * @brief Assembled HDG block system with solve-policy defaults.
 *
 * `HdgSystem` is a thin owner of the assembled `BlockSystem`. Its main job is
 * to inject the HDG condensation policy when the caller did not provide
 * explicit solve-time `_condense`, `_condenser`, or `_local` arguments.
 */
template <typename SpacesType, typename CondensationPolicyType>
class HdgSystem
{
  public:
    using problem_type = HdgProblem<SpacesType, CondensationPolicyType>;
    using block_problem_type = typename problem_type::block_problem_type;
    using block_system_type = BlockSystem<block_problem_type>;
    using value_type = typename problem_type::value_type;
    using backend_type = Backend<value_type>;
    using backend_ptrtype = std::shared_ptr<backend_type>;
    using block_bilinear_form_type = typename problem_type::block_bilinear_form_type;
    using block_linear_form_type = typename problem_type::block_linear_form_type;
    using sparse_matrix_ptrtype = typename block_system_type::sparse_matrix_ptrtype;
    using vector_ptrtype = typename block_system_type::vector_ptrtype;
    using active_system_type = typename block_system_type::active_system_type;

    /**
     * @brief Construct an HDG system from a problem descriptor.
     */
    HdgSystem( problem_type const& problem,
               backend_ptrtype backend = {},
               std::string backendName = {},
               solve::strategy strategy = solve::strategy::monolithic )
        : M_problem( problem ),
          M_blockSystem( M_problem.blockProblem(), std::move( backend ), std::move( backendName ), strategy )
    {
    }

    /**
     * @brief Assemble the underlying block system.
     */
    HdgSystem& assemble()
    {
        M_blockSystem.assemble();
        return *this;
    }

    /**
     * @brief Return true once assembly has run.
     */
    bool assembled() const noexcept { return M_blockSystem.assembled(); }

    /**
     * @brief HDG problem descriptor copied into this system.
     */
    problem_type const& problem() const noexcept { return M_problem; }

    /**
     * @brief Underlying block system for expert workflows.
     */
    block_system_type const& blockSystem() const noexcept { return M_blockSystem; }

    /**
     * @brief Mutable underlying block system for expert workflows.
     */
    block_system_type& blockSystem() noexcept { return M_blockSystem; }

    /**
     * @brief Backend used for block-form allocation.
     */
    backend_ptrtype const& backend() const noexcept { return M_blockSystem.backend(); }

    /**
     * @brief Base backend/solver prefix.
     */
    std::string const& backendName() const noexcept { return M_blockSystem.backendName(); }

    /**
     * @brief Final backend name used by the condensed trace/global solve.
     */
    std::string condensedBackendName() const
    {
        return M_problem.condensationPolicy().backendName( this->backendName() );
    }

    /**
     * @brief Default backend name used by HDG postprocessing solves.
     */
    std::string postBackendName() const
    {
        return M_problem.condensationPolicy().postBackendName( this->backendName() );
    }

    /**
     * @brief Block assembly strategy.
     */
    solve::strategy strategy() const noexcept { return M_blockSystem.strategy(); }

    /**
     * @brief Block bilinear form.
     */
    block_bilinear_form_type const& matrixForm() const noexcept { return M_blockSystem.matrixForm(); }

    /**
     * @brief Block bilinear form.
     */
    block_bilinear_form_type& matrixForm() noexcept { return M_blockSystem.matrixForm(); }

    /**
     * @brief Block linear form.
     */
    block_linear_form_type const& rhsForm() const noexcept { return M_blockSystem.rhsForm(); }

    /**
     * @brief Block linear form.
     */
    block_linear_form_type& rhsForm() noexcept { return M_blockSystem.rhsForm(); }

    /**
     * @brief Base monolithic sparse matrix view when available.
     */
    sparse_matrix_ptrtype baseMatrixPtr() const { return M_blockSystem.baseMatrixPtr(); }

    /**
     * @brief Active matrix/vector pair used by monolithic solves.
     */
    active_system_type activeSystem() const noexcept { return M_blockSystem.activeSystem(); }

    /**
     * @brief Active matrix pointer.
     */
    sparse_matrix_ptrtype activeMatrixPtr() const noexcept { return M_blockSystem.activeMatrixPtr(); }

    /**
     * @brief Active rhs vector pointer.
     */
    vector_ptrtype activeVectorPtr() const noexcept { return M_blockSystem.activeVectorPtr(); }

    /**
     * @brief Return true if constraints materialized a distinct active pair.
     */
    bool constrained() const noexcept { return M_blockSystem.constrained(); }

    /**
     * @brief Solve the assembled HDG system.
     */
    template <typename... Args>
    solve_result solve( Args&&... args )
    {
        if constexpr ( detail::hasNamedArgument<na::condense, Args...> ||
                       detail::hasNamedArgument<na::condenser, Args...> ||
                       detail::hasNamedArgument<na::local, Args...> )
        {
            return M_blockSystem.solve( std::forward<Args>( args )... );
        }
        else if ( M_problem.condensationPolicy().enabled() )
        {
            if constexpr ( detail::hasNamedArgument<na::name, Args...> )
            {
                return M_blockSystem.solve( _condense = true,
                                            _condenser = M_problem.condensationPolicy().condenser(),
                                            std::forward<Args>( args )... );
            }
            else
            {
                auto const solveBackendBaseName =
                    M_problem.condensationPolicy().solveBackendBaseName( this->backendName() );
                if ( solveBackendBaseName == this->backendName() )
                {
                    return M_blockSystem.solve( _condense = true,
                                                _condenser = M_problem.condensationPolicy().condenser(),
                                                std::forward<Args>( args )... );
                }
                return M_blockSystem.solve( _name = solveBackendBaseName,
                                            _condense = true,
                                            _condenser = M_problem.condensationPolicy().condenser(),
                                            std::forward<Args>( args )... );
            }
        }
        else
        {
            return M_blockSystem.solve( std::forward<Args>( args )... );
        }
    }

  private:
    /** @brief HDG problem descriptor copied by value for solve-policy access. */
    problem_type M_problem;
    /** @brief Underlying assembled block system. */
    block_system_type M_blockSystem;
};

} // namespace Feel::systems

namespace Feel
{

/**
 * @brief Build an HDG-aware stationary linear block problem descriptor.
 *
 * Prefer `_spaces = H` for `HdgSpaces`; `_space = H` is accepted as a
 * transitional alias for consistency with `blockProblem`.
 */
template <typename... Ts>
auto hdgProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    std::string name = args.get_else( _name, "hdg_problem" );

    if constexpr ( systems::detail::hasNamedArgument<na::spaces, Ts...> )
    {
        auto spaces = args.get( _spaces );
        return systems::HdgProblem<std::remove_cvref_t<decltype( spaces )>>( std::move( spaces ), std::move( name ) );
    }
    else
    {
        auto spaces = args.get( _space );
        return systems::HdgProblem<std::remove_cvref_t<decltype( spaces )>>( std::move( spaces ), std::move( name ) );
    }
}

} // namespace Feel

#endif // FEELPP_FEELSYSTEM_HDG_HPP
