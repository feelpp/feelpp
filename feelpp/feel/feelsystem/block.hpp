/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_BLOCK_HPP
#define FEELPP_FEELSYSTEM_BLOCK_HPP 1

#include <algorithm>
#include <concepts>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/hana/integral_constant.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/productspaceconcepts.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelsystem/auxiliary.hpp>
#include <feel/feelsystem/backend.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/execution.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/result.hpp>
#include <feel/feelsystem/terms.hpp>
#include <feel/feelvf/blockforms.hpp>

namespace Feel::systems
{

namespace detail
{
/**
 * @brief True for Feel++ DSEL expression wrappers.
 */
template <typename T>
concept ExpressionLike =
    requires {
        typename std::remove_cvref_t<T>::expression_type;
    };

/**
 * @brief True for user callbacks that mutate a form and return nothing.
 */
template <typename Action, typename Form>
concept VoidActionFor =
    ( !ExpressionLike<Action> ) &&
    requires( Action& action, Form& form ) {
        { action( form ) } -> std::same_as<void>;
    };

/**
 * @brief Return a representative communicator for a product space.
 *
 * Static product spaces do not expose `worldCommPtr()` directly today. Their
 * first component space does, which is enough for the default backend factory.
 */
template <typename ProductSpaceType>
worldcomm_ptr_t
productWorldCommPtr( ProductSpaceType const& spaces )
{
    if constexpr ( requires { spaces.worldCommPtr(); } )
    {
        return spaces.worldCommPtr();
    }
    else if constexpr ( requires { spaces[boost::hana::int_c<0>]->worldCommPtr(); } )
    {
        return spaces[boost::hana::int_c<0>]->worldCommPtr();
    }
    else if constexpr ( requires { spaces[0]->worldCommPtr(); } )
    {
        return spaces[0]->worldCommPtr();
    }
    else
    {
        return Environment::worldCommPtr();
    }
}

/**
 * @brief Return an optional integer block index for static or runtime selectors.
 */
template <typename T>
std::optional<int>
blockIndex( T const& value )
{
    using value_type = std::remove_cvref_t<T>;
    if constexpr ( requires { value_type::value; } )
        return static_cast<int>( value_type::value );
    else if constexpr ( requires { static_cast<int>( value ); } )
        return static_cast<int>( value );
    else
        return std::nullopt;
}
} // namespace detail

/**
 * @brief Mathematical descriptor for a linear block or mixed variational system.
 *
 * `BlockProblem` stores replayable callbacks over Feel++ `blockform2` and
 * `blockform1`. The block forms remain the source of truth for monolithic,
 * local, and static-condensation assembly; this class only gives them the same
 * high-level problem/system lifecycle as scalar `Problem`.
 *
 * @tparam ProductSpaceType Static Feel++ product space type returned by
 * `product(Vh, Qh, ...)`.
 * @tparam Tags Problem tags such as `linear_t`, `stationary_t`, `block_t`, or
 * `mixed_t`.
 */
template <typename ProductSpaceType, typename... Tags>
    requires StaticProductSpacesType<ProductSpaceType>
class BlockProblem
{
  public:
    using product_space_type = std::remove_cvref_t<ProductSpaceType>;
    using space_type = product_space_type;
    using value_type = typename product_space_type::value_type;
    using element_type = typename product_space_type::element_type;
    using backend_type = Backend<value_type>;
    using backend_ptrtype = std::shared_ptr<backend_type>;
    using block_bilinear_form_type = decltype( blockform2( std::declval<product_space_type>() ) );
    using block_linear_form_type = decltype( blockform1( std::declval<product_space_type>() ) );
    using block_matrix_action_type = std::function<void( AssemblyContext const& )>;
    using block_vector_action_type = std::function<void( AssemblyContext const& )>;
    using block_constraint_action_type = std::function<void( AssemblyContext const& )>;
    using term_container_type = std::vector<term_descriptor>;

    template <typename Tag>
    static constexpr bool has_tag = ( std::same_as<Tag, Tags> || ... || false );

    /**
     * @brief Lightweight view counting terms by role.
     */
    class role_terms_view
    {
      public:
        role_terms_view( BlockProblem const& problem, term_role role )
            : M_problem( &problem ), M_role( role )
        {
        }

        /**
         * @brief Role selected by this view.
         */
        term_role role() const noexcept { return M_role; }

        /**
         * @brief Number of terms with this role.
         */
        std::size_t size() const noexcept { return M_problem->nTerms( M_role ); }

        /**
         * @brief Return true when the view has no terms.
         */
        bool empty() const noexcept { return this->size() == 0; }

        /**
         * @brief Descriptor at `index` within this role-filtered view.
         */
        term_descriptor const& operator[]( std::size_t index ) const { return M_problem->term( M_role, index ); }

        /**
         * @brief Metadata at `index` within this role-filtered view.
         */
        TermMetadata const& metadata( std::size_t index ) const { return ( *this )[index].metadata(); }

      private:
        BlockProblem const* M_problem = nullptr;
        term_role M_role = term_role::lhs;
    };

    /**
     * @brief Construct a block problem over a product space.
     */
    explicit BlockProblem( product_space_type spaces, std::string name = {} )
        : M_spaces( std::move( spaces ) ),
          M_name( std::move( name ) )
    {
    }

    /**
     * @brief Product space used by the block forms.
     */
    product_space_type const& space() const noexcept { return M_spaces; }

    /**
     * @brief Mutable product space, used when callers need product elements.
     */
    product_space_type& space() noexcept { return M_spaces; }

    /**
     * @brief Problem/backend option prefix.
     */
    std::string const& name() const noexcept { return M_name; }

    /**
     * @brief Default block assembly strategy.
     */
    solve::strategy strategy() const noexcept { return M_strategy; }

    /**
     * @brief Set the default block assembly strategy.
     */
    BlockProblem& strategy( solve::strategy value ) noexcept
    {
        M_strategy = value;
        return *this;
    }

    /**
     * @brief Assemble this descriptor into a `BlockSystem`.
     */
    template <typename... Args>
    auto assemble( Args&&... args ) const;

    /**
     * @brief Assemble and solve this descriptor in one call.
     */
    template <typename... Args>
    auto solve( Args&&... args ) const;

    /**
     * @brief Add a whole-block matrix assembly action.
     */
    template <typename Action>
    BlockProblem& lhs( Action&& action )
    {
        this->addDescriptor( term_role::lhs );
        this->addMatrixAction( std::forward<Action>( action ) );
        return *this;
    }

    /** @brief Declare a floating local operator by replaying all lhs blocks. */
    BlockProblem& localNeumannFromLhs( std::string name = "geneo" )
    {
        this->addAuxiliaryDescriptor( AuxiliaryOperatorDescriptor{
            std::move( name ),
            AuxiliaryOperatorKind::local_neumann,
            AuxiliaryOperatorSource::lhs_terms } );
        return *this;
    }

    /** @brief Add an explicit whole-block contribution to a local operator. */
    template <typename Action>
    BlockProblem& localNeumann( std::string name, Action&& action )
    {
        auto& entry = this->addAuxiliaryDescriptor( AuxiliaryOperatorDescriptor{
            std::move( name ),
            AuxiliaryOperatorKind::local_neumann,
            AuxiliaryOperatorSource::explicit_terms } );
        entry.actions.emplace_back( this->makeMatrixAction( std::forward<Action>( action ) ) );
        return *this;
    }

    /** @brief Declared auxiliary operators in insertion order. */
    std::vector<AuxiliaryOperatorDescriptor> auxiliaryOperatorDescriptors() const
    {
        std::vector<AuxiliaryOperatorDescriptor> descriptors;
        descriptors.reserve( M_auxiliaryOperators.size() );
        for ( auto const& entry : M_auxiliaryOperators )
            descriptors.push_back( entry.descriptor );
        return descriptors;
    }

    bool hasAuxiliaryOperators() const noexcept { return !M_auxiliaryOperators.empty(); }

    /**
     * @brief Add a whole-block vector assembly action.
     */
    template <typename Action>
    BlockProblem& rhs( Action&& action )
    {
        this->addDescriptor( term_role::rhs );
        this->addVectorAction( std::forward<Action>( action ) );
        return *this;
    }

    /**
     * @brief Add a constraint action over live block forms.
     */
    template <typename Constraint>
    BlockProblem& constraint( Constraint&& constraint )
    {
        this->addDescriptor( term_role::constraint );
        this->addConstraintAction( std::forward<Constraint>( constraint ) );
        return *this;
    }

    /**
     * @brief Alias for `.constraint(...)`.
     */
    template <typename Constraint>
    BlockProblem& bc( Constraint&& constraint )
    {
        return this->constraint( std::forward<Constraint>( constraint ) );
    }

    /**
     * @brief Add a contribution to one matrix block.
     *
     * The row and column selectors are the usual blockform selectors, e.g.
     * `0_c`, `1_c`, or runtime integers.
     */
    template <typename Row, typename Column, typename Term>
    BlockProblem& block( Row row, Column column, Term&& term )
    {
        auto storedTerm = std::make_shared<std::decay_t<Term>>( std::forward<Term>( term ) );
        auto metadata = defaultTermMetadata( term_role::lhs );
        metadata.kind = TermKind::Block;
        metadata.name = "block";
        metadata.row = detail::blockIndex( row );
        metadata.column = detail::blockIndex( column );
        this->addDescriptor( term_role::lhs, std::move( metadata ) );
        M_matrixActions.emplace_back(
            [row, column, storedTerm]( AssemblyContext const& context )
            {
                auto termContext = context;
                termContext.row = detail::blockIndex( row );
                termContext.column = detail::blockIndex( column );
                auto& a = termContext.blockBilinearFormAs<block_bilinear_form_type>();
                using term_type = std::decay_t<Term>;
                if constexpr ( detail::VoidActionFor<term_type, block_bilinear_form_type> )
                    ( *storedTerm )( a );
                else
                    a( row, column ) += *storedTerm;
            } );
        return *this;
    }

    /**
     * @brief Add a contribution to one rhs block.
     */
    template <typename Row, typename Term>
    BlockProblem& row( Row row, Term&& term )
    {
        auto storedTerm = std::make_shared<std::decay_t<Term>>( std::forward<Term>( term ) );
        auto metadata = defaultTermMetadata( term_role::rhs );
        metadata.kind = TermKind::Row;
        metadata.name = "row";
        metadata.row = detail::blockIndex( row );
        this->addDescriptor( term_role::rhs, std::move( metadata ) );
        M_vectorActions.emplace_back(
            [row, storedTerm]( AssemblyContext const& context )
            {
                auto termContext = context;
                termContext.row = detail::blockIndex( row );
                auto& l = termContext.blockLinearFormAs<block_linear_form_type>();
                using term_type = std::decay_t<Term>;
                if constexpr ( detail::VoidActionFor<term_type, block_linear_form_type> )
                    ( *storedTerm )( l );
                else
                    l( row ) += *storedTerm;
            } );
        return *this;
    }

    /**
     * @brief Matrix term metadata.
     */
    role_terms_view lhs() const noexcept { return this->terms( term_role::lhs ); }

    /**
     * @brief Matrix term metadata.
     */
    role_terms_view lhsTerms() const noexcept { return this->lhs(); }

    /**
     * @brief Matrix block term metadata.
     */
    role_terms_view blockTerms() const noexcept { return this->lhs(); }

    /**
     * @brief Rhs term metadata.
     */
    role_terms_view rhs() const noexcept { return this->terms( term_role::rhs ); }

    /**
     * @brief Rhs term metadata.
     */
    role_terms_view rhsTerms() const noexcept { return this->rhs(); }

    /**
     * @brief Rhs row term metadata.
     */
    role_terms_view rowTerms() const noexcept { return this->rhs(); }

    /**
     * @brief Constraint term metadata.
     */
    role_terms_view constraints() const noexcept { return this->terms( term_role::constraint ); }

    /**
     * @brief Constraint term metadata.
     */
    role_terms_view constraintTerms() const noexcept { return this->constraints(); }

    /**
     * @brief All stored term descriptors.
     */
    term_container_type const& terms() const noexcept { return M_terms; }

    /**
     * @brief Structured summaries of all stored block terms without assembly.
     */
    std::vector<TermSummary> termSummaries() const { return systems::termSummaries( M_terms ); }

    /**
     * @brief Conservative metadata consistency diagnostics.
     */
    std::vector<TermConsistencyIssue> checkTerms() const { return systems::checkTermConsistency( M_terms ); }

    /**
     * @brief Write a side-effect-free diagnostic description of this block problem.
     */
    std::ostream& describe( std::ostream& os ) const { return systems::describeTerms( os, M_name, M_terms ); }

    /**
     * @brief Return a side-effect-free diagnostic description of this block problem.
     */
    std::string describe() const { return systems::describeTerms( M_name, M_terms ); }

    /**
     * @brief Side-effect-free assembly plan derived from block term metadata.
     */
    AssemblyPlan assemblyPlan() const { return systems::buildAssemblyPlan( M_terms ); }

    /**
     * @brief Write a side-effect-free assembly-plan description.
     */
    std::ostream& describeAssemblyPlan( std::ostream& os ) const
    {
        return systems::describeAssemblyPlan( os, M_name, this->assemblyPlan() );
    }

    /**
     * @brief Return a side-effect-free assembly-plan description.
     */
    std::string describeAssemblyPlan() const
    {
        return systems::describeAssemblyPlan( M_name, this->assemblyPlan() );
    }

    /**
     * @brief Scheduler-ready task graph derived from block term metadata.
     */
    AssemblyTaskGraph assemblyTaskGraph( AssemblyExecutionOptions const& options = {} ) const
    {
        return systems::buildAssemblyTaskGraph( this->assemblyPlan(), options );
    }

    /**
     * @brief Write a side-effect-free scheduler-ready task graph description.
     */
    std::ostream& describeAssemblyTaskGraph( std::ostream& os,
                                             AssemblyExecutionOptions const& options = {} ) const
    {
        return systems::describeAssemblyTaskGraph( os, M_name, this->assemblyTaskGraph( options ) );
    }

    /**
     * @brief Return a side-effect-free scheduler-ready task graph description.
     */
    std::string describeAssemblyTaskGraph( AssemblyExecutionOptions const& options = {} ) const
    {
        return systems::describeAssemblyTaskGraph( M_name, this->assemblyTaskGraph( options ) );
    }

    /**
     * @brief Terms matching a role.
     */
    role_terms_view terms( term_role role ) const noexcept { return role_terms_view{ *this, role }; }

    /**
     * @brief Number of all stored terms.
     */
    std::size_t nTerms() const noexcept { return M_terms.size(); }

    /**
     * @brief Number of stored terms matching one role.
     */
    std::size_t nTerms( term_role role ) const noexcept
    {
        return static_cast<std::size_t>(
            std::count_if( M_terms.begin(), M_terms.end(),
                           [role]( term_descriptor const& term )
                           { return term.role == role; } ) );
    }

    /**
     * @brief Return the `index`-th stored term matching `role`.
     */
    term_descriptor const& term( term_role role, std::size_t index ) const
    {
        std::size_t current = 0;
        for ( auto const& term : M_terms )
        {
            if ( term.role != role )
                continue;
            if ( current == index )
                return term;
            ++current;
        }
        throw std::out_of_range( "Feel::systems::BlockProblem term index out of range" );
    }

    /**
     * @brief Replay all stored actions into live block forms.
     */
    void assembleBlock( block_bilinear_form_type& a, block_linear_form_type& l ) const
    {
        auto context = this->makeAssemblyContext( AssemblyMode::Mixed );
        context.blockBilinearForm = &a;
        context.blockLinearForm = &l;
        for ( auto const& action : M_matrixActions )
            action( context );
        for ( auto const& action : M_vectorActions )
            action( context );
        for ( auto const& action : M_constraintActions )
        {
            auto constraintContext = context;
            constraintContext.mode = AssemblyMode::Constraint;
            constraintContext.role = term_role::constraint;
            action( constraintContext );
        }
    }

    /** @brief Assemble one named auxiliary block operator without constraints. */
    void assembleAuxiliaryOperator( std::string const& name, block_bilinear_form_type& a ) const
    {
        auto const it = std::find_if( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                                      [&name]( auto const& entry )
                                      { return entry.descriptor.name == name; } );
        if ( it == M_auxiliaryOperators.end() )
            throw std::out_of_range( "Feel::systems::BlockProblem auxiliary operator not found: " + name );

        auto context = this->makeAssemblyContext( AssemblyMode::MatrixTerm );
        context.role = term_role::lhs;
        context.blockBilinearForm = &a;
        if ( it->descriptor.source == AuxiliaryOperatorSource::lhs_terms )
        {
            for ( auto const& action : M_matrixActions )
                action( context );
            return;
        }
        for ( auto const& action : it->actions )
            action( context );
    }

  private:
    AssemblyContext makeAssemblyContext( AssemblyMode mode ) const
    {
        AssemblyContext context;
        context.mode = mode;
        context.problem = this;
        context.space = &M_spaces;
        return context;
    }

    void addDescriptor( term_role role )
    {
        this->addDescriptor( role, defaultTermMetadata( role ) );
    }

    void addDescriptor( term_role role, TermMetadata metadata )
    {
        M_terms.push_back( term_descriptor{ role, roleName( role ), std::move( metadata ) } );
    }

    template <typename Action>
    block_matrix_action_type makeMatrixAction( Action&& action )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        return [storedAction]( AssemblyContext const& context )
            {
                auto& a = context.blockBilinearFormAs<block_bilinear_form_type>();
                using action_type = std::decay_t<Action>;
                if constexpr ( detail::VoidActionFor<action_type, block_bilinear_form_type> )
                    ( *storedAction )( a );
                else
                    a += *storedAction;
            };
    }

    template <typename Action>
    void addMatrixAction( Action&& action )
    {
        M_matrixActions.emplace_back( this->makeMatrixAction( std::forward<Action>( action ) ) );
    }

    struct auxiliary_operator_entry
    {
        AuxiliaryOperatorDescriptor descriptor;
        std::vector<block_matrix_action_type> actions;
    };

    auxiliary_operator_entry& addAuxiliaryDescriptor( AuxiliaryOperatorDescriptor descriptor )
    {
        auto const it = std::find_if( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                                      [&descriptor]( auto const& entry )
                                      { return entry.descriptor.name == descriptor.name; } );
        if ( it == M_auxiliaryOperators.end() )
        {
            M_auxiliaryOperators.push_back( auxiliary_operator_entry{ std::move( descriptor ), {} } );
            return M_auxiliaryOperators.back();
        }
        if ( it->descriptor.kind != descriptor.kind || it->descriptor.source != descriptor.source )
            throw std::invalid_argument( "Feel::systems::BlockProblem auxiliary operator declaration conflicts with existing operator: " + descriptor.name );
        return *it;
    }

    template <typename Action>
    void addVectorAction( Action&& action )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        M_vectorActions.emplace_back(
            [storedAction]( AssemblyContext const& context )
            {
                auto& l = context.blockLinearFormAs<block_linear_form_type>();
                using action_type = std::decay_t<Action>;
                if constexpr ( detail::VoidActionFor<action_type, block_linear_form_type> )
                    ( *storedAction )( l );
                else
                    l += *storedAction;
            } );
    }

    template <typename Constraint>
    void addConstraintAction( Constraint&& constraint )
    {
        auto storedConstraint = std::make_shared<std::decay_t<Constraint>>( std::forward<Constraint>( constraint ) );
        M_constraintActions.emplace_back(
            [storedConstraint]( AssemblyContext const& context )
            {
                auto& a = context.blockBilinearFormAs<block_bilinear_form_type>();
                auto& l = context.blockLinearFormAs<block_linear_form_type>();
                ( *storedConstraint )( a, l );
            } );
    }

    /** @brief Product space used by the block forms. */
    product_space_type M_spaces;
    /** @brief Problem/backend option prefix. */
    std::string M_name;
    /** @brief Default block-form assembly/solve strategy. */
    solve::strategy M_strategy = solve::strategy::monolithic;
    /** @brief Stored block term descriptors in insertion order. */
    term_container_type M_terms;
    /** @brief Stored block matrix assembly callbacks. */
    std::vector<block_matrix_action_type> M_matrixActions;
    /** @brief Named auxiliary block operators and explicit actions. */
    std::vector<auxiliary_operator_entry> M_auxiliaryOperators;
    /** @brief Stored block vector assembly callbacks. */
    std::vector<block_vector_action_type> M_vectorActions;
    /** @brief Stored block constraint callbacks. */
    std::vector<block_constraint_action_type> M_constraintActions;
};

/**
 * @brief Assembled block system produced from a `BlockProblem`.
 */
template <BlockProblemLike ProblemType>
class BlockSystem
{
  public:
    using problem_type = std::remove_cvref_t<ProblemType>;
    using value_type = typename problem_type::value_type;
    using backend_type = Backend<value_type>;
    using backend_ptrtype = std::shared_ptr<backend_type>;
    using block_bilinear_form_type = typename problem_type::block_bilinear_form_type;
    using block_linear_form_type = typename problem_type::block_linear_form_type;
    using sparse_matrix_ptrtype = typename block_bilinear_form_type::sparse_matrix_ptrtype;
    using vector_ptrtype = typename backend_type::vector_ptrtype;
    using active_system_type = std::pair<sparse_matrix_ptrtype, vector_ptrtype>;
    using auxiliary_operator_type = AssembledAuxiliaryOperator<value_type>;

    /**
     * @brief Construct a block system from a problem descriptor.
     */
    BlockSystem( problem_type const& problem,
                 backend_ptrtype backend = {},
                 std::string backendName = {},
                 solve::strategy strategy = solve::strategy::monolithic )
        : M_problem( problem ),
          M_backend( std::move( backend ) ),
          M_backendName( detail::defaultBackendName( std::move( backendName ), M_problem.name(), "block_system" ) ),
          M_strategy( strategy )
    {
        if ( M_backend )
            M_backendName = M_backend->prefix();
    }

    /**
     * @brief Assemble all block matrix, vector, and constraint actions.
     */
    BlockSystem& assemble()
    {
        if ( !M_backend )
            M_backend = detail::makeBackend<value_type>( M_backendName, detail::productWorldCommPtr( M_problem.space() ) );

        auto matrixSpaces = M_problem.space();
        auto rhsSpaces = M_problem.space();
        M_matrix.emplace( blockform2( std::move( matrixSpaces ), M_strategy, M_backend ) );
        M_rhs.emplace( blockform1( std::move( rhsSpaces ), M_strategy, M_backend ) );
        auto const distributedMap = M_matrix->baseMatrixPtr()->mapRowPtr();
        M_problem.assembleBlock( *M_matrix, *M_rhs );
        M_matrix->close();
        M_rhs->close();
        M_activeSystem = this->baseActiveSystem();
        if ( M_strategy == solve::strategy::monolithic )
        {
            M_activeSystem = detail::completeActiveSystem( M_matrix->activeSystem( *M_rhs ),
                                                           M_activeSystem.first,
                                                           M_activeSystem.second );
            if ( M_activeSystem.first )
                M_activeSystem.first->close();
            if ( M_activeSystem.second )
                M_activeSystem.second->close();
        }
        this->assembleAuxiliaryOperators( distributedMap );
        M_assembled = true;
        return *this;
    }

    /**
     * @brief Return true once assembly has run.
     */
    bool assembled() const noexcept { return M_assembled; }

    /**
     * @brief Problem descriptor copied into this system.
     */
    problem_type const& problem() const noexcept { return M_problem; }

    /**
     * @brief Backend used for block-form allocation.
     */
    backend_ptrtype const& backend() const noexcept { return M_backend; }

    /**
     * @brief Backend/solver prefix.
     */
    std::string const& backendName() const noexcept { return M_backendName; }

    /**
     * @brief Block assembly strategy.
     */
    solve::strategy strategy() const noexcept { return M_strategy; }

    /**
     * @brief Block bilinear form.
     */
    block_bilinear_form_type const& matrixForm() const noexcept { return *M_matrix; }

    /**
     * @brief Block bilinear form.
     */
    block_bilinear_form_type& matrixForm() noexcept { return *M_matrix; }

    /**
     * @brief Block linear form.
     */
    block_linear_form_type const& rhsForm() const noexcept { return *M_rhs; }

    /**
     * @brief Block linear form.
     */
    block_linear_form_type& rhsForm() noexcept { return *M_rhs; }

    /**
     * @brief Base monolithic sparse matrix view when available.
     */
    sparse_matrix_ptrtype baseMatrixPtr() const { return M_matrix ? M_matrix->baseMatrixPtr() : sparse_matrix_ptrtype{}; }

    /**
     * @brief Active matrix/vector pair used by monolithic solves.
     */
    active_system_type activeSystem() const noexcept { return M_activeSystem; }

    /**
     * @brief Active matrix pointer.
     */
    sparse_matrix_ptrtype activeMatrixPtr() const noexcept { return M_activeSystem.first; }

    /**
     * @brief Active rhs vector pointer.
     */
    vector_ptrtype activeVectorPtr() const noexcept { return M_activeSystem.second; }

    /**
     * @brief Return true when the block form contains Dirichlet constraints.
     */
    bool constrained() const noexcept
    {
        return M_matrix && M_matrix->hasDirichletConstraints();
    }

    std::vector<auxiliary_operator_type> const& auxiliaryOperators() const noexcept
    {
        return M_auxiliaryOperators;
    }

    bool hasAuxiliaryOperator( std::string const& name ) const noexcept
    {
        return std::any_of( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                            [&name]( auto const& op ) { return op.name() == name; } );
    }

    auxiliary_operator_type const& auxiliaryOperator( std::string const& name ) const
    {
        auto const it = std::find_if( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                                      [&name]( auto const& op ) { return op.name() == name; } );
        if ( it == M_auxiliaryOperators.end() )
            throw std::out_of_range( "Feel::systems::BlockSystem auxiliary operator not found: " + name );
        return *it;
    }

    /**
     * @brief Solve the assembled block system.
     */
    template <typename... Args>
    solve_result solve( Args&&... args )
    {
        if ( !M_assembled )
            this->assemble();

        auto rawResult = this->solveBlockForm( std::forward<Args>( args )... );
        return solve_result{
            rawResult.isConverged(),
            static_cast<int>( rawResult.nIterations() ),
            static_cast<double>( rawResult.residual() ),
            rawResult.isConverged()
                ? std::string{ "feelsystem BlockSystem solve converged" }
                : std::string{ "feelsystem BlockSystem solve did not converge" } };
    }

  private:
    void assembleAuxiliaryOperators( std::shared_ptr<DataMap<>> const& distributedMap )
    {
        M_auxiliaryOperators.clear();
        for ( auto const& descriptor : M_problem.auxiliaryOperatorDescriptors() )
        {
            auto layout = detail::makeAuxiliaryDataMap( distributedMap );
            auto localBackend = detail::makeBackendLike<value_type>(
                M_backend,
                M_backendName + ".auxiliary." + descriptor.name,
                Environment::worldCommSeqPtr() );
            auto localMatrix = localBackend->newMatrix( layout.localMap, layout.localMap,
                                                        NON_HERMITIAN, false );
            auto const n = layout.localMap->nDof();
            localMatrix->init( n, n, n, n, 64, 0 );
            localMatrix->zero();

            using condensed_matrix_type = typename block_bilinear_form_type::condensed_matrix_type;
            auto condensed = std::make_shared<condensed_matrix_type>( solve::strategy::monolithic,
                                                                      localMatrix,
                                                                      localBackend );
            auto localSpaces = M_problem.space();
            block_bilinear_form_type localForm( std::move( localSpaces ), condensed );
            M_problem.assembleAuxiliaryOperator( descriptor.name, localForm );
            localForm.close();
            M_auxiliaryOperators.emplace_back( descriptor,
                                               std::move( localBackend ),
                                               std::move( localMatrix ),
                                               std::move( layout.globalDofIds ) );
        }
    }

    template <typename... Args>
    auto solveBlockForm( Args&&... args )
    {
        if ( !M_auxiliaryOperators.empty() )
            return this->solveMonolithicWithAuxiliaryOperators( std::forward<Args>( args )... );

        if constexpr ( detail::hasNamedArgument<na::name, Args...> )
        {
            return M_matrix->solve( _rhs = *M_rhs,
                                    std::forward<Args>( args )... );
        }
        else
        {
            return M_matrix->solve( _rhs = *M_rhs,
                                    _name = M_backendName,
                                    std::forward<Args>( args )... );
        }
    }

    template <typename... Args>
    auto solveMonolithicWithAuxiliaryOperators( Args&&... args )
    {
        CHECK( M_strategy == solve::strategy::monolithic )
            << "auxiliary preconditioner operators currently require a monolithic block solve";
        if constexpr ( detail::hasNamedArgument<na::condense, Args...> )
            CHECK( !detail::namedArgumentValue<na::condense>( args... ) )
                << "auxiliary preconditioner operators cannot be used with solve-time static condensation";
        if constexpr ( detail::hasNamedArgument<na::local, Args...> )
            CHECK( !detail::namedArgumentValue<na::local>( args... ) )
                << "auxiliary preconditioner operators cannot be used with a local-only block solve";

        if constexpr ( detail::hasNamedArgument<na::prec, Args...> )
        {
            auto preconditioner = detail::namedArgumentValue<na::prec>( args... );
            detail::attachAuxiliaryOperators( preconditioner, M_auxiliaryOperators );
            return M_matrix->solve( _rhs = *M_rhs,
                                    std::forward<Args>( args )... );
        }
        else
        {
            std::string pc = M_backend->pcType();
            std::string factorPackage = M_backend->pcFactorMatSolverPackageType();
            if constexpr ( detail::hasNamedArgument<na::pc, Args...> )
                pc = detail::namedArgumentValue<na::pc>( args... );
            if constexpr ( detail::hasNamedArgument<na::pcfactormatsolverpackage, Args...> )
                factorPackage = detail::namedArgumentValue<na::pcfactormatsolverpackage>( args... );

            auto preconditioner = detail::makeAuxiliaryPreconditioner(
                M_backend, M_activeSystem.first, pc, factorPackage );
            detail::attachAuxiliaryOperators( preconditioner, M_auxiliaryOperators );
            if constexpr ( detail::hasNamedArgument<na::name, Args...> )
            {
                return M_matrix->solve( _rhs = *M_rhs,
                                        _prec = preconditioner,
                                        std::forward<Args>( args )... );
            }
            else
            {
                return M_matrix->solve( _rhs = *M_rhs,
                                        _name = M_backendName,
                                        _prec = preconditioner,
                                        std::forward<Args>( args )... );
            }
        }
    }

    active_system_type baseActiveSystem() const
    {
        active_system_type active;
        if ( M_matrix )
            active.first = M_matrix->baseMatrixPtr();
        if ( M_rhs && M_rhs->vectorPtr() )
            active.second = M_rhs->vectorPtr()->getVector();
        return active;
    }

    /** @brief Problem descriptor copied by value for replayable assembly. */
    problem_type M_problem;
    /** @brief Backend used for block-form allocation and solve. */
    backend_ptrtype M_backend;
    /** @brief Backend/solver option prefix. */
    std::string M_backendName;
    /** @brief Block assembly and solve strategy. */
    solve::strategy M_strategy = solve::strategy::monolithic;
    /** @brief Owned block bilinear form, created during assembly. */
    std::optional<block_bilinear_form_type> M_matrix;
    /** @brief Owned block linear form, created during assembly. */
    std::optional<block_linear_form_type> M_rhs;
    /** @brief Active sparse matrix/vector pair for monolithic solves. */
    active_system_type M_activeSystem;
    /** @brief True once assembly has completed. */
    bool M_assembled = false;
    /** @brief Sequential per-rank auxiliary matrices and overlap mappings. */
    std::vector<auxiliary_operator_type> M_auxiliaryOperators;
};

} // namespace Feel::systems

namespace Feel
{

/**
 * @brief Build a stationary block linear problem descriptor.
 */
template <typename... Ts>
auto blockProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto spaces = args.get( _space );
    std::string name = args.get_else( _name, "block_problem" );
    return systems::BlockProblem<std::remove_cvref_t<decltype( spaces )>,
                                 systems::linear_t,
                                 systems::stationary_t,
                                 systems::block_t>( std::move( spaces ), std::move( name ) );
}

/**
 * @brief Build a stationary mixed linear problem descriptor.
 */
template <typename... Ts>
auto mixedProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto spaces = args.get( _space );
    std::string name = args.get_else( _name, "mixed_problem" );
    return systems::BlockProblem<std::remove_cvref_t<decltype( spaces )>,
                                 systems::linear_t,
                                 systems::stationary_t,
                                 systems::block_t,
                                 systems::mixed_t>( std::move( spaces ), std::move( name ) );
}

} // namespace Feel

#endif // FEELPP_FEELSYSTEM_BLOCK_HPP
