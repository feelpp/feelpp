/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_PROBLEM_HPP
#define FEELPP_FEELSYSTEM_PROBLEM_HPP 1

/**
 * @file problem.hpp
 * @brief Fluent scalar problem descriptors for the modern FeelSystem API.
 */

#include <algorithm>
#include <concepts>
#include <functional>
#include <map>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <feel/feelalg/backend.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelsystem/auxiliary.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/execution.hpp>
#include <feel/feelsystem/strategy.hpp>
#include <feel/feelsystem/terms.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
template <typename SpaceType>
class Bdf;
template <typename SpaceType>
class Newmark;
} // namespace Feel

namespace Feel::systems
{

namespace detail
{
/**
 * @brief Forward parameter values to expression-like objects when supported.
 */
template <typename T>
void setParameterValuesIfSupported( T& object, std::map<std::string, double> const& values )
{
    if constexpr ( requires { object.setParameterValues( values ); } )
        object.setParameterValues( values );
}

/**
 * @brief Build a test expression from a first-class test basis or field object.
 */
template <typename TestType>
auto testExpression( TestType const& test )
{
    if constexpr ( Feel::vf::is_test_basis_v<TestType> )
        return Feel::vf::id( test );
    else
        return Feel::vf::id( test );
}

/**
 * @brief Build a trial expression from a first-class trial basis or field object.
 */
template <typename TrialType>
auto trialExpression( TrialType const& trial )
{
    if constexpr ( Feel::vf::is_trial_basis_v<TrialType> )
        return Feel::vf::idt( trial );
    else
        return Feel::vf::idt( trial );
}
} // namespace detail

/**
 * @brief Replayable scalar mathematical problem descriptor.
 *
 * `Problem` stores named variational terms, boundary/constraint callbacks,
 * nonlinear residual/Jacobian callbacks, time-scheme metadata, and optional
 * eigensolver settings. It does not own assembled algebraic state; calls to
 * `assemble()`, `solve()`, or `evolve()` create stateful system wrappers.
 *
 * @tparam SpaceType Scalar Feel++ function-space type.
 * @tparam Tags Semantic tags selecting legal system wrappers and algorithms.
 */
template <typename SpaceType, typename... Tags>
class Problem
{
  public:
    /** @brief Scalar function-space type. */
    using space_type = SpaceType;
    /** @brief Shared function-space pointer type. */
    using space_ptrtype = std::shared_ptr<space_type>;
    /** @brief Scalar value type used by algebraic objects. */
    using value_type = typename space_type::value_type;
    /** @brief Finite-element state type. */
    using element_type = typename space_type::element_type;
    /** @brief Feel++ backend type used by assembled systems. */
    using backend_type = Backend<value_type>;
    /** @brief Sparse matrix pointer type used by assembled systems. */
    using matrix_ptrtype = typename backend_type::sparse_matrix_ptrtype;
    /** @brief Vector pointer type used by assembled systems. */
    using vector_ptrtype = typename backend_type::vector_ptrtype;
    /** @brief Stored term descriptor container. */
    using term_container_type = std::vector<term_descriptor>;
    /** @brief Concrete bilinear form type used to replay matrix actions. */
    using bilinear_form_type = decltype( form2( _test = std::declval<space_ptrtype>(),
                                                _trial = std::declval<space_ptrtype>(),
                                                _matrix = std::declval<matrix_ptrtype>() ) );
    /** @brief Concrete linear form type used to replay vector actions. */
    using linear_form_type = decltype( form1( _test = std::declval<space_ptrtype>(),
                                              _vector = std::declval<vector_ptrtype>() ) );
    /** @brief Matrix/vector pair returned after constraints are applied. */
    using active_system_type = std::pair<matrix_ptrtype, vector_ptrtype>;
    /** @brief Type-erased matrix assembly callback. */
    using matrix_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased vector assembly callback. */
    using vector_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased time-step-aware matrix callback. */
    using timed_matrix_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased time-step-aware vector callback. */
    using timed_vector_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased nonlinear residual callback. */
    using residual_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased nonlinear Jacobian callback. */
    using jacobian_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased stationary constraint callback. */
    using constraint_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased time-step-aware constraint callback. */
    using timed_constraint_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Type-erased Picard iteration callback. */
    using picard_iteration_action_type = std::function<void( AssemblyContext const& )>;
    /** @brief Shared BDF time-scheme pointer type. */
    using bdf_time_scheme_ptrtype = std::shared_ptr<Feel::Bdf<space_type>>;
    /** @brief Shared Newmark time-scheme pointer type. */
    using newmark_time_scheme_ptrtype = std::shared_ptr<Feel::Newmark<space_type>>;
    /** @brief Storage variant for the attached time scheme. */
    using time_scheme_storage_type = std::variant<std::monostate, bdf_time_scheme_ptrtype, newmark_time_scheme_ptrtype>;
    /** @brief Symbolic parameter map type. */
    using parameter_values_type = std::map<std::string, double>;
    /** @brief Shared parameter map pointer type captured by expression callbacks. */
    using parameter_values_ptrtype = std::shared_ptr<parameter_values_type>;

    /** @brief True when `Tag` is present in the problem tag list. */
    template <typename Tag>
    static constexpr bool has_tag = ( std::same_as<Tag, Tags> || ... || false );

    /**
     * @brief Role-filtered view over stored term descriptors.
     */
    class role_terms_view
    {
      public:
        /**
         * @brief Create a view over terms matching `role`.
         */
        role_terms_view( Problem const& problem, term_role role )
            : M_problem( &problem ), M_role( role )
        {
        }

        /** @brief Role selected by this view. */
        term_role role() const noexcept { return M_role; }
        /** @brief Number of terms matching this view role. */
        std::size_t size() const noexcept { return M_problem->nTerms( M_role ); }
        /** @brief Return true when no terms match this view role. */
        bool empty() const noexcept { return this->size() == 0; }
        /** @brief Descriptor at `index` within the role-filtered view. */
        term_descriptor const& operator[]( std::size_t index ) const { return M_problem->term( M_role, index ); }
        /** @brief Metadata at `index` within the role-filtered view. */
        TermMetadata const& metadata( std::size_t index ) const { return ( *this )[index].metadata(); }

      private:
        /** @brief Problem whose terms are inspected. */
        Problem const* M_problem = nullptr;
        /** @brief Role selected by this view. */
        term_role M_role = term_role::lhs;
    };

    /**
     * @brief Construct a problem descriptor from a function space and name.
     */
    explicit Problem( space_ptrtype space, std::string name = {} )
        : M_space( std::move( space ) ),
          M_name( std::move( name ) ),
          M_parameterValues( std::make_shared<parameter_values_type>() )
    {
    }

    /** @brief Function space used by all scalar terms. */
    space_ptrtype const& space() const noexcept { return M_space; }
    /** @brief Problem/backend option prefix. */
    std::string const& name() const noexcept { return M_name; }

    /**
     * @brief Merge symbolic parameter values used by stored expression terms.
     *
     * Values are applied lazily when expression-based form helpers assemble.
     */
    Problem& setParameterValues( parameter_values_type const& values )
    {
        M_parameterValues->insert( values.begin(), values.end() );
        for ( auto const& [key, value] : values )
            ( *M_parameterValues )[key] = value;
        return *this;
    }

    /**
     * @brief Set one symbolic parameter value.
     */
    Problem& setParameterValue( std::string symbol, double value )
    {
        ( *M_parameterValues )[std::move( symbol )] = value;
        return *this;
    }

    /**
     * @brief Access current symbolic parameter values.
     */
    parameter_values_type const& parameterValues() const noexcept { return *M_parameterValues; }

    /**
     * @brief Select the symbolic name updated by evolution time schemes.
     */
    Problem& timeSymbol( std::string symbol )
    {
        M_timeSymbol = std::move( symbol );
        return *this;
    }

    /**
     * @brief Symbolic name updated with the current time before step assembly.
     */
    std::string const& timeSymbol() const noexcept { return M_timeSymbol; }

    /**
     * @brief Update problem parameters from a time-step context.
     */
    Problem& updateForTimeStep( time_step_context const& context )
    {
        return this->setParameterValue( M_timeSymbol, context.time );
    }

    /**
     * @brief Assemble this descriptor through `systems::assemble`.
     */
    template <typename... Args>
    auto assemble( Args&&... args ) const;

    /**
     * @brief Assemble and solve this descriptor through `systems::solve`.
     */
    template <typename... Args>
    auto solve( Args&&... args ) const;

    /**
     * @brief Create a stateful time-step system through `systems::evolve`.
     */
    template <typename... Args>
    auto evolve( Args&&... args ) const;

    /**
     * @brief Apply a user callback to every item of a range while building the problem.
     *
     * This helper is intentionally eager: it iterates immediately and lets the
     * callback add terms or constraints to the current problem. The callback may
     * either accept `(Problem&, item)` when it needs to mutate the problem, or
     * accept only `item` for side-effect-only iteration.
     *
     * @tparam Range Input range type containing model, boundary, or user descriptors.
     * @tparam Function Callback type invoked for each range item.
     * @param range Range of descriptors or values to visit.
     * @param function Callback invoked on every item.
     * @return Reference to this problem for fluent composition.
     */
    template <std::ranges::input_range Range, typename Function>
    Problem& forEach( Range&& range, Function&& function )
    {
        for ( auto&& item : range )
        {
            if constexpr ( std::is_invocable_v<Function&, Problem&, decltype( item )> )
                std::invoke( function, *this, item );
            else
                std::invoke( function, item );
        }
        return *this;
    }

    /**
     * @brief Add a term using a role tag such as `systems::lhs` or `systems::rhs`.
     */
    template <typename Role, typename Term>
        requires requires( Role r ) { systems::role( r ); }
    Problem& add( Role roleTag, Term&& term )
    {
        using role_type = std::decay_t<Role>;
        auto const termRole = systems::role( roleTag );
        this->addDescriptor( termRole );

        if constexpr ( std::same_as<role_type, lhs_t> ||
                       std::same_as<role_type, mass_t> ||
                       std::same_as<role_type, damping_t> ||
                       std::same_as<role_type, stiffness_t> ||
                       std::same_as<role_type, operator_a_t> ||
                       std::same_as<role_type, operator_b_t> )
        {
            this->addMatrixAction( termRole, std::forward<Term>( term ) );
        }
        else if constexpr ( std::same_as<role_type, rhs_t> )
        {
            this->addVectorAction( termRole, std::forward<Term>( term ) );
        }
        else if constexpr ( std::same_as<role_type, constraint_t> )
        {
            this->addConstraintAction( std::forward<Term>( term ) );
        }
        else if constexpr ( std::same_as<role_type, residual_t> )
        {
            this->addResidualAction( std::forward<Term>( term ) );
        }
        else if constexpr ( std::same_as<role_type, jacobian_t> )
        {
            this->addJacobianAction( std::forward<Term>( term ) );
        }
        else if constexpr ( std::same_as<role_type, iteration_t> )
        {
            this->addPicardIterationAction( std::forward<Term>( term ) );
        }
        else if constexpr ( std::same_as<role_type, time_scheme_t> )
        {
            this->setTimeScheme( std::forward<Term>( term ) );
        }

        return *this;
    }

    /**
     * @brief Add a matrix contribution to the lhs role.
     */
    template <typename Term>
    Problem& lhs( Term&& term )
    {
        this->addDescriptor( term_role::lhs );
        this->addMatrixAction( term_role::lhs, std::forward<Term>( term ) );
        return *this;
    }
    /**
     * @brief Add an expression-based lhs contribution with `_range` and `_expr`.
     */
    /**
     * @brief Add Fourier boundary contributions to lhs and rhs.
     *
     * Named arguments are `_range`, `_trial`, `_test`, `_alpha`, and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& lhs( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::lhs );
        this->addMatrixExpressionAction( term_role::lhs, std::move( range ), std::move( expr ) );
        return *this;
    }
    /** @brief View of lhs term descriptors. */
    role_terms_view lhs() const noexcept { return this->terms( term_role::lhs ); }
    /** @brief View of lhs term descriptors. */
    role_terms_view lhsTerms() const noexcept { return this->lhs(); }

    /**
     * @brief Declare a local Neumann auxiliary operator by replaying lhs terms.
     *
     * Constraint actions are intentionally excluded, so the resulting local
     * operator remains floating even when the main system has Dirichlet data.
     */
    Problem& localNeumannFromLhs( std::string name = "geneo" )
    {
        this->addAuxiliaryDescriptor( AuxiliaryOperatorDescriptor{
            std::move( name ),
            AuxiliaryOperatorKind::local_neumann,
            AuxiliaryOperatorSource::lhs_terms } );
        return *this;
    }

    /** @brief Add an explicit contribution to a named local Neumann operator. */
    template <typename Term>
    Problem& localNeumann( std::string name, Term&& term )
    {
        auto& entry = this->addAuxiliaryDescriptor( AuxiliaryOperatorDescriptor{
            std::move( name ),
            AuxiliaryOperatorKind::local_neumann,
            AuxiliaryOperatorSource::explicit_terms } );
        entry.actions.emplace_back( this->makeMatrixAction( std::forward<Term>( term ) ) );
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

    /** @brief Return true when at least one auxiliary operator is declared. */
    bool hasAuxiliaryOperators() const noexcept { return !M_auxiliaryOperators.empty(); }

    /**
     * @brief Add a vector contribution to the rhs role.
     */
    template <typename Term>
    Problem& rhs( Term&& term )
    {
        this->addDescriptor( term_role::rhs );
        this->addVectorAction( term_role::rhs, std::forward<Term>( term ) );
        return *this;
    }
    /**
     * @brief Add an expression-based rhs contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& rhs( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::rhs );
        this->addVectorExpressionAction( term_role::rhs, std::move( range ), std::move( expr ) );
        return *this;
    }
    /** @brief View of rhs term descriptors. */
    role_terms_view rhs() const noexcept { return this->terms( term_role::rhs ); }
    /** @brief View of rhs term descriptors. */
    role_terms_view rhsTerms() const noexcept { return this->rhs(); }

    /**
     * @brief Add a mass matrix contribution used by evolution assembly.
     */
    template <typename Term>
    Problem& mass( Term&& term )
    {
        this->addDescriptor( term_role::mass );
        this->addMatrixAction( term_role::mass, std::forward<Term>( term ) );
        return *this;
    }
    /**
     * @brief Add an expression-based mass contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& mass( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::mass );
        this->addMatrixExpressionAction( term_role::mass, std::move( range ),
                                         std::move( localRange ), std::move( expr ) );
        return *this;
    }
    /** @brief View of mass term descriptors. */
    role_terms_view mass() const noexcept { return this->terms( term_role::mass ); }
    /** @brief View of mass term descriptors. */
    role_terms_view massTerms() const noexcept { return this->mass(); }

    /**
     * @brief Add a damping matrix contribution used by second-order evolution.
     */
    template <typename Term>
    Problem& damping( Term&& term )
    {
        this->addDescriptor( term_role::damping );
        this->addMatrixAction( term_role::damping, std::forward<Term>( term ) );
        return *this;
    }
    /**
     * @brief Add an expression-based damping contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& damping( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::damping );
        this->addMatrixExpressionAction( term_role::damping, std::move( range ),
                                         std::move( localRange ), std::move( expr ) );
        return *this;
    }
    /** @brief View of damping term descriptors. */
    role_terms_view damping() const noexcept { return this->terms( term_role::damping ); }
    /** @brief View of damping term descriptors. */
    role_terms_view dampingTerms() const noexcept { return this->damping(); }

    /**
     * @brief Add a stiffness matrix contribution.
     */
    template <typename Term>
    Problem& stiffness( Term&& term )
    {
        this->addDescriptor( term_role::stiffness );
        this->addMatrixAction( term_role::stiffness, std::forward<Term>( term ) );
        return *this;
    }
    /**
     * @brief Add an expression-based stiffness contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& stiffness( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::stiffness );
        this->addMatrixExpressionAction( term_role::stiffness, std::move( range ),
                                         std::move( localRange ), std::move( expr ) );
        return *this;
    }
    /** @brief View of stiffness term descriptors. */
    role_terms_view stiffness() const noexcept { return this->terms( term_role::stiffness ); }
    /** @brief View of stiffness term descriptors. */
    role_terms_view stiffnessTerms() const noexcept { return this->stiffness(); }

    /**
     * @brief Add a nonlinear residual contribution.
     */
    template <typename Term>
    Problem& residual( Term&& term )
    {
        this->addDescriptor( term_role::residual );
        this->addResidualAction( std::forward<Term>( term ) );
        return *this;
    }
    /** @brief View of residual term descriptors. */
    role_terms_view residual() const noexcept { return this->terms( term_role::residual ); }
    /** @brief View of residual term descriptors. */
    role_terms_view residualTerms() const noexcept { return this->residual(); }

    /**
     * @brief Add a nonlinear Jacobian contribution.
     */
    template <typename Term>
    Problem& jacobian( Term&& term )
    {
        this->addDescriptor( term_role::jacobian );
        this->addJacobianAction( std::forward<Term>( term ) );
        return *this;
    }
    /** @brief View of Jacobian term descriptors. */
    role_terms_view jacobian() const noexcept { return this->terms( term_role::jacobian ); }
    /** @brief View of Jacobian term descriptors. */
    role_terms_view jacobianTerms() const noexcept { return this->jacobian(); }

    /**
     * @brief Add one Picard fixed-point assembly callback.
     *
     * The callback is replayed at every Picard iteration on the current
     * linear forms. It may accept `(state, a, l)` or `(a, l)`.
     */
    template <typename Iteration>
    Problem& iteration( Iteration&& iteration )
    {
        this->addDescriptor( term_role::iteration );
        this->addPicardIterationAction( std::forward<Iteration>( iteration ) );
        return *this;
    }
    /** @brief View of Picard iteration term descriptors. */
    role_terms_view iteration() const noexcept { return this->terms( term_role::iteration ); }
    /** @brief View of Picard iteration term descriptors. */
    role_terms_view iterationTerms() const noexcept { return this->iteration(); }

    /**
     * @brief Alias for `.iteration(...)` in fixed-point problem builders.
     */
    template <typename Iteration>
    Problem& picard( Iteration&& iteration )
    {
        return this->iteration( std::forward<Iteration>( iteration ) );
    }

    /**
     * @brief Add eigen operator A contribution.
     */
    template <typename Term>
    Problem& operatorA( Term&& term )
    {
        this->addDescriptor( term_role::operator_a );
        this->addMatrixAction( term_role::operator_a, std::forward<Term>( term ) );
        return *this;
    }
    /** @brief View of operator A term descriptors. */
    role_terms_view operatorA() const noexcept { return this->terms( term_role::operator_a ); }
    /** @brief View of operator A term descriptors. */
    role_terms_view operatorATerms() const noexcept { return this->operatorA(); }

    /**
     * @brief Add eigen operator B contribution for generalized problems.
     */
    template <typename Term>
    Problem& operatorB( Term&& term )
    {
        this->addDescriptor( term_role::operator_b );
        this->addMatrixAction( term_role::operator_b, std::forward<Term>( term ) );
        return *this;
    }
    /** @brief View of operator B term descriptors. */
    role_terms_view operatorB() const noexcept { return this->terms( term_role::operator_b ); }
    /** @brief View of operator B term descriptors. */
    role_terms_view operatorBTerms() const noexcept { return this->operatorB(); }

    /**
     * @brief Add a deferred constraint callback.
     */
    template <typename Constraint>
    Problem& constraint( Constraint&& c )
    {
        this->addDescriptor( term_role::constraint );
        this->addConstraintAction( std::forward<Constraint>( c ) );
        return *this;
    }
    /** @brief View of constraint term descriptors. */
    role_terms_view constraints() const noexcept { return this->terms( term_role::constraint ); }
    /** @brief View of constraint term descriptors. */
    role_terms_view constraintTerms() const noexcept { return this->constraints(); }

    /**
     * @brief Alias for `constraint(...)` matching boundary-condition terminology.
     */
    template <typename Constraint>
    Problem& bc( Constraint&& c ) { return this->constraint( std::forward<Constraint>( c ) ); }
    /** @brief View of boundary-condition/constraint term descriptors. */
    role_terms_view bcTerms() const noexcept { return this->constraintTerms(); }

    /**
     * @brief Add an essential Dirichlet boundary condition.
     *
     * Supported named arguments include `_range`, `_element`, `_expr`,
     * optional `_type`, and optional `_value_on_diagonal`.
     */
    template <typename... Args>
    Problem& dirichlet( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto element = namedArgs.get( _element );
        auto expr = namedArgs.get( _expr );
        std::string type = namedArgs.get_else( _type, std::string{ "elimination" } );
        double valueOnDiagonal = namedArgs.get_else( _value_on_diagonal, 1.0 );
        auto storedExpr = std::make_shared<std::decay_t<decltype( expr )>>( std::move( expr ) );
        auto parameters = M_parameterValues;

        this->addDescriptor( term_role::constraint );
        this->addConstraintActions(
            [range, element, storedExpr, parameters, type, valueOnDiagonal]( auto& a, auto& l )
            {
                detail::setParameterValuesIfSupported( *storedExpr, *parameters );
                a += on( _range = range,
                         _rhs = l,
                         _element = element,
                         _expr = *storedExpr,
                         _type = type,
                         _value_on_diagonal = valueOnDiagonal );
            },
            [localRange, element, storedExpr, parameters, type, valueOnDiagonal]( auto& a, auto& l )
            {
                detail::setParameterValuesIfSupported( *storedExpr, *parameters );
                a += on( _range = localRange,
                         _rhs = l,
                         _element = element,
                         _expr = *storedExpr,
                         _type = type,
                         _value_on_diagonal = valueOnDiagonal );
            } );
        return *this;
    }

    /**
     * @brief Add a diffusion contribution as stiffness.
     */
    template <typename Term>
    Problem& diffusion( Term&& term ) { return this->stiffness( std::forward<Term>( term ) ); }

    /**
     * @brief Add an expression-based diffusion contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& diffusion( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::stiffness );
        this->addMatrixExpressionAction( term_role::stiffness,
                                         std::move( range ),
                                         std::move( localRange ),
                                         std::move( expr ) );
        return *this;
    }

    /**
     * @brief Add a reaction term to the lhs with reaction metadata.
     */
    template <typename Term>
    Problem& reaction( Term&& term )
    {
        this->addDescriptor( term_role::lhs, namedMetadata( term_role::lhs, TermKind::Reaction, "reaction" ) );
        this->addMatrixAction( term_role::lhs, std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add an expression-based reaction contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& reaction( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::lhs, namedMetadata( term_role::lhs, TermKind::Reaction, "reaction" ) );
        this->addMatrixExpressionAction( term_role::lhs, std::move( range ),
                                         std::move( localRange ), std::move( expr ) );
        return *this;
    }

    /**
     * @brief Add a source term to the rhs with source metadata.
     */
    template <typename Term>
    Problem& source( Term&& term )
    {
        this->addDescriptor( term_role::rhs, namedMetadata( term_role::rhs, TermKind::Source, "source" ) );
        this->addVectorAction( term_role::rhs, std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add an expression-based source contribution with `_range` and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& source( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto expr = namedArgs.get( _expr );
        this->addDescriptor( term_role::rhs, namedMetadata( term_role::rhs, TermKind::Source, "source" ) );
        this->addVectorExpressionAction( term_role::rhs, std::move( range ), std::move( expr ) );
        return *this;
    }

    /**
     * @brief Add a Neumann boundary contribution with `_range`, `_test`, and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& neumann( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto test = namedArgs.get( _test );
        auto expr = namedArgs.get( _expr );
        auto storedExpr = std::make_shared<std::decay_t<decltype( expr )>>( std::move( expr ) );
        auto parameters = M_parameterValues;
        this->addDescriptor( term_role::rhs, boundaryMetadata( term_role::rhs, "neumann" ) );
        this->addVectorActionClosure(
            term_role::rhs,
            [range, test, storedExpr, parameters]( space_ptrtype const& space, vector_ptrtype const& vector )
            {
                detail::setParameterValuesIfSupported( *storedExpr, *parameters );
                auto testExpr = detail::testExpression( test );
                form1( _test = space, _vector = vector ) +=
                    integrate( _range = range,
                               _expr = *storedExpr * testExpr );
            } );
        return *this;
    }

    /**
     * @brief Add an already-built Neumann rhs contribution.
     */
    template <typename Term>
        requires( !NA::is_named_argument_v<std::remove_cvref_t<Term>> )
    Problem& neumann( Term&& term )
    {
        this->addDescriptor( term_role::rhs, boundaryMetadata( term_role::rhs, "neumann" ) );
        this->addVectorAction( term_role::rhs, std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add Robin boundary contributions to lhs and rhs.
     *
     * Named arguments are `_range`, `_trial`, `_test`, `_alpha`, and `_expr`.
     */
    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& robin( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto trial = namedArgs.get( _trial );
        auto test = namedArgs.get( _test );
        auto alpha = namedArgs.get( _alpha );
        auto expr = namedArgs.get( _expr );
        auto storedAlpha = std::make_shared<std::decay_t<decltype( alpha )>>( std::move( alpha ) );
        auto storedExpr = std::make_shared<std::decay_t<decltype( expr )>>( std::move( expr ) );
        auto parameters = M_parameterValues;

        this->addDescriptor( term_role::lhs, boundaryMetadata( term_role::lhs, "robin" ) );
        this->addMatrixActionClosures(
            term_role::lhs,
            [range, trial, test, storedAlpha, parameters]( space_ptrtype const& space, matrix_ptrtype const& matrix )
            {
                detail::setParameterValuesIfSupported( *storedAlpha, *parameters );
                auto trialExpr = detail::trialExpression( trial );
                auto testExpr = detail::testExpression( test );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = range,
                               _expr = *storedAlpha * trialExpr * testExpr );
            },
            [localRange, trial, test, storedAlpha, parameters]( space_ptrtype const& space, matrix_ptrtype const& matrix )
            {
                detail::setParameterValuesIfSupported( *storedAlpha, *parameters );
                auto trialExpr = detail::trialExpression( trial );
                auto testExpr = detail::testExpression( test );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = localRange,
                               _expr = *storedAlpha * trialExpr * testExpr );
            } );
        this->addDescriptor( term_role::rhs, boundaryMetadata( term_role::rhs, "robin" ) );
        this->addVectorActionClosure(
            term_role::rhs,
            [range, test, storedExpr, parameters]( space_ptrtype const& space, vector_ptrtype const& vector )
            {
                detail::setParameterValuesIfSupported( *storedExpr, *parameters );
                auto testExpr = detail::testExpression( test );
                form1( _test = space, _vector = vector ) +=
                    integrate( _range = range,
                               _expr = *storedExpr * testExpr );
            } );
        return *this;
    }

    /**
     * @brief Add an already-built Robin lhs contribution.
     */
    template <typename Term>
        requires( !NA::is_named_argument_v<std::remove_cvref_t<Term>> )
    Problem& robin( Term&& term )
    {
        this->addDescriptor( term_role::lhs, boundaryMetadata( term_role::lhs, "robin" ) );
        this->addMatrixAction( term_role::lhs, std::forward<Term>( term ) );
        return *this;
    }

    template <typename... Args>
        requires( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
    Problem& fourier( Args&&... args )
    {
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto range = namedArgs.get( _range );
        auto localRange = namedArgs.get_else( _range_extended, range );
        auto trial = namedArgs.get( _trial );
        auto test = namedArgs.get( _test );
        auto alpha = namedArgs.get( _alpha );
        auto expr = namedArgs.get( _expr );
        auto storedAlpha = std::make_shared<std::decay_t<decltype( alpha )>>( std::move( alpha ) );
        auto storedExpr = std::make_shared<std::decay_t<decltype( expr )>>( std::move( expr ) );
        auto parameters = M_parameterValues;

        this->addDescriptor( term_role::lhs, boundaryMetadata( term_role::lhs, "fourier" ) );
        this->addMatrixActionClosures(
            term_role::lhs,
            [range, trial, test, storedAlpha, parameters]( space_ptrtype const& space, matrix_ptrtype const& matrix )
            {
                detail::setParameterValuesIfSupported( *storedAlpha, *parameters );
                auto trialExpr = detail::trialExpression( trial );
                auto testExpr = detail::testExpression( test );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = range,
                               _expr = *storedAlpha * trialExpr * testExpr );
            },
            [localRange, trial, test, storedAlpha, parameters]( space_ptrtype const& space, matrix_ptrtype const& matrix )
            {
                detail::setParameterValuesIfSupported( *storedAlpha, *parameters );
                auto trialExpr = detail::trialExpression( trial );
                auto testExpr = detail::testExpression( test );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = localRange,
                               _expr = *storedAlpha * trialExpr * testExpr );
            } );
        this->addDescriptor( term_role::rhs, boundaryMetadata( term_role::rhs, "fourier" ) );
        this->addVectorActionClosure(
            term_role::rhs,
            [range, test, storedAlpha, storedExpr, parameters]( space_ptrtype const& space, vector_ptrtype const& vector )
            {
                detail::setParameterValuesIfSupported( *storedAlpha, *parameters );
                detail::setParameterValuesIfSupported( *storedExpr, *parameters );
                auto testExpr = detail::testExpression( test );
                form1( _test = space, _vector = vector ) +=
                    integrate( _range = range,
                               _expr = *storedAlpha * *storedExpr * testExpr );
            } );
        return *this;
    }

    /**
     * @brief Add an already-built Fourier lhs contribution.
     */
    template <typename Term>
        requires( !NA::is_named_argument_v<std::remove_cvref_t<Term>> )
    Problem& fourier( Term&& term )
    {
        this->addDescriptor( term_role::lhs, boundaryMetadata( term_role::lhs, "fourier" ) );
        this->addMatrixAction( term_role::lhs, std::forward<Term>( term ) );
        return *this;
    }

    /**
     * @brief Add weak Dirichlet lhs and rhs contributions together.
     */
    template <typename LhsTerm, typename RhsTerm>
    Problem& weakDirichlet( LhsTerm&& lhsTerm, RhsTerm&& rhsTerm )
    {
        this->lhs( std::forward<LhsTerm>( lhsTerm ) );
        this->rhs( std::forward<RhsTerm>( rhsTerm ) );
        return *this;
    }

    /**
     * @brief Add weak Dirichlet terms with a distinct owned-plus-ghost lhs.
     *
     * The physical boundary contribution is preserved in the inherited local
     * operator. Only the artificial overlap boundary remains natural.
     */
    template <typename LhsTerm, typename LocalLhsTerm, typename RhsTerm>
    Problem& weakDirichlet( LhsTerm&& lhsTerm, LocalLhsTerm&& localLhsTerm,
                            RhsTerm&& rhsTerm )
    {
        this->addDescriptor( term_role::lhs );
        M_matrixActions.emplace_back(
            term_role::lhs, this->makeMatrixAction( std::forward<LhsTerm>( lhsTerm ) ) );
        M_localMatrixActions.emplace_back(
            term_role::lhs, this->makeMatrixAction( std::forward<LocalLhsTerm>( localLhsTerm ) ) );
        this->rhs( std::forward<RhsTerm>( rhsTerm ) );
        return *this;
    }

    /**
     * @brief Attach a `feelts` time scheme to an evolution problem.
     *
     * Supported schemes are currently `std::shared_ptr<Bdf<space_type>>` for
     * first-order evolution problems and `std::shared_ptr<Newmark<space_type>>`
     * for second-order evolution problems.
     */
    template <typename TimeScheme>
    Problem& timeScheme( TimeScheme&& scheme ) { return this->add( systems::time_scheme, std::forward<TimeScheme>( scheme ) ); }
    /** @brief View of time-scheme term descriptors. */
    role_terms_view timeSchemeTerms() const noexcept { return this->terms( term_role::time_scheme ); }

    /**
     * @brief Return true if the problem has an attached time scheme.
     */
    bool hasTimeScheme() const noexcept { return !std::holds_alternative<std::monostate>( M_timeScheme ); }

    /**
     * @brief Return true when the attached time scheme is BDF.
     */
    bool hasBdfTimeScheme() const noexcept { return std::holds_alternative<bdf_time_scheme_ptrtype>( M_timeScheme ); }

    /**
     * @brief Return true when the attached time scheme is Newmark.
     */
    bool hasNewmarkTimeScheme() const noexcept { return std::holds_alternative<newmark_time_scheme_ptrtype>( M_timeScheme ); }

    /**
     * @brief Access the attached BDF time scheme, or an empty pointer.
     */
    bdf_time_scheme_ptrtype bdfTimeScheme() const
    {
        if ( auto scheme = std::get_if<bdf_time_scheme_ptrtype>( &M_timeScheme ) )
            return *scheme;
        return {};
    }

    /**
     * @brief Access the attached Newmark time scheme, or an empty pointer.
     */
    newmark_time_scheme_ptrtype newmarkTimeScheme() const
    {
        if ( auto scheme = std::get_if<newmark_time_scheme_ptrtype>( &M_timeScheme ) )
            return *scheme;
        return {};
    }

    /** @brief Set requested number of eigenvalues. */
    Problem& nev( int value ) noexcept
    {
        M_nev = value;
        M_hasNev = true;
        return *this;
    }
    /** @brief Requested number of eigenvalues. */
    int nev() const noexcept { return M_nev; }
    /** @brief Return true when `nev` was explicitly set on the problem. */
    bool hasNev() const noexcept { return M_hasNev; }

    /** @brief Set number of eigensolver basis vectors. */
    Problem& ncv( int value ) noexcept
    {
        M_ncv = value;
        return *this;
    }
    /** @brief Number of eigensolver basis vectors. */
    int ncv() const noexcept { return M_ncv; }

    /** @brief Set maximum projected dimension for the eigensolver. */
    Problem& mpd( int value ) noexcept
    {
        M_mpd = value;
        return *this;
    }
    /** @brief Maximum projected dimension for the eigensolver. */
    int mpd() const noexcept { return M_mpd; }

    /** @brief Set maximum number of eigensolver iterations. */
    Problem& maxit( int value ) noexcept
    {
        M_eigenMaxit = value;
        return *this;
    }
    /** @brief Maximum number of eigensolver iterations. */
    int maxit() const noexcept { return M_eigenMaxit; }

    /** @brief Set eigensolver tolerance. */
    Problem& tolerance( double value ) noexcept
    {
        M_eigenTolerance = value;
        return *this;
    }
    /** @brief Eigensolver tolerance. */
    double tolerance() const noexcept { return M_eigenTolerance; }

    /**
     * @brief Set eigensolver spectrum selection.
     */
    template <typename SpectrumSelector>
    Problem& spectrum( SpectrumSelector&& selector )
    {
        M_eigenSpectrum = detail::eigenSpectrum( std::forward<SpectrumSelector>( selector ) );
        M_hasSpectrum = true;
        return *this;
    }
    /** @brief Return true when spectrum selection was explicitly set. */
    bool hasSpectrum() const noexcept { return M_hasSpectrum; }
    /** @brief Selected eigensolver spectrum. */
    PositionOfSpectrum spectrum() const noexcept { return M_eigenSpectrum; }

    /**
     * @brief Set eigensolver spectral transform.
     */
    template <typename TransformSelector>
    Problem& transform( TransformSelector&& selector )
    {
        M_eigenTransform = detail::eigenTransform( std::forward<TransformSelector>( selector ) );
        M_hasTransform = true;
        return *this;
    }
    /** @brief Return true when spectral transform was explicitly set. */
    bool hasTransform() const noexcept { return M_hasTransform; }
    /** @brief Selected eigensolver spectral transform. */
    SpectralTransformType transform() const noexcept { return M_eigenTransform; }

    /** @brief All stored term descriptors. */
    term_container_type const& terms() const noexcept { return M_terms; }
    /** @brief Role-filtered view over stored term descriptors. */
    role_terms_view terms( term_role role ) const noexcept { return role_terms_view{ *this, role }; }

    /**
     * @brief Structured summaries of all stored terms without assembly.
     */
    std::vector<TermSummary> termSummaries() const { return systems::termSummaries( M_terms ); }

    /**
     * @brief Conservative metadata consistency diagnostics.
     */
    std::vector<TermConsistencyIssue> checkTerms() const { return systems::checkTermConsistency( M_terms ); }

    /**
     * @brief Write a side-effect-free diagnostic description of this problem.
     */
    std::ostream& describe( std::ostream& os ) const { return systems::describeTerms( os, M_name, M_terms ); }

    /**
     * @brief Return a side-effect-free diagnostic description of this problem.
     */
    std::string describe() const { return systems::describeTerms( M_name, M_terms ); }

    /**
     * @brief Side-effect-free assembly plan derived from term metadata.
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
     * @brief Scheduler-ready task graph derived from term metadata.
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

    /** @brief Number of all stored term descriptors. */
    std::size_t nTerms() const noexcept { return M_terms.size(); }
    /** @brief Number of stored term descriptors matching `role`. */
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
        throw std::out_of_range( "Feel::systems::Problem term index out of range" );
    }

    /**
     * @brief Assemble affine matrix and vector contributions.
     *
     * Matrix roles `lhs`, `mass`, `damping`, and `stiffness` are assembled into
     * `matrix`; `rhs` roles are assembled into `vector`. Constraints are not
     * applied by this helper.
     */
    void assembleAffine( matrix_ptrtype const& matrix, vector_ptrtype const& vector ) const
    {
        for ( auto const& [role, action] : M_matrixActions )
        {
            if ( role == term_role::lhs || role == term_role::mass ||
                 role == term_role::damping || role == term_role::stiffness )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::Affine, role );
                context.matrix = &matrix;
                action( context );
            }
        }
        for ( auto const& [role, action] : M_vectorActions )
        {
            if ( role == term_role::rhs )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::Affine, role );
                context.vector = &vector;
                action( context );
            }
        }
    }

    /**
     * @brief Assemble matrix contributions matching one role.
     */
    void assembleMatrixTerms( term_role role, matrix_ptrtype const& matrix ) const
    {
        for ( auto const& [actionRole, action] : M_matrixActions )
        {
            if ( actionRole == role )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::MatrixTerm, role );
                context.matrix = &matrix;
                action( context );
            }
        }
    }

    /** @brief Assemble one named auxiliary operator without applying constraints. */
    void assembleAuxiliaryOperator( std::string const& name, matrix_ptrtype const& matrix ) const
    {
        auto const it = std::find_if( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                                      [&name]( auto const& entry )
                                      { return entry.descriptor.name == name; } );
        if ( it == M_auxiliaryOperators.end() )
            throw std::out_of_range( "Feel::systems::Problem auxiliary operator not found: " + name );

        if ( it->descriptor.source == AuxiliaryOperatorSource::lhs_terms )
        {
            this->assembleMatrixTerms( term_role::lhs, matrix );
            return;
        }

        for ( auto const& action : it->actions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::MatrixTerm, term_role::lhs );
            context.matrix = &matrix;
            action( context );
        }
    }

    /**
     * @brief Return true when the stationary problem can be replayed locally.
     *
     * The replay contains every inherited matrix contribution. Physical
     * boundary/PML terms and essential constraints therefore remain part of
     * the local operator; only the overlap boundary introduced by the
     * subdomain decomposition is left natural.
     */
    bool hasInferredLocalOperator() const noexcept
    {
        return !M_localMatrixActions.empty();
    }

    /**
     * @brief Assemble the inherited problem on an owned-plus-ghost local map.
     *
     * Semantic terms may provide an extended range for this replay. The
     * returned pair accounts for constraint implementations that replace the
     * matrix/vector rather than modifying them in place.
     */
    active_system_type assembleInferredLocalOperator( matrix_ptrtype const& matrix,
                                                       vector_ptrtype const& vector ) const
    {
        for ( auto const& [role, action] : M_localMatrixActions )
        {
            if ( role == term_role::lhs || role == term_role::mass ||
                 role == term_role::damping || role == term_role::stiffness )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::MatrixTerm, role );
                context.matrix = &matrix;
                action( context );
            }
        }
        if ( M_localConstraintActions.empty() )
            return { matrix, vector };

        auto a = form2( _test = M_space, _trial = M_space, _matrix = matrix );
        auto l = form1( _test = M_space, _vector = vector );
        for ( auto const& action : M_localConstraintActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Constraint,
                                                      term_role::constraint );
            context.matrix = &matrix;
            context.vector = &vector;
            context.bilinearForm = &a;
            context.linearForm = &l;
            action( context );
        }
        return detail::completeActiveSystem( a.activeSystem( l ), matrix, vector );
    }

    /**
     * @brief Assemble matrix contributions matching one role for a time step.
     */
    void assembleMatrixTerms( time_step_context const& context, term_role role, matrix_ptrtype const& matrix ) const
    {
        for ( auto const& [actionRole, action] : M_matrixActions )
        {
            if ( actionRole == role )
            {
                auto assemblyContext = this->makeAssemblyContext( AssemblyMode::EvolutionStep, role );
                assemblyContext.timeStep = &context;
                assemblyContext.matrix = &matrix;
                action( assemblyContext );
            }
        }
        for ( auto const& [actionRole, action] : M_timedMatrixActions )
        {
            if ( actionRole == role )
            {
                auto assemblyContext = this->makeAssemblyContext( AssemblyMode::EvolutionStep, role );
                assemblyContext.timeStep = &context;
                assemblyContext.matrix = &matrix;
                action( assemblyContext );
            }
        }
    }

    /**
     * @brief Assemble vector contributions matching one role.
     */
    void assembleVectorTerms( term_role role, vector_ptrtype const& vector ) const
    {
        for ( auto const& [actionRole, action] : M_vectorActions )
        {
            if ( actionRole == role )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::VectorTerm, role );
                context.vector = &vector;
                action( context );
            }
        }
    }

    /**
     * @brief Assemble vector contributions matching one role for a time step.
     */
    void assembleVectorTerms( time_step_context const& context, term_role role, vector_ptrtype const& vector ) const
    {
        for ( auto const& [actionRole, action] : M_vectorActions )
        {
            if ( actionRole == role )
            {
                auto assemblyContext = this->makeAssemblyContext( AssemblyMode::EvolutionStep, role );
                assemblyContext.timeStep = &context;
                assemblyContext.vector = &vector;
                action( assemblyContext );
            }
        }
        for ( auto const& [actionRole, action] : M_timedVectorActions )
        {
            if ( actionRole == role )
            {
                auto assemblyContext = this->makeAssemblyContext( AssemblyMode::EvolutionStep, role );
                assemblyContext.timeStep = &context;
                assemblyContext.vector = &vector;
                action( assemblyContext );
            }
        }
    }

    /**
     * @brief Apply stationary constraints to a prepared matrix/vector pair.
     */
    active_system_type applyConstraints( matrix_ptrtype const& matrix, vector_ptrtype const& vector ) const
    {
        if ( M_constraintActions.empty() )
            return { matrix, vector };

        auto a = form2( _test = M_space, _trial = M_space, _matrix = matrix );
        auto l = form1( _test = M_space, _vector = vector );
        for ( auto const& action : M_constraintActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Constraint, term_role::constraint );
            context.matrix = &matrix;
            context.vector = &vector;
            context.bilinearForm = &a;
            context.linearForm = &l;
            action( context );
        }
        return detail::completeActiveSystem( a.activeSystem( l ), matrix, vector );
    }

    /**
     * @brief Apply stationary and time-step constraints to a prepared pair.
     */
    active_system_type applyConstraints( time_step_context const& context,
                                         matrix_ptrtype const& matrix,
                                         vector_ptrtype const& vector ) const
    {
        if ( M_constraintActions.empty() && M_timedConstraintActions.empty() )
            return { matrix, vector };

        auto a = form2( _test = M_space, _trial = M_space, _matrix = matrix );
        auto l = form1( _test = M_space, _vector = vector );
        for ( auto const& action : M_constraintActions )
        {
            auto assemblyContext = this->makeAssemblyContext( AssemblyMode::Constraint, term_role::constraint );
            assemblyContext.timeStep = &context;
            assemblyContext.matrix = &matrix;
            assemblyContext.vector = &vector;
            assemblyContext.bilinearForm = &a;
            assemblyContext.linearForm = &l;
            action( assemblyContext );
        }
        for ( auto const& action : M_timedConstraintActions )
        {
            auto assemblyContext = this->makeAssemblyContext( AssemblyMode::Constraint, term_role::constraint );
            assemblyContext.timeStep = &context;
            assemblyContext.matrix = &matrix;
            assemblyContext.vector = &vector;
            assemblyContext.bilinearForm = &a;
            assemblyContext.linearForm = &l;
            action( assemblyContext );
        }
        return detail::completeActiveSystem( a.activeSystem( l ), matrix, vector );
    }

    /**
     * @brief Assemble the linear active system, including deferred constraints.
     */
    active_system_type assembleLinear( matrix_ptrtype const& matrix, vector_ptrtype const& vector ) const
    {
        this->assembleAffine( matrix, vector );
        return detail::completeActiveSystem( this->applyConstraints( matrix, vector ), matrix, vector );
    }

    /**
     * @brief Assemble all nonlinear residual-only contributions into `vector`.
     */
    void assembleResidualTerms( vector_ptrtype const& vector ) const
    {
        for ( auto const& action : M_residualActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Residual, term_role::residual );
            context.vector = &vector;
            action( context );
        }
    }

    /**
     * @brief Assemble all nonlinear residual-only contributions for `state`.
     */
    void assembleResidualTerms( element_type const& state, vector_ptrtype const& vector ) const
    {
        for ( auto const& action : M_residualActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Residual, term_role::residual );
            context.state = &state;
            context.vector = &vector;
            action( context );
        }
    }

    /**
     * @brief Assemble all nonlinear residual-only contributions for a time step.
     */
    void assembleResidualTerms( time_step_context const& context, element_type const& state, vector_ptrtype const& vector ) const
    {
        for ( auto const& action : M_residualActions )
        {
            auto assemblyContext = this->makeAssemblyContext( AssemblyMode::Residual, term_role::residual );
            assemblyContext.timeStep = &context;
            assemblyContext.state = &state;
            assemblyContext.vector = &vector;
            action( assemblyContext );
        }
    }

    /**
     * @brief Assemble all nonlinear Jacobian-only contributions into `matrix`.
     */
    void assembleJacobianTerms( matrix_ptrtype const& matrix ) const
    {
        for ( auto const& action : M_jacobianActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Jacobian, term_role::jacobian );
            context.matrix = &matrix;
            action( context );
        }
    }

    /**
     * @brief Assemble all nonlinear Jacobian-only contributions for `state`.
     */
    void assembleJacobianTerms( element_type const& state, matrix_ptrtype const& matrix ) const
    {
        for ( auto const& action : M_jacobianActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Jacobian, term_role::jacobian );
            context.state = &state;
            context.matrix = &matrix;
            action( context );
        }
    }

    /**
     * @brief Assemble all nonlinear Jacobian-only contributions for a time step.
     */
    void assembleJacobianTerms( time_step_context const& context, element_type const& state, matrix_ptrtype const& matrix ) const
    {
        for ( auto const& action : M_jacobianActions )
        {
            auto assemblyContext = this->makeAssemblyContext( AssemblyMode::Jacobian, term_role::jacobian );
            assemblyContext.timeStep = &context;
            assemblyContext.state = &state;
            assemblyContext.matrix = &matrix;
            action( assemblyContext );
        }
    }

    /**
     * @brief Assemble Picard iteration contributions for the current state.
     */
    void assemblePicardIteration( element_type const& state,
                                  matrix_ptrtype const& matrix,
                                  vector_ptrtype const& vector ) const
    {
        auto a = form2( _test = M_space, _trial = M_space, _matrix = matrix );
        auto l = form1( _test = M_space, _vector = vector );
        for ( auto const& action : M_picardIterationActions )
        {
            auto context = this->makeAssemblyContext( AssemblyMode::Picard, term_role::iteration );
            context.state = &state;
            context.matrix = &matrix;
            context.vector = &vector;
            context.bilinearForm = &a;
            context.linearForm = &l;
            action( context );
        }
    }

    /**
     * @brief Assemble eigen operator A and optional operator B.
     */
    void assembleEigen( matrix_ptrtype const& matrixA, matrix_ptrtype const& matrixB ) const
    {
        for ( auto const& [role, action] : M_matrixActions )
        {
            if ( role == term_role::operator_a )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::Eigen, role );
                context.matrix = &matrixA;
                action( context );
            }
            else if ( role == term_role::operator_b )
            {
                auto context = this->makeAssemblyContext( AssemblyMode::Eigen, role );
                context.matrix = &matrixB;
                action( context );
            }
        }
    }

  private:
    /**
     * @brief Create the non-owning context passed to stored callbacks.
     */
    AssemblyContext makeAssemblyContext( AssemblyMode mode, term_role role ) const
    {
        AssemblyContext context;
        context.mode = mode;
        context.role = role;
        context.problem = this;
        context.space = &M_space;
        return context;
    }

    /** @brief Add a descriptor with default metadata for `role`. */
    void addDescriptor( term_role role )
    {
        this->addDescriptor( role, defaultTermMetadata( role ) );
    }

    /** @brief Add a descriptor with explicit metadata. */
    void addDescriptor( term_role role, TermMetadata metadata )
    {
        M_terms.push_back( term_descriptor{ role, roleName( role ), std::move( metadata ) } );
    }

    /** @brief Build named metadata for a model term. */
    static TermMetadata namedMetadata( term_role role, TermKind kind, std::string_view name )
    {
        auto metadata = defaultTermMetadata( role );
        metadata.kind = kind;
        metadata.name = std::string( name );
        return metadata;
    }

    /** @brief Build boundary-condition metadata for a role. */
    static TermMetadata boundaryMetadata( term_role role, std::string_view name )
    {
        auto metadata = defaultTermMetadata( role );
        metadata.kind = TermKind::BoundaryCondition;
        metadata.support = TermSupport::BoundaryFacet;
        metadata.name = std::string( name );
        return metadata;
    }

    /** @brief Build a stationary matrix action from a term or callback. */
    template <typename Term>
    matrix_action_type makeMatrixAction( Term&& term )
    {
        auto storedTerm = std::make_shared<std::decay_t<Term>>( std::forward<Term>( term ) );
        using stored_type = std::decay_t<Term>;
        if constexpr ( requires( bilinear_form_type& a, stored_type& t ) { a += t; } )
        {
            return [storedTerm]( AssemblyContext const& context )
                {
                    auto const& space = context.spaceAs<space_ptrtype>();
                    auto const& matrix = context.matrixAs<matrix_ptrtype>();
                    auto a = form2( _test = space, _trial = space, _matrix = matrix );
                    a += *storedTerm;
                };
        }
        else
        {
            return [storedTerm]( AssemblyContext const& context )
                {
                    using stored_type = std::remove_cvref_t<decltype( *storedTerm )>;
                    auto const& space = context.spaceAs<space_ptrtype>();
                    auto const& matrix = context.matrixAs<matrix_ptrtype>();
                    if constexpr ( std::is_invocable_v<stored_type&, space_ptrtype const&, matrix_ptrtype const&> )
                    {
                        ( *storedTerm )( space, matrix );
                    }
                    else
                    {
                        auto a = form2( _test = space, _trial = space, _matrix = matrix );
                        if constexpr ( std::is_invocable_v<stored_type&, bilinear_form_type&> )
                            ( *storedTerm )( a );
                        else
                            a += *storedTerm;
                    }
                };
        }
    }

    /** @brief Store a matrix action, adapting supported term/callback shapes. */
    template <typename Term>
    void addMatrixAction( term_role role, Term&& term )
    {
        using stored_type = std::decay_t<Term>;
        if constexpr ( requires( bilinear_form_type& a, stored_type& t ) { a += t; } )
        {
            auto action = this->makeMatrixAction( std::forward<Term>( term ) );
            M_matrixActions.emplace_back( role, action );
            M_localMatrixActions.emplace_back( role, std::move( action ) );
        }
        else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, space_ptrtype const&, matrix_ptrtype const&> )
        {
            auto storedTerm = std::make_shared<stored_type>( std::forward<Term>( term ) );
            M_timedMatrixActions.emplace_back(
                role,
                [storedTerm]( AssemblyContext const& context )
                {
                    ( *storedTerm )( context.timeStepAs<time_step_context>(),
                                     context.spaceAs<space_ptrtype>(),
                                     context.matrixAs<matrix_ptrtype>() );
                } );
        }
        else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, bilinear_form_type&> )
        {
            auto storedTerm = std::make_shared<stored_type>( std::forward<Term>( term ) );
            M_timedMatrixActions.emplace_back(
                role,
                [storedTerm]( AssemblyContext const& context )
                {
                    auto const& space = context.spaceAs<space_ptrtype>();
                    auto const& matrix = context.matrixAs<matrix_ptrtype>();
                    auto a = form2( _test = space, _trial = space, _matrix = matrix );
                    ( *storedTerm )( context.timeStepAs<time_step_context>(), a );
                } );
        }
        else
        {
            auto action = this->makeMatrixAction( std::forward<Term>( term ) );
            M_matrixActions.emplace_back( role, action );
            M_localMatrixActions.emplace_back( role, std::move( action ) );
        }
    }

    struct auxiliary_operator_entry
    {
        AuxiliaryOperatorDescriptor descriptor;
        std::vector<matrix_action_type> actions;
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
            throw std::invalid_argument( "Feel::systems::Problem auxiliary operator declaration conflicts with existing operator: " + descriptor.name );
        return *it;
    }

    /** @brief Store a vector action, adapting supported term/callback shapes. */
    template <typename Term>
    void addVectorAction( term_role role, Term&& term )
    {
        auto storedTerm = std::make_shared<std::decay_t<Term>>( std::forward<Term>( term ) );
        using stored_type = std::decay_t<Term>;
        if constexpr ( requires( linear_form_type& l, stored_type& t ) { l += t; } )
        {
            M_vectorActions.emplace_back(
                role,
                [storedTerm]( AssemblyContext const& context )
                {
                    auto const& space = context.spaceAs<space_ptrtype>();
                    auto const& vector = context.vectorAs<vector_ptrtype>();
                    auto l = form1( _test = space, _vector = vector );
                    l += *storedTerm;
                } );
        }
        else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, space_ptrtype const&, vector_ptrtype const&> )
        {
            M_timedVectorActions.emplace_back(
                role,
                [storedTerm]( AssemblyContext const& context )
                {
                    ( *storedTerm )( context.timeStepAs<time_step_context>(),
                                     context.spaceAs<space_ptrtype>(),
                                     context.vectorAs<vector_ptrtype>() );
                } );
        }
        else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, linear_form_type&> )
        {
            M_timedVectorActions.emplace_back(
                role,
                [storedTerm]( AssemblyContext const& context )
                {
                    auto const& space = context.spaceAs<space_ptrtype>();
                    auto const& vector = context.vectorAs<vector_ptrtype>();
                    auto l = form1( _test = space, _vector = vector );
                    ( *storedTerm )( context.timeStepAs<time_step_context>(), l );
                } );
        }
        else
        {
            M_vectorActions.emplace_back(
                role,
                [storedTerm]( AssemblyContext const& context )
                {
                    using stored_type = std::remove_cvref_t<decltype( *storedTerm )>;
                    auto const& space = context.spaceAs<space_ptrtype>();
                    auto const& vector = context.vectorAs<vector_ptrtype>();
                    if constexpr ( std::is_invocable_v<stored_type&, space_ptrtype const&, vector_ptrtype const&> )
                    {
                        ( *storedTerm )( space, vector );
                    }
                    else
                    {
                        auto l = form1( _test = space, _vector = vector );
                        if constexpr ( std::is_invocable_v<stored_type&, linear_form_type&> )
                            ( *storedTerm )( l );
                        else
                            l += *storedTerm;
                    }
                } );
        }
    }

    /** @brief Store an expression-based matrix action with shared parameter state. */
    template <typename Range, typename Expr>
    void addMatrixExpressionAction( term_role role, Range&& range, Expr&& expression )
    {
        auto storedRange = std::make_shared<std::decay_t<Range>>( std::forward<Range>( range ) );
        auto storedExpression = std::make_shared<std::decay_t<Expr>>( std::forward<Expr>( expression ) );
        auto parameters = M_parameterValues;
        matrix_action_type action =
            [storedRange, storedExpression, parameters]( AssemblyContext const& context )
            {
                auto const& space = context.spaceAs<space_ptrtype>();
                auto const& matrix = context.matrixAs<matrix_ptrtype>();
                detail::setParameterValuesIfSupported( *storedExpression, *parameters );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = *storedRange,
                               _expr = *storedExpression );
            };
        M_matrixActions.emplace_back( role, action );
        M_localMatrixActions.emplace_back( role, std::move( action ) );
    }

    /** @brief Store distinct distributed and owned-plus-ghost expression actions. */
    template <typename Range, typename LocalRange, typename Expr>
    void addMatrixExpressionAction( term_role role, Range&& range,
                                    LocalRange&& localRange, Expr&& expression )
    {
        auto storedRange = std::make_shared<std::decay_t<Range>>( std::forward<Range>( range ) );
        auto storedLocalRange = std::make_shared<std::decay_t<LocalRange>>( std::forward<LocalRange>( localRange ) );
        auto storedExpression = std::make_shared<std::decay_t<Expr>>( std::forward<Expr>( expression ) );
        auto parameters = M_parameterValues;
        M_matrixActions.emplace_back(
            role,
            [storedRange, storedExpression, parameters]( AssemblyContext const& context )
            {
                auto const& space = context.spaceAs<space_ptrtype>();
                auto const& matrix = context.matrixAs<matrix_ptrtype>();
                detail::setParameterValuesIfSupported( *storedExpression, *parameters );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = *storedRange, _expr = *storedExpression );
            } );
        M_localMatrixActions.emplace_back(
            role,
            [storedLocalRange, storedExpression, parameters]( AssemblyContext const& context )
            {
                auto const& space = context.spaceAs<space_ptrtype>();
                auto const& matrix = context.matrixAs<matrix_ptrtype>();
                detail::setParameterValuesIfSupported( *storedExpression, *parameters );
                form2( _test = space, _trial = space, _matrix = matrix ) +=
                    integrate( _range = *storedLocalRange, _expr = *storedExpression );
            } );
    }

    /** @brief Store an expression-based vector action with shared parameter state. */
    template <typename Range, typename Expr>
    void addVectorExpressionAction( term_role role, Range&& range, Expr&& expression )
    {
        auto storedRange = std::make_shared<std::decay_t<Range>>( std::forward<Range>( range ) );
        auto storedExpression = std::make_shared<std::decay_t<Expr>>( std::forward<Expr>( expression ) );
        auto parameters = M_parameterValues;
        M_vectorActions.emplace_back(
            role,
            [storedRange, storedExpression, parameters]( AssemblyContext const& context )
            {
                auto const& space = context.spaceAs<space_ptrtype>();
                auto const& vector = context.vectorAs<vector_ptrtype>();
                detail::setParameterValuesIfSupported( *storedExpression, *parameters );
                form1( _test = space, _vector = vector ) +=
                    integrate( _range = *storedRange,
                               _expr = *storedExpression );
            } );
    }

    /** @brief Store a low-level matrix callback accepting space and matrix. */
    template <typename Action>
    void addMatrixActionClosure( term_role role, Action&& action )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        matrix_action_type wrapped =
            [storedAction]( AssemblyContext const& context )
            {
                ( *storedAction )( context.spaceAs<space_ptrtype>(),
                                   context.matrixAs<matrix_ptrtype>() );
            };
        M_matrixActions.emplace_back( role, wrapped );
        M_localMatrixActions.emplace_back( role, std::move( wrapped ) );
    }

    /** @brief Store distinct distributed and owned-plus-ghost callbacks. */
    template <typename Action, typename LocalAction>
    void addMatrixActionClosures( term_role role, Action&& action, LocalAction&& localAction )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        auto storedLocalAction = std::make_shared<std::decay_t<LocalAction>>( std::forward<LocalAction>( localAction ) );
        M_matrixActions.emplace_back(
            role,
            [storedAction]( AssemblyContext const& context )
            {
                ( *storedAction )( context.spaceAs<space_ptrtype>(),
                                   context.matrixAs<matrix_ptrtype>() );
            } );
        M_localMatrixActions.emplace_back(
            role,
            [storedLocalAction]( AssemblyContext const& context )
            {
                ( *storedLocalAction )( context.spaceAs<space_ptrtype>(),
                                        context.matrixAs<matrix_ptrtype>() );
            } );
    }

    /** @brief Store a low-level vector callback accepting space and vector. */
    template <typename Action>
    void addVectorActionClosure( term_role role, Action&& action )
    {
        auto storedAction = std::make_shared<std::decay_t<Action>>( std::forward<Action>( action ) );
        M_vectorActions.emplace_back(
            role,
            [storedAction]( AssemblyContext const& context )
            {
                ( *storedAction )( context.spaceAs<space_ptrtype>(),
                                   context.vectorAs<vector_ptrtype>() );
            } );
    }

    /** @brief Store a stationary or time-step-aware constraint action. */
    template <typename Constraint>
    void addConstraintAction( Constraint&& constraint )
    {
        auto storedConstraint = std::make_shared<std::decay_t<Constraint>>( std::forward<Constraint>( constraint ) );
        using stored_type = std::decay_t<Constraint>;
        if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, bilinear_form_type&, linear_form_type&> )
        {
            M_timedConstraintActions.emplace_back(
                [storedConstraint]( AssemblyContext const& context )
                {
                    ( *storedConstraint )( context.timeStepAs<time_step_context>(),
                                           context.bilinearFormAs<bilinear_form_type>(),
                                           context.linearFormAs<linear_form_type>() );
                } );
        }
        else
        {
            constraint_action_type action =
                [storedConstraint]( AssemblyContext const& context )
                {
                    auto& a = context.bilinearFormAs<bilinear_form_type>();
                    auto& l = context.linearFormAs<linear_form_type>();
                    if constexpr ( std::is_invocable_v<std::decay_t<Constraint>&, bilinear_form_type&, linear_form_type&> )
                        ( *storedConstraint )( a, l );
                    else
                        a += *storedConstraint;
                };
            M_constraintActions.emplace_back( action );
            M_localConstraintActions.emplace_back( std::move( action ) );
        }
    }

    /** @brief Store distinct distributed and owned-plus-ghost constraints. */
    template <typename Constraint, typename LocalConstraint>
    void addConstraintActions( Constraint&& constraint, LocalConstraint&& localConstraint )
    {
        auto storedConstraint = std::make_shared<std::decay_t<Constraint>>( std::forward<Constraint>( constraint ) );
        auto storedLocalConstraint = std::make_shared<std::decay_t<LocalConstraint>>( std::forward<LocalConstraint>( localConstraint ) );
        M_constraintActions.emplace_back(
            [storedConstraint]( AssemblyContext const& context )
            {
                ( *storedConstraint )( context.bilinearFormAs<bilinear_form_type>(),
                                       context.linearFormAs<linear_form_type>() );
            } );
        M_localConstraintActions.emplace_back(
            [storedLocalConstraint]( AssemblyContext const& context )
            {
                ( *storedLocalConstraint )( context.bilinearFormAs<bilinear_form_type>(),
                                            context.linearFormAs<linear_form_type>() );
            } );
    }

    /** @brief Store the concrete BDF or Newmark time-scheme pointer. */
    template <typename TimeScheme>
    void setTimeScheme( TimeScheme&& scheme )
    {
        using scheme_type = std::remove_cvref_t<TimeScheme>;
        if constexpr ( std::same_as<scheme_type, bdf_time_scheme_ptrtype> )
        {
            M_timeScheme = std::forward<TimeScheme>( scheme );
        }
        else if constexpr ( std::same_as<scheme_type, newmark_time_scheme_ptrtype> )
        {
            M_timeScheme = std::forward<TimeScheme>( scheme );
        }
        else
        {
            static_assert( sizeof( scheme_type ) == 0,
                           "Feel::systems::Problem::timeScheme supports shared_ptr<Bdf<space_type>> and shared_ptr<Newmark<space_type>>" );
        }
    }

    /** @brief Store a residual callback, adapting supported signatures. */
    template <typename Term>
    void addResidualAction( Term&& term )
    {
        auto storedTerm = std::make_shared<std::decay_t<Term>>( std::forward<Term>( term ) );
        M_residualActions.emplace_back(
            [storedTerm]( AssemblyContext const& context )
            {
                using stored_type = std::decay_t<Term>;
                auto const& space = context.spaceAs<space_ptrtype>();
                auto const& vector = context.vectorAs<vector_ptrtype>();
                auto l = form1( _test = space, _vector = vector );
                if constexpr ( requires( linear_form_type& form, stored_type& t ) { form += t; } )
                {
                    l += *storedTerm;
                }
                else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, element_type const&, space_ptrtype const&, vector_ptrtype const&> )
                {
                    if ( !context.hasTimeStep() || !context.hasState() )
                        throw std::logic_error( "Feel::systems residual action requires a time step and nonlinear state" );
                    ( *storedTerm )( context.timeStepAs<time_step_context>(),
                                     context.stateAs<element_type>(),
                                     space,
                                     vector );
                }
                else if constexpr ( std::is_invocable_v<stored_type&, element_type const&, space_ptrtype const&, vector_ptrtype const&> )
                {
                    if ( !context.hasState() )
                        throw std::logic_error( "Feel::systems residual action requires a nonlinear state" );
                    ( *storedTerm )( context.stateAs<element_type>(), space, vector );
                }
                else
                {
                    if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, element_type const&, linear_form_type&> )
                    {
                        if ( !context.hasTimeStep() || !context.hasState() )
                            throw std::logic_error( "Feel::systems residual action requires a time step and nonlinear state" );
                        ( *storedTerm )( context.timeStepAs<time_step_context>(),
                                         context.stateAs<element_type>(),
                                         l );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, element_type const&, linear_form_type&> )
                    {
                        if ( !context.hasState() )
                            throw std::logic_error( "Feel::systems residual action requires a nonlinear state" );
                        ( *storedTerm )( context.stateAs<element_type>(), l );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, space_ptrtype const&, vector_ptrtype const&> )
                    {
                        if ( !context.hasTimeStep() )
                            throw std::logic_error( "Feel::systems residual action requires a time step" );
                        ( *storedTerm )( context.timeStepAs<time_step_context>(), space, vector );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, space_ptrtype const&, vector_ptrtype const&> )
                    {
                        ( *storedTerm )( space, vector );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, linear_form_type&> )
                    {
                        if ( !context.hasTimeStep() )
                            throw std::logic_error( "Feel::systems residual action requires a time step" );
                        ( *storedTerm )( context.timeStepAs<time_step_context>(), l );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, linear_form_type&> )
                        ( *storedTerm )( l );
                    else
                        l += *storedTerm;
                }
            } );
    }

    /** @brief Store a Jacobian callback, adapting supported signatures. */
    template <typename Term>
    void addJacobianAction( Term&& term )
    {
        auto storedTerm = std::make_shared<std::decay_t<Term>>( std::forward<Term>( term ) );
        M_jacobianActions.emplace_back(
            [storedTerm]( AssemblyContext const& context )
            {
                using stored_type = std::decay_t<Term>;
                auto const& space = context.spaceAs<space_ptrtype>();
                auto const& matrix = context.matrixAs<matrix_ptrtype>();
                auto a = form2( _test = space, _trial = space, _matrix = matrix );
                if constexpr ( requires( bilinear_form_type& form, stored_type& t ) { form += t; } )
                {
                    a += *storedTerm;
                }
                else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, element_type const&, space_ptrtype const&, matrix_ptrtype const&> )
                {
                    if ( !context.hasTimeStep() || !context.hasState() )
                        throw std::logic_error( "Feel::systems jacobian action requires a time step and nonlinear state" );
                    ( *storedTerm )( context.timeStepAs<time_step_context>(),
                                     context.stateAs<element_type>(),
                                     space,
                                     matrix );
                }
                else if constexpr ( std::is_invocable_v<stored_type&, element_type const&, space_ptrtype const&, matrix_ptrtype const&> )
                {
                    if ( !context.hasState() )
                        throw std::logic_error( "Feel::systems jacobian action requires a nonlinear state" );
                    ( *storedTerm )( context.stateAs<element_type>(), space, matrix );
                }
                else
                {
                    if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, element_type const&, bilinear_form_type&> )
                    {
                        if ( !context.hasTimeStep() || !context.hasState() )
                            throw std::logic_error( "Feel::systems jacobian action requires a time step and nonlinear state" );
                        ( *storedTerm )( context.timeStepAs<time_step_context>(),
                                         context.stateAs<element_type>(),
                                         a );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, element_type const&, bilinear_form_type&> )
                    {
                        if ( !context.hasState() )
                            throw std::logic_error( "Feel::systems jacobian action requires a nonlinear state" );
                        ( *storedTerm )( context.stateAs<element_type>(), a );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, space_ptrtype const&, matrix_ptrtype const&> )
                    {
                        if ( !context.hasTimeStep() )
                            throw std::logic_error( "Feel::systems jacobian action requires a time step" );
                        ( *storedTerm )( context.timeStepAs<time_step_context>(), space, matrix );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, space_ptrtype const&, matrix_ptrtype const&> )
                    {
                        ( *storedTerm )( space, matrix );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, time_step_context const&, bilinear_form_type&> )
                    {
                        if ( !context.hasTimeStep() )
                            throw std::logic_error( "Feel::systems jacobian action requires a time step" );
                        ( *storedTerm )( context.timeStepAs<time_step_context>(), a );
                    }
                    else if constexpr ( std::is_invocable_v<stored_type&, bilinear_form_type&> )
                        ( *storedTerm )( a );
                    else
                        a += *storedTerm;
                }
            } );
    }

    /** @brief Store a Picard iteration callback. */
    template <typename Iteration>
    void addPicardIterationAction( Iteration&& iteration )
    {
        auto storedIteration = std::make_shared<std::decay_t<Iteration>>( std::forward<Iteration>( iteration ) );
        M_picardIterationActions.emplace_back(
            [storedIteration]( AssemblyContext const& context )
            {
                using stored_type = std::decay_t<Iteration>;
                auto const& state = context.stateAs<element_type>();
                auto& a = context.bilinearFormAs<bilinear_form_type>();
                auto& l = context.linearFormAs<linear_form_type>();
                if constexpr ( std::is_invocable_v<stored_type&, element_type const&, bilinear_form_type&, linear_form_type&> )
                    ( *storedIteration )( state, a, l );
                else
                    ( *storedIteration )( a, l );
            } );
    }

    /** @brief Function space used by all scalar terms. */
    space_ptrtype M_space;
    /** @brief Problem/backend option prefix. */
    std::string M_name;
    /** @brief Stored term descriptors in insertion order. */
    term_container_type M_terms;
    /** @brief Matrix actions grouped with their storage role. */
    std::vector<std::pair<term_role, matrix_action_type>> M_matrixActions;
    /** @brief Stationary matrix actions replayed on an owned-plus-ghost map. */
    std::vector<std::pair<term_role, matrix_action_type>> M_localMatrixActions;
    /** @brief Named auxiliary operators and their explicit assembly actions. */
    std::vector<auxiliary_operator_entry> M_auxiliaryOperators;
    /** @brief Vector actions grouped with their storage role. */
    std::vector<std::pair<term_role, vector_action_type>> M_vectorActions;
    /** @brief Time-step-aware matrix actions grouped with their storage role. */
    std::vector<std::pair<term_role, timed_matrix_action_type>> M_timedMatrixActions;
    /** @brief Time-step-aware vector actions grouped with their storage role. */
    std::vector<std::pair<term_role, timed_vector_action_type>> M_timedVectorActions;
    /** @brief Nonlinear residual callbacks. */
    std::vector<residual_action_type> M_residualActions;
    /** @brief Nonlinear Jacobian callbacks. */
    std::vector<jacobian_action_type> M_jacobianActions;
    /** @brief Stationary constraint callbacks. */
    std::vector<constraint_action_type> M_constraintActions;
    /** @brief Constraints replayed on the owned-plus-ghost local operator. */
    std::vector<constraint_action_type> M_localConstraintActions;
    /** @brief Time-step-aware constraint callbacks. */
    std::vector<timed_constraint_action_type> M_timedConstraintActions;
    /** @brief Picard iteration callbacks. */
    std::vector<picard_iteration_action_type> M_picardIterationActions;
    /** @brief Attached BDF/Newmark time scheme, when any. */
    time_scheme_storage_type M_timeScheme;
    /** @brief Shared symbolic parameter map captured by expression callbacks. */
    parameter_values_ptrtype M_parameterValues;
    /** @brief Symbol updated from time-step contexts. */
    std::string M_timeSymbol = "t";
    /** @brief Requested number of eigenvalues. */
    int M_nev = 1;
    /** @brief Requested eigensolver basis size. */
    int M_ncv = -1;
    /** @brief Requested eigensolver maximum projected dimension. */
    int M_mpd = -2;
    /** @brief Requested eigensolver maximum iterations. */
    int M_eigenMaxit = -1;
    /** @brief Requested eigensolver tolerance. */
    double M_eigenTolerance = -1.0;
    /** @brief Requested eigensolver spectrum selection. */
    PositionOfSpectrum M_eigenSpectrum = LARGEST_MAGNITUDE;
    /** @brief Requested eigensolver spectral transform. */
    SpectralTransformType M_eigenTransform = SHIFT;
    /** @brief True when spectrum selection was explicitly set. */
    bool M_hasSpectrum = false;
    /** @brief True when spectral transform was explicitly set. */
    bool M_hasTransform = false;
    /** @brief True when `nev` was explicitly set. */
    bool M_hasNev = false;
};

/**
 * @brief Backward-compatible alias for the scalar problem descriptor.
 */
template <typename SpaceType, typename... Tags>
using BasicProblem = Problem<SpaceType, Tags...>;

} // namespace Feel::systems

namespace Feel
{

/**
 * @brief Build a stationary scalar linear problem descriptor.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor tagged with `systems::linear_t` and
 * `systems::stationary_t`.
 */
template <typename... Ts>
auto linearProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "linear_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::linear_t, systems::stationary_t>( space, std::move( name ) );
}

/**
 * @brief Build a stationary elliptic scalar problem descriptor.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor tagged as linear, stationary, and elliptic.
 */
template <typename... Ts>
auto ellipticProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "elliptic_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::linear_t, systems::stationary_t, systems::elliptic_t>( space, std::move( name ) );
}

/**
 * @brief Build a first-order parabolic scalar problem descriptor.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor tagged as linear, first order, and parabolic.
 */
template <typename... Ts>
auto parabolicProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "parabolic_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::linear_t, systems::first_order_t, systems::parabolic_t>( space, std::move( name ) );
}

/**
 * @brief Build a second-order hyperbolic scalar problem descriptor.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor tagged as linear, second order, and hyperbolic.
 */
template <typename... Ts>
auto hyperbolicProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "hyperbolic_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::linear_t, systems::second_order_t, systems::hyperbolic_t>( space, std::move( name ) );
}

/**
 * @brief Build a stationary residual/Jacobian nonlinear problem descriptor.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor tagged with `systems::nonlinear_t`.
 */
template <typename... Ts>
auto nonlinearProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "nonlinear_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::nonlinear_t, systems::stationary_t>( space, std::move( name ) );
}

/**
 * @brief Build a stationary Picard fixed-point problem descriptor.
 *
 * The descriptor stores `.iteration(...)` callbacks replayed by
 * `systems::PicardSystem` as `A(u_k) u_{k+1} = b(u_k)`.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor tagged with `systems::picard_t`.
 */
template <typename... Ts>
auto picardProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "picard_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::picard_t, systems::stationary_t>( space, std::move( name ) );
}

/**
 * @brief Build a stationary generalized eigenvalue problem descriptor.
 *
 * @param ts Feel++ named arguments. Requires `_space`; accepts optional `_name`.
 * @return Problem descriptor carrying `operatorA` and optional `operatorB` terms.
 */
template <typename... Ts>
auto generalizedEigenvalueProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "generalized_eigenvalue_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::stationary_t, systems::generalized_eigen_t>( space, std::move( name ) );
}

/**
 * @brief Build a standard stationary eigenvalue problem descriptor.
 *
 * Omitting `.operatorB(...)` makes `EigenSystem` solve `A x = lambda x`.
 * Adding `.operatorB(...)` upgrades the same descriptor to the generalized
 * form `A x = lambda B x`.
 */
template <typename... Ts>
auto eigenvalueProblem( Ts&&... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>( ts )... );
    auto space = args.get( _space );
    std::string name = args.get_else( _name, "eigenvalue_problem" );
    using space_type = typename std::remove_cvref_t<decltype( space )>::element_type;
    return systems::Problem<space_type, systems::stationary_t, systems::generalized_eigen_t>( space, std::move( name ) );
}

} // namespace Feel

#endif // FEELPP_FEELSYSTEM_PROBLEM_HPP
