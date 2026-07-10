/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
   \file sb9_common.hpp
   \brief Shared SB9 shell operator infrastructure
 */
#ifndef FEELPP_VF_SB9_COMMON_HPP
#define FEELPP_VF_SB9_COMMON_HPP 1

#include <array>
#include <numbers>
#include <utility>

#include <feel/feelvf/basis.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/shellgeometric.hpp>
#include <feel/feelvf/voigt.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Base cache for SB9 element-local geometric data.
 *
 * This class stores the shell cell geometry quantities that are common to all
 * SB9 contributions. Contribution-specific caches derive from it and add only
 * the coefficients required by their operator family.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9KernelBase
{
public:
    /// Scalar type used by the geometry and generated coefficients.
    using value_type = typename GeometryDataType::value_type;
    /// Dense matrix type used by the shell geometry cache.
    using matrix_type = typename GeometryDataType::matrix_type;
    /// Three components of one local frame axis.
    using frame_axes_type = std::array<value_type, 3>;
    /// Three scalar coefficients projected on a local frame.
    using row_type = std::array<value_type, 3>;
    /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = 8;
    /// Number of displacement components handled by SB9 operators.
    static constexpr uint16_type component_count = 3;

    /**
     * \brief Build the shared cache from shell cell geometry data.
     *
     * The constructor stores a reference to the geometry data, computes the
     * inverse-transpose of the mid-surface Jacobian, and copies the local frame
     * components used by all SB9 contribution kernels.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9KernelBase( GeometryDataType const& data )
        :
        M_data( data ),
        M_invJ0( data.jacobian0.fullPivLu().inverse().transpose() )
    {
        for ( uint16_type component = 0; component < component_count; ++component )
            for ( uint16_type axis = 0; axis < component_count; ++axis )
                M_frame[component][axis] = data.frame( component, axis );
    }

protected:
    /**
     * \brief Return the Mandel scaling used for off-diagonal symmetric entries.
     *
     * Feel++ symmetric storage uses Mandel notation for shell strain vectors;
     * shear-like components are therefore scaled by `1/sqrt(2)` when raw
     * tensor coefficients are placed in storage.
     */
    static value_type mandelShearScale()
    {
        return value_type( 1 ) / static_cast<value_type>( std::numbers::sqrt2_v<double> );
    }

    /**
     * \brief Dot product between a coefficient row and a frame vector.
     *
     * \param row Three coefficient values in physical coordinates.
     * \param frame Local frame components for one displacement component.
     * \return Projection of \p row on \p frame.
     */
    static value_type dot( row_type const& row, frame_axes_type const& frame )
    {
        return row[0] * frame[0] + row[1] * frame[1] + row[2] * frame[2];
    }

    /**
     * \brief Fill membrane-like Mandel coefficients for one local dof.
     *
     * The same layout is used by membrane and bending terms: normal components
     * populate storage entries 0 and 1, and the in-plane shear component
     * populates storage entry 3 with Mandel scaling.
     *
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in symmetric storage order.
     * \param component Displacement component index.
     * \param bx First in-plane derivative coefficient.
     * \param by Second in-plane derivative coefficient.
     */
    template <typename VectorType>
    void fillMembraneCoefficients( VectorType& coeff,
                                   uint16_type component,
                                   value_type bx,
                                   value_type by ) const
    {
        auto const& frame = M_frame[component];
        coeff( 0 ) = bx * frame[0];
        coeff( 1 ) = by * frame[1];
        coeff( 3 ) = mandelShearScale() * ( bx * frame[1] + by * frame[0] );
    }

    /**
     * \brief Fill pinching Mandel coefficients for one local dof.
     *
     * The SB9 pinching operator contributes only to the transverse normal
     * component in Mandel storage (entry 2). The coefficient is obtained by
     * projecting the local pinching contribution onto the element frame.
     *
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in symmetric storage order.
     * \param component Displacement component index.
     * \param bz Local pinching coefficient.
     */
    template <typename VectorType>
    void fillPinchingCoefficients( VectorType& coeff,
                                   uint16_type component,
                                   value_type bz ) const
    {
        auto const& frame = M_frame[component];
        coeff( 2 ) = bz * frame[2];
    }

    /**
     * \brief Fill transverse-shearing Mandel coefficients for one local dof.
     *
     * The SB9 transverse-shearing operator contributes only to the shear
     * components in Mandel storage (entries 4 and 5). Each component is obtained
     * by projecting the local transverse-shear derivatives onto the element frame and 
     * applying the Mandel scaling factor.
     * 
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in symmetric storage order.
     * \param component Displacement component index.
     * \param s0 First shear derivative coefficient (correspond to the entry 4).
     * \param s1 Second shear derivative coefficient (correspond to the entry 5).
     */
    template <typename VectorType>
    void fillShearingCoefficients( VectorType& coeff,
                                   uint16_type component,
                                   value_type s0x,
                                   value_type s0y,
                                   value_type s0z,
                                   value_type s1x,
                                   value_type s1y,
                                   value_type s1z ) const
    {
        auto const& frame = M_frame[component];
        coeff( 4 ) = mandelShearScale() * ( s0x * frame[0] + s0y * frame[1] + s0z * frame[2] );
        coeff( 5 ) = mandelShearScale() * ( s1x * frame[0] + s1y * frame[1] + s1z * frame[2] );
    }

    /**
     * \brief Fill mode stabilization coefficients for one local dof.
     *
     * These helpers build the three-component coefficient vectors associated
     * with the four mode stabilization blocks by projecting the coefficients
     * onto the element frame. Bs2 has the same structure as Bs1, so it reuses
     * the same helper.
     * 
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in symmetric storage order.
     * \param component Displacement component index.
     */
    template <typename VectorType>
    void fillBs1Coefficients( VectorType& coeff,
                              uint16_type component,
                              value_type bs1 ) const
    {
        auto const& frame = M_frame[component];
        coeff( 2 ) = bs1 * frame[2];
    }
    template <typename VectorType>
    void fillBs3Coefficients( VectorType& coeff,
                              uint16_type component,
                              value_type bs3 ) const
    {
        auto const& frame = M_frame[component];
        coeff( 0 ) = bs3 * frame[0];
        coeff( 1 ) = bs3 * frame[1];
    }
    template <typename VectorType>
    void fillBs4Coefficients( VectorType& coeff,
                              uint16_type component,
                              value_type bs4 ) const
    {
        auto const& frame = M_frame[component];
        coeff( 0 ) = bs4 * frame[0];
        coeff( 1 ) = bs4 * frame[1];
        coeff( 2 ) = bs4 * frame[2];
    }
    
    /**
     * \brief Fill transverse shear stabilization coefficients for one local dof.
     *
     * 
     * These helpers build the two-component coefficient vectors associated
     * with the two transverse shear stabilization blocks by projecting the
     * coefficients onto the element frame.
     * 
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in symmetric storage order.
     * \param component Displacement component index.
     */
    template <typename VectorType>
    void fillBc1Coefficients( VectorType& coeff,
                              uint16_type component,
                              value_type bc11,
                              value_type bc12,
                              value_type bc13,
                              value_type bc21,
                              value_type bc22,
                              value_type bc23 ) const
    {
        auto const& frame = M_frame[component];
        coeff( 0 ) = bc11 * frame[0] + bc12 * frame[1] + bc13 * frame[2];
        coeff( 1 ) = bc21 * frame[0] + bc22 * frame[1] + bc23 * frame[2];
    }
    template <typename VectorType>
    void fillBc2Coefficients( VectorType& coeff,
                              uint16_type component,
                              value_type bc1,
                              value_type bc2,
                              value_type bc3 ) const
    {
        auto const& frame = M_frame[component];
        coeff( 1 ) = bc1 * frame[0] + bc2 * frame[1] + bc3 * frame[2];
    }

    /// Referenced shell geometry data for the current element.
    GeometryDataType const& M_data;
    /// Inverse-transpose of the mid-surface Jacobian.
    matrix_type M_invJ0;
    /// Local orthonormal frame copied per displacement component.
    std::array<frame_axes_type, component_count> M_frame{};
};

/**
 * \brief Generic Feel++ expression node for an SB9 vector-valued contribution.
 *
 * The operator is parameterized by a contribution cache and a contribution kind.
 * It is reused by bending, pinching, and shear headers. The nested tensor
 * performs the element-local coefficient assembly when the expression is
 * evaluated in a variational form.
 *
 * \tparam ElementType Trial or test basis proxy element type.
 * \tparam Role Feel++ basis role, either __TEST or __TRIAL.
 * \tparam KernelCacheTemplate Contribution-specific cache template.
 * \tparam Kind Contribution selector understood by the cache.
 */
template <typename ElementType, OperatorType Role, template <typename> typename KernelCacheTemplate, auto Kind>
class SB9VectorOperator
{
public:
    static_assert( Role == __TEST || Role == __TRIAL,
                   "SB9VectorOperator requires a trial or test basis role" );

    /// Expression context flags requested from the geometric mapping context.
    static const size_type context = 0;
    /// Marks this expression as a terminal node in the Feel++ expression tree.
    static inline const bool is_terminal = true;

    /// Basis proxy element type carried by the expression.
    using element_type = ElementType;
    /// Function space associated with \ref element_type.
    using functionspace_type = typename element_type::functionspace_type;
    /// Reference finite element type.
    using fe_type = typename functionspace_type::reference_element_type;
    /// Current expression type.
    using this_type = SB9VectorOperator<element_type, Role, KernelCacheTemplate, Kind>;
    /// Expression type returned by Feel++ expression machinery.
    using expression_type = this_type;
    /// Scalar value type of the function space.
    using value_type = typename functionspace_type::value_type;
    /// Size of a 3D symmetric tensor in Feel++ symmetric storage.
    static constexpr int storage_size = symmetric_storage_size_v<3>;
    /// Runtime value type produced by the expression.
    using evaluate_type = Eigen::Matrix<value_type, storage_size, 1>;

    /// Trait reporting whether this expression contains the requested test basis.
    template<typename Func>
    struct HasTestFunction
    {
        /// True when \c Func is the expression test basis finite element.
        static inline const bool result = ( Role == __TEST ) && boost::is_same<Func, fe_type>::value;
    };

    /// Trait reporting whether this expression contains the requested trial basis.
    template<typename Func>
    struct HasTrialFunction
    {
        /// True when \c Func is the expression trial basis finite element.
        static inline const bool result = ( Role == __TRIAL ) && boost::is_same<Func, fe_type>::value;
    };

    /// Compile-time test-basis flag used by Feel++ form assembly.
    template<typename Func>
    static inline const bool has_test_basis = ( Role == __TEST ) && boost::is_same<Func, fe_type>::value;

    /// Compile-time trial-basis flag used by Feel++ form assembly.
    template<typename Func>
    static inline const bool has_trial_basis = ( Role == __TRIAL ) && boost::is_same<Func, fe_type>::value;

    /// Test basis type exposed to the expression system.
    using test_basis = std::conditional_t<Role == __TEST, fe_type, std::nullptr_t>;
    /// Trial basis type exposed to the expression system.
    using trial_basis = std::conditional_t<Role == __TRIAL, fe_type, std::nullptr_t>;

    /// Lambda rebinding hook required by the Feel++ expression system.
    template<typename... TheExpr>
    struct Lambda
    {
        /// Rebound expression type.
        using type = expression_type;
    };

    /**
     * \brief Construct the SB9 vector operator for a trial or test basis proxy.
     * \param element Basis proxy element.
     */
    SB9VectorOperator( element_type element )
        :
        M_element( std::move( element ) )
    {
    }

    /**
     * \brief Return this expression for Feel++ expression expansion.
     * \return A copy of this terminal expression.
     */
    template<typename... TheExpr>
    expression_type
    operator()( TheExpr... ) const
    {
        return *this;
    }

    /// Return the polynomial order advertised to the expression system.
    constexpr uint16_type polynomialOrder() const
    {
        return 0;
    }

    /// SB9 operators are geometry-dependent and are not marked polynomial.
    constexpr bool isPolynomial() const
    {
        return false;
    }

    /**
     * \brief Apply symbolic substitutions.
     *
     * SB9 terminal operators do not contain symbolic sub-expressions, so the
     * operation returns the expression unchanged.
     */
    template <typename SymbolsExprType>
    expression_type applySymbolsExpr( SymbolsExprType const& ) const
    {
        return *this;
    }

    /**
     * \brief Symbolic derivative of an SB9 terminal expression.
     *
     * The operator is linear in the basis function and independent of symbolic
     * scalar symbols, so the symbolic derivative is the zero symmetric vector.
     */
    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&,
               TheSymbolExprType const& ) const
    {
        return vec( cst( value_type( 0 ) ),
                    cst( value_type( 0 ) ),
                    cst( value_type( 0 ) ),
                    cst( value_type( 0 ) ),
                    cst( value_type( 0 ) ),
                    cst( value_type( 0 ) ) );
    }

    /// Return the underlying basis proxy element.
    element_type const& element() const
    {
        return M_element;
    }

    /**
     * \brief Tensor evaluator for an SB9 vector contribution.
     *
     * The tensor owns the coefficient matrix for one geometric context and
     * exposes scalar/vector evaluation functions expected by Feel++ assembly.
     *
     * \tparam Geo_t Geometric mapping context type.
     * \tparam Basis_i_t Test basis tensor type.
     * \tparam Basis_j_t Trial basis tensor type.
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor : public ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>
    {
        /// Base class providing shell geometry cache updates.
        using base_type = ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>;
        /// Parent expression type.
        using expression_type = typename this_type::expression_type;
        /// Geometric mapping context type.
        using gmc_type = typename base_type::gmc_type;
        /// Scalar value type.
        using value_type = typename expression_type::value_type;
        /// Feel++ tensor shape: vector with symmetric-storage length.
        using shape = Shape<storage_size, Vectorial, false, false>;
        /// Vector type returned by vector-valued evaluation overloads.
        using vector_type = Eigen::Matrix<value_type, storage_size, 1>;
        /// Matrix storing one coefficient vector per element dof.
        using coeff_matrix_type = Eigen::Matrix<value_type, storage_size, Eigen::Dynamic>;
        /// Contribution-specific kernel cache for the current geometry data.
        using kernel_cache_type = KernelCacheTemplate<typename base_type::geometry_data_type>;

        /// Zero-trait required by the Feel++ tensor API.
        struct is_zero
        {
            /// SB9 vector contribution tensors are generally nonzero.
            static inline const bool value = false;
        };

        /**
         * \brief Construct from expression, geometry, test basis, and trial basis.
         */
        tensor( expression_type const& expr, Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
            :
            M_expr( expr ),
            M_fe( expr.element().functionSpace()->fe() )
        {
            this->checkFiniteElement();
            this->update( geom );
        }

        /**
         * \brief Construct from expression, geometry, and one basis tensor.
         */
        tensor( expression_type const& expr, Geo_t const& geom, Basis_i_t const& )
            :
            M_expr( expr ),
            M_fe( expr.element().functionSpace()->fe() )
        {
            this->checkFiniteElement();
            this->update( geom );
        }

        /**
         * \brief Construct from expression and geometry only.
         */
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_fe( expr.element().functionSpace()->fe() )
        {
            this->checkFiniteElement();
            this->update( geom );
        }

        /**
         * \brief Construct from an already-expanded expression.
         */
        template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                expression_type const& expr, Geo_t const& geom, TheArgsType const&... args )
            :
            tensor( expr, geom, args... )
        {
        }

        /// Update tensor data from geometry and test/trial basis tensors.
        void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
        {
            this->update( geom );
        }

        /// Update tensor data from geometry and one basis tensor.
        void update( Geo_t const& geom, Basis_i_t const& )
        {
            this->update( geom );
        }

        /// Update geometry data and recompute the coefficient matrix.
        void update( Geo_t const& geom )
        {
            this->updateFromGeom( geom );
            this->computeCoefficients();
        }

        /// Update tensor data from an already-expanded expression.
        template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                     Geo_t const& geom, TheArgsType const&... )
        {
            this->update( geom );
        }

        /// Update from a Feel++ assembly context and recompute coefficients.
        template<typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            this->updateFromContext( ctx... );
            this->computeCoefficients();
        }

        /**
         * \brief Evaluate one scalar component for a test/trial dof pair.
         *
         * \param i Local test dof index.
         * \param j Local trial dof index.
         * \param c1 Symmetric-storage component index.
         * \return Coefficient value for the role-selected dof.
         */
        value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type, uint16_type ) const
        {
            return M_coeff( c1, this->localDofId( i, j ) );
        }

        /// Pattern-context overload forwarding to \ref evalijq.
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return this->evalijq( i, j, c1, c2, q );
        }

        /// Evaluate one scalar component for a single local dof.
        value_type evaliq( uint16_type i, uint16_type c1, uint16_type, uint16_type ) const
        {
            return M_coeff( c1, i );
        }

        /// Evaluate one scalar component for context-independent access.
        value_type evalq( uint16_type c1, uint16_type, uint16_type ) const
        {
            return M_coeff( c1, 0 );
        }

        /// Evaluate the full symmetric-storage vector for a test/trial dof pair.
        Eigen::Map<const vector_type> evalijq( uint16_type i, uint16_type j, uint16_type ) const
        {
            return Eigen::Map<const vector_type>( M_coeff.col( this->localDofId( i, j ) ).data() );
        }

        /// Evaluate the full symmetric-storage vector for a single local dof.
        Eigen::Map<const vector_type> evaliq( uint16_type i, uint16_type ) const
        {
            return Eigen::Map<const vector_type>( M_coeff.col( i ).data() );
        }

        /// Evaluate the first coefficient vector for context-independent access.
        Eigen::Map<const vector_type> evalq( uint16_type ) const
        {
            return Eigen::Map<const vector_type>( M_coeff.col( 0 ).data() );
        }

    private:
        /// Validate that the finite element matches the current SB9 assumptions.
        void checkFiniteElement() const
        {
            CHECK( functionspace_type::nDim == 3 && functionspace_type::nRealDim == 3 )
                << "SB9 shell operators require a 3D displacement space";
            CHECK( M_fe )
                << "SB9 shell operators require a valid reference element";
            CHECK( M_fe->is_vectorial )
                << "SB9 shell operators require a vector-valued displacement basis";
            CHECK( M_fe->nComponents == 3 )
                << "SB9 shell operators require a 3-component displacement basis";
            CHECK( M_fe->nLocalDof == 8 )
                << "SB9 shell operators currently support Q1 eight-node hexahedra only";
        }

        /// Select the active local dof according to the expression role.
        uint16_type localDofId( uint16_type i, uint16_type j ) const
        {
            if constexpr ( Role == __TRIAL )
                return j;
            else
                return i;
        }

        /// Build the coefficient matrix for all vector components and nodes.
        void computeCoefficients()
        {
            uint16_type const nScalarLocalDof = M_fe->nLocalDof;
            uint16_type const nElementDof = M_fe->nComponents * nScalarLocalDof;
            uint16_type const nComponents = M_fe->nComponents;
            kernel_cache_type const cache( this->M_data );

            CHECK( nScalarLocalDof == kernel_cache_type::node_count )
                << "SB9 shell operators currently support 8 hexahedral vertices only";
            CHECK( nComponents == kernel_cache_type::component_count )
                << "SB9 shell operators require a 3-component vector basis";

            M_coeff.resize( storage_size, nElementDof );
            M_coeff.setZero();

            for ( uint16_type component = 0; component < nComponents; ++component )
            {
                for ( uint16_type node = 0; node < nScalarLocalDof; ++node )
                {
                    vector_type coeff = vector_type::Zero();
                    cache.template fillVectorCoefficients<Kind>( coeff, node, component );
                    M_coeff.col( component*nScalarLocalDof + node ) = coeff;
                }
            }
        }

        /// Parent expression instance.
        expression_type const& M_expr;
        /// Reference finite element associated with the basis proxy.
        std::shared_ptr<fe_type> M_fe;
        /// Coefficient matrix with one column per element dof.
        coeff_matrix_type M_coeff;
    };

private:
    /// Basis proxy element stored by the terminal expression.
    element_type M_element;
};

/**
 * \brief Generic Feel++ expression node for SB9 stabilization contributions.
 *
 * Stabilization contributions can have one, two, or three rows rather than the
 * six-component symmetric-storage vector used by the main SB9 strain terms.
 * This operator shares the same basis-role and geometry-update mechanics as
 * \ref SB9VectorOperator, but delegates coefficient generation to
 * `fillStabilizationCoefficients()`.
 *
 * \tparam ElementType Trial or test basis proxy element type.
 * \tparam Role Feel++ basis role, either __TEST or __TRIAL.
 * \tparam KernelCacheTemplate Contribution-specific cache template.
 * \tparam Kind Stabilization selector understood by the cache.
 * \tparam Rows Number of rows in the stabilization contribution.
 */
template <typename ElementType, OperatorType Role, template <typename> typename KernelCacheTemplate, auto Kind, int Rows>
class SB9StabilizationOperator
{
public:
    static_assert( Role == __TEST || Role == __TRIAL,
                   "SB9StabilizationOperator requires a trial or test basis role" );
    static_assert( Rows >= 1, "SB9StabilizationOperator requires at least one row" );

    /// Expression context flags requested from the geometric mapping context.
    static const size_type context = 0;
    /// Marks this expression as a terminal node in the Feel++ expression tree.
    static inline const bool is_terminal = true;

    /// Basis proxy element type carried by the expression.
    using element_type = ElementType;
    /// Function space associated with \ref element_type.
    using functionspace_type = typename element_type::functionspace_type;
    /// Reference finite element type.
    using fe_type = typename functionspace_type::reference_element_type;
    /// Current expression type.
    using this_type = SB9StabilizationOperator<ElementType, Role, KernelCacheTemplate, Kind, Rows>;
    /// Expression type returned by Feel++ expression machinery.
    using expression_type = this_type;
    /// Scalar value type of the function space.
    using value_type = typename functionspace_type::value_type;
    /// Runtime value type produced by the expression.
    using evaluate_type = Eigen::Matrix<value_type, Rows, 1>;

    /// Trait reporting whether this expression contains the requested test basis.
    template<typename Func>
    struct HasTestFunction
    {
        /// True when \c Func is the expression test basis finite element.
        static inline const bool result = ( Role == __TEST ) && boost::is_same<Func, fe_type>::value;
    };

    /// Trait reporting whether this expression contains the requested trial basis.
    template<typename Func>
    struct HasTrialFunction
    {
        /// True when \c Func is the expression trial basis finite element.
        static inline const bool result = ( Role == __TRIAL ) && boost::is_same<Func, fe_type>::value;
    };

    /// Compile-time test-basis flag used by Feel++ form assembly.
    template<typename Func>
    static inline const bool has_test_basis = ( Role == __TEST ) && boost::is_same<Func, fe_type>::value;

    /// Compile-time trial-basis flag used by Feel++ form assembly.
    template<typename Func>
    static inline const bool has_trial_basis = ( Role == __TRIAL ) && boost::is_same<Func, fe_type>::value;

    /// Test basis type exposed to the expression system.
    using test_basis = std::conditional_t<Role == __TEST, fe_type, std::nullptr_t>;
    /// Trial basis type exposed to the expression system.
    using trial_basis = std::conditional_t<Role == __TRIAL, fe_type, std::nullptr_t>;

    /// Lambda rebinding hook required by the Feel++ expression system.
    template<typename... TheExpr>
    struct Lambda
    {
        /// Rebound expression type.
        using type = expression_type;
    };

    /**
     * \brief Construct the SB9 stabilization operator for a basis proxy.
     * \param element Basis proxy element.
     */
    explicit SB9StabilizationOperator( element_type element )
        :
        M_element( std::move( element ) )
    {
    }

    /**
     * \brief Return this expression for Feel++ expression expansion.
     * \return A copy of this terminal expression.
     */
    template<typename... TheExpr>
    expression_type
    operator()( TheExpr... ) const
    {
        return *this;
    }

    /// Return the polynomial order advertised to the expression system.
    constexpr uint16_type polynomialOrder() const
    {
        return 0;
    }

    /// SB9 stabilization operators are geometry-dependent and non-polynomial.
    constexpr bool isPolynomial() const
    {
        return false;
    }

    /**
     * \brief Apply symbolic substitutions.
     *
     * Stabilization terminal operators do not contain symbolic sub-expressions,
     * so the operation returns the expression unchanged.
     */
    template <typename SymbolsExprType>
    expression_type applySymbolsExpr( SymbolsExprType const& ) const
    {
        return *this;
    }

    /**
     * \brief Symbolic derivative of an SB9 stabilization terminal expression.
     *
     * The operator is independent of symbolic scalar symbols, so the symbolic
     * derivative is zero with a shape matching the stabilization row count.
     */
    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&,
               TheSymbolExprType const& ) const
    {
        if constexpr ( Rows == 1 )
            return cst( value_type( 0 ) );
        else
            return vec( cst( value_type( 0 ) ) );
    }

    /// Return the underlying basis proxy element.
    element_type const& element() const
    {
        return M_element;
    }

    /**
     * \brief Tensor evaluator for an SB9 stabilization contribution.
     *
     * \tparam Geo_t Geometric mapping context type.
     * \tparam Basis_i_t Test basis tensor type.
     * \tparam Basis_j_t Trial basis tensor type.
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor : public ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>
    {
        /// Base class providing shell geometry cache updates.
        using base_type = ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>;
        /// Parent expression type.
        using expression_type = typename this_type::expression_type;
        /// Scalar value type.
        using value_type = typename expression_type::value_type;
        /// Feel++ tensor shape: vector with \c Rows entries.
        using shape = Shape<Rows, Vectorial, false, false>;
        /// Vector type returned by vector-valued evaluation overloads.
        using vector_type = Eigen::Matrix<value_type, Rows, 1>;
        /// Matrix storing one coefficient vector per element dof.
        using coeff_matrix_type = Eigen::Matrix<value_type, Rows, Eigen::Dynamic>;
        /// Contribution-specific kernel cache for the current geometry data.
        using kernel_cache_type = KernelCacheTemplate<typename base_type::geometry_data_type>;

        /// Zero-trait required by the Feel++ tensor API.
        struct is_zero
        {
            /// SB9 stabilization tensors are generally nonzero.
            static inline const bool value = false;
        };

        /**
         * \brief Construct from expression, geometry, test basis, and trial basis.
         */
        tensor( expression_type const& expr, Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
            :
            M_expr( expr ),
            M_fe( expr.element().functionSpace()->fe() )
        {
            this->checkFiniteElement();
            this->update( geom );
        }

        /**
         * \brief Construct from expression, geometry, and one basis tensor.
         */
        tensor( expression_type const& expr, Geo_t const& geom, Basis_i_t const& )
            :
            M_expr( expr ),
            M_fe( expr.element().functionSpace()->fe() )
        {
            this->checkFiniteElement();
            this->update( geom );
        }

        /**
         * \brief Construct from expression and geometry only.
         */
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_fe( expr.element().functionSpace()->fe() )
        {
            this->checkFiniteElement();
            this->update( geom );
        }

        /**
         * \brief Construct from an already-expanded expression.
         */
        template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                expression_type const& expr, Geo_t const& geom, TheArgsType const&... args )
            :
            tensor( expr, geom, args... )
        {
        }

        /// Update tensor data from geometry and test/trial basis tensors.
        void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
        {
            this->update( geom );
        }

        /// Update tensor data from geometry and one basis tensor.
        void update( Geo_t const& geom, Basis_i_t const& )
        {
            this->update( geom );
        }

        /// Update geometry data and recompute the coefficient matrix.
        void update( Geo_t const& geom )
        {
            this->updateFromGeom( geom );
            this->computeCoefficients();
        }

        /// Update tensor data from an already-expanded expression.
        template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                     Geo_t const& geom, TheArgsType const&... )
        {
            this->update( geom );
        }

        /// Update from a Feel++ assembly context and recompute coefficients.
        template<typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            this->updateFromContext( ctx... );
            this->computeCoefficients();
        }

        /// Evaluate one scalar component for a test/trial dof pair.
        value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type, uint16_type ) const
        {
            return M_coeff( c1, this->localDofId( i, j ) );
        }

        /// Pattern-context overload forwarding to \ref evalijq.
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return this->evalijq( i, j, c1, c2, q );
        }

        /// Evaluate one scalar component for a single local dof.
        value_type evaliq( uint16_type i, uint16_type c1, uint16_type, uint16_type ) const
        {
            return M_coeff( c1, i );
        }

        /// Evaluate one scalar component for context-independent access.
        value_type evalq( uint16_type c1, uint16_type, uint16_type ) const
        {
            return M_coeff( c1, 0 );
        }

        /// Evaluate the full stabilization vector for a test/trial dof pair.
        Eigen::Map<const vector_type> evalijq( uint16_type i, uint16_type j, uint16_type ) const
        {
            return Eigen::Map<const vector_type>( M_coeff.col( this->localDofId( i, j ) ).data() );
        }

        /// Evaluate the full stabilization vector for a single local dof.
        Eigen::Map<const vector_type> evaliq( uint16_type i, uint16_type ) const
        {
            return Eigen::Map<const vector_type>( M_coeff.col( i ).data() );
        }

        /// Evaluate the first coefficient vector for context-independent access.
        Eigen::Map<const vector_type> evalq( uint16_type ) const
        {
            return Eigen::Map<const vector_type>( M_coeff.col( 0 ).data() );
        }

    private:
        /// Validate that the finite element matches the current SB9 assumptions.
        void checkFiniteElement() const
        {
            CHECK( functionspace_type::nDim == 3 && functionspace_type::nRealDim == 3 )
                << "SB9 stabilization operators require a 3D displacement space";
            CHECK( M_fe )
                << "SB9 stabilization operators require a valid reference element";
            CHECK( M_fe->is_vectorial )
                << "SB9 stabilization operators require a vector-valued displacement basis";
            CHECK( M_fe->nComponents == 3 )
                << "SB9 stabilization operators require a 3-component displacement basis";
            CHECK( M_fe->nLocalDof == 8 )
                << "SB9 stabilization operators currently support Q1 eight-node hexahedra only";
        }

        /// Select the active local dof according to the expression role.
        uint16_type localDofId( uint16_type i, uint16_type j ) const
        {
            if constexpr ( Role == __TRIAL )
                return j;
            else
                return i;
        }

        /// Build the stabilization coefficient matrix for all components/nodes.
        void computeCoefficients()
        {
            uint16_type const nScalarLocalDof = M_fe->nLocalDof;
            uint16_type const nElementDof = M_fe->nComponents * nScalarLocalDof;
            uint16_type const nComponents = M_fe->nComponents;
            kernel_cache_type const cache( this->M_data );

            CHECK( nScalarLocalDof == kernel_cache_type::node_count )
                << "SB9 stabilization operators currently support 8 hexahedral vertices only";
            CHECK( nComponents == kernel_cache_type::component_count )
                << "SB9 stabilization operators require a 3-component vector basis";

            M_coeff.resize( Rows, nElementDof );
            M_coeff.setZero();

            for ( uint16_type component = 0; component < nComponents; ++component )
            {
                for ( uint16_type node = 0; node < nScalarLocalDof; ++node )
                {
                    vector_type coeff = vector_type::Zero();
                    cache.template fillStabilizationCoefficients<Kind>( coeff, node, component );
                    M_coeff.col( component*nScalarLocalDof + node ) = coeff;
                }
            }
        }

        /// Parent expression instance.
        expression_type const& M_expr;
        /// Reference finite element associated with the basis proxy.
        std::shared_ptr<fe_type> M_fe;
        /// Coefficient matrix with one column per element dof.
        coeff_matrix_type M_coeff;
    };

private:
    /// Basis proxy element stored by the terminal expression.
    element_type M_element;
};
} // namespace detail
} // namespace vf
} // namespace Feel

#endif
