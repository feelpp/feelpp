/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file shellgeometric.hpp
   \brief Geometric DSL operators for shell-oriented formulations on hexahedra
 */
#ifndef FEELPP_VF_SHELLGEOMETRIC_HPP
#define FEELPP_VF_SHELLGEOMETRIC_HPP 1

#include <array>
#include <limits>

#include <Eigen/LU>

#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/ones.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
template <typename GMCType>
inline constexpr bool shell_cell_supported_v =
    GMCType::PDim == 3 &&
    GMCType::NDim == 3 &&
    GMCType::element_type::is_hypercube &&
    !GMCType::is_on_face;

template <typename GMCType>
struct ShellCellGeometryData
{
    using value_type = typename GMCType::value_type;
    using vector_type = Eigen::Matrix<value_type, GMCType::NDim, 1>;
    using matrix_type = Eigen::Matrix<value_type, GMCType::NDim, GMCType::NDim>;
    using matrix2_type = Eigen::Matrix<value_type, 2, 2>;
    using node_matrix_type = Eigen::Matrix<value_type, 8, GMCType::NDim>;
    using hallquist_vector_type = Eigen::Matrix<value_type, 8, 1>;
    using gamma_matrix_type = Eigen::Matrix<value_type, 8, 4>;

    value_type area0 = value_type( 0 );
    value_type thickness = value_type( 0 );
    vector_type normal = vector_type::Zero();
    matrix_type frame = matrix_type::Identity();
    matrix_type covariantBasis0 = matrix_type::Identity();
    matrix_type contravariantBasis0 = matrix_type::Identity();
    matrix_type metric0 = matrix_type::Identity();
    matrix_type jacobian0 = matrix_type::Identity();
    value_type invJ0_00 = value_type( 0 );
    value_type invJ0_01 = value_type( 0 );
    value_type invJ0_02 = value_type( 0 );
    value_type invJ0_10 = value_type( 0 );
    value_type invJ0_11 = value_type( 0 );
    value_type invJ0_12 = value_type( 0 );
    value_type invJ0_20 = value_type( 0 );
    value_type invJ0_21 = value_type( 0 );
    value_type invJ0_22 = value_type( 0 );
    node_matrix_type localNodes = node_matrix_type::Zero();
    hallquist_vector_type bx = hallquist_vector_type::Zero();
    hallquist_vector_type by = hallquist_vector_type::Zero();
    hallquist_vector_type bz = hallquist_vector_type::Zero();
    gamma_matrix_type vgamma = gamma_matrix_type::Zero();
    matrix2_type Ja = matrix2_type::Zero();
    matrix2_type Jb = matrix2_type::Zero();
    matrix2_type Jc = matrix2_type::Zero();
    matrix2_type Jd = matrix2_type::Zero();
    matrix_type invJa = matrix_type::Identity();
    matrix_type invJb = matrix_type::Identity();
    matrix_type invJc = matrix_type::Identity();
    matrix_type invJd = matrix_type::Identity();
};

template <typename Matrix2Type, typename Matrix3Type, typename ValueType>
void
fillReducedInverse( Matrix2Type const& J, ValueType h, Matrix3Type& invJ )
{
    ValueType const detJ = J( 0, 0 ) * J( 1, 1 ) - J( 1, 0 ) * J( 0, 1 );
    ValueType const eps = 100 * std::numeric_limits<ValueType>::epsilon();

    CHECK( math::abs( detJ ) > eps )
        << "SB9 reduced Jacobian is singular";
    CHECK( math::abs( h ) > eps )
        << "SB9 reduced Jacobian requires a non-zero shell thickness";

    invJ.setZero();
    invJ( 0, 0 ) =  J( 1, 1 ) / detJ;
    invJ( 0, 1 ) = -J( 0, 1 ) / detJ;
    invJ( 1, 0 ) = -J( 1, 0 ) / detJ;
    invJ( 1, 1 ) =  J( 0, 0 ) / detJ;
    invJ( 2, 2 ) =  ValueType( 2 ) / h;
}

template <typename GMCType>
ShellCellGeometryData<GMCType>
computeShellCellGeometry( GMCType const* gmc )
{
    using value_type = typename GMCType::value_type;
    using data_type = ShellCellGeometryData<GMCType>;

    data_type data;

    CHECK( gmc ) << "invalid geometric mapping context";

    if constexpr ( shell_cell_supported_v<GMCType> )
    {
        using matrix_node_type = typename GMCType::gm_type::matrix_node_t_type;

        matrix_node_type xRefs( GMCType::PDim, 1 );
        xRefs.clear();

        auto centerPc = gmc->geometricMapping()->preCompute( xRefs );
        GMCType centerCtx( gmc->geometricMapping(),
                           gmc->element(),
                           std::integral_constant<int, vm::KB | vm::JACOBIAN>(),
                           centerPc );

        auto const& K = centerCtx.K( 0 );
        auto const xXi = K.col( 0 );
        auto const xEta = K.col( 1 );
        auto const xZeta = K.col( 2 );

        auto crossXiEta = xXi.cross( xEta );
        value_type const crossNorm = crossXiEta.norm();
        value_type const eps = 100 * std::numeric_limits<value_type>::epsilon();

        CHECK( crossNorm > eps )
            << "shell geometric operators require a non-degenerate midsurface Jacobian";

        data.normal = crossXiEta / crossNorm;

        // Align the local shell normal with the positive reference zeta direction.
        if ( data.normal.dot( xZeta ) < value_type( 0 ) )
            data.normal *= value_type( -1 );

        // h_K = V_K / A0_K with A0_K = 4 || X_,xi(0) x X_,eta(0) || on [-1,1]^2.
        data.area0 = value_type( 4 ) * crossNorm;
        data.thickness = gmc->element().measure() / data.area0;
        CHECK( math::abs( data.thickness ) > eps )
            << "shellThickness() requires a non-zero shell thickness";

        value_type const xiNorm = xXi.norm();
        CHECK( xiNorm > eps )
            << "shellFrame() requires a non-degenerate xi tangent at the midsurface";

        auto const t1 = xXi / xiNorm;
        auto t2 = data.normal.cross( t1 );
        value_type const t2Norm = t2.norm();
        CHECK( t2Norm > eps )
            << "shellFrame() requires a non-degenerate eta tangent at the midsurface";
        t2 /= t2Norm;

        data.frame.col( 0 ) = t1;
        data.frame.col( 1 ) = t2;
        data.frame.col( 2 ) = data.normal;
        data.covariantBasis0.col( 0 ) = xXi;
        data.covariantBasis0.col( 1 ) = xEta;
        data.covariantBasis0.col( 2 ) = xZeta;
        data.metric0 = data.covariantBasis0.transpose() * data.covariantBasis0;

        auto const luCovariant = data.covariantBasis0.fullPivLu();
        CHECK( luCovariant.isInvertible() )
            << "shellCovariantBasis0() requires an invertible center Jacobian";
        data.contravariantBasis0 = luCovariant.inverse().transpose();
        data.jacobian0 = data.frame.transpose() * data.covariantBasis0;

        Eigen::Matrix<value_type, GMCType::NDim, 1> center = Eigen::Matrix<value_type, GMCType::NDim, 1>::Zero();
        for ( uint16_type pointId = 0; pointId < gmc->element().nVertices(); ++pointId )
        {
            Eigen::Matrix<value_type, GMCType::NDim, 1> node = Eigen::Matrix<value_type, GMCType::NDim, 1>::Zero();
            for ( uint16_type d = 0; d < GMCType::NDim; ++d )
                node( d ) = gmc->element().point( pointId ).node()[d];
            center += node;
        }
        center /= value_type( gmc->element().nVertices() );

        for ( uint16_type pointId = 0; pointId < gmc->element().nVertices(); ++pointId )
        {
            Eigen::Matrix<value_type, GMCType::NDim, 1> node = Eigen::Matrix<value_type, GMCType::NDim, 1>::Zero();
            for ( uint16_type d = 0; d < GMCType::NDim; ++d )
                node( d ) = gmc->element().point( pointId ).node()[d];
            auto const localNode = data.frame.transpose() * ( node - center );
            data.localNodes.row( pointId ) = localNode.transpose();
        }

        auto const luJacobian0 = data.jacobian0.fullPivLu();
        CHECK( luJacobian0.isInvertible() )
            << "shellJacobian0() requires an invertible local center Jacobian";
        auto const invJacobian0 = luJacobian0.inverse();
        auto const invMatJ0 = invJacobian0.transpose();
        data.invJ0_00 = invMatJ0(0,0);
        data.invJ0_01 = invMatJ0(0,1);
        data.invJ0_02 = invMatJ0(0,2);
        data.invJ0_10 = invMatJ0(1,0);
        data.invJ0_11 = invMatJ0(1,1);
        data.invJ0_12 = invMatJ0(1,2);
        data.invJ0_20 = invMatJ0(2,0);
        data.invJ0_21 = invMatJ0(2,1);
        data.invJ0_22 = invMatJ0(2,2);

        Eigen::Matrix<value_type, 3, 8> BKsi = Eigen::Matrix<value_type, 3, 8>::Zero();
        value_type const u = value_type( 1 ) / value_type( 8 );
        BKsi.row( 0 ) << -u,  u,  u, -u, -u,  u,  u, -u;
        BKsi.row( 1 ) << -u, -u,  u,  u, -u, -u,  u,  u;
        BKsi.row( 2 ) << -u, -u, -u, -u,  u,  u,  u,  u;
        data.bx = ( invMatJ0.row( 0 ) * BKsi ).transpose();
        data.by = ( invMatJ0.row( 1 ) * BKsi ).transpose();
        data.bz = ( invMatJ0.row( 2 ) * BKsi ).transpose();

        Eigen::Matrix<value_type, 8, 4> hallquistH = Eigen::Matrix<value_type, 8, 4>::Zero();
        hallquistH <<
             1.0,  1.0,  1.0, -1.0,
             1.0, -1.0, -1.0,  1.0,
            -1.0, -1.0,  1.0, -1.0,
            -1.0,  1.0, -1.0,  1.0,
            -1.0, -1.0,  1.0,  1.0,
            -1.0,  1.0, -1.0, -1.0,
             1.0,  1.0,  1.0,  1.0,
             1.0, -1.0, -1.0, -1.0;

        for ( int gammaId = 0; gammaId < 4; ++gammaId )
        {
            auto const hVec = data.localNodes.transpose() * hallquistH.col( gammaId );
            for ( int nodeId = 0; nodeId < 8; ++nodeId )
            {
                data.vgamma( nodeId, gammaId ) =
                    ( hallquistH( nodeId, gammaId ) -
                      hVec( 0 ) * data.bx( nodeId ) -
                      hVec( 1 ) * data.by( nodeId ) -
                      hVec( 2 ) * data.bz( nodeId ) ) / value_type( 8 );
            }
        }

        auto const x1 = data.localNodes( 0, 0 );
        auto const y1 = data.localNodes( 0, 1 );
        auto const x2 = data.localNodes( 1, 0 );
        auto const y2 = data.localNodes( 1, 1 );
        auto const x3 = data.localNodes( 2, 0 );
        auto const y3 = data.localNodes( 2, 1 );
        auto const x4 = data.localNodes( 3, 0 );
        auto const y4 = data.localNodes( 3, 1 );
        auto const x5 = data.localNodes( 4, 0 );
        auto const y5 = data.localNodes( 4, 1 );
        auto const x6 = data.localNodes( 5, 0 );
        auto const y6 = data.localNodes( 5, 1 );
        auto const x7 = data.localNodes( 6, 0 );
        auto const y7 = data.localNodes( 6, 1 );
        auto const x8 = data.localNodes( 7, 0 );
        auto const y8 = data.localNodes( 7, 1 );

        data.Ja( 0, 0 ) = (  x2 - x1 + x6 - x5 ) / value_type( 4 );
        data.Ja( 0, 1 ) = (  y2 - y1 + y6 - y5 ) / value_type( 4 );
        data.Ja( 1, 0 ) = ( -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8 ) / value_type( 8 );
        data.Ja( 1, 1 ) = ( -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8 ) / value_type( 8 );

        data.Jb( 0, 0 ) = ( -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8 ) / value_type( 8 );
        data.Jb( 0, 1 ) = ( -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8 ) / value_type( 8 );
        data.Jb( 1, 0 ) = ( -x2 + x3 - x6 + x7 ) / value_type( 4 );
        data.Jb( 1, 1 ) = ( -y2 + y3 - y6 + y7 ) / value_type( 4 );

        data.Jc( 0, 0 ) = (  x3 - x4 + x7 - x8 ) / value_type( 4 );
        data.Jc( 0, 1 ) = (  y3 - y4 + y7 - y8 ) / value_type( 4 );
        data.Jc( 1, 0 ) = ( -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8 ) / value_type( 8 );
        data.Jc( 1, 1 ) = ( -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8 ) / value_type( 8 );

        data.Jd( 0, 0 ) = ( -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8 ) / value_type( 8 );
        data.Jd( 0, 1 ) = ( -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8 ) / value_type( 8 );
        data.Jd( 1, 0 ) = ( -x1 + x4 - x5 + x8 ) / value_type( 4 );
        data.Jd( 1, 1 ) = ( -y1 + y4 - y5 + y8 ) / value_type( 4 );

        fillReducedInverse( data.Ja, data.thickness, data.invJa );
        fillReducedInverse( data.Jb, data.thickness, data.invJb );
        fillReducedInverse( data.Jc, data.thickness, data.invJc );
        fillReducedInverse( data.Jd, data.thickness, data.invJd );
    }
    else
    {
        CHECK( false )
            << "shell geometric operators are defined only on 3D cell integrals "
            << "over hypercube/hexahedron elements";
    }

    return data;
}
} // namespace detail

template <uint16_type Component>
class ReferenceCoordinate
{
public:
    static const size_type context = 0;
    static inline const bool is_terminal = true;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    static inline const bool has_test_basis = false;
    template<typename Func>
    static inline const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    using this_type = ReferenceCoordinate<Component>;
    using expression_type = ReferenceCoordinate<Component>;
    using value_type = double;

    template<typename... TheExpr>
    struct Lambda
    {
        using type = expression_type;
    };

    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... ) const
    {
        return *this;
    }

    constexpr uint16_type polynomialOrder() const
    {
        return 1;
    }

    constexpr bool isPolynomial() const
    {
        return true;
    }

    template <typename SymbolsExprType>
    this_type applySymbolsExpr( SymbolsExprType const& ) const
    {
        return *this;
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&,
               TheSymbolExprType const& ) const
    {
        return cst( value_type( 0 ) );
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expression_type = ReferenceCoordinate<Component>;
        using key_type = key_t<Geo_t>;
        using gmc_ptrtype = gmc_ptr_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        using value_type = typename expression_type::value_type;
        using shape = Shape<gmc_type::NDim, Scalar, false, false>;

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( expression_type const&, Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
        }

        tensor( expression_type const&, Geo_t const& geom, Basis_i_t const& )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
        }

        tensor( expression_type const&, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
        }

        template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                expression_type const& expr, Geo_t const& geom, TheArgsType const&... args )
            :
            tensor( expr, geom, args... )
        {
        }

        void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
        }

        void update( Geo_t const& geom, Basis_i_t const& )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
        }

        void update( Geo_t const& geom )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
        }

        template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                     Geo_t const& geom, TheArgsType const&... )
        {
            this->update( geom );
        }

        template<typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            boost::fusion::vector<CTX...> ctxvec( ctx... );
            M_gmc = boost::fusion::at_c<0>( ctxvec )->gmContext().get();
        }

        value_type evalijq( uint16_type, uint16_type, uint16_type, uint16_type, uint16_type q ) const
        {
            return this->evalq( 0, 0, q );
        }

        template<int PatternContext>
        value_type evalijq( uint16_type, uint16_type, uint16_type, uint16_type, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return this->evalq( 0, 0, q );
        }

        value_type evaliq( uint16_type, uint16_type, uint16_type, uint16_type q ) const
        {
            return this->evalq( 0, 0, q );
        }

        value_type evalq( uint16_type, uint16_type, uint16_type q ) const
        {
            CHECK( M_gmc ) << "invalid geometric mapping context";

            if constexpr ( detail::shell_cell_supported_v<gmc_type> )
            {
                return M_gmc->xRef( q )[Component];
            }
            else
            {
                CHECK( false )
                    << "reference shell coordinates xi(), eta(), zeta() are defined only on "
                    << "3D cell integrals over hypercube/hexahedron elements";
                return value_type( 0 );
            }
        }

        gmc_ptrtype M_gmc = nullptr;
    };
};

using Xi = ReferenceCoordinate<0>;
using Eta = ReferenceCoordinate<1>;
using Zeta = ReferenceCoordinate<2>;

namespace detail
{
template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
struct ShellCellGeometryTensorBase
{
    using key_type = key_t<Geo_t>;
    using gmc_ptrtype = gmc_ptr_t<Geo_t>;
    using gmc_type = gmc_t<Geo_t>;
    using geometry_data_type = ShellCellGeometryData<gmc_type>;

    void updateFromGeom( Geo_t const& geom )
    {
        gmc_ptrtype gmc = fusion::at_key<key_type>( geom ).get();
        M_data = computeShellCellGeometry( gmc );
    }

    template<typename... CTX>
    void updateFromContext( CTX const&... ctx )
    {
        boost::fusion::vector<CTX...> ctxvec( ctx... );
        auto gmc = boost::fusion::at_c<0>( ctxvec )->gmContext().get();
        M_data = computeShellCellGeometry( gmc );
    }

    geometry_data_type M_data;
};
} // namespace detail

namespace detail
{
template <typename Derived>
class ShellCellTerminalBase
{
public:
    static const size_type context = 0;
    static inline const bool is_terminal = true;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    static inline const bool has_test_basis = false;
    template<typename Func>
    static inline const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    using this_type = Derived;
    using expression_type = Derived;
    using value_type = double;

    template<typename... TheExpr>
    struct Lambda
    {
        using type = expression_type;
    };

    template<typename... TheExpr>
    expression_type
    operator()( TheExpr... ) const
    {
        return static_cast<expression_type const&>( *this );
    }

    constexpr uint16_type polynomialOrder() const
    {
        return 0;
    }

    constexpr bool isPolynomial() const
    {
        return false;
    }

    template <typename SymbolsExprType>
    expression_type applySymbolsExpr( SymbolsExprType const& ) const
    {
        return {};
    }
};

template <typename TensorType>
struct ShellThicknessAccessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.thickness ); }
};

template <typename TensorType>
struct ShellArea0Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.area0 ); }
};

template <typename TensorType>
struct ShellNormalAccessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.normal ); }
};

template <typename TensorType>
struct ShellFrameAccessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.frame ); }
};

template <typename TensorType>
struct ShellCovariantBasis0Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.covariantBasis0 ); }
};

template <typename TensorType>
struct ShellContravariantBasis0Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.contravariantBasis0 ); }
};

template <typename TensorType>
struct ShellMetric0Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.metric0 ); }
};

template <typename TensorType>
struct ShellJacobian0Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.jacobian0 ); }
};

template <typename TensorType>
struct ShellInvJ0_00Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_00 ); }
};
template <typename TensorType>
struct ShellInvJ0_01Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_01 ); }
};
template <typename TensorType>
struct ShellInvJ0_02Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_02 ); }
};
template <typename TensorType>
struct ShellInvJ0_10Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_10 ); }
};
template <typename TensorType>
struct ShellInvJ0_11Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_11 ); }
};
template <typename TensorType>
struct ShellInvJ0_12Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_12 ); }
};
template <typename TensorType>
struct ShellInvJ0_20Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_20 ); }
};
template <typename TensorType>
struct ShellInvJ0_21Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_21 ); }
};
template <typename TensorType>
struct ShellInvJ0_22Accessor
{
    static decltype(auto) get( TensorType const& data ) { return ( data.invJ0_22 ); }
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename Derived, template <typename> typename AccessorT>
struct ShellCellScalarTensor : public ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>
{
    using base_type = ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>;
    using expression_type = Derived;
    using gmc_type = typename base_type::gmc_type;
    using value_type = typename expression_type::value_type;
    using shape = Shape<gmc_type::NDim, Scalar, false, false>;

    struct is_zero
    {
        static inline const bool value = false;
    };

    ShellCellScalarTensor( expression_type const&, Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
    {
        this->update( geom );
    }

    ShellCellScalarTensor( expression_type const&, Geo_t const& geom, Basis_i_t const& )
    {
        this->update( geom );
    }

    ShellCellScalarTensor( expression_type const&, Geo_t const& geom )
    {
        this->update( geom );
    }

    template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
    ShellCellScalarTensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                           expression_type const& expr, Geo_t const& geom, TheArgsType const&... args )
        :
        ShellCellScalarTensor( expr, geom, args... )
    {
    }

    void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
    {
        this->update( geom );
    }

    void update( Geo_t const& geom, Basis_i_t const& )
    {
        this->update( geom );
    }

    void update( Geo_t const& geom )
    {
        this->updateFromGeom( geom );
    }

    template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
    void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                 Geo_t const& geom, TheArgsType const&... )
    {
        this->update( geom );
    }

    template<typename... CTX>
    void updateContext( CTX const&... ctx )
    {
        this->updateFromContext( ctx... );
    }

    value_type evalijq( uint16_type, uint16_type, uint16_type, uint16_type, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data );
    }

    template<int PatternContext>
    value_type evalijq( uint16_type, uint16_type, uint16_type, uint16_type, uint16_type,
                        mpl::int_<PatternContext> ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data );
    }

    value_type evaliq( uint16_type, uint16_type, uint16_type, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data );
    }

    value_type evalq( uint16_type, uint16_type, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data );
    }
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename Derived, template <typename> typename AccessorT>
struct ShellCellVectorTensor : public ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>
{
    using base_type = ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>;
    using expression_type = Derived;
    using gmc_type = typename base_type::gmc_type;
    using value_type = typename expression_type::value_type;
    using shape = Shape<gmc_type::NDim, Vectorial, false, false>;
    using vector_type = Eigen::Matrix<value_type, gmc_type::NDim, 1>;

    struct is_zero
    {
        static inline const bool value = false;
    };

    ShellCellVectorTensor( expression_type const&, Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
    {
        this->update( geom );
    }

    ShellCellVectorTensor( expression_type const&, Geo_t const& geom, Basis_i_t const& )
    {
        this->update( geom );
    }

    ShellCellVectorTensor( expression_type const&, Geo_t const& geom )
    {
        this->update( geom );
    }

    template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
    ShellCellVectorTensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                           expression_type const& expr, Geo_t const& geom, TheArgsType const&... args )
        :
        ShellCellVectorTensor( expr, geom, args... )
    {
    }

    void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
    {
        this->update( geom );
    }

    void update( Geo_t const& geom, Basis_i_t const& )
    {
        this->update( geom );
    }

    void update( Geo_t const& geom )
    {
        this->updateFromGeom( geom );
    }

    template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
    void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                 Geo_t const& geom, TheArgsType const&... )
    {
        this->update( geom );
    }

    template<typename... CTX>
    void updateContext( CTX const&... ctx )
    {
        this->updateFromContext( ctx... );
    }

    value_type evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1 );
    }

    Eigen::Map<const vector_type> evalijq( uint16_type, uint16_type, uint16_type ) const
    {
        return Eigen::Map<const vector_type>( AccessorT<typename base_type::geometry_data_type>::get( this->M_data ).data() );
    }

    template<int PatternContext>
    value_type evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type, uint16_type,
                        mpl::int_<PatternContext> ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1 );
    }

    value_type evaliq( uint16_type, uint16_type c1, uint16_type, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1 );
    }

    Eigen::Map<const vector_type> evaliq( uint16_type, uint16_type ) const
    {
        return Eigen::Map<const vector_type>( AccessorT<typename base_type::geometry_data_type>::get( this->M_data ).data() );
    }

    value_type evalq( uint16_type c1, uint16_type, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1 );
    }

    Eigen::Map<const vector_type> evalq( uint16_type ) const
    {
        return Eigen::Map<const vector_type>( AccessorT<typename base_type::geometry_data_type>::get( this->M_data ).data() );
    }
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename Derived, template <typename> typename AccessorT>
struct ShellCellMatrixTensor : public ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>
{
    using base_type = ShellCellGeometryTensorBase<Geo_t, Basis_i_t, Basis_j_t>;
    using expression_type = Derived;
    using gmc_type = typename base_type::gmc_type;
    using value_type = typename expression_type::value_type;
    using shape = Shape<gmc_type::NDim, Tensor2, false, false>;
    using matrix_type = Eigen::Matrix<value_type, gmc_type::NDim, gmc_type::NDim>;

    struct is_zero
    {
        static inline const bool value = false;
    };

    ShellCellMatrixTensor( expression_type const&, Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
    {
        this->update( geom );
    }

    ShellCellMatrixTensor( expression_type const&, Geo_t const& geom, Basis_i_t const& )
    {
        this->update( geom );
    }

    ShellCellMatrixTensor( expression_type const&, Geo_t const& geom )
    {
        this->update( geom );
    }

    template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
    ShellCellMatrixTensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                           expression_type const& expr, Geo_t const& geom, TheArgsType const&... args )
        :
        ShellCellMatrixTensor( expr, geom, args... )
    {
    }

    void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
    {
        this->update( geom );
    }

    void update( Geo_t const& geom, Basis_i_t const& )
    {
        this->update( geom );
    }

    void update( Geo_t const& geom )
    {
        this->updateFromGeom( geom );
    }

    template<typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
    void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&,
                 Geo_t const& geom, TheArgsType const&... )
    {
        this->update( geom );
    }

    template<typename... CTX>
    void updateContext( CTX const&... ctx )
    {
        this->updateFromContext( ctx... );
    }

    value_type evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type c2, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1, c2 );
    }

    template<int PatternContext>
    value_type evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type c2, uint16_type,
                        mpl::int_<PatternContext> ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1, c2 );
    }

    value_type evaliq( uint16_type, uint16_type c1, uint16_type c2, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1, c2 );
    }

    value_type evalq( uint16_type c1, uint16_type c2, uint16_type ) const
    {
        return AccessorT<typename base_type::geometry_data_type>::get( this->M_data )( c1, c2 );
    }

    Eigen::Map<const matrix_type> evalijq( uint16_type, uint16_type, uint16_type ) const
    {
        return Eigen::Map<const matrix_type>( AccessorT<typename base_type::geometry_data_type>::get( this->M_data ).data() );
    }

    Eigen::Map<const matrix_type> evaliq( uint16_type, uint16_type ) const
    {
        return Eigen::Map<const matrix_type>( AccessorT<typename base_type::geometry_data_type>::get( this->M_data ).data() );
    }

    Eigen::Map<const matrix_type> evalq( uint16_type ) const
    {
        return Eigen::Map<const matrix_type>( AccessorT<typename base_type::geometry_data_type>::get( this->M_data ).data() );
    }
};

template <typename Derived, template <typename> typename AccessorT>
class ShellCellScalarTerminal : public ShellCellTerminalBase<Derived>
{
public:
    using value_type = typename ShellCellTerminalBase<Derived>::value_type;

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&,
               TheSymbolExprType const& ) const
    {
        return cst( value_type( 0 ) );
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    using tensor = ShellCellScalarTensor<Geo_t, Basis_i_t, Basis_j_t, Derived, AccessorT>;
};

template <typename Derived, template <typename> typename AccessorT>
class ShellCellVectorTerminal : public ShellCellTerminalBase<Derived>
{
public:
    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&,
               TheSymbolExprType const& ) const
    {
        return vector_zero();
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    using tensor = ShellCellVectorTensor<Geo_t, Basis_i_t, Basis_j_t, Derived, AccessorT>;
};

template <typename Derived, template <typename> typename AccessorT>
class ShellCellMatrixTerminal : public ShellCellTerminalBase<Derived>
{
public:
    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&,
               TheSymbolExprType const& ) const
    {
        return zero<3,3>();
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    using tensor = ShellCellMatrixTensor<Geo_t, Basis_i_t, Basis_j_t, Derived, AccessorT>;
};
} // namespace detail

class ShellThickness : public detail::ShellCellScalarTerminal<ShellThickness, detail::ShellThicknessAccessor>
{
};

class ShellArea0 : public detail::ShellCellScalarTerminal<ShellArea0, detail::ShellArea0Accessor>
{
};

class ShellCovariantBasis0 : public detail::ShellCellMatrixTerminal<ShellCovariantBasis0, detail::ShellCovariantBasis0Accessor>
{
};

class ShellContravariantBasis0 : public detail::ShellCellMatrixTerminal<ShellContravariantBasis0, detail::ShellContravariantBasis0Accessor>
{
};

class ShellMetric0 : public detail::ShellCellMatrixTerminal<ShellMetric0, detail::ShellMetric0Accessor>
{
};

class ShellJacobian0 : public detail::ShellCellMatrixTerminal<ShellJacobian0, detail::ShellJacobian0Accessor>
{
};

class ShellInvJ0_00 : public detail::ShellCellScalarTerminal<ShellInvJ0_00, detail::ShellInvJ0_00Accessor>
{
};
class ShellInvJ0_01 : public detail::ShellCellScalarTerminal<ShellInvJ0_01, detail::ShellInvJ0_01Accessor>
{
};
class ShellInvJ0_02 : public detail::ShellCellScalarTerminal<ShellInvJ0_02, detail::ShellInvJ0_02Accessor>
{
};
class ShellInvJ0_10 : public detail::ShellCellScalarTerminal<ShellInvJ0_10, detail::ShellInvJ0_10Accessor>
{
};
class ShellInvJ0_11 : public detail::ShellCellScalarTerminal<ShellInvJ0_11, detail::ShellInvJ0_11Accessor>
{
};
class ShellInvJ0_12 : public detail::ShellCellScalarTerminal<ShellInvJ0_12, detail::ShellInvJ0_12Accessor>
{
};
class ShellInvJ0_20 : public detail::ShellCellScalarTerminal<ShellInvJ0_20, detail::ShellInvJ0_20Accessor>
{
};
class ShellInvJ0_21 : public detail::ShellCellScalarTerminal<ShellInvJ0_21, detail::ShellInvJ0_21Accessor>
{
};
class ShellInvJ0_22 : public detail::ShellCellScalarTerminal<ShellInvJ0_22, detail::ShellInvJ0_22Accessor>
{
};

class ShellNormal : public detail::ShellCellVectorTerminal<ShellNormal, detail::ShellNormalAccessor>
{
};

class ShellFrame : public detail::ShellCellMatrixTerminal<ShellFrame, detail::ShellFrameAccessor>
{
};

inline
Expr<Xi>
xi()
{
    return Expr<Xi>( Xi() );
}

inline
Expr<Eta>
eta()
{
    return Expr<Eta>( Eta() );
}

inline
Expr<Zeta>
zeta()
{
    return Expr<Zeta>( Zeta() );
}

inline
Expr<ShellThickness>
shellThickness()
{
    return Expr<ShellThickness>( ShellThickness() );
}

inline
Expr<ShellArea0>
shellArea0()
{
    return Expr<ShellArea0>( ShellArea0() );
}

inline
Expr<ShellNormal>
shellNormal()
{
    return Expr<ShellNormal>( ShellNormal() );
}

inline
Expr<ShellCovariantBasis0>
shellCovariantBasis0()
{
    return Expr<ShellCovariantBasis0>( ShellCovariantBasis0() );
}

inline
Expr<ShellContravariantBasis0>
shellContravariantBasis0()
{
    return Expr<ShellContravariantBasis0>( ShellContravariantBasis0() );
}

inline
Expr<ShellMetric0>
shellMetric0()
{
    return Expr<ShellMetric0>( ShellMetric0() );
}

inline
Expr<ShellJacobian0>
shellJacobian0()
{
    return Expr<ShellJacobian0>( ShellJacobian0() );
}

inline
Expr<ShellInvJ0_00>
shellInvJ0_00()
{
    return Expr<ShellInvJ0_00>( ShellInvJ0_00() );
}
inline
Expr<ShellInvJ0_01>
shellInvJ0_01()
{
    return Expr<ShellInvJ0_01>( ShellInvJ0_01() );
}
inline
Expr<ShellInvJ0_02>
shellInvJ0_02()
{
    return Expr<ShellInvJ0_02>( ShellInvJ0_02() );
}
inline
Expr<ShellInvJ0_10>
shellInvJ0_10()
{
    return Expr<ShellInvJ0_10>( ShellInvJ0_10() );
}
inline
Expr<ShellInvJ0_11>
shellInvJ0_11()
{
    return Expr<ShellInvJ0_11>( ShellInvJ0_11() );
}
inline
Expr<ShellInvJ0_12>
shellInvJ0_12()
{
    return Expr<ShellInvJ0_12>( ShellInvJ0_12() );
}
inline
Expr<ShellInvJ0_20>
shellInvJ0_20()
{
    return Expr<ShellInvJ0_20>( ShellInvJ0_20() );
}
inline
Expr<ShellInvJ0_21>
shellInvJ0_21()
{
    return Expr<ShellInvJ0_21>( ShellInvJ0_21() );
}
inline
Expr<ShellInvJ0_22>
shellInvJ0_22()
{
    return Expr<ShellInvJ0_22>( ShellInvJ0_22() );
}


inline
Expr<ShellFrame>
shellFrame()
{
    return Expr<ShellFrame>( ShellFrame() );
}

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_SHELLGEOMETRIC_HPP */
