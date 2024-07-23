/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-04-11

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file ones.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-04-11
 */
#ifndef FEELPP_VF_ONES_H
#define FEELPP_VF_ONES_H 1

// #include <blitz/array.h>
#include <boost/multi_array.hpp>

namespace Feel
{
namespace vf
{
/// \cond DETAIL
namespace detail
{
/**
 * \class Ones
 * \brief Return a matrix expression or N-dimensional array whose elements are all 1
 *
 *If you need to create a matrix whose values are all the same, you
 * should use an expression like
 *
 *\code
 * val * ones<n,m>()
 *\endcode
 *
 * @author Christophe Prud'homme
 */
template <int M, int N>
class Ones
{
  public:
    /** @name Typedefs
     */
    //@{
    static const size_type context = 0;

    static const bool is_terminal = true;

    template <typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    template <typename Func>
    static const bool has_test_basis = false;
    template <typename Func>
    static const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    typedef Ones<M, N> this_type;
    typedef double value_type;
    using evaluate_type = Eigen::Matrix<value_type, M, N>;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    template <typename EigenMatrix>
    Ones( EigenMatrix const& m )
        : M_values( m )
    {
    }

    Ones( Ones const& eig )
        : M_values( eig.M_values )
    {
    }
    ~Ones()
    {
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //! polynomial order
    constexpr uint16_type polynomialOrder() const { return 0; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    // blitz::Array<value_type,2> ones() const { return M_values; }
    Eigen::Matrix<double, M, N> const& ones() const
    {
        return M_values;
    }
    constexpr evaluate_type evaluate( bool ) const
    {
        return M_values;
    }

    template <typename SymbolsExprType>
    this_type applySymbolsExpr( SymbolsExprType const& se ) const
    {
        return *this;
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        return vf::detail::Ones<M, N>( Eigen::Matrix<value_type, M, N>::Zero() );
    }

    //@}
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_type::value_type value_type;
        typedef value_type return_value_type;
        using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t, key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t, key_type>::type::element_type gmc_type;
        
        struct INVALID_SHAPE
        {
        };
        using shape = mp11::mp_if_c<
            ( M == 1 && N == 1 ),
            Shape<gmc_type::nDim, Scalar, false, false>,
            mp11::mp_if_c<
                ( M == gmc_type::nDim && N == 1 ),
                Shape<gmc_type::nDim, Vectorial, false, false>,
                mp11::mp_if_c<
                    ( M == 1 && N == gmc_type::nDim ),
                    Shape<gmc_type::nDim, Vectorial, true, false>,
                    mp11::mp_if_c<
                        ( M == gmc_type::nDim && N == gmc_type::nDim ),
                        Shape<gmc_type::nDim, Tensor2, false, false>,
                        INVALID_SHAPE>>>>;

        using matrix_type = Eigen::Matrix<value_type, shape::M, shape::N>;
        template <class Args>
        struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr, Geo_t const&, Basis_i_t const&, Basis_j_t const& )
            : M_expr( expr ),
              M_values( expr.ones() )
        {
            // std::cout << "tensor::ones = " << M_expr.ones() << "\n";
        }

        tensor( this_type const& expr, Geo_t const&, Basis_i_t const& )
            : M_expr( expr ),
              M_values( expr.ones() )
        {
        }

        tensor( this_type const& expr, Geo_t const& )
            : M_expr( expr ),
              M_values( expr.ones() )
        {
        }
        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                expression_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            : tensor( expr, geom, theInitArgs... )
        {
        }

        void update( Geo_t const&, Basis_i_t const&, Basis_j_t const& )
        {
        }
        void update( Geo_t const&, Basis_i_t const& )
        {
        }
        void update( Geo_t const& )
        {
        }
        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            return evalq( c1, c2, q );
        }

        template <int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            return evalq( c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( q );
            if constexpr ( shape::rank == 0 )
                return M_values( 0, 0 );
            else if constexpr ( shape::rank == 1 )
            {
                if constexpr ( shape::is_transposed )
                    return M_values( 0, c2 );
                else
                    return M_values( c1, 0 );
            }
            else
                return M_values( c1, c2 );
        }
        FEELPP_STRONG_INLINE Eigen::Map<const matrix_type>
        evalijq( uint16_type i, uint16_type j, uint16_type /*q*/ ) const
        {
            return Eigen::Map<const matrix_type>( M_values.data() );
        }
        FEELPP_STRONG_INLINE Eigen::Map<const matrix_type>
        evaliq( uint16_type i, uint16_type /*q*/ ) const
        {
            return Eigen::Map<const matrix_type>( M_values.data() );
        }
        FEELPP_STRONG_INLINE Eigen::Map<const matrix_type>
        evalq( uint16_type /*q*/ ) const
        {
            return Eigen::Map<const matrix_type>( M_values.data() );
        }

      private:
        this_type M_expr;
        Eigen::Matrix<double, M, N> M_values;
    };

  private:
    Eigen::Matrix<double, M, N> M_values;
};
} // namespace detail
/// \endcond

/**
 *
 * \brief Return a matrix expression or N-dimensional array whose elements are all 1
 *
 *If you need to create a matrix whose values are all the same, you
 * should use an expression like
 *
 *\code
 * val * ones<n,m>()
 *\endcode
 *
 * @author Christophe
 */
template <int M, int N = M>
inline Expr<vf::detail::Ones<M, N>>
ones()
{
    return Expr<vf::detail::Ones<M, N>>( vf::detail::Ones<M, N>( Eigen::Matrix<double, M, N>::Ones() ) );
}

template <int M, int N = M>
inline Expr<vf::detail::Ones<M, N>>
zero()
{
    return Expr<vf::detail::Ones<M, N>>( vf::detail::Ones<M, N>( Eigen::Matrix<double, M, N>::Zero() ) );
}

template <int M, int N = M>
inline Expr<vf::detail::Ones<M, N>>
eye()
{
    return Expr<vf::detail::Ones<M, N>>( vf::detail::Ones<M, N>( Eigen::Matrix<double, M, N>::Identity() ) );
}

template <int M, int N = M>
inline Expr<vf::detail::Ones<M, N>>
Id()
{
    return Expr<vf::detail::Ones<M, N>>( vf::detail::Ones<M, N>( Eigen::Matrix<double, M, N>::Identity() ) );
}

template <int M, int N = M>
inline Expr<vf::detail::Ones<M, N>>
constant( double value )
{
    return Expr<vf::detail::Ones<M, N>>( vf::detail::Ones<M, N>( Eigen::Matrix<double, M, N>::Constant( value ) ) );
}

} // namespace vf
} // namespace Feel
#endif /* FEELPP_VF_ONES_H */
