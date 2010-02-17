/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-02

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file val.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-02
 */
#if !defined( LIFE_VF_VAL_HPP )
#define LIFE_VF_VAL_HPP 1

# include <boost/preprocessor/comparison/less.hpp>
# include <boost/preprocessor/logical/and.hpp>
# include <boost/preprocessor/control/if.hpp>
# include <boost/preprocessor/list/at.hpp>
# include <boost/preprocessor/list/cat.hpp>
# include <boost/preprocessor/list/for_each_product.hpp>
# include <boost/preprocessor/logical/or.hpp>
# include <boost/preprocessor/tuple/to_list.hpp>
# include <boost/preprocessor/tuple/eat.hpp>
# include <boost/preprocessor/facilities/empty.hpp>
# include <boost/preprocessor/punctuation/comma.hpp>
# include <boost/preprocessor/facilities/identity.hpp>

/// \cond detail
#include <life/lifecore/traits.hpp>
#include <life/lifevf/unaryfunctor.hpp>


#if defined( HAVE_QD_H ) && defined(HAVE_MPFR)
# define VF_CHECK_ARITHMETIC_TYPE()                                     \
    BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value || \
                          ::boost::is_same<value_1_type, std::complex<float> >::value || \
                          ::boost::is_same<value_1_type, std::complex<double> >::value || \
                          ::boost::is_same<value_1_type,mp_type>::value || \
                          ::boost::is_same<value_1_type,dd_real>::value || \
                          ::boost::is_same<value_1_type,qd_real>::value) ); \
    /**/
#elif defined( HAVE_QD_H )
# define VF_CHECK_ARITHMETIC_TYPE()                                     \
    BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value || \
                          ::boost::is_same<value_1_type, std::complex<float> >::value || \
                          ::boost::is_same<value_1_type, std::complex<double> >::value || \
                          ::boost::is_same<value_1_type,dd_real>::value || \
                          ::boost::is_same<value_1_type,qd_real>::value) ); \
    /**/
#elif defined( HAVE_MPFR )
# define VF_CHECK_ARITHMETIC_TYPE()                                     \
    BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value || \
                          ::boost::is_same<value_1_type, std::complex<float> >::value || \
                          ::boost::is_same<value_1_type, std::complex<double> >::value || \
                          ::boost::is_same<value_1_type,mp_type>::value) ); \
    /**/
#else
# define VF_CHECK_ARITHMETIC_TYPE()                                     \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<value_1_type>::value || \
                           ::boost::is_same<value_1_type, std::complex<float> >::value || \
                           ::boost::is_same<value_1_type, std::complex<double> >::value ) \
                         );                                             \
    /**/
#endif

namespace Life
{
namespace vf
{

template < typename ExprT1 >
class Val
    :
        public UnaryFunctor<typename ExprT1::value_type>
{
public:

    static const size_type context = ExprT1::context;

    static const uint16_type imorder = ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    typedef UnaryFunctor<typename ExprT1::value_type> super;
    typedef typename super::functordomain_type functordomain_type;
    typedef typename super::functordomain_ptrtype functordomain_ptrtype;
    typedef ExprT1 expression_1_type;
    typedef Val<ExprT1> this_type;
    typedef typename expression_1_type::value_type value_1_type;
    typedef value_1_type value_type;

    VF_CHECK_ARITHMETIC_TYPE()

    explicit Val( expression_1_type const& __expr1  )
    :
    super( "value", functordomain_ptrtype(new UnboundedDomain<value_type>() )),
        _M_expr_1( __expr1 )
    {
        Debug( 5051 ) << "Val::Val default constructorn";
    }

    Val( Val const& __vfp  )
        :
        super( "value", functordomain_ptrtype(new UnboundedDomain<value_type>() )),
        _M_expr_1( __vfp._M_expr_1 )
    {
        Debug( 5051 ) << "Val::Val copy constructorn";
    }

    bool isSymetric() const { return false; }

    void eval( int nx, value_type const* x, value_type* f ) const
    {
        for( int i = 0; i < nx; ++i )
            f[i] = x[i];
    }

    expression_1_type const& expression() const { return _M_expr_1; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t,Basis_j_t> tensor2_expr_type;
        typedef typename tensor2_expr_type::value_type value_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                                  mpl::identity<detail::gmc<0> >,
                                  mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef typename tensor2_expr_type::shape shape;

        struct is_zero { static const bool value = tensor2_expr_type::is_zero::value; };

        template<typename ExprT>
        tensor( ExprT const& expr, Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            _M_expr( expr.expression(), geom ),
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_loc( boost::extents[shape::M][shape::N][_M_gmc->nPoints()] )
        {
            update( geom );
        }
        template<typename ExprT>
        tensor( ExprT const& expr,Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            _M_expr( expr.expression(), geom ),
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_loc( boost::extents[shape::M][shape::N][_M_gmc->nPoints()] )
        {
            update( geom );
        }
        template<typename ExprT>
        tensor( ExprT const& expr, Geo_t const& geom )
            :
            _M_expr( expr.expression(), geom ),
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_loc( boost::extents[shape::M][shape::N][_M_gmc->nPoints()] )
        {
            update( geom );
        }
        template<typename IM>
        void init( IM const& im )
        {
            _M_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom )
        {
            _M_expr.update( geom );
            for( int c1 = 0; c1 < shape::M; ++c1 )
                for( int c2 = 0; c2 < shape::N; ++c2 )
                    for ( int q = 0; q < _M_gmc->nPoints(); ++q )
                        {
                            _M_loc[c1][c2][q] = _M_expr.evalq( c1, c2, q );
                        }
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            _M_expr.update( geom, face );
            for( int c1 = 0; c1 < shape::M; ++c1 )
                for( int c2 = 0; c2 < shape::N; ++c2 )
                    for ( int q = 0; q < _M_gmc->nPoints(); ++q )
                        {
                            _M_loc[c1][c2][q] = _M_expr.evalq( c1, c2, q );
                        }
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalq( c1, c2, q );
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q, mpl::int_<shape::rank>() );
        }
    private:
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const
        {
            Life::detail::ignore_unused_variable_warning(c1);
            Life::detail::ignore_unused_variable_warning(c2);
            return _M_loc[0][0][q];
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const
        {
            if ( shape::M > shape::N )
                return _M_loc[c1][0][q];
            return _M_loc[0][c2][q];
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const
        {
            return _M_loc[c1][c2][q];
        }
    private:
        tensor2_expr_type _M_expr;
        gmc_ptrtype _M_gmc;
        boost::multi_array<value_type,3> _M_loc;
    };

protected:
    Val() {}

    expression_1_type _M_expr_1;
};
/// \endcond

/**
 * \brief precompute expression tensor
 *
 * This allows for more efficient  bi/linear form assembly
 */
template<typename ExprT1>
inline
Expr< Val<typename mpl::if_<boost::is_arithmetic<ExprT1>,
                            mpl::identity<Cst<ExprT1> >,
                            mpl::identity<ExprT1> >::type::type > >
val( ExprT1 const& __e1 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
        mpl::identity<Cst<ExprT1> >,
        mpl::identity<ExprT1> >::type::type t1;
    typedef Val<t1> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ) ) );
}
} // vf
} // Life
#endif /* LIFE_VF_VAL_HPP */

