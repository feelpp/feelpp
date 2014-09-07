/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-20

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
   \file matvec.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-20
 */
#ifndef __VFVec_H
#define __VFVec_H 1

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/mpl/bitor.hpp>
#include <boost/mpl/bitwise.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/mpl/plus.hpp>
#include <boost/mpl/arithmetic.hpp>


namespace Feel
{
namespace fusion = boost::fusion;
namespace mpl = boost::mpl;
namespace vf
{
/// \cond detail
namespace detail
{
struct GetContext
{
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<GetContext( Lhs,Rhs )>
#if BOOST_VERSION < 104200
: boost::remove_reference<Lhs>
#else
: boost::remove_reference<Rhs>
#endif
    {

        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
#if BOOST_VERSION < 104200
        typedef typename mpl::bitor_<mpl::size_t<lhs_noref_type::context>, rhs_noref_type >::type type;
#else
        typedef typename mpl::bitor_<mpl::size_t<rhs_noref_type::context>, lhs_noref_type >::type type;
#endif
    };
#if 0
    template<typename Expr1, typename Expr2>
    struct result
    {
        typedef typename mpl::bitor_<mpl::size_t<Expr1::context>, Expr2 >::type type;
    };
#endif
    template<typename Lhs, typename Rhs>
#if BOOST_VERSION < 104200
    Lhs
#else
    Rhs
#endif
    operator()( const Lhs& lhs, const Rhs& rhs ) const
    {
        return lhs | rhs;
    }
};
struct GetImOrder
{
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<GetImOrder( Lhs,Rhs )>
#if BOOST_VERSION < 104200
: boost::remove_reference<Lhs>
#else
: boost::remove_reference<Rhs>
#endif
    {

        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
#if BOOST_VERSION < 104200
        typedef typename boost::mpl::max< boost::mpl::size_t<lhs_noref_type::imorder>, rhs_noref_type >::type type;
#else
        typedef typename boost::mpl::max< boost::mpl::size_t<rhs_noref_type::imorder>, lhs_noref_type >::type type;
#endif
    };

    template<typename Lhs, typename Rhs>
#if BOOST_VERSION < 104200
    Lhs
#else
    Rhs
#endif
    operator()( const Lhs& lhs, const Rhs& rhs ) const
    {
    }
};
struct GetImIsPoly
{
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<GetImIsPoly( Lhs,Rhs )>
#if BOOST_VERSION < 104200
: boost::remove_reference<Lhs>
#else
: boost::remove_reference<Rhs>
#endif
    {

        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
#if BOOST_VERSION < 104200
        typedef typename boost::mpl::and_< boost::mpl::bool_<lhs_noref_type::imIsPoly>, rhs_noref_type >::type type;
#else
        typedef typename boost::mpl::and_< boost::mpl::bool_<rhs_noref_type::imIsPoly>, lhs_noref_type >::type type;
#endif
    };

    template<typename Lhs, typename Rhs>
#if BOOST_VERSION < 104200
    Lhs
#else
    Rhs
#endif
    operator()( const Lhs& lhs, const Rhs& rhs ) const
    {
    }
};
template<typename Func>
struct ExprHasTestFunction
{
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<ExprHasTestFunction( Lhs,Rhs )>
#if BOOST_VERSION < 104200
: boost::remove_reference<Lhs>
#else
: boost::remove_reference<Rhs>
#endif
    {

        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
#if BOOST_VERSION < 104200
        typedef typename mpl::or_<mpl::bool_<lhs_noref_type::template HasTestFunction<Func>::result>, rhs_noref_type >::type type;
#else
        typedef typename mpl::or_<mpl::bool_<rhs_noref_type::template HasTestFunction<Func>::result>, lhs_noref_type >::type type;
#endif
    };
    template<typename Lhs, typename Rhs>
#if BOOST_VERSION < 104200
    Lhs
#else
    Rhs
#endif
    operator()( const Lhs& lhs, const Rhs& rhs ) const
    {
#if 0
        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
        return ( lhs_noref_type::template ExprHasTestFunction<Func>::result ||
                 lhs_noref_type::template ExprHasTestFunction<Func>::result );
#endif
    }
};
template<typename Func>
struct ExprHasTrialFunction
{
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<ExprHasTrialFunction( Lhs,Rhs )>
#if BOOST_VERSION < 104200
: boost::remove_reference<Lhs>
#else
: boost::remove_reference<Rhs>
#endif
    {

        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
#if BOOST_VERSION < 104200
        typedef typename mpl::or_<mpl::bool_<lhs_noref_type::template HasTrialFunction<Func>::result>, rhs_noref_type >::type type;
#else
        typedef typename mpl::or_<mpl::bool_<rhs_noref_type::template HasTrialFunction<Func>::result>, lhs_noref_type >::type type;
#endif
    };
    template<typename Lhs, typename Rhs>
#if BOOST_VERSION < 104200
    Lhs
#else
    Rhs
#endif
    operator()( const Lhs& lhs, const Rhs& rhs ) const
    {
#if 0
        return ( Lhs::template HasTrialFunction<Func>::result ||
                 Rhs::template HasTrialFunction<Func>::result );
#endif
    }
};
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
struct initialize_expression_gij
{
    template<typename Sig>
    struct result;

    template<typename ExprT>
    struct result<initialize_expression_gij( ExprT )>
    {
        typedef typename boost::remove_reference<ExprT>::type::template tensor<Geo_t, Basis_i_t, Basis_j_t> type;
    };

    initialize_expression_gij( Geo_t const& geom , Basis_i_t const& fev, Basis_j_t const& feu )
        :
        M_geom( geom ),
        M_fev( fev ),
        M_feu( feu )
    {}

    template <typename ExprT>
    typename ExprT::template tensor<Geo_t, Basis_i_t, Basis_j_t>
    operator()( ExprT& expr ) const
    {
        return typename ExprT::template tensor<Geo_t, Basis_i_t, Basis_j_t>( expr, M_geom, M_fev, M_feu );
    }
    const Geo_t& M_geom;
    const Basis_i_t& M_fev;
    const Basis_j_t& M_feu;
};
template<typename Geo_t, typename Basis_i_t>
struct initialize_expression_gi
{
    template<typename Sig>
    struct result;

    template<typename ExprT>
    struct result<initialize_expression_gi( ExprT )>
    {
        typedef typename boost::remove_reference<ExprT>::type::template tensor<Geo_t, Basis_i_t> type;
    };
    initialize_expression_gi( Geo_t const& geom , Basis_i_t const& fev )
        :
        M_geom( geom ),
        M_fev( fev )
    {}

    template <typename ExprT>
    typename ExprT::template tensor<Geo_t, Basis_i_t>
    operator()( ExprT& expr ) const
    {
        return typename ExprT::template tensor<Geo_t, Basis_i_t>( expr, M_geom, M_fev );
    }

    const Geo_t& M_geom;
    const Basis_i_t& M_fev;
};

template<typename Geo_t>
struct initialize_expression_g
{
    template<typename Sig>
    struct result;

    template<typename ExprT>
    struct result<initialize_expression_g( ExprT )>
    {
        typedef typename boost::remove_reference<ExprT>::type::template tensor<Geo_t> type;
    };
    initialize_expression_g( Geo_t const& geom )
        :
        M_geom( geom )
    {}

    template <typename ExprT>
    typename ExprT::template tensor<Geo_t>
    operator()( ExprT& expr ) const
    {
        return typename ExprT::template tensor<Geo_t>( expr, M_geom );
        //return typename ExprT::template tensor<Geo_t>( expr, M_geom );
    }

    const Geo_t& M_geom;
};

template<typename CTX>
struct update_context
{
    update_context( CTX const& ctx )
        :
        M_ctx( ctx )
    {}

    template <typename ExprT>
    void operator()( ExprT& expr ) const
    {
        expr.updateContext( M_ctx );
    }

    const CTX & M_ctx;
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
struct update_expression_gij
{
    update_expression_gij( Geo_t const& geom , Basis_i_t const& fev, Basis_j_t const& feu )
        :
        M_geom( geom ),
        M_fev( fev ),
        M_feu( feu )
    {}

    template <typename ExprT>
    void operator()( ExprT& expr ) const
    {
        expr.update( M_geom, M_fev, M_feu );
    }
    const Geo_t& M_geom;
    const Basis_i_t& M_fev;
    const Basis_j_t& M_feu;
};
template<typename Geo_t, typename Basis_i_t>
struct update_expression_gi
{
    update_expression_gi( Geo_t const& geom , Basis_i_t const& fev )
        :
        M_geom( geom ),
        M_fev( fev )
    {}

    template <typename ExprT>
    void operator()( ExprT& expr ) const
    {
        expr.update( M_geom, M_fev );
    }
    const Geo_t& M_geom;
    const Basis_i_t& M_fev;
};

template<typename Geo_t>
struct update_expression_g
{
    update_expression_g( Geo_t const& geom )
        :
        M_geom( geom )
    {}

    template <typename ExprT>
    void operator()( ExprT& expr ) const
    {
        expr.update( M_geom );
    }
    const Geo_t& M_geom;
};


template<typename Geo_t>
struct update_expression_face_g
{
    update_expression_face_g( Geo_t const& geom, uint16_type face )
        :
        M_geom( geom ),
        M_face( face )
    {}

    template <typename ExprT>
    void operator()( ExprT& expr ) const
    {
        expr.update( M_geom, M_face );
    }
    const Geo_t& M_geom;
    const uint16_type M_face;
};

template<typename IM_t>
struct init_expression
{
    init_expression( IM_t const& im )
        :
        M_im( im )
    {}

    template <typename ExprT>
    void operator()( ExprT& expr ) const
    {
        expr.init( M_im );
    }
    const IM_t& M_im;
};

template<typename T>
struct evaluate_expression_gij
{
    template<typename Sig>
    struct result;

    template<typename ExprT, typename V>
#if BOOST_VERSION < 104200
    struct result<evaluate_expression_gij<T>( ExprT,V )>
#else
    struct result<evaluate_expression_gij<T>( V,ExprT )>
#endif
:
    boost::remove_reference<V>
    {};
    evaluate_expression_gij( uint16_type indexi , uint16_type indexj, uint16_type c1, uint16_type q )
        :
        M_indexi( indexi ),
        M_indexj( indexj ),
        M_c1( c1 ),
        M_c2( 0 ),
        M_q( q ),
        M_index( M_c1 ),
        M_current( 0 )
    {}

    evaluate_expression_gij( int n, uint16_type indexi, uint16_type indexj, uint16_type c1, uint16_type c2, uint16_type q )
        :
        M_indexi( indexi ),
        M_indexj( indexj ),
        M_c1( c1 ),
        M_c2( c2 ),
        M_q( q ),
        M_index( M_c1*n+M_c2 ),
        M_current( 0 )

    {}

    template <typename ExprT>
    T
#if BOOST_VERSION < 104200
    operator()( ExprT const& expr, T const& res ) const
#else
    operator()( T const& res, ExprT const& expr ) const
#endif
    {
        T ret = res;

        if ( M_current == M_index )
            ret = expr.evalijq( M_indexi, M_indexj, 0, 0, M_q );

        ++M_current;
        return ret;
    }
    const uint16_type& M_indexi;
    const uint16_type& M_indexj;
    const uint16_type M_c1;
    const uint16_type M_c2;
    const uint16_type M_q;
    const uint16_type M_index;
    mutable int M_current;

};
template<typename T=double>
struct evaluate_expression_gi
{
    template<typename Sig>
    struct result;

    template<typename ExprT,typename V>
#if BOOST_VERSION < 104200
    struct result<evaluate_expression_gi<T>( ExprT,V )>
#else
    struct result<evaluate_expression_gi<T>( V,ExprT )>
#endif
:
    boost::remove_reference<V>
    {};

    evaluate_expression_gi( uint16_type indexi, uint16_type c1, uint16_type q )
        :
        M_indexi( indexi ),
        M_c1( c1 ),
        M_c2( 0 ),
        M_q( q ),
        M_index( M_c1 ),
        M_current( 0 )
    {}

    evaluate_expression_gi( int n, uint16_type indexi, uint16_type c1, uint16_type c2, uint16_type q )
        :
        M_indexi( indexi ),
        M_c1( c1 ),
        M_c2( c2 ),
        M_q( q ),
        M_index( M_c1*n+M_c2 ),
        M_current( 0 )

    {}
    template <typename ExprT>
    T
#if BOOST_VERSION < 104200
    operator()( ExprT const& expr, T const& res ) const
#else
    operator()( T const& res, ExprT const& expr ) const
#endif
    {
        T ret = res;

        if ( M_current == M_index )
            ret = expr.evaliq( M_indexi, 0, 0, M_q );

        ++M_current;
        return ret;

    }
    const uint16_type& M_indexi;
    const uint16_type M_c1;
    const uint16_type M_c2;
    const uint16_type M_q;
    const uint16_type M_index;
    mutable int M_current;
};

template<typename T=double>
struct evaluate_expression_g
{
    template<typename Sig>
    struct result;

    template<typename ExprT,typename V>
#if BOOST_VERSION < 104200
    struct result<evaluate_expression_g<T>( ExprT,V )>
#else
    struct result<evaluate_expression_g<T>( V,ExprT )>
#endif
:
    boost::remove_reference<V>
    {};

    evaluate_expression_g( uint16_type c1, uint16_type q )
        :
        M_c1( c1 ),
        M_c2( 0 ),
        M_q( q ),
        M_index( M_c1 ),
        M_current( 0 )
    {}

    evaluate_expression_g( int n, uint16_type c1, uint16_type c2, uint16_type q )
        :
        M_c1( c1 ),
        M_c2( c2 ),
        M_q( q ),
        M_index( M_c1*n+M_c2 ),
        M_current( 0 )

    {}

    template <typename ExprT>
    T
#if BOOST_VERSION < 104200
    operator()( ExprT const& expr, T const& res ) const
#else
    operator()( T const& res, ExprT const& expr ) const
#endif
    {
        T ret = res;

        if ( M_current == M_index )
            ret = expr.evalq( 0, 0, M_q );

        ++M_current;
        return ret;
    }
    const uint16_type M_c1;
    const uint16_type M_c2;
    const uint16_type M_q;
    const uint16_type M_index;
    mutable int M_current;
};
/**
 * \class Vec
 * \brief class that represents a matrix in the language
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename VectorExpr>
class Vec
{
public:


    /** @name Typedefs
     */
    //@{
    static const size_type context = fusion::result_of::accumulate<VectorExpr,mpl::size_t<0>,GetContext>::type::value;
    //static const size_type context = fusion::accumulate( VectorExpr(),size_type(0),GetContext() );
    static const bool is_terminal = false;
    static const uint16_type imorder = fusion::result_of::accumulate<VectorExpr,mpl::size_t<0>,GetImOrder>::type::value;
    static const bool imIsPoly = fusion::result_of::accumulate<VectorExpr,mpl::bool_<true>,GetImIsPoly>::type::value;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = fusion::result_of::accumulate<VectorExpr,mpl::bool_<false>,ExprHasTestFunction<Func> >::type::value;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = fusion::result_of::accumulate<VectorExpr,mpl::bool_<false>,ExprHasTrialFunction<Func> >::type::value;
    };

    typedef VectorExpr expression_vector_type;
    typedef Vec<expression_vector_type> this_type;

    typedef double value_type;
    typedef value_type evaluate_type;

    static const uint16_type vector_size =  fusion::result_of::size<expression_vector_type>::type::value;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Vec( VectorExpr const& expr )
        :
        M_expr( expr )
    {}
    Vec( Vec const & expr )
        :
        M_expr( expr.M_expr )
    {}
    ~Vec()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    expression_vector_type const&  expression() const
    {
        return M_expr;
    }
    expression_vector_type      &  expression()
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}

    template<typename Geo_t, typename Basis_i_t = boost::none_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_type::expression_vector_type expression_vector_type;
        typedef Shape<expression_type::vector_size, Vectorial, false, false> shape;
        static const bool theshape = ( shape::M == expression_type::vector_size && shape::N == 1 );
        BOOST_MPL_ASSERT_MSG( theshape,
                              INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_1,
                              ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );

        typedef typename boost::remove_reference<typename fusion::result_of::at<expression_vector_type, boost::mpl::int_<0> >::type>::type::value_type value_type;
        //typedef typename expression_type::value_type value_type;

        template<typename ExprT>
        struct ExpressionToTensor
        {
            typedef typename boost::remove_reference<ExprT>::type expr_type;
            typedef typename expr_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> type;
        };
        typedef typename mpl::transform<expression_vector_type, ExpressionToTensor<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type tensor_vector_type;


        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_expr( fusion::transform( expr.expression(), initialize_expression_gij<Geo_t,Basis_i_t,Basis_j_t>( geom, fev, feu ) ) )
        {

            update( geom, fev, feu );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_expr( fusion::transform( expr.expression(), initialize_expression_gi<Geo_t,Basis_i_t>( geom, fev ) ) )
        {
            update( geom, fev );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            M_expr( fusion::transform( expr.expression(), initialize_expression_g<Geo_t>( geom ) ) )
        {
            update( geom );
        }
        template<typename IM>
        void init( IM const& im )
        {
            fusion::for_each( M_expr,vf::detail::init_expression<IM>( im ) );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_gij<Geo_t, Basis_i_t, Basis_j_t>( geom, fev, feu ) );

        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_gi<Geo_t, Basis_i_t>( geom, fev ) );
        }
        void update( Geo_t const& geom )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_g<Geo_t>( geom ) );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_face_g<Geo_t>( geom, face ) );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            fusion::for_each( M_expr,vf::detail::update_context<CTX>( ctx ) );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type /*c2*/, uint16_type q ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),vf::detail::evaluate_expression_gij<value_type>( i, j, c1, q  ) );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type /*c2*/, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),
                                      vf::detail::evaluate_expression_gij<value_type>( i, j, c1, q  ) );

        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type /*c2*/, uint16_type q ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),vf::detail::evaluate_expression_gi<value_type>( i, c1, q ) );
        }
        value_type
        evalq( uint16_type c1, uint16_type /*c2*/, uint16_type q ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),vf::detail::evaluate_expression_g<value_type>( c1, q ) );
        }
        tensor_vector_type M_expr;
    };



protected:

private:
    Vec();

    expression_vector_type M_expr;
};
} // detail
/// \endcond

/**
 * \brief vector definition
 */
template<typename Expr1>
inline
Expr<vf::detail::Vec<fusion::vector<Expr1> > >
vec( Expr1  expr1 )
{
    typedef vf::detail::Vec<fusion::vector<Expr1> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1>( expr1 ) ) );
}
/**
 * \brief vector definition
 */
template<typename Expr1, typename Expr2>
inline
Expr<vf::detail::Vec<fusion::vector<Expr1, Expr2> > >
vec( Expr1  expr1, Expr2  expr2 )
{
    typedef vf::detail::Vec<fusion::vector<Expr1, Expr2> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2>( expr1, expr2 ) ) );
}
/**
 * \brief vector definition
 */
template<typename Expr1, typename Expr2, typename Expr3>
inline
Expr<vf::detail::Vec<fusion::vector<Expr1, Expr2, Expr3> > >
vec( Expr1  expr1, Expr2  expr2, Expr3  expr3 )
{
    typedef vf::detail::Vec<fusion::vector<Expr1, Expr2, Expr3> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2, Expr3>( expr1, expr2, expr3 ) ) );
}

/// \cond detail
namespace detail
{

/**
 * \class Mat
 * \brief class that represents a matrix in the language
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int M, int N, typename MatrixExpr>
class Mat
{
public:


    /** @name Typedefs
     */
    //@{

    static const size_type context = fusion::result_of::accumulate<MatrixExpr,mpl::size_t<0>,GetContext>::type::value;
    //static const size_type context = fusion::accumulate( MatrixExpr(),size_type(0),GetContext() );
    static const bool is_terminal = false;
    static const uint16_type imorder = fusion::result_of::accumulate<MatrixExpr,mpl::size_t<0>,GetImOrder>::type::value;
    static const bool imIsPoly = fusion::result_of::accumulate<MatrixExpr,mpl::bool_<true>,GetImIsPoly>::type::value;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = fusion::result_of::accumulate<MatrixExpr,mpl::bool_<false>,ExprHasTestFunction<Func> >::type::value;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = fusion::result_of::accumulate<MatrixExpr,mpl::bool_<false>,ExprHasTrialFunction<Func> >::type::value;
    };

    typedef MatrixExpr expression_matrix_type;
    typedef Mat<M, N, expression_matrix_type> this_type;

    typedef double value_type;
    typedef value_type evaluate_type;

    static const uint16_type matrix_size1 = M;
    static const uint16_type matrix_size2 = N;
    static const uint16_type matrix_size  = M*N;

#if 0
    BOOST_MPL_ASSERT_MSG( ( M*N == fusion::result_of::size<expression_matrix_type>::type::value ),
                          ( INVALID_MATRIX_SIZE ),
                          ( mpl::int_<M>, mpl::int_<N>, mpl::int_<M*N>,
                            mpl::int_<fusion::result_of::size<expression_matrix_type>::type::value> ) );
#endif
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Mat( MatrixExpr const& expr )
        :
        M_expr( expr )
    {
    }
    Mat( Mat const & expr )
        :
        M_expr( expr.M_expr )
    {
    }
    ~Mat()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    expression_matrix_type const&  expression() const
    {
        return M_expr;
    }
    expression_matrix_type      &  expression()
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_type::expression_matrix_type expression_matrix_type;
        typedef Shape<expression_type::matrix_size1, Tensor2, false, false> shape;

#if 0
        static const bool theshape = ( shape::M == expression_type::matrix_size && shape::N == 1 );
        BOOST_MPL_ASSERT_MSG( theshape,
                              INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_1,
                              ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );
#endif

        typedef typename boost::remove_reference<typename fusion::result_of::at<expression_matrix_type, boost::mpl::int_<0> >::type>::type::value_type value_type;
        typedef typename boost::remove_reference<typename fusion::result_of::at<expression_matrix_type, boost::mpl::int_<0> >::type>::type::value_type first_value_type;
        //typedef fusion::result_of::accumulate<VectorExpr,first_value_type,GetValueType>::type value_type;
        //typedef typename expression_type::value_type value_type;

        template<typename ExprT>
        struct ExpressionToTensor
        {
            typedef typename boost::remove_reference<ExprT>::type expr_type;
            typedef typename expr_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> type;
        };
        typedef typename mpl::transform<expression_matrix_type, ExpressionToTensor<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type tensor_matrix_type;


        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_expr( fusion::transform( expr.expression(), initialize_expression_gij<Geo_t,Basis_i_t,Basis_j_t>( geom, fev, feu ) ) )
        {

            update( geom, fev, feu );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_expr( fusion::transform( expr.expression(), initialize_expression_gi<Geo_t,Basis_i_t>( geom, fev ) ) )
        {
            update( geom, fev );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            M_expr( fusion::transform( expr.expression(), initialize_expression_g<Geo_t>( geom ) ) )
        {
            update( geom );
        }
        template<typename IM>
        void init( IM const& im )
        {
            fusion::for_each( M_expr,vf::detail::init_expression<IM>( im ) );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_gij<Geo_t, Basis_i_t, Basis_j_t>( geom, fev, feu ) );

        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_gi<Geo_t, Basis_i_t>( geom, fev ) );
        }
        void update( Geo_t const& geom )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_g<Geo_t>( geom ) );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            fusion::for_each( M_expr,vf::detail::update_expression_face_g<Geo_t>( geom, face ) );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            fusion::for_each( M_expr,vf::detail::update_context<CTX>( ctx ) );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),
                                      vf::detail::evaluate_expression_gij<value_type>( expression_type::matrix_size2,
                                               i, j, c1, c2, q ) );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),
                                      vf::detail::evaluate_expression_gij<value_type>( expression_type::matrix_size2,
                                               i, j, c1, c2, q ) );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),
                                      vf::detail::evaluate_expression_gi<value_type>( expression_type::matrix_size2, i, c1, c2, q ) );

        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return fusion::accumulate( M_expr, value_type( 0 ),vf::detail::evaluate_expression_g<value_type>( expression_type::matrix_size2,
                                       c1, c2, q ) );
        }
        tensor_matrix_type M_expr;
    };



protected:

private:
    Mat();

    expression_matrix_type M_expr;
};
} // detail
/// \endcond

/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1> > >
mat( Expr1  expr1 )
{
    BOOST_MPL_ASSERT_MSG( ( M == 1 && N == 1 ),  INVALID_MATRIX_SIZE_SHOULD_BE_1_1, ( mpl::int_<M>, mpl::int_<N> ) );
    typedef vf::detail::Mat<M, N, fusion::vector<Expr1> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1>( expr1 ) ) );
}

/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1, typename Expr2>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2> > >
mat( Expr1  expr1, Expr2  expr2 )
{
    BOOST_MPL_ASSERT_MSG( ( M == 2 && N == 1 || M == 1 && N == 2 ), INVALID_MATRIX_SIZE, ( mpl::int_<M>, mpl::int_<N> ) );
    typedef vf::detail::Mat<M,N,fusion::vector<Expr1, Expr2> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2>( expr1, expr2 ) ) );
}

/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1, typename Expr2, typename Expr3>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3> > >
mat( Expr1  expr1, Expr2  expr2, Expr3  expr3 )
{
#if 0
    BOOST_MPL_ASSERT_MSG( ( M == 3 && N == 1 ||
                            M == 1 && N == 3 ),
                          ( INVALID_MATRIX_SIZE ),
                          ( mpl::int_<M>, mpl::int_<N> ) );
#endif
    typedef vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2, Expr3>( expr1, expr2, expr3 ) ) );
}

/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1, typename Expr2, typename Expr3, typename Expr4>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4> > >
mat( Expr1  expr1, Expr2  expr2, Expr3  expr3, Expr4 expr4 )
{
#if 0
    BOOST_MPL_ASSERT_MSG( ( M == 4 && N == 1 ||
                            M == 1 && N == 4 ||
                            M == 2 && N == 2 ),
                          ( INVALID_MATRIX_SIZE ),
                          ( mpl::int_<M>, mpl::int_<N> ) );
#endif
    typedef vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2, Expr3, Expr4>( expr1, expr2, expr3, expr4 ) ) );
}

/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1, typename Expr2, typename Expr3, typename Expr4, typename Expr5>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5> > >
mat( Expr1  expr1, Expr2  expr2, Expr3  expr3, Expr4 expr4, Expr5 expr5 )
{
#if 0
    BOOST_MPL_ASSERT_MSG( ( M == 4 && N == 1 ||
                            M == 1 && N == 4 ||
                            M == 2 && N == 2 ),
                          ( INVALID_MATRIX_SIZE ),
                          ( mpl::int_<M>, mpl::int_<N> ) );
#endif
    typedef vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5>( expr1, expr2, expr3, expr4, expr5 ) ) );
}

/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1, typename Expr2, typename Expr3, typename Expr4, typename Expr5, typename Expr6>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5, Expr6> > >
mat( Expr1  expr1, Expr2  expr2, Expr3  expr3, Expr4 expr4, Expr5 expr5, Expr6 expr6 )
{
#if 0
    BOOST_MPL_ASSERT_MSG( ( M == 4 && N == 1 ||
                            M == 1 && N == 4 ||
                            M == 2 && N == 2 ),
                          ( INVALID_MATRIX_SIZE ),
                          ( mpl::int_<M>, mpl::int_<N> ) );
#endif
    typedef vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5, Expr6> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5, Expr6>( expr1, expr2, expr3, expr4, expr5, expr6 ) ) );
}

/**
 * \brief matrix definition
 */
template<int M, int N,
         typename Expr1, typename Expr2, typename Expr3, typename Expr4,
         typename Expr5, typename Expr6, typename Expr7, typename Expr8,
         typename Expr9>
inline
Expr<vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5, Expr6, Expr7, Expr8, Expr9 > > >
mat( Expr1  expr1, Expr2  expr2, Expr3  expr3, Expr4 expr4,
     Expr5  expr5, Expr6  expr6, Expr7  expr7, Expr8 expr8, Expr9 expr9 )
{
#if 0
    BOOST_MPL_ASSERT_MSG( ( M == 4 && N == 1 ||
                            M == 1 && N == 4 ||
                            M == 2 && N == 2 ),
                          ( INVALID_MATRIX_SIZE ),
                          ( mpl::int_<M>, mpl::int_<N> ) );
#endif
    typedef vf::detail::Mat<M, N, fusion::vector<Expr1, Expr2, Expr3, Expr4,Expr5, Expr6, Expr7, Expr8, Expr9> > expr_t;
    return Expr<expr_t>( expr_t( fusion::vector<Expr1, Expr2, Expr3, Expr4, Expr5, Expr6, Expr7, Expr8, Expr9>( expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9 ) ) );
}

} // vf
} // Feel
#endif /* __VFVec_H */
