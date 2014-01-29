/* -*- Mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

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
   \file expr.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#ifndef FEELPP_EXPR_HPP
#define FEELPP_EXPR_HPP 1

#undef max
#include <boost/version.hpp>
#if (BOOST_VERSION >= 103400)
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <algorithm>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_class.hpp>
#include <boost/static_assert.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/multi_array.hpp>

#include <Eigen/Core>

#include <feel/feelcore/environment.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/context.hpp>

#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/shape.hpp>

namespace Feel
{
namespace vf
{
class GiNaCBase {};

/// \cond detail
typedef node<double>::type node_type;

enum
{
    CONTEXT_1 = ( 1<<0 ), /**< identifier 1 for the context */
    CONTEXT_2 = ( 1<<1 )  /**< identifier 2 for the context */
};

template<typename ExprT>
class ComponentsExpr
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;
    //integration order
    static const uint16_type imorder = ExprT::imorder;
    //the expression is a polynomial type?
    static const bool imIsPoly = ExprT::imIsPoly;


    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };




    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef ComponentsExpr<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    ComponentsExpr()
        :
        M_expr()
    {}

    explicit ComponentsExpr( expression_type const & __expr, int c1, int c2 )
        :
        M_expr( __expr ),
        M_c1( c1 ),
        M_c2( c2 )
    {}
    ~ComponentsExpr()
    {}

    //@}

    expression_type const& expression() const
    {
        return M_expr;
    }

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {

        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::NDim, Scalar, false> shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu ),
            M_c1( expr.M_c1 ),
            M_c2( expr.M_c2 )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev ),
            M_c1( expr.M_c1 ),
            M_c2( expr.M_c2 )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom ),
            M_c1( expr.M_c1 ),
            M_c2( expr.M_c2 )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, M_c1, M_c2, q );
        }

        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr.evalijq( i, j, M_c1, M_c2, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evaliq( i, M_c1, M_c2, q );
        }

        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evalq( M_c1, M_c2, q );
        }

        tensor_expr_type M_tensor_expr;
        const int M_c1, M_c2;
    };
    expression_type M_expr;
    int M_c1, M_c2;
};
class CstBase {};
class IntegratorBase {};
class LambdaExprBase {};
class LambdaExpr1 : public LambdaExprBase
{
public:

    static const size_type context = 0;
    static const bool is_terminal = false;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = false;

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

    typedef double value_type;

    template<typename TheExpr>
    struct Lambda
    {
        typedef typename TheExpr::expression_type type;
    };

    template<typename ExprT>
    typename Lambda<ExprT>::type
    operator()( ExprT const& e ) {
        return e.expression();
    }

    template<typename ExprT>
    typename Lambda<ExprT>::type
    operator()( ExprT const& e ) const { return e.expression(); }

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > >,fusion::pair<vf::detail::gmc<1>,boost::shared_ptr<vf::detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef LambdaExpr1 expression_type;
        typedef typename LambdaExpr1::value_type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;


        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
        }
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
        {
        }
        tensor( expression_type const& expr, Geo_t const& /*geom*/ )
        {
        }
        template<typename IM>
        void init( IM const& /*im*/ )
        {
        }
        void update( Geo_t const&, Basis_i_t const& , Basis_j_t const&  )
        {
        }
        void update( Geo_t const& , Basis_i_t const&  )
        {
        }
        void update( Geo_t const& )
        {
        }
        void update( Geo_t const&, uint16_type )
        {
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
        }


        value_type
        evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return 0;
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return 0;
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return 0;
        }
        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return 0;
        }
    };

};



/*!
  \class Expr
  \brief Variational Formulation Expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class Expr//: public boost::enable_shared_from_this<Expr<ExprT> >
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = ExprT::is_terminal;

    //integration order
    static const uint16_type imorder = ExprT::imorder;
    //the expression is a polynomial type?
    static const bool imIsPoly = ExprT::imIsPoly;


    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };




    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef Expr<ExprT> this_type;
    typedef boost::shared_ptr<this_type> this_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Expr()
        :
        M_expr()
    {}

    explicit Expr( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~Expr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{
    Expr<ComponentsExpr<Expr<ExprT> > >
    operator()( int c1 = 0, int c2 = 0 ) const
    {
        auto ex = ComponentsExpr<Expr<ExprT> >( Expr<ExprT>( M_expr ), c1, c2 );
        return Expr<ComponentsExpr<Expr<ExprT> > >( ex );
    }

    template<typename TheExpr>
    struct Lambda
    {
        typedef typename ExprT::template Lambda<TheExpr>::type expr_type;
        typedef Expr<expr_type> type;
        //typedef expr_type type;
    };



    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  )
        {
            //typename Lambda<TheExpr>::expr_type e1( M_expr(e) );
            //typename Lambda<TheExpr>::type r( Expr(e1 ) );
            //return r;
            return expr( M_expr( e ) );
        }

    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  ) const { return expr(M_expr(e)); }

    void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            this->setParameterValues( mp, boost::is_base_of<Feel::vf::GiNaCBase,expression_type>() );
        }
    void setParameterValues( std::map<std::string,value_type> const& mp, mpl::bool_<true> )
        {
            M_expr.setParameterValues( mp );
        }
    void setParameterValues( std::map<std::string,value_type> const& mp, mpl::bool_<false> )
        {
        }

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > >,fusion::pair<vf::detail::gmc<1>,boost::shared_ptr<vf::detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {

        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            M_tensor_expr.updateContext( ctx );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }

        Eigen::Matrix<value_type, shape::M, shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, q );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, c1, c2, q );
        }

        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evaliq( i, c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evalq( c1, c2, q );
        }

        tensor_expr_type M_tensor_expr;
    };
#if 0
    class Diff
    {
    public:

        value_type value() const
        {
            return __expression.value();
        }

        value_type grad( int __ith ) const
        {
            return __expression.grad( __ith );
        }

        value_type hessian( int __i, int __j ) const
        {
            return __expression.hessian( __i, __j );
        }

    };
#endif /**/
    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    //this_ptrtype ptr() { return boost::shared_from_this(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    expression_type& expression()
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

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __u,
                   boost::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const
    {
        DVLOG(2) << "calling assemble(u,v)\n";
        M_expr.assemble( __u, __v, __f );
        DVLOG(2) << "calling assemble(u,v) done\n";
    }

    template<typename Elem1, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __v,
                   FormType& __f ) const
    {
        DVLOG(2) << "calling assemble(v)\n";
        M_expr.assemble( __v, __f );
        DVLOG(2) << "calling assemble(v) done\n";
    }

    template<typename P0hType>
    typename P0hType::element_type
    broken( boost::shared_ptr<P0hType>& P0h ) const
    {
        return M_expr.broken( P0h );
    }
    //__typeof__( M_expr.evaluate() )
    //ublas::matrix<typename expression_type::value_type>

    typename expression_type::value_type
    evaluate( bool parallel = true, WorldComm const& worldcomm = Environment::worldComm() ) const
    {
        return M_expr.evaluate( parallel,worldcomm );
    }

    typename expression_type::value_type
    evaluateAndSum() const
    {
        return M_expr.evaluateAndSum();
    }
    std::string expressionStr() const
    {
        return std::string();
        //return M_expr.expressionStr();
    }


    //@}

protected:

private:

    mutable expression_type  M_expr;
};

template <typename ExprT>
Expr<ExprT>
expr( ExprT const& exprt )
{
    return Expr<ExprT>( exprt );
}

template <typename ExprT>
boost::shared_ptr<Expr<ExprT> >
exprPtr( ExprT const& exprt )
{
    return boost::shared_ptr<Expr<ExprT> >( new Expr<ExprT>( exprt ) );
}

extern Expr<LambdaExpr1> _e1;

/**
 * \class ExpressionOrder
 *
 * Class that compute the expression polynomial order of \p ExprT
 *
 * \tparam ExprT expression whose approximate order must be computed
 *
 * Note that if the expression is polynomial then the order is exact, however if
 * analytic functions such as exp, cos, sin ... then the order is only an
 * approximation.
 */
template<typename IntElts,typename ExprT>
struct ExpressionOrder
{

    typedef typename boost::tuples::template element<1, IntElts>::type element_iterator_type;
    typedef typename boost::remove_reference<typename element_iterator_type::reference>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_face_element_type;
    typedef typename the_face_element_type::super2::template Element<the_face_element_type>::type the_element_type;

    static const uint16_type nOrderGeo = the_element_type::nOrder;

    static const bool is_polynomial = ExprT::imIsPoly;
#if 0
    static const int value = boost::mpl::if_< boost::mpl::bool_< ExprT::imIsPoly > ,
                     typename boost::mpl::if_< boost::mpl::greater< boost::mpl::int_<ExprT::imorder>,
                     boost::mpl::int_<19> > ,
                     boost::mpl::int_<19>,
                     boost::mpl::int_<ExprT::imorder> >::type,
                     boost::mpl::int_<10> >::type::value;
#else
    // this is a very rough approximation
    static const int value = ( ExprT::imorder )?( ExprT::imorder*nOrderGeo ):( nOrderGeo );
    static const int value_1 = ExprT::imorder+(the_element_type::is_hypercube?nOrderGeo:0);
#endif


};


template<typename PrintExprT>
class PrintExpr
{
public:

    static const size_type context = PrintExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = PrintExprT::imorder;
    static const bool imIsPoly = PrintExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = PrintExprT::template HasTestFunction<Func>::result;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = PrintExprT::template HasTrialFunction<Func>::result;
    };

    /** @name Typedefs
     */
    //@{

    typedef PrintExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef PrintExpr<PrintExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit PrintExpr( expression_type const & __expr,
                        std::string const & __tag )
        :
        M_expr( __expr ),
        M_tag( __tag )
    {}
    ~PrintExpr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };
        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu ),
            M_tag( expr.tag() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev ),
            M_tag( expr.tag() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom ),
            M_tag( expr.tag() )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res= M_tensor_expr.evalijq( i, j, c1, c2, q );
            std::cout << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ") evalijq( " << i << "," << j << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            LOG(INFO) << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ") evalijq( " << i << "," << j << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            value_type res= M_tensor_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
            std::cout << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ") evalijq( " << i << "," << j << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            LOG(INFO) << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ") evalijq( " << i << "," << j << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }



        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res= M_tensor_expr.evaliq( i, c1, c2, q );
            std::cout << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ")  evaliq( " << i  << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            LOG(INFO) << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ")  evaliq( " << i  << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res= M_tensor_expr.evalq( c1, c2, q );
            std::cout << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ")  evalq( " << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            LOG(INFO) << "[print] " << M_tag << " shape(" << shape::M << "," << shape::N << ")  evalq( " << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }

        tensor_expr_type M_tensor_expr;
        std::string M_tag;
    };
#if 0
    class Diff
    {
    public:

        value_type value() const
        {
            return __expression.value();
        }

        value_type grad( int __ith ) const
        {
            return __expression.grad( __ith );
        }

        value_type hessian( int __i, int __j ) const
        {
            return __expression.hessian( __i, __j );
        }

    };
#endif /**/
    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    const std::string& tag() const
    {
        return M_tag;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __u,
                   boost::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const
    {
        DVLOG(2) << "calling assemble(u,v)\n";
        M_expr.assemble( __u, __v, __f );
        DVLOG(2) << "calling assemble(u,v) done\n";
    }

    template<typename Elem1, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __v,
                   FormType& __f ) const
    {
        DVLOG(2) << "calling assemble(v)\n";
        M_expr.assemble( __v, __f );
        DVLOG(2) << "calling assemble(v) done\n";
    }
#if 0
    //__typeof__( M_expr.evaluate() )
    ublas::matrix<typename expression_type::value_type>
    evaluate() const
    {
        return M_expr.evaluate();
    }
#endif


    //@}

protected:

private:

    mutable expression_type  M_expr;
    const std::string M_tag;
};
template<typename ExprT>
inline
Expr< PrintExpr<ExprT> >
print( ExprT v, std::string t="" )
{
    typedef PrintExpr<ExprT> print_t;
    return Expr< print_t >(  print_t( v, t ) );
}

/*!
  \class Trans
  \brief Transpose expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class Trans
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef Trans<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Trans( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~Trans()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{
    template<typename TheExpr>
    struct Lambda
    {
        typedef Trans<typename expression_type::template Lambda<TheExpr>::type> type;
    };
    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  ) { return trans(M_expr(e)); }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename Transpose<typename tensor_expr_type::shape>::type shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            M_tensor_expr.updateContext( ctx );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, c2, c1, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr.evalijq( i, j, c2, c1, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evaliq( i, c2, c1, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evalq( c2, c1, q );
        }

        tensor_expr_type M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
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

protected:

private:

    mutable expression_type  M_expr;
};

template<typename ExprT>
inline
Expr< Trans<ExprT> >
trans( ExprT v )
{
    typedef Trans<ExprT> trans_t;
    return Expr< trans_t >(  trans_t( v ) );
}


template < class T>
class Cst : public CstBase
{
public:

    //BOOST_STATIC_ASSERT( ::boost::is_arithmetic<T>::value );

    static const size_type context = vm::JACOBIAN;
    static const bool is_terminal = false;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = true;

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

    typedef typename mpl::if_<boost::is_reference_wrapper<T>,
            mpl::identity<T>,
            mpl::identity<mpl::identity<T> > >::type::type::type value_type;

    typedef Cst<T> expression_type;

    constexpr explicit Cst( const T& value )
        :
        M_constant( value )
    {
    }

    Cst( Cst const& __cst )
        :
        M_constant( __cst.M_constant )
    {
    }

    Cst&
    operator=( Cst const& c )
        {
            if ( this != &c )
                M_constant = c.M_constant;
            return *this;
        }

    constexpr value_type value() const
    {
        return M_constant;
    }

    constexpr value_type evaluate() const
    {
        return M_constant;
    }

    constexpr value_type evaluate( bool ) const
    {
        return M_constant;
    }

    constexpr value_type evaluate( bool, WorldComm const& ) const
    {
        return M_constant;
    }

    template<typename TheExpr>
    struct Lambda
    {
        typedef expression_type type;
    };
    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  ) { return typename Lambda<TheExpr>::type(M_constant); }

    template<typename Geo_t, typename Basis_i_t=mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename Cst<T>::expression_type expression_type;
        typedef typename Cst<T>::value_type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;


        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_constant( expr.value() )
        {
        }
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
            :
            M_constant( expr.value() )
        {
        }
        tensor( expression_type const& expr, Geo_t const& /*geom*/ )
            :
            M_constant( expr.value() )
        {
        }
        template<typename IM>
        void init( IM const& /*im*/ )
        {
        }
        void update( Geo_t const&, Basis_i_t const& , Basis_j_t const&  )
        {
        }
        void update( Geo_t const& , Basis_i_t const&  )
        {
        }
        void update( Geo_t const& )
        {
        }
        void update( Geo_t const&, uint16_type )
        {
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
        }

        constexpr value_type
        evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return M_constant;
        }


        constexpr value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return M_constant;
        }
        template<int PatternContext>
        constexpr value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return M_constant;
        }

        constexpr value_type
        evaliq( uint16_type /*i*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return M_constant;
        }
        constexpr value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return M_constant;
        }
        const value_type M_constant;
    };

protected:
    Cst() : M_constant( 0 )
    {
        //DVLOG(2) << "Cst::Cst( default ) : constant value: " << M_constant << "\n";
    }

    const T M_constant;
};

template<typename T>
inline
Expr< Cst<T> >
constant( T v )
{
    typedef Cst<T> cst_t;
    return Expr< cst_t >(  cst_t( v ) );
}

template<typename T>
inline
Expr< Cst<T> >
cst( T v )
{
    typedef Cst<T> cst_t;
    return Expr< cst_t >(  cst_t( v ) );
}
template<typename T>
inline
Expr< Cst<boost::reference_wrapper<T> > >
cst_ref( T& v )
{
    typedef Cst<boost::reference_wrapper<T> > cst_t;
    return Expr< cst_t >(  cst_t( boost::ref( v ) ) );
}
template<typename T>
inline
Expr< Cst<boost::reference_wrapper<T> > >
constant_ref( T& v )
{
    return cst_ref( v );
}



template<int CType>
class One
{
public:
    static const size_type context = 0;
    static const bool is_terminal = false;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = true;


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


    typedef One<CType> this_type;

    typedef double value_type;

    One() {}
    One( One const& /*__vff*/ ) {}

    template<typename TheExpr>
    struct Lambda
    {
        typedef this_type type;
    };
    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  ) { return this_type(); }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Vectorial, false, false> shape;
        static const bool theshape = ( shape::M == gmc_type::nDim && shape::N == 1 );
        BOOST_MPL_ASSERT_MSG( theshape,
                              INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_1,
                              ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );

        typedef typename expression_type::value_type value_type;

        static const uint16_type nComponents = gmc_type::nDim;
        static const int16_type vector_comp = ( CType==-1 )?1:CType;

        typedef typename mpl::if_<mpl::equal_to<mpl::int_<CType>,mpl::int_<-1> >,
                mpl::identity<ublas::scalar_vector<scalar_type> >,
                mpl::identity<ublas::unit_vector<scalar_type> > >::type::type vector_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_one( nComponents, vector_comp )
        {
            //std::cout << "one = " << M_one << "\n";
        }
        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/,
                Basis_i_t const& /*fev*/ )
            :
            M_one( nComponents, vector_comp )
        {
        }
        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/ )
            :
            M_one( nComponents, vector_comp )
        {
            //                 std::cout << "one = " << M_one << "\n"
            //                           << "M=" << shape::M << "\n"
            //                           << "N=" << shape::N << "\n";
        }
        template<typename IM>
        void init( IM const& /*im*/ )
        {
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
        {
        }
        void update( Geo_t const& /*geom*/ )
        {
        }
        void update( Geo_t const& /*geom*/, uint16_type /*face*/ )
        {
        }

        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }
        template<int PatternContext>
        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }

        FEELPP_STRONG_INLINE value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }
        FEELPP_STRONG_INLINE value_type
        evalq( uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }
        vector_type M_one;
    };

};

inline
Expr<One<-1> >
one()
{
    return Expr< One<-1> >(  One<-1>() );
}

inline
Expr<One<0> >
oneX()
{
    return Expr< One<0> >(  One<0>() );
}

inline
Expr<One<1> >
oneY()
{
    return Expr< One<1> >(  One<1>() );
}

inline
Expr<One<2> >
oneZ()
{
    return Expr< One<2> >(  One<2>() );
}
inline
Expr<One<0> >
unitX()
{
    return Expr< One<0> >(  One<0>() );
}

inline
Expr<One<1> >
unitY()
{
    return Expr< One<1> >(  One<1>() );
}

inline
Expr<One<2> >
unitZ()
{
    return Expr< One<2> >(  One<2>() );
}



/**
   \class UnaryPlus
   \brief handler for unary plus expression
*/
template < class ExprT >
class UnaryPlus
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    typedef ExprT expression_type;
    typedef typename ExprT::value_type value_type;
    typedef UnaryPlus<ExprT> this_type;

    UnaryPlus( const ExprT& expr )
        :
        M_expr( expr )
    {
    }
    UnaryPlus( UnaryPlus const& e )
        :
        M_expr( e.M_expr )
    {
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape shape;

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_t_expr( expr.expression(), geom, fev, feu )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_t_expr( expr.expression(), geom, fev )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_t_expr( expr.expression(), geom )
        {}
        template<typename IM>
        void init( IM const& im )
        {
            M_t_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev , Basis_j_t const& feu )
        {
            M_t_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev  )
        {
            M_t_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_t_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_t_expr.update( geom, face );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            M_t_expr.updateContext( ctx );
        }

        value_type
        evalij( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2 ) const
        {
            return M_t_expr.evalij( i, j, c1, c2 );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_t_expr.evalijq( i, j, c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_t_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_t_expr.evaliq( i, c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_t_expr.evalq( c1, c2, q );
        }
        tensor_expr_type M_t_expr;
    };
protected:
    UnaryPlus() {}

    expression_type M_expr;
};
template <class T> inline
Expr< UnaryPlus< Expr<T> > >
operator + ( const Expr<T>& expr )
{
    typedef UnaryPlus< Expr<T> > expr_t;

    return Expr< expr_t >( expr_t( expr ) );
}

/**
   \class UnaryMinus
   \brief handler for unary minus expression
*/
template < class ExprT >
class UnaryMinus
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    typedef ExprT expression_type;
    typedef typename ExprT::value_type value_type;
    typedef UnaryMinus<ExprT> this_type;

    UnaryMinus( const ExprT& expr )
        :
        M_expr( expr )
    {
        ;
    }

    UnaryMinus( UnaryMinus const& u )
        :
        M_expr( u.M_expr )
    {
        ;
    }

    UnaryMinus( UnaryMinus&& u )
        :
        M_expr( u.M_expr )
    {
    }

    UnaryMinus&
    operator=( UnaryMinus const& u )
    {
        if ( this != &u )
            M_expr = u.M_expr;
        return *this;
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape shape;

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_t_expr( expr.expression(), geom, fev, feu )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_t_expr( expr.expression(), geom, fev )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_t_expr( expr.expression(), geom )
        {}
        template<typename IM>
        void init( IM const& im )
        {
            M_t_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev , Basis_j_t const& feu )
        {
            M_t_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev  )
        {
            M_t_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_t_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_t_expr.update( geom, face );
        }

        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            M_t_expr.updateContext( ctx );
        }

        value_type
        evalij( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2 ) const
        {
            return -M_t_expr.evalij( i, j, c1, c2 );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -M_t_expr.evalijq( i, j, c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return -M_t_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -M_t_expr.evaliq( i, c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -M_t_expr.evalq( c1, c2, q );
        }
        tensor_expr_type M_t_expr;
    };
protected:
    UnaryMinus() {}

    expression_type M_expr;
};
template <class T> inline
Expr< UnaryMinus< Expr<T> > >
operator - ( const Expr<T>& expr )
{
    typedef UnaryMinus< Expr<T> > expr_t;

    return Expr< expr_t >( expr_t( expr ) );
}

template < typename ExprT1, typename ExprT2 >
class OpMax
{
public:
    static const size_type context = ExprT1::context | ExprT2::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ( ExprT1::imorder<ExprT2::imorder )*ExprT2::imorder + ( ExprT1::imorder>=ExprT2::imorder )*ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly && ExprT2::imIsPoly;

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

    typedef OpMax<ExprT1, ExprT2> this_type;
    typedef ExprT1 expression_1_type;
    typedef ExprT2 expression_2_type;
    typedef typename strongest_numeric_type<typename expression_1_type::value_type,
            typename expression_2_type::value_type>::type value_type;
    explicit OpMax( expression_1_type const& __expr1, expression_2_type const& __expr2  )
        :
        M_expr_1( __expr1 ),
        M_expr_2( __expr2 )
    {
        DVLOG(2) << "OpMax::OpMax default constructor\n";
    }

    OpMax( OpMax const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 ),
        M_expr_2( __vfp.M_expr_2 )
    {
        DVLOG(2) << "OpMax::OpMax copy constructor\n";
    }

    expression_1_type const& left() const
    {
        return M_expr_1;
    }
    expression_2_type const& right() const
    {
        return M_expr_2;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                typename r_type::value_type>::type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        BOOST_MPL_ASSERT_MSG( ( boost::is_same<typename l_type::shape, typename r_type::shape>::value ),
                              INVALID_SHAPES_FOR_MIN,
                              ( mpl::int_<l_type::shape::M>,mpl::int_<l_type::shape::N>,mpl::int_<r_type::shape::M>,mpl::int_<r_type::shape::N> ) );
        typedef typename l_type::shape shape;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev ),
            M_right( expr.right(), geom, fev )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom ),
            M_right( expr.right(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_left.init( im );
            M_right.init( im );
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
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom );
            M_right.update( geom );

        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom, face );
            M_right.update( geom, face );

        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );

        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type left = M_left.evalq( c1, c2, q );
            value_type right = M_right.evalq( c1, c2, q );
            return std::max( left, right );
        }

    private:
        gmc_ptrtype M_gmc;
        l_type M_left;
        r_type M_right;
    };

protected:
    OpMax() {}

    expression_1_type M_expr_1;
    expression_2_type M_expr_2;
};

template<typename ExprT1, typename ExprT2>
inline
Expr< OpMax<typename mpl::if_<boost::is_arithmetic<ExprT1>,
      mpl::identity<Cst<ExprT1> >,
      mpl::identity<ExprT1> >::type::type,
      typename mpl::if_<boost::is_arithmetic<ExprT2>,
      mpl::identity<Cst<ExprT2> >,
      mpl::identity<ExprT2> >::type::type
      > >
      max( ExprT1 const& __e1, ExprT2 const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef OpMax<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}

template < typename ExprT1, typename ExprT2 >
class OpMin
{
public:
    static const size_type context = ExprT1::context | ExprT2::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ( ExprT1::imorder<ExprT2::imorder )*ExprT2::imorder + ( ExprT1::imorder>=ExprT2::imorder )*ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly && ExprT2::imIsPoly;

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

    typedef OpMin<ExprT1, ExprT2> this_type;
    typedef ExprT1 expression_1_type;
    typedef ExprT2 expression_2_type;
    typedef typename strongest_numeric_type<typename expression_1_type::value_type,
            typename expression_2_type::value_type>::type value_type;

    explicit OpMin( expression_1_type const& __expr1, expression_2_type const& __expr2  )
        :
        M_expr_1( __expr1 ),
        M_expr_2( __expr2 )
    {
        DVLOG(2) << "OpMin::OpMin default constructor\n";
    }

    OpMin( OpMin const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 ),
        M_expr_2( __vfp.M_expr_2 )
    {
        DVLOG(2) << "OpMin::OpMin copy constructor\n";
    }

    expression_1_type const& left() const
    {
        return M_expr_1;
    }
    expression_2_type const& right() const
    {
        return M_expr_2;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                typename r_type::value_type>::type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        BOOST_MPL_ASSERT_MSG( ( boost::is_same<typename l_type::shape, typename r_type::shape>::value ),
                              INVALID_SHAPES_FOR_MIN,
                              ( mpl::int_<l_type::shape::M>,mpl::int_<l_type::shape::N>,mpl::int_<r_type::shape::M>,mpl::int_<r_type::shape::N> ) );
        typedef typename l_type::shape shape;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev ),
            M_right( expr.right(), geom, fev )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom ),
            M_right( expr.right(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_left.init( im );
            M_right.init( im );
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
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom );
            M_right.update( geom );

        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom, face );
            M_right.update( geom, face );

        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );

        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type left = M_left.evalq( c1, c2, q );
            value_type right = M_right.evalq( c1, c2, q );
            return std::min( left, right );
        }

    private:
        gmc_ptrtype M_gmc;
        l_type M_left;
        r_type M_right;
    };

protected:
    OpMin() {}

    expression_1_type M_expr_1;
    expression_2_type M_expr_2;
};

template<typename ExprT1, typename ExprT2>
inline
Expr< OpMin<typename mpl::if_<boost::is_arithmetic<ExprT1>,
      mpl::identity<Cst<ExprT1> >,
      mpl::identity<ExprT1> >::type::type,
      typename mpl::if_<boost::is_arithmetic<ExprT2>,
      mpl::identity<Cst<ExprT2> >,
      mpl::identity<ExprT2> >::type::type> >
      min( ExprT1 const& __e1, ExprT2 const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef OpMin<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}

template < typename ExprT1, typename ExprT2 >
class Pow
{
public:

    static const size_type context = ExprT1::context|ExprT2::context;
    static const bool is_terminal = false;

    /**
     * \warning the Pow order computation is wrong here, we actually need the
     * ExprT2 value (and not imorder) to multiply by ExprT1::imorder.
     */
    static const uint16_type imorder = ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly && ExprT2::imIsPoly;

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

    typedef Pow<ExprT1, ExprT2> this_type;
    typedef ExprT1 expression_1_type;
    typedef ExprT2 expression_2_type;
    typedef typename expression_1_type::value_type value_1_type;
    typedef typename expression_2_type::value_type value_2_type;
    typedef value_1_type value_type;

    // verify that all returning types are integral or floating types
    BOOST_STATIC_ASSERT( ::boost::is_arithmetic<value_1_type>::value  &&
                         ::boost::is_arithmetic<value_2_type>::value );

    explicit Pow( expression_1_type const& __expr1, expression_2_type const& __expr2  )
        :
        M_expr_1( __expr1 ),
        M_expr_2( __expr2 )
    {
        DVLOG(2) << "Pow::Pow default constructor\n";
    }

    Pow( Pow const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 ),
        M_expr_2( __vfp.M_expr_2 )
    {
        DVLOG(2) << "Pow::Pow copy constructor\n";
    }

    bool isSymetric() const
    {
        return false;
    }

    expression_1_type const& left() const
    {
        return M_expr_1;
    }
    expression_2_type const& right() const
    {
        return M_expr_2;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                typename r_type::value_type>::type value_type;


        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef typename l_type::shape shape;

        typedef typename Eigen::Matrix<value_type,shape::M,shape::N> loc_type;

        struct is_zero
        {
            static const bool value = l_type::is_zero::value;
        };

        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu ),
            M_loc( boost::extents[M_gmc->nPoints()]  )
        {
            update( geom );
        }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev ),
            M_right( expr.right(), geom, fev ),
            M_loc( boost::extents[M_gmc->nPoints()] )
        {
            update( geom );
        }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom ),
            M_right( expr.right(), geom ),
            M_loc(  boost::extents[M_gmc->nPoints()] )
        {
            update( geom );
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_left.init( im );
            M_right.init( im );
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
            const int npts =  fusion::at_key<key_type>( geom ).get()->nPoints();
            M_left.update( geom );
            M_right.update( geom );

            for ( int q = 0; q < npts; ++q )
                for ( int c1 = 0; c1 < shape::M; ++c1 )
                    for ( int c2 = 0; c2 < shape::N; ++c2 )
                    {
                        value_type left = M_left.evalq( c1, c2, q );
                        value_type right = M_right.evalq( c1, c2, q );
                        M_loc[q]( c1,c2 ) = std::pow( left, right );
                    }
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom, face );
            M_right.update( geom, face );

            for ( int q = 0; q < M_gmc->nPoints(); ++q )
                for ( int c1 = 0; c1 < shape::M; ++c1 )
                    for ( int c2 = 0; c2 < shape::N; ++c2 )
                    {
                        value_type left = M_left.evalq( c1, c2, q );
                        value_type right = M_right.evalq( c1, c2, q );
                        M_loc[q]( c1,c2 ) = std::pow( left, right );
                    }
        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q  ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q  ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q, mpl::int_<shape::rank>() );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return M_loc[q]( 0,0 );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const
        {
            if ( shape::M > shape::N )
                return M_loc[q]( c1,0 );

            return M_loc[q]( 0,c2 );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const
        {
            return M_loc[q]( c1,c2 );
        }

    private:
        gmc_ptrtype M_gmc;
        l_type M_left;
        r_type M_right;
        //ublas::vector<double> M_loc;
        boost::multi_array<loc_type,1> M_loc;
    };

protected:
    Pow() {}

    expression_1_type M_expr_1;
    expression_2_type M_expr_2;
};

template<typename ExprT1,  typename ExprT2>
inline
Expr< Pow<typename mpl::if_<boost::is_arithmetic<ExprT1>,
      mpl::identity<Cst<ExprT1> >,
      mpl::identity<ExprT1> >::type::type,
      typename mpl::if_<boost::is_arithmetic<ExprT2>,
      mpl::identity<Cst<ExprT2> >,
      mpl::identity<ExprT2> >::type::type> >
      pow( ExprT1 const& __e1, ExprT2 const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef Pow<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}


template<typename ExprT1,  typename ExprT2>
inline
Expr< Pow<typename mpl::if_<boost::is_arithmetic<ExprT1>,
                            mpl::identity<Cst<ExprT1> >,
                            mpl::identity<Expr<ExprT1> > >::type::type,
          typename mpl::if_<boost::is_arithmetic<ExprT2>,
                            mpl::identity<Cst<ExprT2> >,
                            mpl::identity<Expr<ExprT2> > >::type::type> >
operator^( typename mpl::if_<boost::is_arithmetic<ExprT1>,
                             mpl::identity<ExprT1>,
                             mpl::identity<Expr<ExprT1> > >::type::type const& __e1,
           typename mpl::if_<boost::is_arithmetic<ExprT2>,
                             mpl::identity<ExprT2>,
                             mpl::identity<Expr<ExprT2> > >::type::type const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef Pow<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}

/**
 * \class EvalFace
 * \brief Variational Formulation Expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int GeoId, typename ExprT>
class EvalFace
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;
    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef EvalFace<GeoId,ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit EvalFace( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~EvalFace()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename VecGeo_t, typename Basis_i_t = boost::none_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename fusion::result_of::at_c<VecGeo_t,GeoId>::type Geo_t;

        typedef typename expression_type::template tensor<Geo_t,
                Basis_i_t,
                Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        tensor( this_type const& expr,
                VecGeo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), fusion::at_c<GeoId>(  geom ), fev, feu )
        {}

        tensor( this_type const& expr,
                VecGeo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), fusion::at_c<GeoId>(  geom ), fev )
        {}

        tensor( this_type const& expr,
                VecGeo_t const& geom )
            :
            M_tensor_expr( expr.expression(), fusion::at_c<GeoId>(  geom ) )
        {}
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }
        void update( VecGeo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( fusion::at_c<GeoId>(  geom ), fev, feu );
        }
        void update( VecGeo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( fusion::at_c<GeoId>(  geom ), fev );
        }
        void update( VecGeo_t const& geom )
        {
            M_tensor_expr.update( fusion::at_c<GeoId>(  geom ) );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, int q ) const
        {
            return M_tensor_expr.evalijq( i, j, q );
        }


        value_type
        evaliq( uint16_type i, int q ) const
        {
            return M_tensor_expr.evaliq( i, q );
        }

        value_type
        evalq( int q, int c ) const
        {
            return M_tensor_expr.evalq( q, c );
        }

        tensor_expr_type M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    expression_type const& expression() const
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

protected:

private:

    mutable expression_type  M_expr;
};

template<int GeoId, typename ExprT>
inline
Expr< EvalFace<GeoId, ExprT> >
evalface( ExprT const& v )
{
    typedef EvalFace<GeoId, ExprT> eval_t;
    return Expr< eval_t >(  eval_t( v ) );
}

template < class Element, int Type>
class GElem
{
public:

    static const size_type context = vm::JACOBIAN |vm::POINT;
    static const bool is_terminal = false;

    typedef Element element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef GElem<element_type, Type> this_type;
    typedef this_type self_type;

    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename functionspace_type::reference_element_type* fe_ptrtype;
    typedef typename functionspace_type::reference_element_type fe_type;
    typedef typename functionspace_type::geoelement_type geoelement_type;
    typedef typename functionspace_type::gm_type gm_type;
    typedef typename functionspace_type::value_type value_type;
    static const uint16_type rank = fe_type::rank;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;
    typedef std::map<size_type,std::vector<element_ptrtype> > basis_type;

    static const uint16_type imorder = element_type::functionspace_type::basis_type::nOrder;
    static const bool imIsPoly = true;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ( Type==0 );
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ( Type==1 );
    };


    GElem ( std::map<size_type,std::vector<element_ptrtype> > const& v )
        :
        M_basis ( v )
    {
        typename basis_type::iterator it = M_basis.begin();
        typename basis_type::iterator en = M_basis.end();

        for ( ; it != en; ++it )
            for ( uint16_type i = 0; i < it->second.size(); ++i )
                it->second[i]->updateGlobalValues();
    }
    GElem( GElem const& op )
        :
        M_basis ( op.M_basis )
    {

    }

    basis_type const&  basis() const
    {
        return M_basis;
    }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        typedef typename mpl::if_<mpl::equal_to<mpl::int_<rank>,
                mpl::int_<0> >,
                mpl::identity<Shape<gmc_type::NDim, Scalar, false> >,
                typename mpl::if_<mpl::equal_to<mpl::int_<rank>,
                mpl::int_<1> >,
                mpl::identity<Shape<gmc_type::NDim, Vectorial, false> >,
                mpl::identity<Shape<gmc_type::NDim, Tensor2, false> > >::type>::type::type shape;
        typedef typename fe_type::PreCompute pc_type;
        typedef boost::shared_ptr<pc_type> pc_ptrtype;
        typedef typename fe_type::template Context<context, fe_type, gm_type,geoelement_type,gmc_type::context> ctx_type;
        typedef boost::shared_ptr<ctx_type> ctx_ptrtype;

        typedef typename expression_type::value_type value_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_expr( expr ),
            M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            M_ctx( new ctx_type( expr.basis().begin()->second[0]->functionSpace()->fe(),
                                 fusion::at_key<key_type>( geom ), ( pc_ptrtype const& )M_pc ) ),
            M_loc( expr.basis().begin()->second.size() )
        {
            for ( uint16_type i = 0; i < M_loc.size(); ++i )
                M_loc[i].resize( boost::extents[M_pc.nPoints()][nComponents1][nComponents2] );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            M_expr( expr ),
            M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            M_ctx( new ctx_type( expr.basis().begin()->second[0]->functionSpace()->fe(),
                                 fusion::at_key<key_type>( geom ), ( pc_ptrtype const& )M_pc ) ),
            M_loc( expr.basis().begin()->second.size() )
        {
            for ( uint16_type i = 0; i < M_loc.size(); ++i )
                M_loc[i].resize( boost::extents[M_pc.nPoints()] );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            M_expr( expr ),
            M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            M_ctx( new ctx_type( expr.basis().begin()->second[0]->functionSpace()->fe(),
                                 fusion::at_key<key_type>( geom ), ( pc_ptrtype const& )M_pc ) ),
            M_loc( expr.basis().begin()->second.size() )
        {
            for ( uint16_type i = 0; i < M_loc.size(); ++i )
                M_loc[i].resize( boost::extents[M_pc.nPoints()] );
        }
        template<typename IM>
        void init( IM const& /*im*/ )
        {

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
            //VLOG(1) << "[GElem] updating element " << fusion::at_key<key_type>( geom )->id() << "\n";
            typename basis_type::iterator it = const_cast<basis_type&>( M_expr.basis() ).find( fusion::at_key<key_type>( geom )->id() );
            typename basis_type::iterator en = const_cast<basis_type&>( M_expr.basis() ).end();

            FEELPP_ASSERT( it != en )( fusion::at_key<key_type>( geom )->id() ).error ( "invalid basis function to integrate" );

            for ( uint16_type i = 0; i < M_loc.size(); ++i )
            {
                //M_loc[i] = it->second[i]->id( *M_ctx, M_pc, M_loc[i] );
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( Type == 0 )
                return M_loc[i]( c1,c2,q );

            return M_loc[j]( c1,c2,q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext> ) const
        {
            if ( Type == 0 )
                return M_loc[i]( c1,c2,q );

            return M_loc[j]( c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_loc[i]( c1,c2,q );
        }
        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {

        }
    private:
        this_type const& M_expr;
        pc_type M_pc;
        ctx_ptrtype M_ctx;
        std::vector<typename element_type::id_type> M_loc;
    };
private:
    basis_type M_basis;

};

template<typename Elem>
inline
Expr< GElem<Elem,1> >
basist( std::map<size_type,std::vector<boost::shared_ptr<Elem> > > const& v )
{
    typedef GElem<Elem,1> expr_t;
    return Expr< expr_t >(  expr_t( v ) );
}
template<typename Elem>
inline
Expr< GElem<Elem,0> >
basis( std::map<size_type,std::vector<boost::shared_ptr<Elem> > > const& v )
{
    typedef GElem<Elem,0> expr_t;
    return Expr< expr_t >(  expr_t( v ) );
}

/// \endcond
} // vf


using namespace vf;

} // feel
#endif /* FEELPP_EXPR_HPP */
