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
#include <boost/none.hpp>
#include <algorithm>
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

#include <feel/feelvf/exprbase.hpp>
#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/shape.hpp>
#include <feel/feelvf/lambda.hpp>

namespace Feel
{
namespace vf
{
class GiNaCBase;

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
    typedef value_type evaluate_type;
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

class IntegratorBase {};



/*!
  \class Expr
  \brief Variational Formulation Expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class Expr : public ExprBase //: public boost::enable_shared_from_this<Expr<ExprT> >
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
    typedef typename expression_type::evaluate_type evaluate_type;
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
    virtual ~Expr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename... TheExprs>
    struct Lambda
    {
        typedef typename ExprT::template Lambda<TheExprs...>::type expr_type;
        typedef Expr<expr_type> type;
        //typedef expr_type type;
    };


#if 0
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr...  e  )
        {
            return expr( M_expr( e... ) );
        }
#else
    template<typename TheExpr1>
    typename Lambda<TheExpr1>::type
    operator()( TheExpr1 const&  e   )
        {
            return expr( M_expr( e ) );
        }
    template<typename TheExpr1,typename TheExpr2>
    typename Lambda<Expr<TheExpr1>,Expr<TheExpr2>>::type
    operator()( Expr<TheExpr1> const&  e1, Expr<TheExpr2> const& e2 )
        {
            return expr( M_expr( e1, e2 ) );
        }

    template<typename TheExpr1,typename TheExpr2,typename TheExpr3>
    typename Lambda<TheExpr1,TheExpr2,TheExpr3>::type
    operator()( TheExpr1 const&  e1, TheExpr2 const& e2, TheExpr3 const& e3  )
        {
            return expr( M_expr( e1, e2, e3 ) );
        }
#endif

#if 0
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) const { return expr(M_expr(e...)); }
#else
    template<typename TheExpr1>
    typename Lambda<TheExpr1>::type
    operator()( TheExpr1 const&  e   ) const
        {
            return expr( M_expr( e ) );
        }
    template<typename TheExpr1,typename TheExpr2>
    typename Lambda<Expr<TheExpr1>,Expr<TheExpr2>>::type
        operator()( Expr<TheExpr1> const&  e1, Expr<TheExpr2> const& e2 ) const
        {
            return expr( M_expr( e1, e2 ) );
        }
    template<typename TheExpr1,typename TheExpr2,typename TheExpr3>
    typename Lambda<TheExpr1,TheExpr2,TheExpr3>::type
    operator()( TheExpr1 const&  e1, TheExpr2 const& e2, TheExpr3 const& e3  ) const
        {
            return expr( M_expr( e1, e2, e3 ) );
        }

#endif

    Expr<ComponentsExpr<Expr<ExprT> > >
    operator()( int c1, int c2 )
    {
        auto ex = ComponentsExpr<Expr<ExprT> >( Expr<ExprT>( M_expr ), c1, c2 );
        return Expr<ComponentsExpr<Expr<ExprT> > >( ex );
    }

    Expr<ComponentsExpr<Expr<ExprT> > >
    operator()( int c1, int c2 ) const
    {
        auto ex = ComponentsExpr<Expr<ExprT> >( Expr<ExprT>( M_expr ), c1, c2 );
        return Expr<ComponentsExpr<Expr<ExprT> > >( ex );
    }

    void setParameterValues( std::pair<std::string,value_type> const& mp )
        {
            this->setParameterValues( { { mp.first, mp.second } } );
        }
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

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

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
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom )
        {
        }

        gmc_ptrtype geom() const { return M_geo; }

        int nPoints() const { return M_geo->nPoints(); }

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

        gmc_ptrtype M_geo;
        //Geo_t const& M_geo;
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

    evaluate_type
    evaluate( std::pair<std::string,value_type> const& mp  )
    {
        return M_expr.evaluate( { { mp.first, mp.second } } );
    }
    evaluate_type
    evaluate( std::map<std::string,value_type> const& mp  )
    {
        return M_expr.evaluate( mp );
    }
    evaluate_type
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

template <typename ExprT>
std::ostream&
operator<<( std::ostream& os, Expr<ExprT> const& exprt )
{
    os << exprt.expression();
    return os;
}


extern Expr<LambdaExpr1> _e1;
extern Expr<LambdaExpr2> _e2;
extern Expr<LambdaExpr3> _e3;

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
