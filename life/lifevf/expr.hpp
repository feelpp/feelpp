/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-17
 */
#ifndef __Expr_H
#define __Expr_H 1

#undef max
#include <boost/version.hpp>
#if (BOOST_VERSION >= 103400)
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <algorithm>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/multi_array.hpp>

#include <life/lifepoly/policy.hpp>
#include <life/lifepoly/context.hpp>
#include <life/lifevf/shape.hpp>

namespace Life
{
namespace vf
{
/// \cond detail
typedef node<double>::type node_type;

enum
{
    CONTEXT_1 = ( 1<<0 ), /**< identifier 1 for the context */
    CONTEXT_2 = ( 1<<1 )  /**< identifier 2 for the context */
};

/*!
  \class Expr
  \brief Variational Formulation Expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class Expr
{
public:

    static const size_type context = ExprT::context;

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

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Expr()
        :
        _M_expr()
    {}

    explicit Expr( expression_type const & __expr )
        :
        _M_expr( __expr )
    {}
    ~Expr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<detail::gmc<0>,boost::shared_ptr<detail::gmc<0> > >,fusion::pair<detail::gmc<1>,boost::shared_ptr<detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {

        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;

        template <class Args> struct sig { typedef value_type type; };

        struct is_zero { static const bool value = tensor_expr_type::is_zero::value; };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), geom )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            _M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
            {
                _M_tensor_expr.update( geom, face );
            }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q );
        }

        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }

        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evaliq( i, c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evalq( c1, c2, q );
        }

        tensor_expr_type _M_tensor_expr;
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

    bool isSymetric() const { return _M_expr.isSymetric(); }

    expression_type const& expression() const { return _M_expr; }

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
        Debug( 5051 ) << "calling assemble(u,v)\n";
        _M_expr.assemble( __u, __v, __f );
        Debug( 5051 ) << "calling assemble(u,v) done\n";
    }

    template<typename Elem1, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __v,
                   FormType& __f ) const
    {
        Debug( 5051 ) << "calling assemble(v)\n";
        _M_expr.assemble( __v, __f );
        Debug( 5051 ) << "calling assemble(v) done\n";
    }

    //__typeof__( _M_expr.evaluate() )
    ublas::matrix<typename expression_type::value_type>
    evaluate() const
    {
        return _M_expr.evaluate();
    }

    ublas::matrix<typename expression_type::value_type>
    evaluateAndSum() const
    {
        return _M_expr.evaluateAndSum();
    }



    //@}

protected:

private:

    mutable expression_type  _M_expr;
};

template<typename PrintExprT>
class PrintExpr
{
public:

    static const size_type context = PrintExprT::context;

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
        _M_expr( __expr ),
        _M_tag( __tag )
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

        template <class Args> struct sig { typedef value_type type; };
        struct is_zero { static const bool value = tensor_expr_type::is_zero::value; };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), geom, fev, feu ),
            _M_tag( expr.tag() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), geom, fev ),
            _M_tag( expr.tag() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), geom ),
            _M_tag( expr.tag() )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            _M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
            {
                _M_tensor_expr.update( geom, face );
            }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res= _M_tensor_expr.evalijq( i, j, c1, c2, q );
            std::cout << "[print] " << _M_tag << " shape(" << shape::M << "," << shape::N << ") evalijq( " << i << "," << j << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            value_type res= _M_tensor_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
            std::cout << "[print] " << _M_tag << " shape(" << shape::M << "," << shape::N << ") evalijq( " << i << "," << j << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }


        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res= _M_tensor_expr.evaliq( i, c1, c2, q );
            std::cout << "[print] " << _M_tag << " shape(" << shape::M << "," << shape::N << ")  evaliq( " << i  << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res= _M_tensor_expr.evalq( c1, c2, q );
            std::cout << "[print] " << _M_tag << " shape(" << shape::M << "," << shape::N << ")  evalq( " << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }

        tensor_expr_type _M_tensor_expr;
        std::string _M_tag;
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

    bool isSymetric() const { return _M_expr.isSymetric(); }

    expression_type const& expression() const { return _M_expr; }

    const std::string& tag() const { return _M_tag; }

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
        Debug( 5051 ) << "calling assemble(u,v)\n";
        _M_expr.assemble( __u, __v, __f );
        Debug( 5051 ) << "calling assemble(u,v) done\n";
    }

    template<typename Elem1, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __v,
                   FormType& __f ) const
    {
        Debug( 5051 ) << "calling assemble(v)\n";
        _M_expr.assemble( __v, __f );
        Debug( 5051 ) << "calling assemble(v) done\n";
    }

    //__typeof__( _M_expr.evaluate() )
    ublas::matrix<typename expression_type::value_type>
    evaluate() const
    {
        return _M_expr.evaluate();
    }



    //@}

protected:

private:

    mutable expression_type  _M_expr;
    const std::string _M_tag;
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
        _M_expr( __expr )
    {}
    ~Trans()
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

        typedef typename Transpose<typename tensor_expr_type::shape>::type shape;

        template <class Args> struct sig { typedef value_type type; };

        struct is_zero { static const bool value = tensor_expr_type::is_zero::value; };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            _M_tensor_expr.init( im );
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
            {
                _M_tensor_expr.update( geom, face );
            }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evalijq( i, j, c2, c1, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return _M_tensor_expr.evalijq( i, j, c2, c1, q, mpl::int_<PatternContext>() );
        }

        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evaliq( i, c2, c1, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evalq( c2, c1, q );
        }

        tensor_expr_type _M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const { return _M_expr.isSymetric(); }

    expression_type const& expression() const { return _M_expr; }

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

    mutable expression_type  _M_expr;
};

template<typename ExprT>
inline
Expr< Trans<ExprT> >
trans( ExprT v )
{
    typedef Trans<ExprT> trans_t;
    return Expr< trans_t >(  trans_t( v ) );
}

/*!
  \class DiagExpr
  \brief Diag expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class DiagExpr
{
public:

    static const size_type context = ExprT::context;

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
    typedef DiagExpr<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit DiagExpr( expression_type const & __expr )
        :
        _M_expr( __expr )
    {}
    ~DiagExpr()
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

        typedef typename Diag<typename tensor_expr_type::shape>::type shape;

        template <class Args> struct sig { typedef value_type type; };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), geom, fev, feu )
        {
//             std::cout << "diag ijq rank : " << shape::rank << "\n";
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), geom, fev )
        {
//             std::cout << "diag iq rank : " << shape::rank << "\n";
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), geom )
        {
//             std::cout << "diag q rank : " << shape::rank << "\n";
        }

        template<typename IM>
        void init( IM const& im )
        {
            _M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_tensor_expr.update( geom );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalijq( i, j, c1, c2, q, mpl::int_<shape::rank>() );
        }

    private:
        /* evalijq */
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const
        {
            if ( i.component() == c1 && j.component() == c1 )
                return _M_tensor_expr.evalijq( i, j, c1, c2, q );
            return value_type( 0 );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const
        {
            if ( i.component() == c1 && j.component() == c1 )
                return _M_tensor_expr.evalijq( i, j, c1, c2, q );
            return value_type( 0 );
        }
        tensor_expr_type _M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const { return _M_expr.isSymetric(); }

    expression_type const& expression() const { return _M_expr; }

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

    mutable expression_type  _M_expr;
};

template<typename ExprT>
inline
Expr< DiagExpr<ExprT> >
diag( ExprT v )
{
    typedef DiagExpr<ExprT> diag_t;
    return Expr< diag_t >(  diag_t( v ) );
}


template < class T>
class Cst
{
public:

    //BOOST_STATIC_ASSERT( ::boost::is_arithmetic<T>::value );

    static const size_type context = vm::JACOBIAN;

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

    explicit Cst(const T& value)
        :
        _M_constant(value)
    {
    }

    Cst( Cst const& __cst )
        :
        _M_constant( __cst._M_constant )
    {
    }


    value_type value() const
    {
        return _M_constant;
    }

    value_type evaluate() const { return _M_constant; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename Cst<T>::expression_type expression_type;
        typedef typename Cst<T>::value_type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >, mpl::identity<detail::gmc<0> >, mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;


        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            typedef value_type type;
        };

        struct is_zero { static const bool value = false; };

        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            _M_constant( expr.value() )
        {
        }
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
            :
            _M_constant( expr.value() )
        {
        }
        tensor( expression_type const& expr, Geo_t const& /*geom*/ )
            :
            _M_constant( expr.value() )
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
        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& /*i*/, IndexJ const& /*j*/ ) const
        {
            return _M_constant;
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return _M_constant;
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return _M_constant;
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& /*i*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return _M_constant;
        }
        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return _M_constant;
        }
        const value_type _M_constant;
    };

protected:
    Cst() : _M_constant( 0 )
    {
        //Debug( 5051 ) << "Cst::Cst( default ) : constant value: " << _M_constant << "\n";
    }

    const T _M_constant;
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
    One( One const& /*__vff*/) {}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                                  mpl::identity<detail::gmc<0> >,
                                  mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Vectorial, false, false> shape;
        static const bool theshape = (shape::M == gmc_type::nDim && shape::N == 1);
        BOOST_MPL_ASSERT_MSG( theshape,
                              INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_1,
                              (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

        typedef typename expression_type::value_type value_type;

        static const uint16_type nComponents = gmc_type::nDim;
        static const int16_type vector_comp = (CType==-1)?1:CType;

        typedef typename mpl::if_<mpl::equal_to<mpl::int_<CType>,mpl::int_<-1> >,
                                  mpl::identity<ublas::scalar_vector<scalar_type> >,
                                  mpl::identity<ublas::unit_vector<scalar_type> > >::type::type vector_type;

        struct is_zero { static const bool value = false; };

        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            _M_one( nComponents, vector_comp )
            {
                //std::cout << "one = " << _M_one << "\n";
            }
        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/,
                Basis_i_t const& /*fev*/ )
            :
            _M_one( nComponents, vector_comp )
            {
            }
        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/ )
            :
            _M_one( nComponents, vector_comp )
            {
//                 std::cout << "one = " << _M_one << "\n"
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
        template<typename IndexI, typename IndexJ>
        value_type const&
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return _M_one[c1];
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return _M_one[c1];
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& /*i*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return _M_one[c1];
        }
        value_type
        evalq( uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return _M_one[c1];
        }
        vector_type _M_one;
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


/**
   \class UnaryPlus
   \brief handler for unary plus expression
*/
template < class ExprT >
class UnaryPlus
{
public:

    static const size_type context = ExprT::context;

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

    UnaryPlus(const ExprT& expr)
        :
        _M_expr(expr)
    {;}

    expression_type const& expression() const { return _M_expr; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape shape;

        struct is_zero { static const bool value = tensor_expr_type::is_zero::value; };

        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            _M_t_expr( expr.expression(), geom, fev, feu )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            _M_t_expr( expr.expression(), geom, fev )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            _M_t_expr( expr.expression(), geom )
        {}
        template<typename IM>
        void init( IM const& im )
        {
            _M_t_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev , Basis_j_t const& feu )
        {
            _M_t_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev  )
        {
            _M_t_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_t_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            _M_t_expr.update( geom, face );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2 ) const
        {
            return _M_t_expr.evalij( i, j, c1, c2 );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_t_expr.evalijq( i, j, c1, c2, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return _M_t_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_t_expr.evaliq( i, c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_t_expr.evalq( c1, c2, q );
        }
        tensor_expr_type _M_t_expr;
    };
protected:
    UnaryPlus() {}

    const expression_type& _M_expr;
};
template <class T> inline
Expr< UnaryPlus< Expr<T> > >
operator + (const Expr<T>& expr)
{
    typedef UnaryPlus< Expr<T> > expr_t;

    return Expr< expr_t >( expr_t(expr) );
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

    UnaryMinus(const ExprT& expr)
        :
        _M_expr(expr)
    {;}

    expression_type const& expression() const { return _M_expr; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape shape;

        struct is_zero { static const bool value = tensor_expr_type::is_zero::value; };

        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            _M_t_expr( expr.expression(), geom, fev, feu )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            _M_t_expr( expr.expression(), geom, fev )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            _M_t_expr( expr.expression(), geom )
        {}
        template<typename IM>
        void init( IM const& im )
        {
            _M_t_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev , Basis_j_t const& feu )
        {
            _M_t_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev  )
        {
            _M_t_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_t_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            _M_t_expr.update( geom, face );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2 ) const
        {
            return -_M_t_expr.evalij( i, j, c1, c2 );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -_M_t_expr.evalijq( i, j, c1, c2, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return -_M_t_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -_M_t_expr.evaliq( i, c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -_M_t_expr.evalq( c1, c2, q );
        }
        tensor_expr_type _M_t_expr;
    };
protected:
    UnaryMinus() {}

    const expression_type& _M_expr;
};
template <class T> inline
Expr< UnaryMinus< Expr<T> > >
operator - (const Expr<T>& expr)
{
    typedef UnaryMinus< Expr<T> > expr_t;

    return Expr< expr_t >( expr_t(expr) );
}

template < typename ExprT1, typename ExprT2 >
class OpMax
{
public:
    static const size_type context = ExprT1::context | ExprT2::context;

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
        _M_expr_1( __expr1 ),
        _M_expr_2( __expr2 )
        {
            Debug( 5051 ) << "OpMax::OpMax default constructor\n";
        }

    OpMax( OpMax const& __vfp  )
        :
        _M_expr_1( __vfp._M_expr_1 ),
        _M_expr_2( __vfp._M_expr_2 )
        {
            Debug( 5051 ) << "OpMax::OpMax copy constructor\n";
        }

    expression_1_type const& left() const { return _M_expr_1; }
    expression_2_type const& right() const { return _M_expr_2; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                                                typename r_type::value_type>::type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                                  mpl::identity<detail::gmc<0> >,
                                  mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        BOOST_MPL_ASSERT_MSG( (boost::is_same<typename l_type::shape, typename r_type::shape>::value ),
                              INVALID_SHAPES_FOR_MIN,
                              ( mpl::int_<l_type::shape::M>,mpl::int_<l_type::shape::N>,mpl::int_<r_type::shape::M>,mpl::int_<r_type::shape::N> ) );
        typedef typename l_type::shape shape;

        struct is_zero { static const bool value = false; };

        tensor( expression_type const& expr, Geo_t const& geom,
                 Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom, fev, feu ),
            _M_right( expr.right(), geom, fev, feu )
            {
            }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom, fev ),
            _M_right( expr.right(), geom, fev )
            {
            }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom ),
            _M_right( expr.right(), geom )
            {
            }
        template<typename IM>
        void init( IM const& im )
        {
            _M_left.init( im );
            _M_right.init( im );
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
            _M_gmc = fusion::at_key<key_type>( geom ).get();
            _M_left.update( geom );
            _M_right.update( geom );

        }
        void update( Geo_t const& geom, uint16_type face )
        {
            _M_gmc = fusion::at_key<key_type>( geom ).get();
            _M_left.update( geom, face );
            _M_right.update( geom, face );

        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Life::detail::ignore_unused_variable_warning(i);
            Life::detail::ignore_unused_variable_warning(j);
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
            value_type left = _M_left.evalq( c1, c2, q );
            value_type right = _M_right.evalq( c1, c2, q );
            return std::max( left, right );
        }

    private:
        gmc_ptrtype _M_gmc;
        l_type _M_left;
        r_type _M_right;
    };

protected:
    OpMax() {}

    expression_1_type _M_expr_1;
    expression_2_type _M_expr_2;
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
        _M_expr_1( __expr1 ),
        _M_expr_2( __expr2 )
        {
            Debug( 5051 ) << "OpMin::OpMin default constructor\n";
        }

    OpMin( OpMin const& __vfp  )
        :
        _M_expr_1( __vfp._M_expr_1 ),
        _M_expr_2( __vfp._M_expr_2 )
        {
            Debug( 5051 ) << "OpMin::OpMin copy constructor\n";
        }

    expression_1_type const& left() const { return _M_expr_1; }
    expression_2_type const& right() const { return _M_expr_2; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                                                typename r_type::value_type>::type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                                  mpl::identity<detail::gmc<0> >,
                                  mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        BOOST_MPL_ASSERT_MSG( (boost::is_same<typename l_type::shape, typename r_type::shape>::value ),
                              INVALID_SHAPES_FOR_MIN,
                              ( mpl::int_<l_type::shape::M>,mpl::int_<l_type::shape::N>,mpl::int_<r_type::shape::M>,mpl::int_<r_type::shape::N> ) );
        typedef typename l_type::shape shape;

        struct is_zero { static const bool value = false; };

        tensor( expression_type const& expr, Geo_t const& geom,
                 Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom, fev, feu ),
            _M_right( expr.right(), geom, fev, feu )
            {
            }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom, fev ),
            _M_right( expr.right(), geom, fev )
            {
            }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom ),
            _M_right( expr.right(), geom )
            {
            }
        template<typename IM>
        void init( IM const& im )
        {
            _M_left.init( im );
            _M_right.init( im );
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
            _M_gmc = fusion::at_key<key_type>( geom ).get();
            _M_left.update( geom );
            _M_right.update( geom );

        }
        void update( Geo_t const& geom, uint16_type face )
        {
            _M_gmc = fusion::at_key<key_type>( geom ).get();
            _M_left.update( geom, face );
            _M_right.update( geom, face );

        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
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
            value_type left = _M_left.evalq( c1, c2, q );
            value_type right = _M_right.evalq( c1, c2, q );
            return std::min( left, right );
        }

    private:
        gmc_ptrtype _M_gmc;
        l_type _M_left;
        r_type _M_right;
    };

protected:
    OpMin() {}

    expression_1_type _M_expr_1;
    expression_2_type _M_expr_2;
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
        _M_expr_1( __expr1 ),
        _M_expr_2( __expr2 )
        {
            Debug( 5051 ) << "Pow::Pow default constructor\n";
        }

    Pow( Pow const& __vfp  )
        :
        _M_expr_1( __vfp._M_expr_1 ),
        _M_expr_2( __vfp._M_expr_2 )
        {
            Debug( 5051 ) << "Pow::Pow copy constructor\n";
        }

    bool isSymetric() const { return false; }

    expression_1_type const& left() const { return _M_expr_1; }
    expression_2_type const& right() const { return _M_expr_2; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                                                typename r_type::value_type>::type value_type;


        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                                  mpl::identity<detail::gmc<0> >,
                                  mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        typedef typename l_type::shape shape;

        struct is_zero { static const bool value = l_type::is_zero::value; };

        tensor( expression_type const& expr, Geo_t const& geom,
                 Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom, fev, feu ),
            _M_right( expr.right(), geom, fev, feu ),
            _M_loc( boost::extents[shape::M][shape::N][_M_gmc->nPoints()]  )
            {
                update( geom );
            }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom, fev ),
            _M_right( expr.right(), geom, fev ),
            _M_loc( boost::extents[shape::M][shape::N][_M_gmc->nPoints()] )
            {
                update( geom );
            }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            _M_gmc( fusion::at_key<key_type>( geom ).get() ),
            _M_left( expr.left(),  geom ),
            _M_right( expr.right(), geom ),
            _M_loc(  boost::extents[shape::M][shape::N][_M_gmc->nPoints()] )
            {
                update( geom );
            }
        template<typename IM>
        void init( IM const& im )
        {
            _M_left.init( im );
            _M_right.init( im );
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
            _M_left.update( geom );
            _M_right.update( geom );

            for( int c1 = 0; c1 < shape::M; ++c1 )
                for( int c2 = 0; c2 < shape::N; ++c2 )
                    for ( int q = 0; q < npts; ++q )
                        {
                            value_type left = _M_left.evalq( c1, c2, q );
                            value_type right = _M_right.evalq( c1, c2, q );
                            _M_loc[c1][c2][q] = std::pow( left, right );
                        }
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            _M_gmc = fusion::at_key<key_type>( geom ).get();
            _M_left.update( geom, face );
            _M_right.update( geom, face );

            for( int c1 = 0; c1 < shape::M; ++c1 )
                for( int c2 = 0; c2 < shape::N; ++c2 )
                    for ( int q = 0; q < _M_gmc->nPoints(); ++q )
                        {
                            value_type left = _M_left.evalq( c1, c2, q );
                            value_type right = _M_right.evalq( c1, c2, q );
                            _M_loc[c1][c2][q] = std::pow( left, right );
                        }
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q  ) const
        {
            return evalq( c1, c2, q );
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Life::detail::ignore_unused_variable_warning(i);
            Life::detail::ignore_unused_variable_warning(j);
            return evalq( c1, c2, q );
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& /*i*/, uint16_type c1, uint16_type c2, uint16_type q  ) const
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
        gmc_ptrtype _M_gmc;
        l_type _M_left;
        r_type _M_right;
        //ublas::vector<double> _M_loc;
        boost::multi_array<value_type,3> _M_loc;
    };

protected:
    Pow() {}

    expression_1_type _M_expr_1;
    expression_2_type _M_expr_2;
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
                            mpl::identity<ExprT1> >::type::type,
          typename mpl::if_<boost::is_arithmetic<ExprT2>,
                            mpl::identity<Cst<ExprT2> >,
                            mpl::identity<ExprT2> >::type::type> >
operator^( ExprT1 const& __e1, ExprT2 const& __e2 )
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
        _M_expr( __expr )
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

        template <class Args> struct sig { typedef value_type type; };

        tensor( this_type const& expr,
                VecGeo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), fusion::at_c<GeoId>(  geom ), fev, feu )
        {}

        tensor( this_type const& expr,
                VecGeo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), fusion::at_c<GeoId>(  geom ), fev )
        {}

        tensor( this_type const& expr,
                VecGeo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), fusion::at_c<GeoId>(  geom ) )
        {}
        template<typename IM>
        void init( IM const& im )
        {
            _M_tensor_expr.init( im );
        }
        void update( VecGeo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( fusion::at_c<GeoId>(  geom ), fev, feu );
        }
        void update( VecGeo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( fusion::at_c<GeoId>(  geom ), fev );
        }
        void update( VecGeo_t const& geom )
        {
            _M_tensor_expr.update( fusion::at_c<GeoId>(  geom ) );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, int q ) const
        {
            return _M_tensor_expr.evalijq( i, j, q );
        }

        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, int q ) const
        {
            return _M_tensor_expr.evaliq( i, q );
        }

        value_type
        evalq( int q, int c ) const
        {
            return _M_tensor_expr.evalq( q, c );
        }

        tensor_expr_type _M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    expression_type const& expression() const { return _M_expr; }

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

    mutable expression_type  _M_expr;
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

    typedef Element element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef GElem<element_type, Type> this_type;
    typedef this_type self_type;

    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename functionspace_type::reference_element_type* fe_ptrtype;
    typedef typename functionspace_type::reference_element_type fe_type;
    typedef typename functionspace_type::value_type value_type;
    static const uint16_type rank = fe_type::rank;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;
    typedef std::map<size_type,std::vector<element_ptrtype> > basis_type;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = (Type==0);
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = (Type==1);
    };


    GElem ( std::map<size_type,std::vector<element_ptrtype> > const& v )
        :
        _M_basis ( v )
    {
        typename basis_type::iterator it = _M_basis.begin();
        typename basis_type::iterator en = _M_basis.end();
        for( ; it != en; ++it )
            for( uint16_type i = 0; i < it->second.size(); ++i )
                it->second[i]->updateGlobalValues();
    }
    GElem( GElem const& op )
        :
        _M_basis ( op._M_basis )
    {

    }

    basis_type const&  basis() const { return _M_basis; }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                                  mpl::identity<detail::gmc<0> >,
                                  mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        typedef typename mpl::if_<mpl::equal_to<mpl::int_<rank>,
                                                mpl::int_<0> >,
                                  mpl::identity<Shape<gmc_type::NDim, Scalar, false> >,
                                  typename mpl::if_<mpl::equal_to<mpl::int_<rank>,
                                                                  mpl::int_<1> >,
                                                    mpl::identity<Shape<gmc_type::NDim, Vectorial, false> >,
                                                    mpl::identity<Shape<gmc_type::NDim, Tensor2, false> > >::type>::type::type shape;
        typedef typename fe_type::PreCompute pc_type;

        typedef typename expression_type::value_type value_type;

        struct is_zero { static const bool value = false; };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            _M_expr( expr ),
            _M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            _M_loc( expr.basis().begin()->second.size() )
            {
                for( uint16_type i = 0; i < _M_loc.size(); ++i )
                    _M_loc[i].resize( boost::extents[nComponents1][nComponents2][_M_pc.nPoints()] );
            }
        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            _M_expr( expr ),
            _M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            _M_loc( expr.basis().begin()->second.size() )
            {
                for( uint16_type i = 0; i < _M_loc.size(); ++i )
                    _M_loc[i].resize( boost::extents[nComponents1][nComponents2][_M_pc.nPoints()] );
            }
        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            _M_expr( expr ),
            _M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            _M_loc( expr.basis().begin()->second.size() )
        {
            for( uint16_type i = 0; i < _M_loc.size(); ++i )
                _M_loc[i].resize( boost::extents[nComponents1][nComponents2][_M_pc.nPoints()] );
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
            //Debug() << "[GElem] updating element " << fusion::at_key<key_type>( geom )->id() << "\n";
            typename basis_type::iterator it = const_cast<basis_type&>(_M_expr.basis()).find( fusion::at_key<key_type>( geom )->id() );
            typename basis_type::iterator en = const_cast<basis_type&>(_M_expr.basis()).end();

            LIFE_ASSERT( it != en )( fusion::at_key<key_type>( geom )->id() ).error ("invalid basis function to integrate" );
            for( uint16_type i = 0; i < _M_loc.size(); ++i )
                {
                    _M_loc[i] = it->second[i]->id( *fusion::at_key<key_type>( geom ), _M_pc );
                }
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( Type == 0 )
                return _M_loc[i](c1,c2,q);
            return _M_loc[j](c1,c2,q);
        }
        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext> ) const
        {
            if ( Type == 0 )
                return _M_loc[i](c1,c2,q);
            return _M_loc[j](c1,c2,q);
        }
        template<typename IndexI>
        value_type
        evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_loc[i](c1,c2,q);
        }
        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {

        }
    private:
        this_type const& _M_expr;
        pc_type _M_pc;
        std::vector<typename element_type::id_type> _M_loc;
    };
private:
    basis_type _M_basis;

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
} // life
#endif /* __Expr_H */
