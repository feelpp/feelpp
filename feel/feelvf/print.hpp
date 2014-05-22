/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014 Feel++ Consortium

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
   \file print.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_PRINT_HPP
#define FEELPP_VF_PRINT_HPP 1

namespace Feel
{
namespace vf
{
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
    typedef value_type evaluate_type;
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

/**
   \ingroup vf

   print an expression during evaluation. This function is typically used for
   debugging purposes and it allows to check the numerical results of
   expressions and sub-expressions.

   The following code prints the intermediary results of the evaluation of
   \f$\int_\Omega x y\f$

   \code
   auto e = Px()*Py();
   std::cout << integrate( _range=elements(mesh), _expr=print(e,"e=") ).evaluate();
   \endcode

   \return the expression that was printed.
 */
template<typename ExprT>
inline
Expr< PrintExpr<ExprT> >
print( ExprT v, std::string t="" )
{
    typedef PrintExpr<ExprT> print_t;
    return Expr< print_t >(  print_t( v, t ) );
}

}
}
#endif /* FEELPP_VF_PRINT_HPP */
