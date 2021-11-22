/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-03-05

  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file sym.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-03-05
 */
#ifndef __Sym_H
#define __Sym_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Sym
 * \brief symmetric  or antisymmetric part  of a matricial expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprT, int Part = 1>
class Sym
{
public:

    static const size_type context = ExprT::context;
    static const bool is_symetric = Part;
    static const bool is_terminal = false;

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
    template<typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    using evaluate_type = typename expression_type::evaluate_type;
    typedef Sym<ExprT,Part> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Sym( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    Sym( Sym const & te )
        :
        M_expr( te.M_expr )
    {}
    ~Sym()
    {}

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
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            auto a = M_expr.evaluate(p,worldcomm);
            auto at = a.transpose();
            if constexpr ( Part == 1 )
                return 0.5*( a+at );

            else
                return 0.5*( a-at );
        }

    void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            M_expr.setParameterValues( mp );
        }
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
             M_expr.updateParameterValues( pv );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newExpr = M_expr.applySymbolsExpr( se );
            using new_expr_type = std::decay_t<decltype(newExpr)>;
            return Sym<new_expr_type,Part>( newExpr );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr.hasSymbolDependency( symb, se );
        }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_expr.dependentSymbols( symb,res,se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        CHECK( false ) << "TODO";
        return *this;
    }


    //@}

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
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
            {}

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
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type a = M_tensor_expr.evalijq( i, j, c1, c2, q );
            value_type at = M_tensor_expr.evalijq( i, j, c2, c1, q );

            if ( Part == 1 )
                return 0.5*( a+at );

            else
                return 0.5*( a-at );
        }
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type a = M_tensor_expr.evaliq( i, c1, c2, q );
            value_type at = M_tensor_expr.evaliq( i, c2, c1, q );

            if ( Part == 1 )
                return 0.5*( a+at );

            else
                return 0.5*( a-at );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type a = M_tensor_expr.evalq( c1, c2, q );
            value_type at = M_tensor_expr.evalq( c2, c1, q );

            if ( Part == 1 )
                return 0.5*( a+at );

            else
                return 0.5*( a-at );
        }

    private:
        tensor_expr_type M_tensor_expr;
    };

private:
    mutable expression_type  M_expr;
};
/// \endcond

/**
 * \brief symetric part of a matricial expression
 */
template<typename ExprT>
inline
Expr< Sym<ExprT,1> >
sym( ExprT v )
{
    typedef Sym<ExprT,1> sym_t;
    return Expr< sym_t >(  sym_t( v ) );
}

/**
 * \brief symetric part of a matricial expression
 */
template<typename ExprT>
inline
Expr< Sym<ExprT,0> >
antisym( ExprT v )
{
    typedef Sym<ExprT,0> sym_t;
    return Expr< sym_t >(  sym_t( v ) );
}

}
}
#endif /* __Unsym_H */
