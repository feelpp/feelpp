/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library
 
  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014-2016 Feel++ Consortium

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
   \file trans.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_TRANS_HPP
#define FEELPP_VF_TRANS_HPP 1

namespace Feel
{
namespace vf
{

/*!
  \ingroup vf
  \brief Transpose an scalar, vectorial or matricial expression


  @author Christophe Prud'homme
*/
template <typename ExprT>
class Trans : public ExprDynamicBase
{
  public:
    using super = ExprDynamicBase;
    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };
    template <typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template <typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    using evaluate_type = Eigen::Matrix<value_type,
                                        expression_type::evaluate_type::ColsAtCompileTime,
                                        expression_type::evaluate_type::RowsAtCompileTime >;
    typedef Trans<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Trans( expression_type const& __expr )
        : super( Feel::vf::dynamicContext( __expr ) ),
          M_expr( __expr )
    {
    }
    ~Trans() = default;

    //@}

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            return M_expr.evaluate(p,worldcomm).transpose();
        }
    /** @name Operator overloads
     */
    //@{
    template <typename... TheExpr>
    struct Lambda
    {
        typedef Trans<typename expression_type::template Lambda<TheExpr...>::type> type;
    };
    template <typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_expr( e... ) ); }

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
            return Trans<new_expr_type>( newExpr );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr.hasSymbolDependency( symb, se );
        }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_expr.dependentSymbols( symb, res, se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            auto theDiffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
            using new_expr_type = std::decay_t<decltype(theDiffExpr)>;
            return Trans<new_expr_type>( theDiffExpr );
        }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        static constexpr size_type context = this_type::context;
        using expr_type = typename this_type::expression_type;
        typedef typename Transpose<typename tensor_expr_type::shape>::type shape;

        template <class Args>
        struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensor_expr( expr.expression(), geom )
        {
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse,  expr.expression(), geom, theInitArgs... )
            {}

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
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
        }

        template <typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_tensor_expr.updateContext( ctx... );
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
        template <int PatternContext>
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
    mutable expression_type M_expr;
};

/**
   \ingroup vf

   Transpose an expression \p v. The expression can be scalar (no-op), vectorial
   or matricial.

   \code
   // define a column vector of two element
   auto e = vec(Px(),Py());
   // transpose it and get a row vector expression of two elements
   auto eT = trans(e);
   \endcode

   \return the transposed expression
 */
template <typename ExprT>
inline Expr<Trans<ExprT>>
trans( ExprT v )
{
    typedef Trans<ExprT> trans_t;
    return Expr<trans_t>( trans_t( v ) );
}

} // namespace vf
} // namespace Feel
#endif /* FEELPP_VF_TRANS_HPP */
