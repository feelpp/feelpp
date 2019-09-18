/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 22 Feb 2018

 Copyright (C) 2018 Feel++ Consortium

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
#ifndef FEELPP_MODELEXPRESSION_HPP
#define FEELPP_MODELEXPRESSION_HPP 1

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelmodels/modelindexes.hpp>

namespace Feel {

class FEELPP_EXPORT ModelExpression
{
public :

    static const uint16_type expr_order = 2;
    typedef scalar_field_expression<expr_order> expr_scalar_type;
    typedef vector_field_expression<2,1,expr_order> expr_vectorial2_type;
    typedef vector_field_expression<3,1,expr_order> expr_vectorial3_type;
    typedef matrix_field_expression<2,2,expr_order> expr_matrix22_type;
    typedef matrix_field_expression<3,3,expr_order> expr_matrix33_type;

    ModelExpression() = default;
    ModelExpression( ModelExpression const& ) = default;
    ModelExpression( ModelExpression && ) = default;
    ModelExpression& operator=( ModelExpression const& ) = default;
    ModelExpression& operator=( ModelExpression && ) = default;

    bool hasExprScalar() const { return (M_exprScalar)? true : false; }
    bool hasExprVectorial2() const { return (M_exprVectorial2)? true : false; }
    bool hasExprVectorial3() const { return (M_exprVectorial3)? true : false; }
    bool hasExprMatrix22() const { return (M_exprMatrix22)? true : false; }
    bool hasExprMatrix33() const { return (M_exprMatrix33)? true : false; }
    template <int M,int N>
    bool hasExprMatrix() const
    {
        if ( M==2 && N==2 ) return this->hasExprMatrix22();
        else if ( M==3 && N==3 ) return this->hasExprMatrix33();
        else return false;
    }
    template <int M,int N>
    bool hasExpr() const
    {
        if ( M==1 && N == 1 ) return this->hasExprScalar();
        else if ( M == 2 && N == 1 ) return this->hasExprVectorial2();
        else if ( M == 3 && N == 1 ) return this->hasExprVectorial3();
        else if ( M==2 && N==2 ) return this->hasExprMatrix22();
        else if ( M==3 && N==3 ) return this->hasExprMatrix33();
        else return false;
    }

    bool hasAtLeastOneExpr() const { return  this->hasExprScalar() || this->hasExprVectorial2() || this->hasExprVectorial3() || this->hasExprMatrix22() || this->hasExprMatrix33(); }

    bool isScalar() const { return hasExprScalar(); }
    bool isVector() const { return hasExprVectorial2() || hasExprVectorial3(); }
    bool isMatrix() const { return hasExprMatrix22() || hasExprMatrix33(); }

    // TODO : not necessary specific to scalar expression
    bool isConstant() const { return this->hasExprScalar() && this->exprScalar().expression().isConstant(); }
    double value() const { CHECK( this->isConstant() ) << "expression is not constant";return this->exprScalar().evaluate(); }

    expr_scalar_type const& exprScalar() const { CHECK( this->hasExprScalar() ) << "no Scalar expression"; return *M_exprScalar; }
    expr_vectorial2_type const& exprVectorial2() const { CHECK( this->hasExprVectorial2() ) << "no Vectorial2 expression"; return *M_exprVectorial2; }
    expr_vectorial3_type const& exprVectorial3() const { CHECK( this->hasExprVectorial3() ) << "no Vectorial3 expression"; return *M_exprVectorial3; }
    expr_matrix22_type const& exprMatrix22() const { CHECK( this->hasExprMatrix22() ) << "no Matrix22 expression"; return *M_exprMatrix22; }
    expr_matrix33_type const& exprMatrix33() const { CHECK( this->hasExprMatrix33() ) << "no Matrix33 expression"; return *M_exprMatrix33; }

    template <int M=1,int N=1>
        expr_scalar_type const& expr( typename std::enable_if< M==1 && N == 1>::type* = nullptr ) const { return this->exprScalar(); }
    template <int M=1,int N=1>
        expr_vectorial2_type const& expr( typename std::enable_if< M==2 && N == 1>::type* = nullptr ) const { return this->exprVectorial2(); }
    template <int M=1,int N=1>
        expr_vectorial3_type const& expr( typename std::enable_if< M==3 && N == 1>::type* = nullptr ) const { return this->exprVectorial3(); }
    template <int M=1,int N=1>
        expr_matrix22_type const& expr( typename std::enable_if< M==2 && N == 2>::type* = nullptr ) const { return this->exprMatrix22(); }
    template <int M=1,int N=1>
        expr_matrix33_type const& expr( typename std::enable_if< M==3 && N == 3>::type* = nullptr ) const { return this->exprMatrix33(); }

    template <int M,int N>
        expr_matrix22_type const& exprMatrix( typename std::enable_if< M==2 && N == 2>::type* = nullptr ) const { return this->exprMatrix22(); }
    template <int M,int N>
        expr_matrix33_type const& exprMatrix( typename std::enable_if< M==3 && N == 3>::type* = nullptr ) const { return this->exprMatrix33(); }

    void setExprScalar( expr_scalar_type const& expr ) { M_exprScalar = boost::optional<expr_scalar_type>( expr ); }
    void setExprVectorial2( expr_vectorial2_type const& expr ) { M_exprVectorial2 = boost::optional<expr_vectorial2_type>( expr ); }
    void setExprVectorial3( expr_vectorial3_type const& expr ) { M_exprVectorial3 = boost::optional<expr_vectorial3_type>( expr ); }
    void setExprMatrix22( expr_matrix22_type const& expr ) { M_exprMatrix22 = boost::optional<expr_matrix22_type>( expr ); }
    void setExprMatrix33( expr_matrix33_type const& expr ) { M_exprMatrix33 = boost::optional<expr_matrix33_type>( expr ); }

    //! set an expression from a key of a ptree p
    void setExpr( std::string const& key, pt::ptree const& p, WorldComm const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes = ModelIndexes() );
    void setExpr( std::string const& expr, WorldComm const& worldComm = Environment::worldComm(), std::string const& directoryLibExpr = "" );

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        if ( this->hasExprScalar() )
            M_exprScalar->setParameterValues( mp );
        if ( this->hasExprVectorial2() )
            M_exprVectorial2->setParameterValues( mp );
        if ( this->hasExprVectorial3() )
            M_exprVectorial3->setParameterValues( mp );
        if ( this->hasExprMatrix22() )
            M_exprMatrix22->setParameterValues( mp );
        if ( this->hasExprMatrix33() )
            M_exprMatrix33->setParameterValues( mp );
    }

    template<typename ExprT>
    constexpr auto //Expr< GinacExVF<ExprT,expr_order> >
    exprScalar( std::string const& symb, ExprT const& e ) const
    {
        return Feel::vf::expr( this->exprScalar(), symbolExpr(symb, e) );
    }

private :
    boost::optional<expr_scalar_type> M_exprScalar;
    boost::optional<expr_vectorial2_type> M_exprVectorial2;
    boost::optional<expr_vectorial3_type> M_exprVectorial3;
    boost::optional<expr_matrix22_type> M_exprMatrix22;
    boost::optional<expr_matrix33_type> M_exprMatrix33;
};

class FEELPP_EXPORT ModelExpressionScalar : private ModelExpression
{
    typedef ModelExpression super_type;
public :
    static const uint16_type expr_order = super_type::expr_order;
    typedef typename super_type::expr_scalar_type expr_scalar_type;

    ModelExpressionScalar() = default;
    ModelExpressionScalar( ModelExpressionScalar const& ) = default;
    ModelExpressionScalar( ModelExpressionScalar && ) = default;
    ModelExpressionScalar& operator=( ModelExpressionScalar const& ) = default;
    ModelExpressionScalar& operator=( ModelExpressionScalar && ) = default;

    bool isConstant() const { return super_type::isConstant(); }
    double value() const { return super_type::value(); }
    expr_scalar_type const& expr() const { return this->exprScalar(); }
    bool hasExpr() const { return this->hasExprScalar(); }

    void setExpr( expr_scalar_type const& expr ) { this->setExprScalar( expr ); }

    void setParameterValues( std::map<std::string,double> const& mp ) { super_type::setParameterValues( mp ); }

    template<typename ExprT>
    constexpr auto //Expr< GinacExVF<ExprT,expr_order> >
    expr( std::string const& symb, ExprT const& e ) const
    {
        return this->exprScalar( symb,e );
    }

};

#if 0
template <int M,int N>
class ModelExpressionMatrix : private ModelExpression
{
    typedef ModelExpression super_type;
public :
    static const uint16_type expr_order = super_type::expr_order;
    typedef typename mpl::if_<mpl::bool_< M==2 && N==2 >,
                              typename super_type::expr_matrix22_type,
                              typename super_type::expr_matrix33_type>::type expr_matrix_type;

    ModelExpressionMatrix() = default;
    ModelExpressionMatrix( ModelExpressionMatrix const& ) = default;
    ModelExpressionMatrix( ModelExpressionMatrix && ) = default;
    ModelExpressionMatrix& operator=( ModelExpressionMatrix const& ) = default;
    ModelExpressionMatrix& operator=( ModelExpressionMatrix && ) = default;

    //bool isConstant() const { return !this->hasExprScalar(); }
    //double value() const { return super_type::value(); }
    expr_matrix_type const& expr() const { return this->exprMatrix<M,N>(); }

    void setValue( double d ) { super_type::setValue( d ); }
    void setExpr( expr_matrix_type const& expr ) { this->setExprScalar( expr ); }

    void setParameterValues( std::map<std::string,double> const& mp ) { super_type::setParameterValues( mp ); }
#if 0
    template<typename ExprT>
    Expr< GinacExVF<ExprT,expr_order> >
    expr( std::string const& symb, ExprT const& e ) const
    {
        return this->exprScalar( symb,e );
    }
#endif
};
#endif
} // namespace Feel

#endif
