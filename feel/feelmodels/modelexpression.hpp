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

namespace Feel {

class FEELPP_EXPORT ModelExpression
{
public :

    static const uint16_type expr_order = 2;
    typedef scalar_field_expression<expr_order> expr_scalar_type;
    typedef vector_field_expression<2,1,expr_order> expr_vectorial2_type;
    typedef vector_field_expression<3,1,expr_order> expr_vectorial3_type;

    ModelExpression() = default;
    ModelExpression( ModelExpression const& ) = default;
    ModelExpression( ModelExpression && ) = default;
    ModelExpression& operator=( ModelExpression const& ) = default;
    ModelExpression& operator=( ModelExpression && ) = default;

    bool hasValue(int c=0) const { return M_values.find(c) != M_values.end(); }
    bool hasExprScalar() const { return (M_exprScalar)? true : false; }
    bool hasExprVectorial2() const { return (M_exprVectorial2)? true : false; }
    bool hasExprVectorial3() const { return (M_exprVectorial3)? true : false; }

    double value( int c=0 ) const { CHECK( this->hasValue(c) ) << "invalid comp"; return M_values.find(c)->second; }
    expr_scalar_type const& exprScalar() const { CHECK( this->hasExprScalar() ) << "no Scalar expression"; return *M_exprScalar; }
    expr_vectorial2_type const& exprVectorial2() const { CHECK( this->hasExprVectorial2() ) << "no Vectorial2 expression"; return *M_exprVectorial2; }
    expr_vectorial3_type const& exprVectorial3() const { CHECK( this->hasExprVectorial3() ) << "no Vectorial3 expression"; return *M_exprVectorial3; }

    void setValue( double val, int c = 0 ) { M_values[c] = val; }
    void setExprScalar( expr_scalar_type const& expr ) { M_exprScalar = boost::optional<expr_scalar_type>( expr ); }
    void setExprVectorial2( expr_vectorial2_type const& expr ) { M_exprVectorial2 = boost::optional<expr_vectorial2_type>( expr ); }
    void setExprVectorial3( expr_vectorial3_type const& expr ) { M_exprVectorial3 = boost::optional<expr_vectorial3_type>( expr ); }

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        if ( this->hasExprScalar() )
            M_exprScalar->setParameterValues( mp );
        if ( this->hasExprVectorial2() )
            M_exprVectorial2->setParameterValues( mp );
        if ( this->hasExprVectorial3() )
            M_exprVectorial3->setParameterValues( mp );
    }

private :
    std::map<int,double> M_values;
    boost::optional<expr_scalar_type> M_exprScalar;
    boost::optional<expr_vectorial2_type> M_exprVectorial2;
    boost::optional<expr_vectorial3_type> M_exprVectorial3;
};

class FEELPP_EXPORT ModelExpressionScalar : private ModelExpression
{
    typedef ModelExpression super_type;
public :
    typedef typename super_type::expr_scalar_type expr_scalar_type;

    ModelExpressionScalar() = default;
    ModelExpressionScalar( ModelExpressionScalar const& ) = default;
    ModelExpressionScalar( ModelExpressionScalar && ) = default;
    ModelExpressionScalar& operator=( ModelExpressionScalar const& ) = default;
    ModelExpressionScalar& operator=( ModelExpressionScalar && ) = default;

    bool isConstant() const { return !this->hasExprScalar(); }
    double value() const { return super_type::value(); }
    expr_scalar_type const& expr() const { return this->exprScalar(); }

    void setValue( double d ) { super_type::setValue( d ); }
    void setExpr( expr_scalar_type const& expr ) { this->setExprScalar( expr ); }

    void setParameterValues( std::map<std::string,double> const& mp ) { super_type::setParameterValues( mp ); }
};

} // namespace Feel

#endif
