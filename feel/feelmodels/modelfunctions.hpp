/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 28/01/2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_MODELFUNCTIONS_HPP
#define FEELPP_MODELFUNCTIONS_HPP 1

#include <map>

#include <boost/property_tree/ptree.hpp>

#include <feel/feelvf/ginac.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

struct FEELPP_EXPORT ModelFunction
{
    static const uint16_type expr_order = 2;
    typedef scalar_field_expression<expr_order> expr_scalar_type;
    typedef vector_field_expression<2,1,expr_order> expr_vectorial2_type;
    typedef vector_field_expression<3,1,expr_order> expr_vectorial3_type;

    ModelFunction() = default;
    ModelFunction( ModelFunction const& ) = default;
    ModelFunction( ModelFunction&& ) = default;
    ModelFunction& operator=( ModelFunction const& ) = default;
    ModelFunction& operator=( ModelFunction && ) = default;
    ModelFunction( std::string const& name, std::string const& expression,
                   std::string const& dirLibExpr = "", WorldComm const& world = Environment::worldComm() );

    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }
    std::string const& expressionString() const { return M_exprString; }

    bool isScalar() const { return M_exprScalar.get_ptr() != 0; }
    bool isVectorial2() const { return M_exprVectorial2.get_ptr() != 0; }
    bool isVectorial3() const { return M_exprVectorial3.get_ptr() != 0; }
    bool hasSymbol( std::string const& symb ) const;

    void setParameterValues( std::map<std::string,double> const& mp );

    expr_scalar_type const& expressionScalar() const { CHECK( this->isScalar() ) << "no scalar expression"; return *M_exprScalar; }
    expr_vectorial2_type const& expressionVectorial2() const { CHECK( this->isVectorial2() ) << "no vectorial2 expression"; return *M_exprVectorial2; }
    expr_vectorial3_type const& expressionVectorial3() const { CHECK( this->isVectorial3() ) << "no vectorial3 expression"; return *M_exprVectorial3; }

private:
    std::string M_name;
    std::string M_exprString;
    boost::optional<expr_scalar_type> M_exprScalar;
    boost::optional<expr_vectorial2_type> M_exprVectorial2;
    boost::optional<expr_vectorial3_type> M_exprVectorial3;

};
class FEELPP_EXPORT ModelFunctions: public std::map<std::string,ModelFunction>
{
public:
    ModelFunctions( WorldComm const& world = Environment::worldComm() );
    ModelFunctions( ModelFunctions const& ) = default;
    virtual ~ModelFunctions();
    void setPTree( pt::ptree const& _p );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setParameterValues( std::map<std::string,double> const& mp );

private:
    void setup();
private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
};


}
#endif
