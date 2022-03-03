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

#include <feel/feelmodels/modelexpression.hpp>

namespace Feel {

void
ModelExpression::setExpr( std::string const& key, pt::ptree const& p, WorldComm const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    if( boost::optional<double> itvald = p.get_optional<double>( key ) )
    {
        double val = *itvald;
        VLOG(1) << "expr at key " << key << " is constant : " << val;
        this->setExprScalar( Feel::vf::expr( val ) );
    }
    else if( boost::optional<std::string> itvals = p.get_optional<std::string>( key ) )
        this->setExpr( indexes.replace( *itvals ), worldComm, directoryLibExpr );
}
void
ModelExpression::setExpr( nl::json const& jarg, WorldComm const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    if ( jarg.is_number() )
    {
        double val = jarg.get<double>();
        this->setExprScalar( Feel::vf::expr( val ) );
    }
    else if ( jarg.is_string() )
    {
        this->setExpr( indexes.replace( jarg.get<std::string>() ), worldComm, directoryLibExpr );
    }
    else
        CHECK( false ) << "invalid json arg";
}

void
ModelExpression::setExpr( std::string const& feelExprString, WorldComm const& worldComm, std::string const& directoryLibExpr )
{
    auto parseExpr = GiNaC::parse( feelExprString );
    auto const& exprSymbols = parseExpr.second;
    auto ginacEvalm = parseExpr.first.evalm();

    bool isLst = GiNaC::is_a<GiNaC::lst>( ginacEvalm );
    int nComp = 1;
    if ( isLst )
        nComp = ginacEvalm.nops();

    VLOG(1) << "expr " << feelExprString << " build symbolic expr with nComp=" << nComp;
    if ( nComp == 1 )
        this->setExprScalar( Feel::vf::expr<expr_order>( feelExprString,"",worldComm,directoryLibExpr ) );
    else if ( nComp == 2 )
        this->setExprVectorial2( Feel::vf::expr<2,1,expr_order>( feelExprString,"",worldComm,directoryLibExpr ) );
    else if ( nComp == 3 )
        this->setExprVectorial3( Feel::vf::expr<3,1,expr_order>( feelExprString,"",worldComm,directoryLibExpr ) );
    else if ( nComp == 4 )
        this->setExprMatrix22( Feel::vf::expr<2,2,expr_order>( feelExprString,"",worldComm,directoryLibExpr ) );
    else if ( nComp == 9 )
        this->setExprMatrix33( Feel::vf::expr<3,3,expr_order>( feelExprString,"",worldComm,directoryLibExpr ) );
}


} // namespace Feel
