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

#include <feel/feelmodels/modelfunctions.hpp>

namespace Feel {


ModelFunction::ModelFunction( std::string const& name, std::string const& expression,
                              std::string const& dirLibExpr, WorldComm const& world )
    :
    M_name( name ),
    M_exprString( expression )
{
    auto parseExpr = GiNaC::parse(expression);

    auto ginacEvalm = parseExpr.first.evalm();
    bool isLst = GiNaC::is_a<GiNaC::lst>(  ginacEvalm );
    int nComp = 1;
    if ( isLst )
        nComp = ginacEvalm.nops();

    //M_exprScalar.emplace( expr<expr_order>( expression,"",world,dirLibExpr ) );
    if ( nComp == 1 )
        M_exprScalar =  boost::optional<expr_scalar_type>( expr<expr_order>( expression,"",world,dirLibExpr ) );
    else if ( nComp == 2 )
        M_exprVectorial2 = boost::optional<expr_vectorial2_type>( expr<2,1,expr_order>( expression,"",world,dirLibExpr ) );
    else if ( nComp == 3 )
        M_exprVectorial3 = boost::optional<expr_vectorial3_type>( expr<3,1,expr_order>( expression,"",world,dirLibExpr ) );
}

bool
ModelFunction::hasSymbol( std::string const& symb ) const
{
    if ( this->isScalar() )
        return M_exprScalar->expression().hasSymbol( symb );
    else if ( this->isVectorial2() )
        return M_exprVectorial2->expression().hasSymbol( symb );
    else if ( this->isVectorial3() )
        return M_exprVectorial3->expression().hasSymbol( symb );
    return false;
}

void
ModelFunction::setParameterValues( std::map<std::string,double> const& mp )
{
    if ( this->isScalar() )
        M_exprScalar->setParameterValues( mp );
    else if ( this->isVectorial2() )
        M_exprVectorial2->setParameterValues( mp );
    else if ( this->isVectorial3() )
        M_exprVectorial3->setParameterValues( mp );
}


ModelFunctions::ModelFunctions( WorldComm const& world )
    :
    M_worldComm( world )
{}

ModelFunctions::~ModelFunctions()
{}

void
ModelFunctions::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}


void
ModelFunctions::setup()
{
    for( auto const& v : M_p )
    {
        std::string funcName = v.first;
        if ( funcName.empty() )
        {
            LOG(WARNING) << "ignore function because no name given";
            continue;
        }

        try
        {
            std::string exprStr = v.second.get<std::string>("expr");
            this->operator[](funcName) = ModelFunction( funcName, exprStr, M_directoryLibExpr, M_worldComm );
        }
        catch( ... )
        {
            LOG(WARNING) << "ignore function " << funcName << " because no expression given";
        }
    }
}


void
ModelFunctions::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto & p : *this )
        p.second.setParameterValues( mp );
}




}
