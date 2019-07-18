/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Mar 2015

 Copyright (C) 2015 Feel++ Consortium

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
#include <iostream>
#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>

#include <feel/feelmodels/modelmaterials.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>

namespace Feel {

ModelMaterial::ModelMaterial( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}
ModelMaterial::ModelMaterial( std::string const& name, pt::ptree const& p, worldcomm_ptr_t const& worldComm, std::string const& directoryLibExpr )
    :
    super( worldComm ),
    M_name( name ),
    M_p( p ),
    M_directoryLibExpr( directoryLibExpr ),
    M_meshMarkers( name )
{
    if ( auto fnameOpt = p.get_optional<std::string>("filename") )
    {
        std::string fname = Environment::expand( fnameOpt.get() );
        LOG(INFO) << "  - filename = " << fname << std::endl;
        pt::ptree pf;
        auto json_str_wo_comments = removeComments(readFromFile(fname));
        std::istringstream istr( json_str_wo_comments );
        pt::read_json( istr, pf );
        for ( auto const& subpf : pf )
        {
            if ( M_p.count( subpf.first ) == 0 )
            {
                M_p.add_child(subpf.first,subpf.second );
            }
        }
    }


    if ( auto markers = M_p.get_child_optional("markers") )
        M_meshMarkers.setPTree(*markers);

    if ( auto physics = M_p.get_child_optional("physics") )
    {
        for( auto const& item : M_p.get_child("physics") )
            M_physics.insert(item.second.template get_value<std::string>());
        if( M_physics.empty() )
            M_physics.insert(M_p.get<std::string>("physics") );
    }

    for( auto const& [k,v] : M_p )
    {
        if ( (k!= "markers") &&  (k!= "physics") && (k!= "name") && (k!= "filename")
             && v.empty() && !v.data().empty() )
        {
            this->setProperty( k,M_p );
        }
    }
#if 0
    std::set<std::string> matProperties = { "rho","mu","Cp","Cv","Tref","beta",
                                            "k","k11","k12","k13","k22","k23","k33",
                                            "E","nu","sigma","C","Cs","Cl","L",
                                            "Ks","Kl","Tsol","Tliq" };
    for ( std::string const& prop : matProperties )
        this->setProperty( prop,M_p );
#endif
}

bool
ModelMaterial::hasProperty( std::string const& prop ) const
{
    auto p = M_p.get_child_optional( prop );
    if( !p ) // child is missing
        return false;
    else
        return true;
}

bool
ModelMaterial::hasPropertyConstant( std::string const& prop ) const
{
    auto itFindProp = M_materialProperties.find( prop );
    if ( itFindProp == M_materialProperties.end() )
        return false;
    auto const& matProp = itFindProp->second;
    return matProp.isConstant();//hasValue();
}
bool
ModelMaterial::hasPropertyExprScalar( std::string const& prop ) const
{
    auto itFindProp = M_materialProperties.find( prop );
    if ( itFindProp == M_materialProperties.end() )
        return false;
    auto const& matProp = itFindProp->second;
    return matProp.hasExprScalar();
}
bool
ModelMaterial::hasPropertyExprVectorial2( std::string const& prop ) const
{
    auto itFindProp = M_materialProperties.find( prop );
    if ( itFindProp == M_materialProperties.end() )
        return false;
    auto const& matProp = itFindProp->second;
    return matProp.hasExprVectorial2();
}
bool
ModelMaterial::hasPropertyExprVectorial3( std::string const& prop ) const
{
    auto itFindProp = M_materialProperties.find( prop );
    if ( itFindProp == M_materialProperties.end() )
        return false;
    auto const& matProp = itFindProp->second;
    return matProp.hasExprVectorial3();
}
ModelMaterial::mat_property_expr_type const&
ModelMaterial::property( std::string const& prop ) const
{
    CHECK( this->hasProperty( prop ) ) << "no prop";
    return M_materialProperties.find( prop )->second;
}
double
ModelMaterial::propertyConstant( std::string const& prop ) const
{
    if ( this->hasPropertyConstant( prop ) )
        return M_materialProperties.find( prop )->second.value();
    else
        return 0;
}
ModelMaterial::expr_scalar_type const&
ModelMaterial::propertyExprScalar( std::string const& prop ) const
{
    CHECK( this->hasPropertyExprScalar( prop ) ) << "no scalar expr";
    return M_materialProperties.find( prop )->second.exprScalar();
}
ModelMaterial::expr_vectorial2_type const&
ModelMaterial::propertyExprVectorial2( std::string const& prop ) const
{
    CHECK( this->hasPropertyExprVectorial2( prop ) ) << "no vectorial2 expr";
    return M_materialProperties.find( prop )->second.exprVectorial2();
}
ModelMaterial::expr_vectorial3_type const&
ModelMaterial::propertyExprVectorial3( std::string const& prop ) const
{
    CHECK( this->hasPropertyExprVectorial3( prop ) ) << "no vectorial3 expr";
    return M_materialProperties.find( prop )->second.exprVectorial3();
}

void
ModelMaterial::setProperty( std::string const& property, pt::ptree const& p )
{
    M_materialProperties[property] = mat_property_expr_type();
    M_materialProperties[property].setExpr( property,p,this->worldComm(),M_directoryLibExpr );
}

void
ModelMaterial::setProperty( std::string const& property, std::string const& e )
{
    M_materialProperties[property] = mat_property_expr_type();
    M_materialProperties[property].setExpr( e,this->worldComm(),M_directoryLibExpr );
}

void
ModelMaterial::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & matPropPair : M_materialProperties )
    {
        matPropPair.second.setParameterValues( mp );
    }
}


std::ostream& operator<<( std::ostream& os, ModelMaterial const& m )
{
    os << "Material " << m.name()
       << "[ rho: " << m.rho()
       << ", mu: " << m.mu()
       << ", Cp: " << m.Cp()
       << ", Cv: " << m.Cv()
       << ", Tref: " << m.Tref()
       << ", beta: " << m.beta()
       << ", k11: " << m.k11()
       << ", k12: " << m.k12()
       << ", k13: " << m.k13()
       << ", k22: " << m.k22()
       << ", k23: " << m.k23()
       << ", k33: " << m.k33()
       << ", E: " << m.E()
       << ", nu: " << m.nu()
       << ", sigma: " << m.sigma()
       << ", Cs: " <<  m.Cs()
       << ", Cl: " <<  m.Cl()
       << ", L: " <<  m.L()
       << ", Ks: " <<  m.Ks()
       << ", Kl: " <<  m.Kl()
       << ", Tsol: " <<  m.Tsol()
       << ", Tliq: " <<  m.Tliq()
       << "]";
    return os;
}

ModelMaterials::ModelMaterials( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}

ModelMaterials::ModelMaterials( pt::ptree const& p, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    M_p( p )
{
    setup();
}

void
ModelMaterials::setup()
{
    for( auto const& v : M_p )
    {
        LOG(INFO) << "Material Physical/Region :" << v.first  << "\n";
        std::string name = v.first;
        this->insert( std::make_pair( name, ModelMaterial( name, v.second, this->worldCommPtr(), M_directoryLibExpr ) ) );
    }
}

void
ModelMaterials::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto & mat : *this )
        mat.second.setParameterValues( mp );
}

void
ModelMaterials::saveMD(std::ostream &os)
{
  os << "### Materials\n";
  os << "|Material Physical Region Name|Rho|mu|Cp|Cv|k11|k12|k13|k22|k23|k33|Tref|beta|C|YoungModulus|nu|Sigma|Cs|Cl|L|Ks|Kl|Tsol|Tliq|\n";
  os << "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n";
  for(auto it = this->begin(); it!= this->end(); it++ )
    os << "|" << it->first
       << "|" << it->second.name()
       << "|" << it->second.rho()
       << "|" << it->second.mu()
       << "|" << it->second.Cp()
       << "|" << it->second.Cv()
       << "|" << it->second.Tref()
       << "|" << it->second.beta()
       << "|" << it->second.k11()
       << "|" << it->second.k12()
       << "|" << it->second.k13()
       << "|" << it->second.k22()
       << "|" << it->second.k23()
       << "|" << it->second.k33()
       << "|" << it->second.E()
       << "|" << it->second.nu()
       << "|" << it->second.sigma()
       << "|" << it->second.Cs()
       << "|" << it->second.Cl()
       << "|" << it->second.L()
       << "|" << it->second.Ks()
       << "|" << it->second.Kl()
       << "|" << it->second.Tsol()
       << "|" << it->second.Tliq()
       << "|\n";
  os << "\n";
}

}
