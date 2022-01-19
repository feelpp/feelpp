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

//#include <iostream>
#include <regex>
//#include <boost/property_tree/json_parser.hpp>
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
ModelMaterial::ModelMaterial( std::string const& name, nl::json const& jarg, worldcomm_ptr_t const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes )
    :
    super( worldComm ),
    M_name( name ),
    M_directoryLibExpr( directoryLibExpr ),
    M_meshMarkers( name )
{

    std::optional<nl::json> j_fileLoaded;
    if ( jarg.contains("filename") )
    {
        auto const& j_filename = jarg.at("filename");
        if ( j_filename.is_string() )
        {
            std::string fname = Environment::expand( indexes.replace( j_filename.get<std::string>() ) );
            LOG(INFO) << "  - filename = " << fname << std::endl;
            pt::ptree pf;
            auto json_str_wo_comments = removeComments(readFromFile(fname));
            std::istringstream istr( json_str_wo_comments );
            j_fileLoaded = nl::json{};
            istr >> *j_fileLoaded;
        }
    }

    if ( jarg.contains("markers") )
        M_meshMarkers.setup( jarg.at("markers"), indexes);

    if ( jarg.contains("physics") )
    {
        auto const& j_physics = jarg.at("physics");
        if ( j_physics.is_string() )
            M_physics.insert( indexes.replace( j_physics.get<std::string>() ) );
        else if ( j_physics.is_array() )
        {
            for ( auto const& [j_physicskey,j_physicsval] : j_physics.items() )
                if ( j_physicsval.is_string() )
                    M_physics.insert( indexes.replace( j_physicsval.get<std::string>() ) );
        }
    }


    std::map<std::string,double> constantMaterialProperty;
    auto addMatPropLambda = [&constantMaterialProperty,&indexes,this]( nl::json const& jarg ) {
                                std::regex index_pattern { "index[1-9][0-9]*" }; // motif
                                for ( auto const& [jargkey,jargval] : jarg.items() )
                                {
                                    std::string const& propName = indexes.replace( jargkey );
                                    if ( propName == "markers" || propName == "physics" || propName == "name" || propName == "filename" )
                                        continue;
                                    if ( jargval.empty() )
                                        continue;
                                    if ( std::regex_match( propName, index_pattern ) )
                                        continue;

                                    VLOG(1) << "propName: " << propName;

                                    this->setProperty( propName,jargval,indexes );
                                    if ( this->hasPropertyConstant( propName ) && this->hasPropertyExprScalar( propName ) )
                                    {
                                        VLOG(1) << "key: " << propName << " value: " << this->property( propName ).value() << "  constant prop";
                                        constantMaterialProperty[propName] = this->property( propName ).value();
                                    }
                                }
                            };

    if ( j_fileLoaded )
        addMatPropLambda( *j_fileLoaded );
    addMatPropLambda( jarg );
    this->setParameterValues( constantMaterialProperty );
}

bool
ModelMaterial::hasProperty( std::string const& prop ) const
{
    auto itFindProp = M_materialProperties.find( prop );
    if ( itFindProp == M_materialProperties.end() )
        return false;
    auto const& matProp = itFindProp->second;
    return matProp.hasAtLeastOneExpr();
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
ModelMaterial::setProperty( std::string const& property, nl::json const& jarg, ModelIndexes const& indexes )
{
    try
    {
        mat_property_expr_type mexpr;
        mexpr.setExpr( jarg,this->worldComm(),M_directoryLibExpr,indexes );
        M_materialProperties.emplace( std::make_pair( property, std::move( mexpr ) ) );
    } catch (std::exception &p) {
        LOG(WARNING) << p.what() << std::endl;
        return;
    }
}

void
ModelMaterial::setProperty( std::string const& property, std::string const& e )
{
    M_materialProperties[property] = mat_property_expr_type();
    try
    {
        M_materialProperties[property].setExpr( e,this->worldComm(),M_directoryLibExpr );
    } catch (std::exception &p) {
        LOG(WARNING) << p.what() << std::endl;
        M_materialProperties.erase(property);
        return;
    }

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
#if 0
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
#endif
    return os;
}

ModelMaterials::ModelMaterials( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}

void
ModelMaterials::setup()
{
    for ( auto const& [jargkey,jargval] : M_p.items() )
    {
        auto indexesAllCases = ModelIndexes::generateAllCases( jargval );
        for ( auto const& indexes : indexesAllCases )
        {
            std::string const& name = indexes.replace( jargkey );
            LOG(INFO) << "Material Physical/Region :" << name  << "\n";
            this->insert( std::make_pair( name, ModelMaterial( name, jargval, this->worldCommPtr(), M_directoryLibExpr, indexes ) ) );
        }
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
#if 0
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
#endif
}

ModelMaterial const&
ModelMaterials::material( std::string const& m ) const
{
    auto it = this->find( m );
    if ( it == this->end() )
        throw std::invalid_argument( std::string("ModelMaterial: Invalid material name ") + m );
    return it->second;
}

std::map<std::string,ModelMaterial> ModelMaterials::materialWithPhysic(std::string const& physic) const
{
    std::map<std::string,ModelMaterial> mat;
    std::copy_if(this->begin(),this->end(),std::inserter(mat,mat.begin()),
                 [physic](std::pair<std::string,ModelMaterial> const& mp)
                 { return mp.second.hasPhysics(physic); } );
    return mat;
}

std::map<std::string,ModelMaterial> ModelMaterials::materialWithPhysic(std::vector<std::string> const& physics) const
{
    std::map<std::string,ModelMaterial> mat;
    std::copy_if(this->begin(),this->end(),std::inserter(mat,mat.begin()),
                 [physics](std::pair<std::string,ModelMaterial> const& mp)
                 {
                     bool b = false;
                     for( auto const& p : physics )
                         b = b || mp.second.hasPhysics(p);
                     return b;
                 });
    return mat;
}

std::set<std::string> ModelMaterials::markersWithPhysic(std::string const& physic) const
{
    std::set<std::string> markers;
    auto mat = this->materialWithPhysic(physic);
    for( auto const& m : mat )
        markers.insert(m.second.meshMarkers().begin(), m.second.meshMarkers().end());
    return markers;
}

std::set<std::string> ModelMaterials::markersWithPhysic(std::vector<std::string> const& physic) const
{
    std::set<std::string> markers;
    auto mat = this->materialWithPhysic(physic);
    for( auto const& m : mat )
        markers.insert(m.second.meshMarkers().begin(), m.second.meshMarkers().end());
    return markers;
}

}
