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

#include <regex>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>

#include <feel/feelmodels/modelmaterials.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>

namespace Feel {


void
ModelMaterialProperty::setup( nl::json const& jarg, worldcomm_t const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    if ( jarg.is_string() || jarg.is_number() )
        M_mexpr.setExpr( jarg,worldComm,directoryLibExpr,indexes );
    else if ( jarg.is_object() )
    {
        if ( jarg.contains("expr") )
            M_mexpr.setExpr( jarg.at("expr"),worldComm,directoryLibExpr,indexes );

        if ( jarg.contains("data") )
            M_data = jarg.at("data");
    }
}

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
    std::regex index_pattern { "index[1-9][0-9]*" }; // motif

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

    if ( jarg.contains("laws") )
    {
        auto const& j_laws = jarg.at("laws");
        if ( j_laws.is_object() )
        {
            for ( auto const& [jlawkey, jlawval]: j_laws.items() )
            {
                this->setLaw( jlawkey, jlawval );
            }
        }
    }

    if ( jarg.contains("submaterials") )
    {
        auto const& j_submat = jarg.at("submaterials");
        if ( j_submat.is_object() )
        {
            for ( auto const& [jargkey,jargval]: j_submat.items() )
            {
                CHECK( jargval.is_object() ) << "submaterial properties must be an object";

                std::string subMatName = indexes.replace( jargkey );
                for ( auto const& [jpropkey,jpropval]: jargval.items() )
                {
                    std::string const& propName = indexes.replace( jpropkey );
                    if ( jpropval.empty() || std::regex_match( propName, index_pattern ) )
                        continue;
                    VLOG(1) << "subPropName: " << propName;
                    this->setSubMaterialProperty( subMatName, propName, jpropval, indexes );
                }
            }
        }
    }


    std::map<std::string,double> constantMaterialProperty;
    auto addMatPropLambda = [&constantMaterialProperty,&indexes,&index_pattern,this]( nl::json const& jarg ) {
        for ( auto const& [jargkey,jargval] : jarg.items() )
        {
            std::string const& propName = indexes.replace( jargkey );
            if ( propName == "markers" || propName == "physics" || propName == "name" || propName == "filename" || propName == "submaterials" )
                continue;
            if ( jargval.empty() )
                continue;
            if ( std::regex_match( propName, index_pattern ) )
                continue;

            VLOG(1) << "propName: " << propName;

            this->setProperty( propName,jargval,indexes );
#if 0 // VINCENT CHECK IF REALLY NECESSARY
            if ( this->hasPropertyConstant( propName ) && this->hasPropertyExprScalar( propName ) )
            {
                VLOG(1) << "key: " << propName << " value: " << this->property( propName ).value() << "  constant prop";
                constantMaterialProperty[propName] = this->property( propName ).value();
            }
#endif
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
    return M_materialProperties.find( prop ) != M_materialProperties.end();
}


ModelMaterial::material_property_type const&
ModelMaterial::property( std::string const& prop ) const
{
    if ( !this->hasProperty( prop ) )
        throw std::out_of_range( fmt::format("property {} is not registered",prop) );
    return M_materialProperties.at( prop );
}

void
ModelMaterial::setProperty( std::string const& property, nl::json const& jarg, ModelIndexes const& indexes )
{
    M_materialProperties[property].setup( jarg,this->worldComm(),M_directoryLibExpr,indexes );
}

void
ModelMaterial::setProperty( std::string const& property, std::string const& e )
{
    nl::json j_matprop( { { "expr",e } } );
    this->setProperty( property, j_matprop );
}

bool
ModelMaterial::hasLaw( std::string const& law ) const
{
    return M_materialLaws.find( law ) != M_materialLaws.end();
}

ModelMaterial::material_law_ptrtype const&
ModelMaterial::law( std::string const& law ) const
{
    if ( !this->hasLaw( law ) )
        throw std::out_of_range( fmt::format("law {} is not registered", law) );
    return M_materialLaws.at( law );
}

void
ModelMaterial::setLaw( std::string const& law, nl::json const& jarg )
{
    auto [it, inserted] = M_materialLaws.try_emplace( law, ModelMaterialLawFactory::instance().createObject( law ) );
    it->second->setup( jarg );
}

void
ModelMaterial::setLaw( std::string const& law, std::string const& e )
{
    nl::json j_matlaw( { e } );
    this->setLaw( law, j_matlaw );
}

void
ModelMaterial::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & [matName,matProp] : M_materialProperties )
        matProp.setParameterValues( mp );
}

bool
ModelMaterial::hasSubMaterial( std::string const& subMat ) const
{
    return M_subMaterialProperties.find( subMat ) != M_subMaterialProperties.end();
}

ModelMaterial::material_properties_type const&
ModelMaterial::subMaterialProperties( std::string const& subMat ) const
{
    if ( !this->hasSubMaterial( subMat ) )
        throw std::out_of_range( fmt::format("sub-material {} is not registered", subMat) );
    return M_subMaterialProperties.at( subMat );
}

void
ModelMaterial::setSubMaterialProperty( std::string const& submat, std::string const& property, nl::json const& jarg, ModelIndexes const& indexes )
{
    M_subMaterialProperties[submat][property].setup( jarg,this->worldComm(),M_directoryLibExpr,indexes );
}

void
ModelMaterial::setSubMaterialProperty( std::string const& submat, std::string const& property, std::string const& e )
{
    nl::json j_matprop( { { "expr",e } } );
    this->setSubMaterialProperty( submat, property, j_matprop );
}

bool
ModelMaterial::hasSubMaterialProperty( std::string const& submat, std::string const& prop ) const
{
    return this->hasSubMaterial(submat) 
        && M_subMaterialProperties.at(submat).find( prop ) != M_subMaterialProperties.at(submat).end();
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
