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

namespace Feel {

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
       << "]";
    return os;
}

ModelMaterials::ModelMaterials( pt::ptree const& p )  : M_p( p )
{
    setup();
}
ModelMaterial
ModelMaterials::loadMaterial( std::string const& s )
{
    pt::ptree p;
    pt::read_json( s, p );
    return this->getMaterial( p );
}
void
ModelMaterials::setup()
{
    for( auto const& v : M_p )
    {
        LOG(INFO) << "Material Physical/Region :" << v.first  << "\n";
        
        if ( auto fname = v.second.get_optional<std::string>("filename") )
        {
            LOG(INFO) << "  - filename = " << Environment::expand( fname.get() ) << std::endl;
            
            this->insert( std::make_pair( v.first, this->loadMaterial( Environment::expand( fname.get() ) ) ) );
        }
        else
        {
            this->insert( std::make_pair( v.first, this->getMaterial( v.second ) ) );
        }
    }
}
ModelMaterial
ModelMaterials::getMaterial( pt::ptree const& v )
{
    std::string t = v.get<std::string>( "name" );
    LOG(INFO) << "loading material name: " << t << std::endl;
    ModelMaterial m(t, v);
    m.setRho( v.get( "rho", 1.f ) );
    m.setMu( v.get( "mu", 1.f ) );
    m.setCp( v.get( "Cp", 1.f ) );
    m.setCv( v.get( "Cv", 1.f ) );
    m.setTref( v.get( "Tref", 0.f ) );
    m.setBeta( v.get( "beta", 0.f ) );
    m.setK11( v.get( "k11", 1.f ) );
    m.setK12( v.get( "k12", 0.f ) );
    m.setK13( v.get( "k13", 0.f ) );
    m.setK22( v.get( "k22", 1.f ) );
    m.setK23( v.get( "k23", 0.f ) );
    m.setK33( v.get( "k33", 1.f ) );
    m.setE( v.get( "E", 1.f ) );
    m.setNu( v.get( "nu", 1.f ) );
    m.setSigma( v.get( "sigma", 1.f ) );
    m.setC( v.get( "C", 1.f ) );

    LOG(INFO) << "adding material " << m;
    return m;
}
void
ModelMaterials::saveMD(std::ostream &os)
{
  os << "### Materials\n";
  os << "|Material Physical Region Name|Rho|mu|Cp|Cv|k11|k12|k13|k22|k23|k33|Tref|beta|C|YoungModulus|nu|Sigma|\n";
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
       << "|\n";
  os << "\n";
}

}


