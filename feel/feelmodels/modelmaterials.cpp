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
#include <feel/feelcore/feel.hpp>


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
       << "]";
    return os;
}
ModelMaterials::ModelMaterials( pt::ptree const& p )  : M_p( p )
{
    setup();
}
void
ModelMaterials::setup()
{
    for( auto const& v : M_p )
    {
        
        std::cout << "Material :" << v.first  << "\n";
        std::string t = v.first; // parameter name
        auto f= v.second;
        {
            ModelMaterial m(t);
            m.setRho( f.get( "rho", 1.f ) );
            m.setMu( f.get( "mu", 1.f ) );
            m.setCp( f.get( "Cp", 1.f ) );
            m.setCv( f.get( "Cv", 1.f ) );
            m.setTref( f.get( "Tref", 0.f ) );
            m.setBeta( f.get( "beta", 0.f ) );
            m.setK11( f.get( "k11", 1.f ) );
            m.setK12( f.get( "k12", 0.f ) );
            m.setK13( f.get( "k13", 0.f ) );
            m.setK22( f.get( "k22", 1.f ) );
            m.setK23( f.get( "k23", 0.f ) );
            m.setK33( f.get( "k33", 1.f ) );
            
            LOG(INFO) << "adding material " << m;
            this->push_back( std::move(m) );
        }
    }
}    
}


