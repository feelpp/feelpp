/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 Mar 2015
 
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
#include <feel/feelcore/environment.hpp>

#include <feel/feelmodels/modelparameters.hpp>

namespace Feel {

ModelParameters::ModelParameters()
{}

ModelParameters::ModelParameters(pt::ptree const& p)
{}

ModelParameters::~ModelParameters()
{}

void
ModelParameters::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}


void
ModelParameters::setup()
{
    for( auto const& v : M_p )
    {

        if ( Environment::isMasterRank() )
            std::cout << "reading parameter " << v.first  << "\n";
        std::string t = v.first; // parameter name
        auto f= v.second;
        {
            try
            {
                auto val= M_p.get<double>(v.first);
                LOG(INFO) << "adding parameter " << t << " with value " << val;
                this->operator[](t) = ModelParameter( t, val, val, val );
            }
            catch( ... )
            {
                try
                {
                    auto val= f.get<double>("value");
                    LOG(INFO) << "adding parameter " << t << " with value " << val;
                    this->operator[](t) = ModelParameter( t, val, val, val );
                }
                catch( ... )
                {
                    this->operator[](t) = ModelParameter( t, 0., 0., 0. );
                }
            }
        }
        
    }
}


}
