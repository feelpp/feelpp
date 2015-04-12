/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 11 Apr 2015
 
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

#include <feel/feelmodels/modelpostprocess.hpp>

namespace Feel {

ModelPostprocess::ModelPostprocess()
{}

ModelPostprocess::ModelPostprocess(pt::ptree const& p)
{}

ModelPostprocess::~ModelPostprocess()
{}

void
ModelPostprocess::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}

template <typename T>
std::vector<T> as_vector(pt::ptree const& pt, pt::ptree::key_type const& key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
        r.push_back(item.second.template get_value<T>());
    return r;
}

void
ModelPostprocess::setup()
{
    auto fields = M_p.get_child_optional("Fields");
    if ( fields )
    {
        
        for (auto i : as_vector<std::string>(M_p, "Fields"))
        {
            this->operator[]("Fields").push_back( i );
            LOG(INFO) << "add to postprocess field  " << i;
        }
        
    }
    auto forces = M_p.get_child_optional("Force");
    if ( forces )
    {
        
        for (auto i : as_vector<std::string>(M_p, "Force"))
        {
            this->operator[]("Force").push_back( i );
            LOG(INFO) << "add to postprocess force  " << i;
            
        }
        
    }
}


}
