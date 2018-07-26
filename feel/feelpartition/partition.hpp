/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 May 2015
 
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
#ifndef FEELPP_PARTITION_HPP
#define FEELPP_PARTITION_HPP 1

#include <functional>

#include <feel/feelpartition/partitioner.hpp>
#include <feel/feelpartition/partitionermetis.hpp>
#include <feel/feelpartition/partitionerlinear.hpp>



namespace Feel {



template<typename MeshType>
std::map<std::string, std::unique_ptr<Partitioner<MeshType> > > 
registerPartitioner()
{
    typedef std::function< Partitioner<MeshType>*() > factory_partitioner;
    typedef boost::function< std::unique_ptr<Partitioner<MeshType>>() > a_factory;

    std::map<std::string, std::unique_ptr<Partitioner<MeshType> > > m;
    
    return m;
    
}
template<typename MeshType>
void partition( std::string partitioner, std::shared_ptr<MeshType> mesh, int n_parts )
{
    

}


}
#endif
