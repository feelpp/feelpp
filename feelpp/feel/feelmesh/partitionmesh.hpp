/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 janv. 2016
 
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
#ifndef FEELPP_PARTITIONMESH_HPP
#define FEELPP_PARTITIONMESH_HPP 1

#include <boost/iterator/counting_iterator.hpp>
#include <feel/feelmesh/meshpartitionset.hpp>
#include <feel/feelpartition/partitioner.hpp>
#if defined(FEELPP_HAS_METIS)
#include <feel/feelpartition/partitionermetis.hpp>
#endif




namespace Feel {

/**
 * partition a mesh and provide a helper data struct navigate in the newly
 * partitioned mesh
 *
 * @return a std::unique_ptr to a MeshPartitionSet
 */
template<typename MeshType>
std::unique_ptr<MeshPartitionSet<MeshType>>
partitionMesh( std::shared_ptr<MeshType> mesh,
               rank_type nGlobalParts,
               std::vector<Range<MeshType,MESH_ELEMENTS>> partitionByRange = std::vector<Range<MeshType,MESH_ELEMENTS>>(),
               json const& partconfig = json() )
{
#if defined(FEELPP_HAS_METIS)
    // metis is hard coded for now, this will be customizable with different
    // partitioners
    PartitionerMetis<MeshType> metis( partconfig );
    metis.partition( mesh, nGlobalParts, partitionByRange );
    std::set<rank_type> localPartitionIds (boost::counting_iterator<int
    >(0), boost::counting_iterator<int>(nGlobalParts));
    return std::make_unique<MeshPartitionSet<MeshType>>( mesh, nGlobalParts, localPartitionIds );
#else
    CHECK(false) << "no partitioner implementation";
    return NULL;
#endif
}



}
#endif
