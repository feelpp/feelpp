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
#ifndef FEELPP_PARTITIONER_IMPL_HPP
#define FEELPP_PARTITIONER_IMPL_HPP 1

#include <feel/feelmesh/filters.hpp>



namespace Feel {

// Partitioner static data
template<typename MeshType>
const dof_id_type 
Partitioner<MeshType>::communication_blocksize = 1000000;

// Partitioner implementation
template<typename MeshType>
void 
Partitioner<MeshType>::partition (mesh_ptrtype mesh)
{
    this->partition( mesh, Environment::numberOfProcessors() );
}


template<typename MeshType>
void 
Partitioner<MeshType>::partition ( mesh_ptrtype mesh, rank_type n, std::vector<range_element_type> const& partitionByRange )
{
    LOG(INFO) << "Partitioner::partition starts";
    // we cannot partition into more pieces than we have
    // elements!
    size_type n_parts = std::min(mesh->numElements(), static_cast<dof_id_type>(n));

    // Set the number of partitions in the mesh
    mesh->setNumberOfPartitions(n_parts);
    LOG(INFO) << "Partitioner::partition:: number of partitions=" << n_parts;

    if (n_parts == 1)
    {
        this->singlePartition (mesh);
        return;
    }

    // Call the partitioning function
    this->partitionImpl(mesh,n_parts,partitionByRange);

    // Set the node's processor ids
    Partitioner::setNodeProcessorIds(mesh);
    LOG(INFO) << "Partitioner::partition done";
}


template<typename MeshType>
void 
Partitioner<MeshType>::repartition (mesh_ptrtype mesh)
{
    this->repartition( mesh, Environment::numberOfProcessors() );
}


template<typename MeshType>
void 
Partitioner<MeshType>::repartition ( mesh_ptrtype mesh, rank_type n )
{
  // we cannot partition into more pieces than we have
  // active elements!
  size_type n_parts = std::min(mesh->numElements(), static_cast<dof_id_type>(n));

  // Set the number of partitions in the mesh
  mesh.setNumberOfPartitions(n_parts);
  LOG(INFO) << "Partitioner::partition:: number of partitions=" << n_parts;

  if (n_parts == 1)
    {
      this->singlePartition (mesh);
      return;
    }

  // Call the partitioning function
  this->repartitionImpl(mesh,n_parts);

  // Set the node's processor ids
  this->setNodeProcessorIds(mesh);
}




template<typename MeshType>
void 
Partitioner<MeshType>::singlePartition ( mesh_ptrtype mesh )
{
    for( auto const& elt : elements(mesh) )
        mesh->elementIterator( boost::unwrap_ref( elt ).id() )->second.setProcessId( 0 );
    for( auto const& elt : faces(mesh) )
        mesh->faceIterator( boost::unwrap_ref( elt ).id() )->second.setProcessId( 0 );
    for ( auto itp = mesh->beginPoint(), enp = mesh->endPoint(); itp != enp; ++itp )
        itp->second.setProcessId( 0 );
}


template<typename MeshType>
void 
Partitioner<MeshType>::setNodeProcessorIds(mesh_ptrtype mesh)
{
    LOG(INFO) << "setNodeProcessIds starts...";
    LOG(INFO) << "setNodeProcessIds done.";
}



} // Feel

#endif
