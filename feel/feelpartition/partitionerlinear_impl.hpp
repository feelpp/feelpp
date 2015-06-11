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
#ifndef FEELPP_PARTITIONERLINEAR_IMPL_HPP
#define FEELPP_PARTITIONERLINEAR_IMPL_HPP 1


#include <feel/feelpartition/partitionerlinear.hpp>


namespace Feel {

template<typename MeshType>
void PartitionerLinear<MeshType>::partitionImpl ( mesh_ptrtype& mesh, 
                                                  rank_type n)
{
    CHECK( n != invalid_rank_type_value ) 
        << "Invalid number of partitions : " 
        << n << " must be != " << invalid_rank_type_value;
  

  // Check for an easy return
  if (n == 1)
    {
      this->singlePartition (mesh);
      return;
    }

  LOG(INFO) << "linear partition mesh starts...";
  tic();

  dof_id_type blksize       = mesh->numElements()/n;
  
  dof_id_type e = 0;
  
  for( auto elt : elements( mesh ) )
  {
      auto pid = ((e/blksize) < n)?:(e/blksize):0;
      auto it = mesh->elementIterator( elt );
      mesh->modify( it, []( element_type& e ){ e.setProcessId(pid); });
      e++;
  }
  auto t = toc("linear partition",FLAGS_v > 0);
  LOG(INFO) << "linear partition mesh done in " << t << "s";
}


} // Feel

#endif
