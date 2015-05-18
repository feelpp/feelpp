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
#ifndef FEELPP_PARTITIONERMETIS_IMPL_HPP
#define FEELPP_PARTITIONERMETIS_IMPL_HPP 1

namespace Metis {
extern "C" {
#include <metis.h>
} //"C"
}

namespace Feel {

template<typename MeshType>
void 
PartitionerMetis<MeshType>::partitionImpl ( mesh_ptrtype mesh, rank_type np )
{
    LOG(INFO) << "PartitionerMetis::partitionImpl starts...";
    tic();

    // Check for an easy return
    if (np == 1)
    {
        this->singlePartition (mesh);
        return;
    }
    const dof_id_type n_elems = mesh->numElements();

    // build the graph
    // std::vector<Metis::idx_t> options(5);
    std::vector<Metis::idx_t> vwgt(n_elems);
    std::vector<Metis::idx_t> part(n_elems);

    // number of "nodes" (elements) in the graph
    Metis::idx_t n = static_cast<Metis::idx_t>(n_elems);
    // number of subdomains to create
    Metis::idx_t nparts  = static_cast<Metis::idx_t>(np); 
    // number of edges cut by the resulting partition
    Metis::idx_t edgecut = 0;
    
    auto t = toc("PartitionerMetis::partitionImpl", FLAGS_v > 0 );
    LOG(INFO) << "PartitionerMetis::partitionImpl done in " << t << "s";
}


} // Feel


#endif
