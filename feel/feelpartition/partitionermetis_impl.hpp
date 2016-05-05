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

namespace Metis
{
extern "C" {
//#include <metis.h>
#include <feelmetis.h>
} //"C"
}
#include <feel/feelpartition/csrgraphmetis.hpp>

namespace Feel
{

template <typename MeshType>
void PartitionerMetis<MeshType>::partitionImpl( mesh_ptrtype mesh, rank_type np )
{
    LOG( INFO ) << "PartitionerMetis::partitionImpl starts...";
    tic();

    // Check for an easy return
    if ( np == 1 )
    {
        this->singlePartition( mesh );
        return;
    }
    const dof_id_type n_elems = mesh->numElements();
    if ( n_elems == 0 )
        return;

    // build the graph
    // std::vector<Metis::idx_t> options(5);
    std::vector<Metis::idx_t> vwgt( n_elems );
    std::vector<Metis::idx_t> part( n_elems );

    // number of "nodes" (elements) in the graph
    Metis::idx_t n = static_cast<Metis::idx_t>( n_elems );
    // number of subdomains to create
    Metis::idx_t nparts = static_cast<Metis::idx_t>( np );
    // number of edges cut by the resulting partition
    Metis::idx_t edgecut = 0;

    std::map<dof_id_type, std::pair<dof_id_type, rank_type>> global_index_map;

    {
        std::vector<dof_id_type> global_index( n_elems, 0 );
        std::iota( global_index.begin(), global_index.end(), 0 );

        size_type cnt = 0;

        for ( auto const& elt : allelements( mesh ) )
        {
            global_index_map.insert( std::make_pair( elt.id(), std::make_pair( global_index[cnt++], elt.processId() ) ) );
        }
    }

    // Invoke METIS, but only on processor 0.
    // Then broadcast the resulting decomposition
    if ( mesh->worldComm().isMasterRank() )
    {
        CSRGraphMetis<Metis::idx_t> csr_graph;
        csr_graph.offsets.resize( n_elems + 1, 0 );

        // Local scope for these
        {
#ifndef NDEBUG
            std::size_t graph_size = 0;
#endif
            // build the graph in CSR format.  Note that
            // the edges in the graph will correspond to
            // face neighbors
            for ( auto const& elt : allelements( mesh ) )
            {

                // (1) first pass - get the row sizes for each element by counting the number
                // of face neighbors.  Also populate the vwght array if necessary
                const dof_id_type gid = global_index_map[elt.id()].first;

                CHECK( gid < vwgt.size() ) << "Invalid gid " << gid << " greater or equal than " << vwgt.size();

                // maybe there is a better weight?
                // The weight is used to define what a balanced graph is
                //if(!_weights)
                vwgt[gid] = elt.numPoints;
                //else
                //vwgt[gid] = static_cast<Metis::idx_t>((*_weights)[elem->id()]);

                unsigned int num_neighbors = 0;

                // Loop over the element's neighbors.  An element
                // adjacency corresponds to a face neighbor
                for ( uint16_type ms = 0; ms < elt.nNeighbors(); ms++ )
                {
                    element_type const* neighbor = NULL;
                    size_type neighbor_id = elt.neighbor( ms ).first;
                    if ( neighbor_id != invalid_size_type_value )
                    {
                        num_neighbors++;
                    }
                }
#if 0
                std::cout << "element id " << elt.id() << " gid: " << gid << " w: " << vwgt[gid] 
                          << " neigh: " << num_neighbors << std::endl;
#endif
                csr_graph.prepareNumberNonZeros( gid, num_neighbors );
#ifndef NDEBUG
                graph_size += num_neighbors;
#endif
            }

            csr_graph.prepareForUse();

            // (2) second pass - fill the compressed adjacency array
            for ( auto const& elt : allelements( mesh ) )
            {
                dof_id_type gid = global_index_map[elt.id()].first;

                unsigned int connection = 0;

                // Loop over the element's neighbors.  An element
                // adjacency corresponds to a face neighbor
                for ( uint16_type ms = 0; ms < elt.nNeighbors(); ms++ )
                {
                    element_type const* neighbor = NULL;
                    size_type neighbor_id = elt.neighbor( ms ).first;
                    if ( neighbor_id != invalid_size_type_value )
                    {
                        csr_graph( gid, connection++ ) = global_index_map[neighbor_id].first;
                    }
                }
            }
#ifndef NDEBUG
            // We create a non-empty vals for a disconnected graph, to
            // work around a segfault from METIS.
            DCHECK( csr_graph.vals.size() == std::max( graph_size, std::size_t( 1 ) ) )
                << "Invalid graph";
#endif
        } // done building the graph

        Metis::idx_t ncon = 1;

        // Select which type of partitioning to create

        // Use recursive if the number of partitions is less than or equal to 8
        if ( np <= 8 )
            Metis::Feel_METIS_PartGraphRecursive( &n, &ncon, &csr_graph.offsets[0], &csr_graph.vals[0], &vwgt[0], NULL,
                                                  NULL, &nparts, NULL, NULL, NULL,
                                                  &edgecut, &part[0] );

        // Otherwise  use kway
        else
            Metis::Feel_METIS_PartGraphKway( &n, &ncon, &csr_graph.offsets[0], &csr_graph.vals[0], &vwgt[0], NULL,
                                             NULL, &nparts, NULL, NULL, NULL,
                                             &edgecut, &part[0] );

    } // end processor 0 part

    // Assign the returned processor ids.  The part array contains the processor
    // id for each element, but in terms of the contiguous indexing we defined
    // above
    LOG( INFO ) << "PartitionerMetis::partitionImpl nelements : " << n_elems;

    for ( auto const& pairElt : global_index_map )
    {
        dof_id_type eltId = pairElt.first;
        dof_id_type gid = pairElt.second.first;
        rank_type initialPid = pairElt.second.second;
        rank_type newPid = static_cast<rank_type>( part[gid] );
        auto eltToUpdate = mesh->elementIterator( eltId, initialPid );
        mesh->elements().modify( eltToUpdate, Feel::detail::UpdateProcessId( newPid ) );
    }

    auto t = toc( "PartitionerMetis::partitionImpl", FLAGS_v > 0 );
    LOG( INFO ) << "PartitionerMetis::partitionImpl done in " << t << "s";
}

} // Feel

#endif
