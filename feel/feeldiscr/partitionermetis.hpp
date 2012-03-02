/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-09

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foun dation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#ifndef __metis_partitioner_h__
#define __metis_partitioner_h__

#include <feel/feelconfig.h>


namespace Metis
{
extern "C" {
#if defined( HAVE_METIS_H )
# include <metis.h>
#else
#if defined( HAVE_METIS_METIS_H )
# include <metis/metis.h>
#else
#if defined( HAVE_METIS )
# include <metis/metis.h>
#endif // HAVE_METIS
#endif // HAVE_METIS_METIS_H
#endif // HAVE_METIS_H
} // "C"
} // Metis


#include <feel/feeldiscr/partitioner.hpp>

namespace Feel
{

/**
 * The \p PartitionerMetis uses the Metis graph partitioner
 * to partition the elements.
 */
template<typename Mesh>
class PartitionerMetis
    :
    public Partitioner<Mesh>
{
    typedef Partitioner<Mesh> super;
public:

    typedef PartitionerMetis<Mesh> self_type;

    typedef typename super::mesh_type mesh_type;
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::element_iterator element_iterator;
    typedef typename mesh_type::element_const_iterator element_const_iterator;

    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::face_iterator face_iterator;
    typedef typename mesh_type::pid_face_iterator pid_face_iterator;
    typedef typename mesh_type::location_face_iterator location_face_iterator;
    /**
     * Constructor.
     */
    PartitionerMetis () {}

protected:
    /**
     * Partition the \p mesh_type into \p n subdomains.
     */
    virtual void doPartition ( mesh_type & mesh,
                               size_type n);

private:
};

template<typename Mesh>
void
PartitionerMetis<Mesh>::doPartition ( mesh_type& mesh,
                                      size_type n_pieces )
{
    FEELPP_ASSERT (n_pieces > 0)( n_pieces ).error( "the number of partitions should be >0" );

    // Check for an easy return
    if (n_pieces == 1)
        {
            this->singlePartition (mesh);
            return;
        }

    // What to do if the Metis library IS NOT present
#if !defined( HAVE_METIS_H ) && !defined( HAVE_METIS ) && !defined( HAVE_METIS_METIS_H )

    std::cerr << "ERROR: The library has been built without"    << std::endl
              << "Metis support. "  << std::endl;
    this->singlePartition (mesh);
    return;

    // What to do if the Metis library IS present
#else

    const unsigned int n_elem        = mesh.numElements();

    // build the graph
    // the forward_map maps each active element id
    // into a contiguous block of indices for Metis
    std::vector<size_type> forward_map (n_elem,  invalid_size_type_value );

    std::vector<int> xadj;
    std::vector<int> adjncy;
    std::vector<int> options(5);
    std::vector<int> vwgt(n_elem);
    std::vector<int> part(n_elem);

    xadj.reserve(n_elem+1);

    int
        n = static_cast<int>(n_elem),  // number of "nodes" (elements)
        //   in the graph
        wgtflag = 2,                          // weights on vertices only,
        //   none on edges
        numflag = 0,                          // C-style 0-based numbering
        nparts  = static_cast<int>(n_pieces), // number of subdomains to create
        edgecut = 0;                          // the numbers of edges cut by the
    //   resulting partition

    // Set the options
    options[0] = 0; // use default options


    // Metis will only consider the active elements.
    // We need to map the active element ids into a
    // contiguous range.
    {

        element_iterator       elem_it  = mesh.beginElement();
        const element_iterator elem_end = mesh.endElement();

        unsigned int el_num = 0;

        for (; elem_it != elem_end; ++elem_it)
            {
                FEELPP_ASSERT (elem_it->id() < forward_map.size())( elem_it->id() )( forward_map.size() ).error( "invalid dimensions" );

                forward_map[elem_it->id()] = el_num++;
            }

        FEELPP_ASSERT (el_num == n_elem)( el_num )( n_elem ).error( "incompatible number of elements" );
    }


    // build the graph in CSR format.  Note that
    // the edges in the graph will correspond to
    // face neighbors
    {
        std::vector<element_type const*> neighbors_offspring;


        element_iterator       elem_it  = mesh.beginElement();
        const element_iterator elem_end = mesh.endElement();

        // This will be exact when there is no refinement and all the
        // elements are of the same type.
        adjncy.reserve (n_elem * mesh.numLocalFaces() );

        for (; elem_it != elem_end; ++elem_it)
            {
                const element_type* elem = boost::addressof( *elem_it );

                FEELPP_ASSERT (elem->id() < forward_map.size()).error( "element id and forward_map incompatible" );
                FEELPP_ASSERT (forward_map[elem->id()] != invalid_size_type_value ).error( "forward_map problem" );

                // maybe there is a better weight?
                // The weight is used to define what a balanced graph is
                vwgt[forward_map[elem->id()]] = elem->nPoints();

                // The beginning of the adjacency array for this elem
                xadj.push_back(adjncy.size());

                size_type counter = 0;
                for (uint16_type ms=0; ms < mesh.numLocalFaces(); ms++)
                    {
                        Debug( 4021 ) << "** element " << elem_it->id() << " face " << ms << " neighbor " << mesh.localFaceId( elem->id(), ms ).template get<1>() << "\n";
                        if ( mesh.localFaceId( elem->id(), ms ).template get<1>() != invalid_size_type_value )
                            ++counter;

                    }
                Debug( 4021 ) << "[PartitionerMetis] element " << elem->id() << " number of neighbors: " << counter << "\n";
                FEELPP_ASSERT( counter >= 1 )( elem->id() )( counter ).error( "invalid neighboring data" );


                // Loop over the element's neighbors.  An element
                // adjacency corresponds to a face neighbor
                for (uint16_type ms=0; ms < mesh.numLocalFaces(); ms++)
                    {
                        size_type neighbor_id = mesh.localFaceId( elem->id(), ms ).template get<1>();
                        if ( neighbor_id != invalid_size_type_value )
                            {
                                FEELPP_ASSERT (neighbor_id < forward_map.size())
                                ( neighbor_id )( forward_map.size() ).error( "problem with neighbor id and forward_map" );
                                FEELPP_ASSERT (forward_map[neighbor_id] != invalid_size_type_value )
                                ( (forward_map[neighbor_id] ) ).error( "invalid forward_map" );

                                adjncy.push_back (forward_map[neighbor_id]);
                            }
                    }
            }

        // The end of the adjacency array for the last elem
        xadj.push_back(adjncy.size());

    } // done building the graph

    // Select which type of partitioning to create

    // Use recursive if the number of partitions is less than or equal to 8
    if (n_pieces <= 8)
        Metis::METIS_PartGraphRecursive(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
                                        &wgtflag, &numflag, &nparts, &options[0],
                                        &edgecut, &part[0]);

    // Otherwise  use kway
    else
        Metis::METIS_PartGraphKway(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL,
                                   &wgtflag, &numflag, &nparts, &options[0],
                                   &edgecut, &part[0]);


    // Assign the returned processor ids.  The part array contains
    // the processor id for each active element, but in terms of
    // the contiguous indexing we defined above
    {

        element_iterator elem_it  = mesh.beginElement();
        element_iterator elem_end = mesh.endElement();

        for (size_type i = 0; i < mesh.numElements(); ++i )
            {
                element_iterator  elem_it = mesh.elementIterator( i, 0 );
                element_type elem = *elem_it;

                FEELPP_ASSERT ( elem.id() < forward_map.size() )( elem.id() )( forward_map.size() ).error( "invalid size" );
                FEELPP_ASSERT ( forward_map[elem.id()] != invalid_size_type_value )( forward_map[elem.id()] ).error( "invalid forward map" );;
                elem.setProcessId( static_cast<short int>(part[forward_map[elem.id()]]) );


                // update neighbor process id
                for (uint16_type ms=0; ms < mesh.numLocalFaces(); ms++)
                    {
                        size_type eid = mesh.localFaceId( elem, ms ).template get<1>();
                        mesh.localFaceId( elem, ms ).template get<1>() = static_cast<short int>( part[forward_map[eid]] );
#if 0
                        size_type neighbor_id = elem.neighbor(ms).first;
                        if ( neighbor_id != invalid_size_type_value )
                            {
                                FEELPP_ASSERT ( neighbor_id < forward_map.size() )( neighbor_id )( forward_map.size() ).error( "invalid size" );
                                FEELPP_ASSERT ( forward_map[neighbor_id] != invalid_size_type_value )( forward_map[neighbor_id] ).error( "invalid forward map" );

                                elem.setNeighbor( ms, neighbor_id, static_cast<short int>(part[forward_map[neighbor_id]]) );

                                Debug( 4021 ) << "[PartitionerMetis] neighbor element " << neighbor_id << " of element " << elem_it->id() << " is on proc " << elem.neighbor(ms).second << "\n";

                            }
#endif
                    }


                // here we invalidate any pointers to *elem_it that
                // are stored in other data sutrctures such as faces
                mesh.elements().replace( elem_it, elem );

#if 0
                // go through the faces and update the pointers to the
                // element
                element_type const& __element = *elem_it;
                for ( size_type j = 0; j < numLocalFaces(); j++ )
                    {
                        detail::UpdateFaceConnection1<typename face_type::element_connectivity_type> update1( boost::make_tuple( boost::addressof( __element ), __element.id(), j ) );
                        update1( __f );

                    }
#endif
                Debug( 4021 ) << "[PartitionerMetis] element " << elem_it->id() << " will be on proc " << elem_it->processId() << "\n";

            }
        for ( int i = 0; i < n_pieces; ++i )
            {
                size_type dist = std::distance( mesh.beginElementWithProcessId( i ),
                                                mesh.endElementWithProcessId( i ) );
                Debug( 4020 ) << "[PartitionerMetis] " << dist << " elts on proc " << i << "\n";
            }
    }
#if 0
    // Assign the returned processor ids.  The part array contains
    // the processor id for each active1 boundary faces, but in terms of
    // the contiguous indexing we defined above
    {

        face_iterator face_it = mesh.beginFace();
        face_iterator face_end = mesh.endFace();
        Debug( 4020 ) << "[PartitionerMetis] " << " boundary faces  " << std::distance( face_it,face_end ) << "\n";


        for (; face_it != face_end; ++face_it)
            {
                face_type  face = *face_it;

                // boundary faces can belong to only one process id
                if ( face.isOnBoundary() )
                    {
                        FEELPP_ASSERT( face.isConnectedTo0() &&
                                     !face.isConnectedTo1() )
                            ( face.ad_first() )( face.pos_first() )
                            ( face.ad_second() )( face.pos_second() ).error( "invalid boundary face" );

                        Debug( 4021 ) << "[PartitionerMetis] face " << face.id() << " will be on proc " << face.element0().processId() << "\n";
                        face.setProcessId( face.element0().processId() );
                        mesh.faces().replace( face_it, face );
                    }
                else
                    {
                        Debug( 4021 ) << "[PartitionerMetis] face " << face.id() << " will be on proc " << face.element0().processId() << "\n";

                        FEELPP_ASSERT( face.isConnectedTo0() &&
                                     face.isConnectedTo1() )
                            ( face.ad_first() )( face.pos_first() )
                            ( face.ad_second() )( face.pos_second() ).error( "invalid internal face" );

                        //
                        // Assign process ids to faces. we don't
                        // really about faces that don't belong to the
                        // current processors as they won't be used by
                        // the system.
                        //
                        uint16_type proc0 = face.element0().processId();
                        uint16_type proc1 = face.element1().processId();
                        if ( proc0 == proc1 )
                            face.setProcessId( proc0 );
                        else if ( proc0 == this->comm().rank() )
                            face.setProcessId( proc0 );
                        else if ( proc1 == this->comm().rank() )
                            face.setProcessId( proc1 );
                        else
                            face.setProcessId( proc0 );

                    }

                FEELPP_ASSERT( face.element(0).facePtr( face.pos_first() ) )
                    ( face.ad_first() )( face.pos_first() )( face.element(0).id() ).error( "invalid face in element" );
                if ( face.isConnectedTo1() )
                    FEELPP_ASSERT( face.element(1).facePtr( face.pos_second() ) )
                        ( face.ad_second() )( face.pos_second() )( face.element(1).id() ).error( "invalid internal face in element" );

                mesh.faces().replace( face_it, face );
            }
        for ( int i = 0; i < n_pieces; ++i )
            {
                std::pair<location_face_iterator, location_face_iterator> p = mesh.facesOnBoundary( i );
                size_type dist = std::distance( p.first, p.second );
                Debug( 4020 ) << "[PartitionerMetis] " << dist << " boundary faces on proc " << i << "\n";
            }
    }
#endif // 0
    for ( int i = 0; i < n_pieces; ++i )
        {
            typename mesh_type::element_iterator pid_it = mesh.beginElementWithProcessId( i );
            typename mesh_type::element_iterator pid_en = mesh.endElementWithProcessId( i );
            while( pid_it != pid_en )
                {
                    Debug( 4021 ) << "piece " << i << " element " << pid_it->id() << " on proc " << pid_it->processId() << "\n";
                    ++pid_it;
                }
        }

#endif


}


} // Feel


#endif
