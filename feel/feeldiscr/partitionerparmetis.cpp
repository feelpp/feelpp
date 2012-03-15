/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-18

  Copyright (C) 2007 Universit� Joseph Fourier (Grenoble I)

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
/**
   \file partitionerparmetis.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-18
 */
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/partitionerparmetis.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{


template<typename Mesh>
void
PartitionerParmetis<Mesh>::doPartition ( mesh_type& mesh,
                                         const uint16_type n_pieces )
{
    FEELPP_ASSERT( mesh.isUpdatedForUse() ).error( "invalid call, must call mesh.updateForUse first" );
    FEELPP_ASSERT (n_pieces > 0)( n_pieces ).error( "the number of partitions should be >0" );

    // Check for an easy return
    if (n_pieces == 1)
        {
            this->singlePartition (mesh);
            return;
        }
    Debug() << "[PartitionerParmetis::doPartition] n_pieces = " << n_pieces << "\n";
    // What to do if the Parmetis library IS NOT present
#if !defined( FEELPP_HAS_PARMETIS_H ) && !defined(FEELPP_HAS_PARMETIS_PARMETIS_H)

    std::cerr << "ERROR: The library has been built without"    << std::endl
              << "Parmetis support.   Using a Metis"           << std::endl
              << "partitioner instead!"                       << std::endl;
    PartitionerParmetis<Mesh> metis;
    metis.partition( mesh, n_pieces );
    return;

    // What to do if the Parmetis library IS present
#else

    // Initialize the data structures required by ParMETIS
    this->initialize (mesh, n_pieces);

    // build the graph corresponding to the mesh
    this->buildGraph (mesh);


    // Partition the graph
    std::vector<int>  local_part(M_part);
    MPI_Comm mpi_comm = this->comm();


    // Call the ParMETIS k-way partitioning algorithm.
    //Parmetis::ParMETIS_V3_PartKway(&M_vtxdist[0], &M_xadj[0], &M_adjncy[0], &M_vwgt[0], NULL,
    Metis::ParMETIS_V3_PartKway(&M_vtxdist[0], &M_xadj[0], &M_adjncy[0], &M_vwgt[0], NULL,
                         &M_wgtflag, &M_numflag, &M_ncon, &M_nparts, &M_tpwgts[0],
                         &M_ubvec[0], &M_options[0], &M_edgecut,
                         &local_part[M_first_local_elem],
                         &mpi_comm);

    // Collect the partioning information from all the processors.
    FEELPP_ASSERT (M_part.size() == local_part.size())
        (M_part.size())(local_part.size()).error( "invalid partitioning info" );
    MPI_Allreduce (&local_part[0], &M_part[0], M_part.size(), MPI_INT, MPI_SUM,
                   this->comm());

    // Assign the returned processor ids
    this->assignPartitioning (mesh);

#endif // FEELPP_HAS_PARMETIS_H
}


// Only need to compile these methods if ParMETIS is present
#if defined( FEELPP_HAS_PARMETIS_H ) || defined( FEELPP_HAS_PARMETIS_PARMETIS_H )
template<typename Mesh>
void
PartitionerParmetis<Mesh>::initialize (const mesh_type& mesh,
                                       const unsigned int n_sbdmns)
{
    const unsigned int n_elem               = mesh.elements().size();
    const unsigned int n_active_local_elem  = std::distance( mesh.beginElementWithProcessId( this->comm().rank() ),
                                                             mesh.endElementWithProcessId( this->comm().rank() ) );
    const unsigned int n_active_elem        = n_elem;//mesh.n_active_elem();
    const unsigned int n_procs              = this->comm().size();

    // Set parameters.
    M_wgtflag = 2;                          // weights on vertices only
    M_ncon    = 1;                          // one weight per vertex
    M_numflag = 0;                          // C-style 0-based numbering
    M_nparts  = static_cast<int>(n_sbdmns); // number of subdomains to create
    M_edgecut = 0;                          // the numbers of edges cut by the
    //   partition

    // Initialize data structures for ParMETIS
    M_vtxdist.resize (n_procs+1);     std::fill (M_vtxdist.begin(), M_vtxdist.end(), 0);
    M_tpwgts.resize  (M_nparts);      std::fill (M_tpwgts.begin(),  M_tpwgts.end(),  1./M_nparts);
    M_ubvec.resize   (M_ncon);        std::fill (M_ubvec.begin(),   M_ubvec.end(),   1.);
    M_part.resize    (n_active_elem); std::fill (M_part.begin(),    M_part.end(),    0);
    M_options.resize (5);
    M_vwgt.resize    (n_active_local_elem);


    // Set the options
    M_options[0] = 0; // use default options


    // Set up the vtxdist array.  This will be the same on each processor.
    // Consult the Parmetis documentation.
    {
        FEELPP_ASSERT (M_vtxdist.size() == n_procs+1)
            ( M_vtxdist.size() )( n_procs+1 ).error( "invalid vtxdist array size" );
        FEELPP_ASSERT (M_vtxdist[0] == 0)( M_vtxdist[0] ).error( "invalid vtxdist entry 0" );

        for (uint16_type proc_id=0; proc_id<this->comm().size(); proc_id++)
            M_vtxdist[proc_id+1] = M_vtxdist[proc_id] + std::distance( mesh.beginElementWithProcessId( proc_id ),
                                                                       mesh.endElementWithProcessId( proc_id ) );

        assert (M_vtxdist[this->comm().size()] == static_cast<int>(n_active_elem));
    }

    // Metis will only consider the active elements.
    // We need to map the active element ids into a
    // contiguous range.
    M_forward_map.resize (n_elem); std::fill (M_forward_map.begin(),
                                              M_forward_map.end(),
                                              invalid_size_type_value );
    M_first_local_elem = 0;
    size_type el_num = 0;
    size_type local_el_num = 0;

    for (unsigned int proc_id=0; proc_id<this->comm().size(); proc_id++)
        {
            if (proc_id == this->comm().rank() )
                M_first_local_elem = el_num;

            element_const_iterator elem_it  = mesh.beginElementWithProcessId(proc_id);
            element_const_iterator elem_end = mesh.endElementWithProcessId(proc_id);

            for (; elem_it != elem_end; ++elem_it)
                {
                    FEELPP_ASSERT ( elem_it->id() < M_forward_map.size() ).error( "invalid element id" );
                    FEELPP_ASSERT ( M_forward_map[elem_it->id()] == invalid_size_type_value );

                    M_forward_map[elem_it->id()] = el_num;
                    el_num++;

                    // maybe there is a better weight?
                    if ( proc_id == this->comm().rank() )
                        M_vwgt[local_el_num++] = elem_it->nPoints();
                }
        }

    FEELPP_ASSERT (el_num       == n_active_elem)( el_num )( n_active_elem ).error( "invalid element number" );
    FEELPP_ASSERT (local_el_num == n_active_local_elem)( local_el_num )( n_active_local_elem ).error( "invalid local element number" );
}


template<typename Mesh>
void
PartitionerParmetis<Mesh>::buildGraph (const mesh_type& mesh)
{
    // build the graph in distributed CSR format.  Note that
    // the edges in the graph will correspond to
    // face neighbors

    // Reserve space in the adjacency array
    const unsigned int n_active_local_elem  = mesh.numElements();
    M_xadj.reserve (n_active_local_elem + 1);

    std::vector<element_type const*> neighbors_offspring;

    element_iterator       elem_it  = mesh.beginElementWithProcessId( this->comm().rank() );
    const element_iterator elem_end = mesh.endElementWithProcessId( this->comm().rank() );

    for (; elem_it != elem_end; ++elem_it)
        {
            const element_type* elem = boost::addressof( *elem_it );

            FEELPP_ASSERT (elem->id() < M_forward_map.size()).error( "element id and forward_map incompatible" );
            FEELPP_ASSERT (M_forward_map[elem->id()] != invalid_size_type_value ).error( "forward_map problem" );

            // The beginning of the adjacency array for this elem
            M_xadj.push_back(M_adjncy.size());

            // Loop over the element's neighbors.  An element
            // adjacency corresponds to a face neighbor
            for (unsigned int ms=0; ms<elem->nNeighbors(); ms++)
                {
                    element_type const* neighbor = NULL;
                    size_type neighbor_id = elem->neighbor(ms).first;

                    if ( neighbor_id != invalid_size_type_value )
                        {
                            neighbor = boost::addressof( mesh.element( neighbor_id, 0 ) );
                            Debug( 4021 ) << " neighbor id = " << neighbor_id << " elem = " << neighbor << "\n";
                            if ( neighbor )
                                Debug( 4021 ) << " elem id = " << neighbor->id() << " proc = " << neighbor->processId() << "\n";
                        }
                    if (neighbor != NULL)
                        {
                            // If the neighbor is active treat it as a
                            // connection
                            if (neighbor->active())
                                {
                                    FEELPP_ASSERT (neighbor->id() < M_forward_map.size())
                                        ( neighbor->id() )( M_forward_map.size() ).error( "problem with neighbor id and forward_map" );
                                    FEELPP_ASSERT (M_forward_map[neighbor->id()] != invalid_size_type_value )
                                        ( (M_forward_map[neighbor->id()] ) ).error( "invalid forward_map" );

                                    M_adjncy.push_back (M_forward_map[neighbor->id()]);
                                }
                        }
                }
        }

    // The end of the adjacency array for this elem
    M_xadj.push_back (M_adjncy.size());
}


template<typename Mesh>
void
PartitionerParmetis<Mesh>::assignPartitioning (mesh_type& mesh)
{
    // TODO: refactoring: move it to partitioner base class, it is
    // _exactely_ the same as in partitionermetis

    // Assign the returned processor ids
    element_iterator elem_it  = mesh.beginElement();
    element_iterator elem_end = mesh.endElement();

    for (; elem_it != elem_end; ++elem_it)
        {
            element_type elem = *elem_it;

            FEELPP_ASSERT ( elem.id() < M_forward_map.size() )( elem.id() )( M_forward_map.size() ).error( "invalid size" );
            FEELPP_ASSERT ( M_forward_map[elem.id()] != invalid_size_type_value )( M_forward_map[elem.id()] ).error( "invalid forward map" );
            FEELPP_ASSERT (M_forward_map[elem.id()] < M_part.size())
                ( M_forward_map[elem.id()] )( M_part.size() ).error( "invalid forward map entry" );

            elem.setProcessId( static_cast<short int>(M_part[M_forward_map[elem.id()]]) );
            mesh.elements().replace( elem_it, elem );

        }

    for ( int i = 0; i < M_nparts; ++i )
        {
            size_type dist = std::distance( mesh.beginElementWithProcessId( i ),
                                            mesh.endElementWithProcessId( i ) );
            Debug( 4020 ) << "[PartitionerParmetis] " << dist << " elts on proc " << i << "\n";
        }

        // Assign the returned processor ids.  The part array contains
    // the processor id for each active1 boundary faces, but in terms of
    // the contiguous indexing we defined above
    {

        face_iterator face_it = mesh.beginFace();
        face_iterator face_end = mesh.endFace();
        Debug( 4020 ) << "[PartitionerParmetis] " << " boundary faces  " << std::distance( face_it,face_end ) << "\n";

        // boundary faces can belong to only one process id
        for (; face_it != face_end; ++face_it)
        {

            if ( face_it->isOnBoundary() )
            {
                face_type  face = *face_it;

                Debug( 4021 ) << "[PartitionerParmetis] face " << face.id() << " will be on proc " << face.element0().processId() << "\n";
                face.setProcessId( face.element0().processId() );
                mesh.faces().replace( face_it, face );
            }
        }
        for ( int i = 0; i < M_nparts; ++i )
        {
            std::pair<location_face_iterator, location_face_iterator> p = mesh.facesOnBoundary( i );
            size_type dist = std::distance( p.first, p.second );
            Debug( 4020 ) << "[PartitionerParmetis] " << dist << " boundary faces on proc " << i << "\n";
        }
    }
    for ( int i = 0; i < M_nparts; ++i )
        {
            typename mesh_type::element_iterator pid_it = mesh.beginElementWithProcessId( i );
            typename mesh_type::element_iterator pid_en = mesh.endElementWithProcessId( i );
            while( pid_it != pid_en )
                {
                    Debug( 4021 ) << "piece " << i << " element " << pid_it->id() << " on proc " << pid_it->processId() << "\n";
                    ++pid_it;
                }
        }

}

#endif // #ifdef FEELPP_HAS_PARMETIS



#if defined( FEELPP_INSTANTIATION_MODE )
//
// explicit instantiations
//
template class PartitionerParmetis<Mesh<Simplex<1,1> > >;
template class PartitionerParmetis<Mesh<Simplex<2,1> > >;
template class PartitionerParmetis<Mesh<Simplex<3,1> > >;

#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template class PartitionerParmetis<Mesh<Simplex<1,2> > >;
template class PartitionerParmetis<Mesh<Simplex<2,2> > >;
template class PartitionerParmetis<Mesh<Simplex<3,2> > >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template class PartitionerParmetis<Mesh<Simplex<1,3> > >;
template class PartitionerParmetis<Mesh<Simplex<2,3> > >;
template class PartitionerParmetis<Mesh<Simplex<3,3> > >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template class PartitionerParmetis<Mesh<Simplex<1,4> > >;
template class PartitionerParmetis<Mesh<Simplex<2,4> > >;
template class PartitionerParmetis<Mesh<Simplex<3,4> > >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 5 )
template class PartitionerParmetis<Mesh<Simplex<1,5> > >;
template class PartitionerParmetis<Mesh<Simplex<2,5> > >;
template class PartitionerParmetis<Mesh<Simplex<3,5> > >;
#endif

template class PartitionerParmetis<Mesh<Hypercube<1,1> > >;
template class PartitionerParmetis<Mesh<Hypercube<2,1> > >;
template class PartitionerParmetis<Mesh<Hypercube<3,1> > >;
template class PartitionerParmetis<Mesh<Hypercube<3,2> > >;

#endif // FEELPP_INSTANTIATION_MODE
}


