/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-07-21

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file createsubmesh.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-07-21
 */

#ifndef __createsubmesh_H
#define __createsubmesh_H 1

#include <boost/mpl/if.hpp>
#include <boost/mpl/identity.hpp>

namespace Feel
{

template <typename C, typename V, int T> class Mesh;

template <typename MeshType,typename IteratorRange>
class createSubmeshTool
{
public :

    typedef IteratorRange range_type;
    typedef typename boost::tuples::template element<0, range_type>::type idim_type;
    typedef typename boost::tuples::template element<1, range_type>::type iterator_type;

    static const uint16_type tag = MeshType::tag;
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;

    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > >,
                              mpl::identity< Mesh< Hypercube<mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > > >::type::type mesh_faces_type;
    typedef boost::shared_ptr<mesh_faces_type> mesh_faces_ptrtype;

    typedef typename mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_ELEMENTS> >,
                               mesh_type,
                               mesh_faces_type>::type mesh_build_type;

    typedef boost::shared_ptr<mesh_build_type> mesh_build_ptrtype;

    createSubmeshTool( boost::shared_ptr<MeshType> inputMesh,IteratorRange const& range )
        :
        M_mesh( inputMesh ),
        M_range( range )
    {}

    mesh_build_ptrtype
    build()
    {
        return build( mpl::int_<idim_type::value>() );
    }

private:

    mesh_ptrtype build( mpl::int_<MESH_ELEMENTS> /**/ );

    mesh_faces_ptrtype build( mpl::int_<MESH_FACES> /**/ );

    mesh_ptrtype M_mesh;
    range_type M_range;


};

template <typename MeshType,typename IteratorRange>
typename createSubmeshTool<MeshType,IteratorRange>::mesh_ptrtype
createSubmeshTool<MeshType,IteratorRange>::build( mpl::int_<MESH_ELEMENTS> /**/ )
{
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::face_type face_type;

    mesh_ptrtype newMesh( new mesh_type );

    //-----------------------------------------------------------//

    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        newMesh->addMarkerName( itMark.first,itMark.second.template get<0>(),itMark.second.template get<1>() );
    }

    //-----------------------------------------------------------//

    // How the nodes on this mesh will be renumbered to nodes
    // on the new_mesh.
    std::vector<size_type> new_node_numbers ( M_mesh->numPoints() );
    std::vector<size_type> new_vertex ( M_mesh->numPoints() );

    std::fill ( new_node_numbers.begin(),
                new_node_numbers.end(),
                invalid_size_type_value );

    std::fill ( new_vertex.begin(),
                new_vertex.end(),
                0 );


    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_faces = 0;

    //-----------------------------------------------------------//

    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = M_range;

    for ( ; it != en; ++ it )
    {
        element_type const& old_elem = *it;

        // copy element so that we can modify it
        element_type new_elem = old_elem;

        /*
        // get element markers
        new_elem.setMarker(old_elem.marker().value());
        new_elem.setMarker2(old_elem.marker2().value());
        new_elem.setMarker3(old_elem.marker3().value());
        */
        // Loop over the nodes on this element.
        for ( unsigned int n=0; n < old_elem.nPoints(); n++ )
        {
            //FEELPP_ASSERT (old_elem.point( n ).id() < new_node_numbers.size()).error( "invalid point id()" );

            if ( new_node_numbers[old_elem.point( n ).id()] == invalid_size_type_value )
            {
                new_node_numbers[old_elem.point( n ).id()] = n_new_nodes;

                Debug( 4015 ) << "[Mesh<Shape,T>::createSubmesh] insert point " << old_elem.point( n ) << "\n";

                point_type pt( old_elem.point( n ) );
                pt.setId( n_new_nodes );

                // Add this node to the new mesh
                newMesh->addPoint ( pt );

                Debug( 4015 ) << "[Mesh<Shape,T>::createSubmesh] number of  points " << newMesh->numPoints() << "\n";

                // Increment the new node counter
                n_new_nodes++;

                if ( n < element_type::numVertices )
                {
                    FEELPP_ASSERT( new_vertex[old_elem.point( n ).id()] == 0 ).error( "already seen this point?" );
                    new_vertex[old_elem.point( n ).id()]=1;
                }
            }

            // Define this element's connectivity on the new mesh
            FEELPP_ASSERT ( new_node_numbers[old_elem.point( n ).id()] < newMesh->numPoints() ).error( "invalid connectivity" );

            Debug( 4015 ) << "[Mesh<Shape,T>::createSubmesh] adding point old(" << old_elem.point( n ).id()
                          << ") as point new(" << new_node_numbers[old_elem.point( n ).id()]
                          << ") in element " << new_elem.id() << "\n";

            new_elem.setPoint( n, newMesh->point( new_node_numbers[old_elem.point( n ).id()] ) );

        } // for (unsigned int n=0 ... )

        // set id of element
        new_elem.setId ( n_new_elem );

        // increment the new element counter
        n_new_elem++;


        // Add an equivalent element type to the new_mesh
        newMesh->addElement( new_elem );

        // Maybe add faces for this element
        for ( unsigned int s=0; s<old_elem.numTopologicalFaces; s++ )
        {
            if ( !old_elem.facePtr( s ) ) continue;

            // only add face on the boundary: they have some data
            // (boundary ids) which cannot be retrieved otherwise
            //if ( old_elem.neighbor(s) == invalid_size_type_value )
            size_type global_face_id = old_elem.face( s ).id();

            if ( M_mesh->hasFace( global_face_id ) )
            {

                // get the corresponding face
                face_type const& old_face = old_elem.face( s );
                face_type new_face = old_face;

                // disconnect from elements of old mesh,
                // the connection will be redone in
                // \c updateForUse()
                new_face.disconnect();
                //std::cout << "disconnect face\n";
                // update points info

                //this line is very important!!!!!!!!!!
                //updateForUse put false for internalfaces
                new_face.setOnBoundary( true );

                /*
                if (old_face.isOnBoundary()) new_face.setOnBoundary( true );
                else new_face.setOnBoundary( false );

                new_face.setMarker(old_face.marker().value());
                new_face.setMarker2(old_face.marker2().value());
                new_face.setMarker2(old_face.marker3().value());
                */
                for ( uint16_type p = 0; p < new_face.nPoints(); ++p )
                {
                    new_face.setPoint( p, newMesh->point( new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] ) );
                }

                new_face.setId( n_new_faces++ );

                // add it to the list of faces
                newMesh->addFace( new_face );
            }

        } // for (unsigned int s=0 ... )

    } //  for( ; it != en; ++ it )


    newMesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0 ) );

    Debug( 4015 ) << "[Mesh<Shape,T>::createSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    newMesh->updateForUse();

    return newMesh;
}

template <typename MeshType,typename IteratorRange>
typename createSubmeshTool<MeshType,IteratorRange>::mesh_faces_ptrtype
createSubmeshTool<MeshType,IteratorRange>::build( mpl::int_<MESH_FACES> /**/ )
{
    mesh_faces_ptrtype newMesh( new mesh_faces_type );

    //-----------------------------------------------------------//

    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        newMesh->addMarkerName( itMark.first,itMark.second.template get<0>(),itMark.second.template get<1>() );
    }

    //-----------------------------------------------------------//

    typedef typename mesh_faces_type::element_type new_element_type;

    std::vector<size_type> new_node_numbers ( M_mesh->numPoints() );
    std::vector<size_type> new_vertex ( M_mesh->numPoints() );

    std::fill ( new_node_numbers.begin(),
                new_node_numbers.end(),
                invalid_size_type_value );

    std::fill ( new_vertex.begin(),
                new_vertex.end(),
                0 );

    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_faces = 0;


    //-----------------------------------------------------------//

    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = M_range;

    for ( ; it != en; ++ it )
    {
        // create a new element
        //element_type const& old_elem = *it;
        auto oldElem = *it;

        // copy element so that we can modify it
        new_element_type newElem;// = oldElem;

        // get element markers
        newElem.setMarker( oldElem.marker().value() );
        newElem.setMarker2( oldElem.marker2().value() );
        newElem.setMarker3( oldElem.marker3().value() );

        //std::cout << "\n oldElem.nPoints " << oldElem.nPoints();
        // Loop over the nodes on this element.
        for ( unsigned int n=0; n < oldElem.nPoints(); n++ )
        {

            if ( new_node_numbers[oldElem.point( n ).id()] == invalid_size_type_value )
            {
                new_node_numbers[oldElem.point( n ).id()] = n_new_nodes;

                //Debug( 4015 ) << "[Mesh<Shape,T>::createP1mesh] insert point " << old_elem.point(n) << "\n";

                typename mesh_faces_type::point_type pt( oldElem.point( n ) );
                pt.setId( n_new_nodes );

                // Add this node to the new mesh
                newMesh->addPoint( pt );

                //Debug( 4015 ) << "[Mesh<Shape,T>::createSubmesh] number of  points " << new_mesh->numPoints() << "\n";

                // Increment the new node counter
                n_new_nodes++;

                if ( n < new_element_type::numVertices )
                {
                    FEELPP_ASSERT( new_vertex[oldElem.point( n ).id()] == 0 ).error( "already seen this point?" );
                    new_vertex[oldElem.point( n ).id()]=1;
                }

            }

            newElem.setPoint( n, newMesh->point( new_node_numbers[oldElem.point( n ).id()] ) );

        } // end for n

        // set id of element
        newElem.setId ( n_new_elem );

        // increment the new element counter
        n_new_elem++;

        // Add an equivalent element type to the new_mesh
        newMesh->addElement( newElem );

    } // end for it

    newMesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0 ) );

    //Debug( 4015 ) << "[Mesh<Shape,T>::createP1mesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    newMesh->updateForUse();



    return newMesh;
}





template <typename MeshType,typename IteratorRange>
typename createSubmeshTool<MeshType,IteratorRange>::mesh_build_ptrtype
createSubmesh( boost::shared_ptr<MeshType> inputMesh,IteratorRange const& range )
{
    createSubmeshTool<MeshType,IteratorRange> cSmT( inputMesh,range );
    return cSmT.build();
}


} // namespace Feel

#endif // createsubmesh
