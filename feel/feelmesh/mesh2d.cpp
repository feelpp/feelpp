/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-09

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file mesh2d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-09
 */
#include <iostream>

#include <feel/feelcore/feel.hpp>
#include <feel/feelmesh/mesh2d.hpp>


namespace Feel
{

template <typename GEOSHAPE>
Mesh2D<GEOSHAPE>::Mesh2D()
    :
    super_visitable(),
    super(),
    super_elements(),
    super_points(),
    super_faces()
{}

template <typename GEOSHAPE>
Mesh2D<GEOSHAPE>::Mesh2D( Mesh2D const & m )
    :
    super_visitable(),
    super( m ),
    super_elements( m ),
    super_points( m ),
    super_faces( m )
{}

template <typename GEOSHAPE>
Mesh2D<GEOSHAPE>::~Mesh2D()
{}

template <typename GEOSHAPE>
Mesh2D<GEOSHAPE>&
Mesh2D<GEOSHAPE>::operator=( Mesh2D const& m )
{
    if ( this != &m )
        {
            super::operator=( m );
            super_elements::operator=( m );
            super_points::operator=( m );
            super_faces::operator=( m );
        }
    return *this;
}

template <typename GEOSHAPE>
void
Mesh2D<GEOSHAPE>::clear()
{
    this->elements().clear();
    this->points().clear();
    this->faces().clear();
    FEEL_ASSERT( isEmpty() ).error( "all mesh containers should be empty after a clear." );
}

template <typename GEOSHAPE>
void
Mesh2D<GEOSHAPE>::updateEntitiesCoDimensionTwo()
{
    // no-op
}

#if 0

template <typename GEOSHAPE>
void
Mesh2D<GEOSHAPE>::updateFaces()
{
    boost::timer ti;

    BareItemsHandler<BareEdge> _be;
    std::pair<size_type, bool> e;

    std::pair<BareEdge, bool> _face;
    GEOSHAPE ele;
    // First We check if we have already Faces stored
    if ( ! this->faces().empty() )
    {
        // dump first faces, to maintain the correct numbering
        // if everything is correct the numbering in the bareface
        // structure will reflect the actual face numbering
        std::pair<size_type, bool> _check;
        for ( typename super_faces::face_iterator j = this->beginFace(); j != this->endFace();++j )
        {
            Debug( 4015 )  << "[mesh2d::updateFaces] get face " << j->id()
                           << " marker : " << j->marker() << "\n";
            size_type i1 = ( j->point( 0 ) ).id();
            size_type i2 = ( j->point( 1 ) ).id();
            Debug( 4015 )  << "[mesh2d::updateFaces] point index in face " << i1 << "\n";
            Debug( 4015 )  << "[mesh2d::updateFaces] point index in face " << i2 << "\n";
            _face = makeBareEdge( i1, i2 );
            _check = _be.addIfNotThere( _face.first, j->id() );
            Debug( 4015 )  << "[mesh2d::updateFaces] adding face ?  " << _check.second << "\n";
            Debug( 4015 )  << "[mesh2d::updateFaces] face id in container:  " << _check.first << "\n";
            _check = _be.addIfNotThere( _face.first, j->id() );
            Debug( 4015 )  << "[mesh2d::updateFaces] 2 adding face ?  " << _check.second << "\n";
            Debug( 4015 )  << "[mesh2d::updateFaces] 2 face id in container:  " << _check.first << "\n";

            //FEEL_ASSERT( !_check.second )( i1 )( i2 ).error( "Two identical Faces stored" );
            //FEEL_ASSERT( _check.first != j->id() )( _check.first )( j->id() ).error( "face has incorrect id" );
        }
    }
    Debug( 4015 ) << "[Mesh2D::updateFaces] adding edges : " << ti.elapsed() << "s\n";
    ti.restart();

    // --verbose
#if 0
    for ( typename super_elements::element_const_iterator elt = this->beginElement();
          elt != this->endElement(); ++elt )
    {
        element_type const& __element = *elt;
        Debug( 4015 ) << "element " << __element.id() << " proc " << __element.processId() << "\n";
    }
#endif

    for ( typename super_elements::element_iterator elt = this->beginElement();
          elt != this->endElement(); ++elt )
    {
        element_type const& __element = *elt;
        size_type __element_id = __element.id();
        Debug( 4015 ) << "[updateFaces] element: " << __element_id << "\n";

        for ( size_type j = 0;j < size_type(element_type::numEdges);j++ )
        {
            size_type i1 = __element.point( __element.eToP( j, 0 ) ).id();
            size_type i2 = __element.point( __element.eToP( j, 1 ) ).id();
            Debug( 4015 )  << "[mesh2d::updateFaces] point index in face " << j << ": " << i1 << "\n";
            Debug( 4015 )  << "[mesh2d::updateFaces] point index in face " << j << ": " << i2 << "\n";
            _face = makeBareEdge( i1, i2 );
            e = _be.addIfNotThere( _face.first );

            Debug( 4015 ) << "[updateFaces] (element: " << __element_id
                          << ") (local face: " << j << ") (face id in cont: " << e.first << ")\n";
            Debug( 4015 )  << "[mesh2d::updateFaces] face in container ?  " << e.second << "\n";
            // the edge was not in the constainer
            if ( e.second )
            {
                Debug( 4015 ) << "creating the face:" << e.first << "\n";

                face_type face;

                // set face id
                face.setId( e.first );

                // the face is not on the boundary, should already been in there
                face.setOnBoundary( true );

                // a new face It must be internal.
                for ( size_type k = 0;k < face_type::numPoints;++k )
                    face.setPoint( k, __element.point( __element.fToP( j, k ) ) );

                // set the connection with the element
                face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );

                // no-op anyway
#if 0
                this->elements().modify( elt, lambda::bind( &element_type::setEdgePermutation,
                                                            lambda::_1,
                                                            face.pos_first(),
                                                            edge_permutation_type(edge_permutation_type::IDENTITY) ) );
#endif

                Debug( 4015 ) << "[updateFaces][face not in container] adding face info : \n";
                Debug( 4015 ) << "[updateFaces][face not in container] id: " << face.id() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] bdy: " << face.isOnBoundary() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] ad_first: " << face.ad_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] pos_first: " << face.pos_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] ad_second: " << face.ad_second() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] pos_second: " << face.pos_second() << "\n";

                // insert face in container
                face_type const& new_face = this->addFace( face );

                FEEL_ASSERT( new_face.ad_first() == __element_id &&
                              new_face.pos_first() == j &&
                              elt->edgePermutation( j ).value() == edge_permutation_type::IDENTITY )
                    ( new_face.ad_first() )( new_face.pos_first() )
                    ( __element_id )( j )
                    ( elt->edgePermutation( j ).value() )( edge_permutation_type::IDENTITY ).error( "invalid face" );


                Debug( 4015 ) << "[updateFaces][face not in container] adding face info : \n";
                Debug( 4015 ) << "[updateFaces][face not in container] id: " << new_face.id() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] bdy: " << new_face.isOnBoundary() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] ad_first: " << new_face.ad_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] pos_first: " << new_face.pos_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] ad_second: " << new_face.ad_second() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] pos_second: " << new_face.pos_second() << "\n";

                this->elements().modify( this->elementIterator( new_face.ad_first(), 0 ),
                                         lambda::bind( &element_type::setEdge,
                                                       lambda::_1,
                                                       new_face.pos_first(),
                                                       boost::cref( new_face ) ) );

            }
            else
            {
                // the edge was already in the container, so it
                // must be an internal face

                Debug( 4015 ) << "[updateFaces][face] index= : " << e.first << "\n";

                typename super_faces::face_iterator __fit = this->faces().find( face_type( e.first ) );

                FEEL_ASSERT( __fit != this->endFace() ).error( "invalid face iterator" );

                //face_type __f = *__fit;

                FEEL_ASSERT( __fit->id() == __fit->id() )
                    ( __fit->id() )
                    ( __fit->id() ).error( "error in face construction" );

                if ( !__fit->isConnectedTo0() )
                {
                    Debug( 4015 ) << "[updateFaces][boundary] element: " << __element_id
                                  << " face: " << j << " id: " << e.first
                                  << " face_id: " << __fit->id() << " marker:  " << __fit->marker()
                                  << "\n";


                    //__fit->setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );
                    this->faces().modify( __fit,
                                          detail::UpdateFaceConnection0<typename edge_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j ) ) );


#if 0
                    this->elements().modify( elt, lambda::bind( &element_type::setEdgePermutation,
                                                                lambda::_1,
                                                                __fit->pos_first(),
                                                                edge_permutation_type(edge_permutation_type::IDENTITY) ) );
#endif

                }
                else
                {
                    FEEL_ASSERT( __fit->isConnectedTo0() )( __element_id )( j ).error( "invalid data : ad_first not filled" );

                    edge_type __edge = *__fit;

                    // make sure that the face is not on the boundary
                    //__fit->setConnection1( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );
                    __edge.setConnection1( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );
#if 0
                    this->faces().modify( __fit,
                                          detail::UpdateFaceConnection1<typename edge_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j ) ) );
#else
                    __edge.setOnBoundary( false );
                    this->faces().replace( __fit,__edge);
#endif

                    FEEL_ASSERT( __element_id != __fit->ad_first() )
                        ( __element_id )( __fit->id() )
                        ( __fit->ad_first() )
                        ( __fit->pos_first() )
                        ( __fit->ad_second() )
                        ( __fit->pos_second() )
                        ( j )
                        .error( "identical element id" );
                    FEEL_ASSERT( ublas::norm_2( __element.point( __element.fToP( j, 0 ) ).node() -
                                                 this->element(  __fit->ad_first(), 0 ).point( __element.fToP( __fit->pos_first(), 1 ) ).node() ) < 1e-10 )
                        ( j )
                        ( __fit->ad_first() )( __fit->pos_first() )
                        ( __fit->ad_second() )( __fit->pos_second() )
                        ( __element.point( __element.fToP( j, 0 ) ).node() )
                        ( this->element(  __fit->ad_first(), 0 ).point( __element.fToP( __fit->pos_first(), 1 ) ).node() )
                        ( __element.point( __element.fToP( j, 1 ) ).node() )
                        ( this->element(  __fit->ad_first(), 0 ).point( __element.fToP( __fit->pos_first(), 0 ) ).node() )
                        ( this->element(  __fit->ad_first(), 0 ).G() )
                        ( this->element(  __fit->ad_second(), 0 ).G() )
                        .warn( "inconsistent points" );

                    FEEL_ASSERT( __fit->isConnectedTo0() && __fit->isConnectedTo1() )
                        ( __fit->isConnectedTo0() )( __fit->isConnectedTo1() ).error( "invalid face connection" );

                    // now the edge is reversed
#if 0
                    this->elements().modify( elt, lambda::bind( &element_type::setEdgePermutation,
                                                                lambda::_1,
                                                                __fit->pos_second(),
                                                                edge_permutation_type( edge_permutation_type::REVERSE_PERMUTATION ) ) );
#else
                    this->elements().modify( elt,
                                             detail::UpdateEdgePermutation<edge_permutation_type>( __fit->pos_second(),
                                                                                                   edge_permutation_type( edge_permutation_type::REVERSE_PERMUTATION ) ) );
#endif
                    // update neighbors for each element and replace in element container
                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), 0 );
                    FEEL_ASSERT( elt1 != this->endElement() )( __fit->ad_first() ).error( "invalid element" );
                    this->elements().modify( elt1, update_element_neighbor_type( __fit->pos_first(),
                                                                                 __fit->ad_second() ) );
                    this->elements().modify( elt, update_element_neighbor_type( __fit->pos_second(),
                                                                                __fit->ad_first() ) );
                }

                edge_type const& __edge = *__fit;

                Debug( 4015 ) << "[updateFaces][face in container] adding face info : \n";
                Debug( 4015 ) << "[updateFaces][face in container] id: " << __edge.id() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] bdy: " << __edge.isOnBoundary() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] ad_first: " << __edge.ad_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] pos_first: " << __edge.pos_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] ad_second: " << __edge.ad_second() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] pos_second: " << __edge.pos_second() << "\n";



                {
                    this->elements().modify( this->elementIterator( __edge.ad_first(), 0 ),
                                             lambda::bind( &element_type::setEdge,
                                                           lambda::_1,
                                                           __edge.pos_first(),
                                                           boost::cref(__edge) ) );
                }

                if ( __edge.isConnectedTo1() )
                    {
                        this->elements().modify( this->elementIterator( __edge.ad_second(), 0 ),
                                                 lambda::bind( &element_type::setEdge,
                                                               lambda::_1,
                                                               __edge.pos_second(),
                                                               boost::cref(__edge) ) );
                    }

            }
        }
    }
    Debug( 4015 ) << "[Mesh2D::updateFaces] element/face connectivity : " << ti.elapsed() << "s\n";


}


template <typename GEOSHAPE>
void
Mesh2D<GEOSHAPE>::check() const
{
}
#endif // 0

template <typename GEOSHAPE>
void
Mesh2D<GEOSHAPE>::updateEntitiesCoDimensionOnePermutation()
{
    for ( typename super_elements::element_iterator elt = this->beginElement();
          elt != this->endElement(); ++elt )
        {
            for ( uint16_type j = 0;j < element_type::numEdges;j++ )
                {
                    if ( elt->face( j ).isConnectedTo1() &&
                         elt->face( j ).ad_second() == elt->id() )
                        {
                            this->elements().modify( elt,
                                                     detail::UpdateEdgePermutation<edge_permutation_type>( elt->face(j).pos_second(),
                                                                                                           edge_permutation_type( edge_permutation_type::REVERSE_PERMUTATION ) ) );
                        }
                }
        }
}

//
// explicit instantiation
//

template class Mesh2D<GeoEntity<Simplex<2, 1, 2> > >;
template class Mesh2D<GeoEntity<Simplex<2, 2, 2> > >;
template class Mesh2D<GeoEntity<Simplex<2, 3, 2> > >;
template class Mesh2D<GeoEntity<Simplex<2, 4, 2> > >;
template class Mesh2D<GeoEntity<Simplex<2, 5, 2> > >;
template class Mesh2D<GeoEntity<SimplexProduct<2, 1, 2> > >;
template class Mesh2D<GeoEntity<SimplexProduct<2, 2, 2> > >;

}
