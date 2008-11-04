/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-12

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
/**
   \file mesh1d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-12
 */
#include <life/lifecore/context.hpp>
#include <life/lifemesh/mesh1d.hpp>

namespace Life
{
template <typename GEOSHAPE>
Mesh1D<GEOSHAPE>::Mesh1D()
    :
    super_visitable(),
    super(),
    super_elements(),
    super_points(),
    super_faces()
{}

template <typename GEOSHAPE>
Mesh1D<GEOSHAPE>::Mesh1D( Mesh1D const & m )
    :
    super_visitable(),
    super( m ),
    super_elements( m ),
    super_points( m ),
    super_faces( m )
{}

template <typename GEOSHAPE>
Mesh1D<GEOSHAPE>::~Mesh1D()
{}

template <typename GEOSHAPE>
Mesh1D<GEOSHAPE>&
Mesh1D<GEOSHAPE>::operator=( Mesh1D const& m )
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
Mesh1D<GEOSHAPE>::clear()
    {
        this->elements().clear();
        this->points().clear();
        this->faces().clear();
        LIFE_ASSERT( isEmpty() ).error( "all mesh containers should be empty after a clear." );
    }

template <typename GEOSHAPE>
void
Mesh1D<GEOSHAPE>::updateEntitiesCoDimensionOnePermutation()
{
    // no-op
}
template <typename GEOSHAPE>
void
Mesh1D<GEOSHAPE>::updateEntitiesCoDimensionTwo()
{
    // no-op
}
#if 0


template <typename GEOSHAPE>
void
Mesh1D<GEOSHAPE>::updateFaces()
{
    BareItemsHandler<BarePoint> _be;
    std::pair<size_type, bool> e;

    std::pair<BarePoint, bool> _face;
    GEOSHAPE ele;
    // First We check if we have already Faces stored
    if ( ! this->faces().empty() )
    {
        // dump first faces, to maintain the correct numbering
        // if everything is correct the numbering in the bareface
        // structure will reflect the actual face numbering
        std::pair<size_type, bool> _check;
        for ( typename super_faces::face_iterator j = this->beginFace(); j != this->endFace(); ++j )
        {
            Debug( 4015 )  << "[mesh1d::updateFaces] get face " << j->id() << "\n";
            size_type i1 = ( j->point( 0 ) ).id();
            Debug( 4015 )  << "[mesh1d::updateFaces] point index in face " << i1 << "\n";
            _face = makeBarePoint( i1 );
            _check = _be.addIfNotThere( _face.first, j->id() );
            Debug( 4015 )  << "[mesh1d::updateFaces] adding face ?  " << _check.second << "\n";
            Debug( 4015 )  << "[mesh1d::updateFaces] face id in container:  " << _check.first << "\n";
            //LIFE_ASSERT( !_check.second )( i1 )( i2 ).error( "Two identical Faces stored" );
            //LIFE_ASSERT( _check.first != j->id() )( _check.first )( j->id() ).error( "face has incorrect id" );
        }
    }
    for ( typename super_elements::element_const_iterator elt = this->beginElement();
          elt != this->endElement(); ++elt )
    {
        element_type const& __element = *elt;
        Debug( 4015 ) << "element " << __element.id() << " proc " << __element.processId() << "\n";
    }
    for ( typename super_elements::element_iterator elt = this->beginElement();
          elt != this->endElement(); ++elt )
    {
        element_type const& __element = *elt;

        size_type __element_id = __element.id();
        for ( size_type j = 0;j < size_type(element_type::numVertices);j++ )
        {
            //size_type i1 = __element.point( __element.eToP( j, 0 ) ).id();

            size_type i1 = __element.point( j ).id();
            _face = makeBarePoint( i1 );
            e = _be.addIfNotThere( _face.first, i1 );
            Debug( 4015 ) << "[updateFaces] (element: " << __element_id << ") (local face: " << j << ") (face id in cont: " << e.first << ")\n";
            Debug( 4015 )  << "[mesh1d::updateFaces] face in container ?  " << e.second << "\n";

            // the edge was not in the constainer
            if ( e.second )
            {
                //face_type face( e.first, true, 0 );

                // set face id
                //face.setId( e.first );
                //face.setMarker( 0 );

                // the face is not on the boundary, should already been in there
                //face.setOnBoundary( false );

                // a new face It must be internal.
                //for ( size_type k = 0;k < face_type::numPoints;++k )
                //Debug( 4015 ) << "[updateFaces][face not in container before] id: " << face.id() << "\n";
                //face.setPoint( 0, __element.point( j ) );
                face_type face( __element.point( j ) );
                face.disconnect();
                face.setOnBoundary( true );

                // set the connection with the element
                face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );

                // insert face in container
                face_type const& new_face = this->addFace( face );

                Debug( 4015 ) << "[updateFaces][face not in container] adding face info : \n";
                Debug( 4015 ) << "[updateFaces][face not in container] id: " << new_face.id() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] marker: " << new_face.marker() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] bdy: " << new_face.isOnBoundary() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] ad_first: " << new_face.ad_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] pos_first: " << new_face.pos_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] ad_second: " << new_face.ad_second() << "\n";
                Debug( 4015 ) << "[updateFaces][face not in container] pos_second: " << new_face.pos_second() << "\n";

            }
            else
            {
                // the edge was already in the container, so it
                // must be an internal face

                Debug( 4015 ) << "[updateFaces][face] index= : " << e.first << "\n";

                typename super_faces::face_iterator __fit = this->faces().find( face_type( e.first ) );

                LIFE_ASSERT( __fit != this->endFace() )( e.first )( j )( i1 ).error( "invalid face iterator" );

                face_type __f = *__fit;

                LIFE_ASSERT( __f.id() == __fit->id() )
                    ( __f.id() )
                    ( __fit->id() ).error( "error in face construction" );

                if ( !__f.isConnectedTo0() )
                {

                    __f.setOnBoundary( true );
                    __f.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );

                }
                else
                {
                    LIFE_ASSERT( __f.isConnectedTo0() )( __element_id )( j ).error( "invalid data : ad_first not filled" );



                    // make sure that the face is not on the boundary
                    __f.setOnBoundary( false );
                    __f.setConnection1( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );

                    LIFE_ASSERT( __element_id != __f.ad_first() )
                        ( __element_id )( __f.ad_first() ).error( "identical element id" );
#if 0
                    LIFE_ASSERT( ublas::norm_2( __element.point( __element.fToP( j, 0 ) ).node() -
                                                 this->element(  __f.ad_first(), 0 ).point( __element.fToP( __f.pos_first(), 1 ) ).node() ) < 1e-10 )
                        ( __f.ad_first() )( __f.pos_first() )
                        ( __f.ad_second() )( __f.pos_second() )
                        ( __element.point( __element.fToP( j, 0 ) ).node() )
                        ( this->element(  __f.ad_first(), 0 ).point( __element.fToP( __f.pos_first(), 1 ) ).node() ) .error( "inconsistent points" );

#endif
                    LIFE_ASSERT( __f.isConnectedTo0() && __f.isConnectedTo1() )
                        ( __f.isConnectedTo0() )( __f.isConnectedTo1() ).error( "invalid face connection" );


                    // update neighbors for each element and replace in element container
                    element_iterator elt1 = this->elementIterator( __f.ad_first(), 0 );
                    LIFE_ASSERT( elt1 != this->endElement() )( __f.ad_first() ).error( "invalid element" );
                    this->elements().modify( elt1, update_element_neighbor_type( __f.pos_first(),
                                                                                 __f.ad_second() ) );
                    this->elements().modify( elt, update_element_neighbor_type( __f.pos_second(),
                                                                                __f.ad_first() ) );
                }

                this->faces().replace( __fit, __f );
                face_type const& __point = *__fit;

                Debug( 4015 ) << "[updateFaces][face in container] adding face info : \n";
                Debug( 4015 ) << "[updateFaces][face in container] id: " << __point.id() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] marker: " << __point.marker() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] bdy: " << __point.isOnBoundary() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] ad_first: " << __point.ad_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] pos_first: " << __point.pos_first() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] ad_second: " << __point.ad_second() << "\n";
                Debug( 4015 ) << "[updateFaces][face in container] pos_second: " << __point.pos_second() << "\n";

            }
        }
    }

} // Mesh1D::updateFaces

template <typename GEOSHAPE>
void
Mesh1D<GEOSHAPE>::check() const
{
}
#endif // 0

// Explicit instantiation
template class Mesh1D<GeoEntity<Simplex<1,1,1> > >;
template class Mesh1D<GeoEntity<SimplexProduct<1,1,1> > >;

}
