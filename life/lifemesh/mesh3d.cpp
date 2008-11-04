/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-14

  Copyright (C) 2005,2006 EPFL
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
   \file mesh3d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-14
 */
#include <life/lifemesh/mesh3d.hpp>
#include <life/lifemesh/meshutil.hpp>

namespace Life
{

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>::Mesh3D()
    :
    super_visitable(),
    super(),
    super_elements(),
    super_points(),
    super_faces(),
    super_edges(),
    _M_e2e()
{}

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>::Mesh3D( Mesh3D const & m )
    :
    super_visitable(),
    super( m ),
    super_elements( m ),
    super_points( m ),
    super_faces( m ),
    super_edges( m ),
    _M_e2e( m._M_e2e )
{}

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>::~Mesh3D()
{}

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>&
Mesh3D<GEOSHAPE>::operator=( Mesh3D const& m )
{
    if ( this != &m )
        {
            super::operator=( m );
            super_elements::operator=( m );
            super_points::operator=( m );
            super_faces::operator=( m );
            super_edges::operator=( m );

            _M_e2e = m._M_e2e;
        }
    return *this;
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::clear()
{
    this->elements().clear();
    this->points().clear();
    this->faces().clear();
    this->edges().clear();

    _M_e2e.resize( boost::extents[0][0] );
    LIFE_ASSERT( isEmpty() ).error( "all mesh containers should be empty after a clear." );
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                                            std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                                            face_permutation_type& permutation, mpl::bool_<true>)
{
    if ( numZeros == 0 )
        {
            for (uint16_type i = 0; i < def.size(); ++i)
                diff[i] = def[i] - cur[2-i];
        }

    std::vector<uint32_type>::iterator _id_it = find( diff.begin(),
                                                      diff.end(),
                                                      uint32_type(0) );

    uint16_type pos = distance( diff.begin(), _id_it);

    if ( numZeros == 0 )
        {
            if ( pos == 0 )
                permutation = face_permutation_type::ROTATION_CLOCKWISE;
            else
                permutation = face_permutation_type::ROTATION_ANTICLOCK;
        }
    else if ( numZeros == 1 )
        {
            if ( pos == 0 )
                permutation = face_permutation_type::REVERSE_HYPOTENUSE;
            else if ( pos == 1 )
                permutation = face_permutation_type::REVERSE_HEIGHT;
            else
                permutation = face_permutation_type::REVERSE_BASE;
        }
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                                            std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                                            face_permutation_type& permutation, mpl::bool_<false>)
{
    std::vector<uint32_type>::iterator _id_it = find( diff.begin(),
                                                      diff.end(),
                                                      uint32_type(0) );

    uint16_type pos = distance( diff.begin(), _id_it);

    if ( numZeros == 2 )
        {
            if ( pos == 0 )
                permutation = face_permutation_type::SECOND_DIAGONAL;
            else
                permutation = face_permutation_type::PRINCIPAL_DIAGONAL;
        }
    else if ( numZeros == 0 )
        {
            if ( cur[0] == def[1])
                if (cur[2] == def[3] )
                    permutation = face_permutation_type::REVERSE_BASE;
                else
                    permutation = face_permutation_type::ROTATION_CLOCKWISE;
            else if ( cur[0] == def[2] )
                if (cur[2] == def[0] )
                    permutation = face_permutation_type::ROTATION_ANTICLOCK;
                else
                    permutation = face_permutation_type::REVERSE_HEIGHT;
            else
                permutation = face_permutation_type::ROTATION_TWICE_CLOCKWISE;
        }
}

#if 0

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::updateFaces()
{
    boost::timer ti;
    face_type face;

    BareItemsHandler<BareFace> _be;
    std::pair<size_type, bool> e;

    size_type i1, i2, i3, i4;
    std::pair<BareFace, bool>_face;
    element_type ele;

    // First We check if we have already Faces stored
    if ( ! this->faces().empty() )
        {
            // dump all faces in the container, to maintain the correct numbering
            // if everything is correct the numbering in the bareface structure
            // will reflect the actual face numbering However, if I want to create
            // the internal faces I need to make sure that I am processing only the
            // boundary ones. So I resize the container!
            //if ( cf ) this->faces().resize( _numBFaces );

            std::pair<size_type, bool> _check;
            face_iterator __it = this->beginFace();
            face_iterator __en = this->endFace();
            for ( ;__it!=__en; )
                {
#if 0

                    i1 = ( __it->point( 0 ) ).id();
                    i2 = ( __it->point( 1 ) ).id();
                    i3 = ( __it->point( 2 ) ).id();
                    if ( face_type::numVertices == 4 )
                        {
                            i4 = ( __it->point( 3 ) ).id();
                            _face = makeBareFace( i1, i2, i3, i4 );
                        }
                    else
                        {
                            _face = makeBareFace( i1, i2, i3 );
                        }
                    _check = _be.addIfNotThere( _face.first );

#else
                    MakeBareEntity<face_type,nDim> baremaker( *__it );
                    _check = _be.addIfNotThere( baremaker() );

#endif
                    if ( _check.second )
                        Debug( 4015 ) << "added face with id " << __it->id () << "\n";
                    else
                        Debug( 4015 ) << "not added face with id " << __it->id ()
                                      << " was alread face with id = " << _check.first << "\n";
                    LIFE_ASSERT( _check.second )
                        (_check.first )
                        ( __it->id() )
                        ( i1 )( i2 )( i3 ).warn( "duplicated face" );
                    if ( !_check.second )
                        {
                            // here we get the next face or \c end()
                            size_type theid = __it->id();
                            __it = this->eraseFace( __it );
                            face_iterator __other = this->faces().find( face_type( _check.first ) );
                            LIFE_ASSERT( __other->id() != theid )
                                ( __other->id() )
                                ( theid ).error( "faces should have different ids " );


                        }
                    else
                        {
                            // ensure that item handler ids are in sync with
                            // faces ids
                            face_type __f = *__it;
                            __f.setId( _check.first );
                            this->faces().replace( __it, __f );
                            ++__it;
                        }
                }
        }
    Debug( 4015 ) << "[Mesh3D::updateFaces] adding faces : " << ti.elapsed() << "\n";
    ti.restart();


    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();
    for ( ;iv != en; ++iv )
        {
            element_type const& __element = *iv;
            size_type __element_id = __element.id();

            MakeBareEntity<element_type,nDim> baremaker( __element );
            for ( size_type j = 0; j < numLocalFaces(); j++ )
                {
                    Debug( 4015 ) << "------------------------------------------------------------\n";
                    Debug( 4015 ) << "Element id: " << iv->id() << " local face id: " << j << "\n";
#if 0
                    i1 = ele.fToP( j, 0 );
                    i2 = ele.fToP( j, 1 );
                    i3 = ele.fToP( j, 2 );
                    // go to global
                    i1 = ( iv->point( i1 ) ).id();
                    i2 = ( iv->point( i2 ) ).id();
                    i3 = ( iv->point( i3 ) ).id();
                    if ( face_type::numVertices == 4 )
                        {
                            i4 = ele.fToP( j, 3 );
                            i4 = ( iv->point( i4 ) ).id();
                            _face = makeBareFace( i1, i2, i3, i4 );
                        }
                    else
                        {
                            _face = makeBareFace( i1, i2, i3 );
                        }
                    e = _be.addIfNotThere( _face.first );
#else

                    e = _be.addIfNotThere( baremaker( j ) );

#endif


                    if ( e.second )
                        {
                            Debug( 4015 ) << "creating the face:" << e.first << "\n";

                            // set face id
                            face.setId( e.first );

                            // set the vertices of the face
                            for ( size_type k = 0;k < face_type::numPoints;++k )
                                face.setPoint( k, iv->point( ele.fToP( j, k ) ) );

                            // set the connection with the element
                            face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );
                            face.setOnBoundary( true );

                            // adding the face
                            face = this->addFace( face );
                            face_iterator __fit = this->faces().find( face_type( e.first ) );
                            LIFE_ASSERT( __fit != this->endFace() )( e.first )( iv->id() ).error( "invalid face iterator" );
                            this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref(*__fit) ) );

                            Debug( 4015 ) << "Adding [new] face info : \n";
                            Debug( 4015 ) << "element id: " << __element_id << "\n";
                            Debug( 4015 ) << "id: " << face.id() << "\n";
                            Debug( 4015 ) << "bdy: " << face.isOnBoundary() << "\n";
                            Debug( 4015 ) << "marker: " << face.marker() << "\n";
                            Debug( 4015 ) << "ad_first: " << face.ad_first() << "\n";
                            Debug( 4015 ) << "pos_first: " << face.pos_first() << "\n";
                            Debug( 4015 ) << "ad_second: " << face.ad_second() << "\n";
                            Debug( 4015 ) << "pos_second: " << face.pos_second() << "\n";

                        }
                    else
                        {
                            Debug( 4015 ) << "found the face:" << e.first << " in element " << __element_id << " and local face: " << j << "\n";

                            // look in the face table for the face
                            face_iterator __fit = this->faces().find( face_type( e.first ) );
                            LIFE_ASSERT( __fit != this->endFace() )( e.first ).error( "face is not in face container" );


                            //face_type __f = *__fit;

                            // the face could have been entered apriori given by
                            // the mesh generator, so just set the connection0
                            // properly .
                            if ( !__fit->isConnectedTo0() )
                                {
                                    Debug( 4015 ) << "[updateFaces][boundary] element: " << __element_id
                                                  << " face: " << j << " id: " << e.first << "\n";


                                    this->faces().modify( __fit,
                                                          detail::UpdateFaceConnection0<typename face_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j ) ) );

                                    this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref(*__fit) ) );

                                    Debug( 4015 ) << "adding [!isConnectedTo0] face info : \n";
                                    Debug( 4015 ) << "id: " << __fit->id() << "\n";
                                    Debug( 4015 ) << "bdy: " << __fit->isOnBoundary() << "\n";
                                    Debug( 4015 ) << "marker: " << __fit->marker() << "\n";
                                    Debug( 4015 ) << "ad_first: " << __fit->ad_first() << "\n";
                                    Debug( 4015 ) << "pos_first: " << __fit->pos_first() << "\n";
                                    Debug( 4015 ) << "ad_second: " << __fit->ad_second() << "\n";
                                    Debug( 4015 ) << "pos_second: " << __fit->pos_second() << "\n";

                                    // we must use replace here since we modify the
                                    // boundary property which is used for indexing
                                    //this->faces().replace( __fit, __f );

                                }
                            // we found an internal face
                            else
                                {

                                    face_type __f = *__fit;
#if 0
                                    this->faces().modify( __fit,
                                                          detail::UpdateFaceConnection1<typename face_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j ) ) );


                                    LIFE_ASSERT( __fit->isConnectedTo0() && __fit->isConnectedTo1() )
                                        ( __fit->isConnectedTo0() )( __fit->isConnectedTo1() ).error( "invalid face connection" );
#endif
                                    detail::UpdateFaceConnection1<typename face_type::element_connectivity_type> update1( boost::make_tuple( boost::addressof( __element ), __element_id, j ) );
                                    update1( __f );
                                    __f.setOnBoundary( false );
                                    this->faces().replace( __fit, __f );
#if 0
                                    // update neighbors for each element and replace in element container
                                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), __fit->proc_first() );
                                    this->elements().modify( elt1, update_element_neighbor_type( __fit->pos_first(),
                                                                                                 __fit->ad_second() ) );
                                    this->elements().modify( iv, update_element_neighbor_type( __fit->pos_second(),
                                                                                               __fit->ad_first() ) );
#endif
                                    // be careful: this step is crucial to set proper neighbor
                                    // connectivity
                                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), __fit->proc_first() );
                                    this->elements().modify( elt1, detail::UpdateFace<face_type>( boost::cref(*__fit) ) );
                                    this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref(*__fit) ) );

                                    Debug( 4015 ) << "adding face info : \n";
                                    Debug( 4015 ) << "id: " << __fit->id() << "\n";
                                    Debug( 4015 ) << "bdy: " << __fit->isOnBoundary() << "\n";
                                    Debug( 4015 ) << "marker: " << __fit->marker() << "\n";
                                    Debug( 4015 ) << "ad_first: " << __fit->ad_first() << "\n";
                                    Debug( 4015 ) << "pos_first: " << __fit->pos_first() << "\n";
                                    Debug( 4015 ) << "ad_second: " << __fit->ad_second() << "\n";
                                    Debug( 4015 ) << "pos_second: " << __fit->pos_second() << "\n";


                                }


                        }

                    LIFE_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).error( "invalid element face error" );
                }

        }

    Debug( 4015 ) << "[Mesh3D::updateFaces] element/face connectivity : " << ti.elapsed() << "\n";
    ti.restart();


    std::vector<size_type> _left(face_type::numVertices);
    std::vector<size_type> _right(face_type::numVertices);
    std::vector<uint32_type> _diff(face_type::numVertices);

    //determine permutation for the faces
    for ( face_iterator elt_it = this->beginFace();
          elt_it != this->endFace(); ++elt_it )
        {
            face_permutation_type permutation( face_permutation_type::IDENTITY );

            // if on boundary don't do anything
            if ( elt_it->isOnBoundary() )
                continue;

            for ( uint16_type i = 0; i < face_type::numVertices; ++i )
                {
                    _left[i] = elt_it->element0().point( elt_it->element0().fToP( elt_it->pos_first(), i ) ).id();
                    _right[i] = elt_it->element1().point( elt_it->element1().fToP( elt_it->pos_second(), i ) ).id();
                    _diff[i] = _left[i] - _right[i];
                }

            uint16_type _numZeros = count(_diff.begin(), _diff.end(), uint32_type(0));

            determineFacePermutation( _numZeros, _left, _right, _diff,
                                      permutation, mpl::bool_<( SHAPE == SHAPE_TETRA )>());

            if ( permutation.value() != face_permutation_type::IDENTITY )
                this->elements().modify( this->elementIterator( elt_it->ad_second(), elt_it->proc_second() ),
                                         detail::UpdateFacePermutation<face_permutation_type>( elt_it->pos_second(),
                                                                                               permutation ) );
        }

    boost::tie( iv, en ) = this->elementsRange();
    for ( ;iv != en; ++iv )
        {
            for ( size_type j = 0; j < numLocalFaces(); j++ )
                {
                    LIFE_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).warn( "invalid element face check" );
                }
        }
    Debug( 4015 ) << "[Mesh3D::updateFaces] element/face permutation : " << ti.elapsed() << "\n";

}
#endif // 0

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::updateEntitiesCoDimensionOnePermutation()
{
    boost::timer ti;
    std::vector<size_type> _left(face_type::numVertices);
    std::vector<size_type> _right(face_type::numVertices);
    std::vector<uint32_type> _diff(face_type::numVertices);

    //determine permutation for the faces
    for ( face_iterator elt_it = this->beginFace();
          elt_it != this->endFace(); ++elt_it )
        {
            face_permutation_type permutation( face_permutation_type::IDENTITY );

            // if on boundary don't do anything
            if ( elt_it->isOnBoundary() )
                continue;

            for ( uint16_type i = 0; i < face_type::numVertices; ++i )
                {
                    _left[i] = elt_it->element0().point( elt_it->element0().fToP( elt_it->pos_first(), i ) ).id();
                    _right[i] = elt_it->element1().point( elt_it->element1().fToP( elt_it->pos_second(), i ) ).id();
                    _diff[i] = _left[i] - _right[i];
                }

            uint16_type _numZeros = count(_diff.begin(), _diff.end(), uint32_type(0));

            determineFacePermutation( _numZeros, _left, _right, _diff,
                                      permutation, mpl::bool_<( SHAPE == SHAPE_TETRA )>());

            if ( permutation.value() != face_permutation_type::IDENTITY )
                this->elements().modify( this->elementIterator( elt_it->ad_second(), elt_it->proc_second() ),
                                         detail::UpdateFacePermutation<face_permutation_type>( elt_it->pos_second(),
                                                                                               permutation ) );
        }
    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();
    for ( ;iv != en; ++iv )
        {
            for ( size_type j = 0; j < numLocalFaces(); j++ )
                {
                    LIFE_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).warn( "invalid element face check" );
                }
        }
    Debug( 4015 ) << "[Mesh3D::updateFaces] element/face permutation : " << ti.elapsed() << "\n";
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::updateEntitiesCoDimensionTwo()
{
    boost::timer ti;
    BareItemsHandler<BareEdge> _be;
    std::pair<size_type, bool> e;
    _M_e2e.resize( boost::extents[this->numElements()][this->numLocalEdges()] );

    size_type vid, i1, i2;
    std::pair<BareEdge, bool> _edge;
    element_type ele;
    face_type bele;
    // First We check if we have already Edges stored
    if ( ! this->edges().empty() )
        {
            // dump first the existing edges, to maintain the correct numbering
            // if everything is correct the numbering in the bareedge
            // structure will reflect the actual edge numbering
            std::pair<size_type, bool> _check;
            for ( size_type j = 0; j < this->edges().size();++j )
                {
                    i1 = ( this->edge( j ).point( 0 ) ).id();
                    i2 = ( this->edge( j ).point( 1 ) ).id();
                    _edge = makeBareEdge( i1, i2 );
                    _check = _be.addIfNotThere( _edge.first );

                    LIFE_ASSERT( _check.second )( i1 )( i2 ).error( "Two identical Edges stored in EdgeList" );
                    LIFE_ASSERT( _check.first == this->edge( j ).id() )( _check.first )( this->edge( j ).id() ).error( "Edges in EdgeList have inconsistent id" );

                }
        }
    Debug( 4015 ) << "[Mesh3D::updateEdges] adding edges : " << ti.elapsed() << "\n";
    ti.restart();

    edge_type edg;

    if ( this->edges().empty() )
        {
            // We want that the first edges be those on the boundary, in order to obey the paradigm for
            // a Mesh3D
            const location_faces& location_index=this->faces().template get<detail::by_location>();
            location_face_iterator ifa;
            location_face_iterator efa;
            boost::tie( ifa, efa ) = location_index.equal_range( boost::make_tuple( true ) );
            for ( ; ifa!=efa; ++ifa )
                {
                    for ( uint16_type j = 0; j < face_type::numEdges; j++ )
                        {
                            i1 = bele.eToP( j, 0 );
                            i2 = bele.eToP( j, 1 );
                            // go to global
                            i1 = ( ifa->point( i1 ) ).id();
                            i2 = ( ifa->point( i2 ) ).id();
                            _edge = makeBareEdge( i1, i2 );
                            e = _be.addIfNotThere( _edge.first );
                            if ( e.second )
                                {
                                    // set edge id
                                    edg.setId( e.first );
                                    edg.setOnBoundary( true );

                                    for ( uint16_type k = 0;k < 2 + face_type::nbPtsPerEdge;k++ )
                                        edg.setPoint( k, ifa->point( k ) );
                                    inheritWeakerMarker( edg );
                                    this->addEdge( edg );
                                }
                        }
                }
        }
    Debug( 4015 ) << "[Mesh3D::updateEdges] adding edges : " << ti.elapsed() << "\n";
    ti.restart();

    std::map<size_type, edge_pair_type> _oriented_edges;
    typedef typename std::map<size_type, edge_pair_type>::iterator oe_iterator;

    for ( element_iterator elt_it = this->beginElement();
          elt_it != this->endElement(); ++elt_it )
        {
            vid = elt_it->id();
            for ( uint16_type j = 0; j < element_type::numEdges; ++j )
                {
                    uint16_type j1 = ele.eToP( j, 0 );
                    uint16_type j2 = ele.eToP( j, 1 );

                    // go to global
                    i1 = ( elt_it->point( j1 ) ).id();
                    i2 = ( elt_it->point( j2 ) ).id();
                    _edge = makeBareEdge( i1, i2 );
                    e = _be.addIfNotThere( _edge.first );
                    _M_e2e[ vid ][ j] = boost::make_tuple( e.first, 1 );

                    if ( e.second )
                        {
                            LIFE_ASSERT( e.first == this->numEdges() )( e.first )( this->numEdges() ).error( "invalid edge index" );
                            // set edge id
                            edg.setId( e.first );
                            edg.setOnBoundary( false );
                            for ( uint16_type k = 0; k < 2 + element_type::nbPtsPerEdge; k++ )
                                edg.setPoint( k, elt_it->point( k ) );
                            inheritWeakerMarker( edg );
                            this->addEdge( edg );
                        }

                    this->elements().modify( elt_it,
                                             detail::UpdateEdge<edge_type>( j, boost::cref( this->edge(e.first) ) ) );
                }
        }
    Debug( 4015 ) << "[Mesh3D::updateEdges] updating element/edges : " << ti.elapsed() << "\n";
    ti.restart();
    for ( element_iterator elt_it = this->beginElement();
          elt_it != this->endElement(); ++elt_it )
        {
            for ( size_type j = 0; j < (size_type)element_type::numEdges; ++j )
                {
                    size_type j1 = ele.eToP( j, 0 );
                    size_type j2 = ele.eToP( j, 1 );

                    // go to global
                    i1 = ( elt_it->point( j1 ) ).id();
                    i2 = ( elt_it->point( j2 ) ).id();

                    _edge = makeBareEdge( i1, i2 );
                    e = _be.addIfNotThere( _edge.first );
                    edge_pair_type _current = std::make_pair( i1, i2 );

                    edge_permutation_type permutation( edge_permutation_type::IDENTITY );

                    oe_iterator _edge_it = _oriented_edges.find( e.first );
                    if (  _edge_it != _oriented_edges.end() )
                        {
                            edge_pair_type _default = _edge_it->second;

                            LIFE_ASSERT( _default.first == _current.first ||
                                         _default.first == _current.second ).error("invalid edge index");

                            if ( _default.first != _current.first )
                                {
                                    permutation = edge_permutation_type::REVERSE_PERMUTATION;
                                    this->elements().modify( elt_it,
                                                             detail::UpdateEdgePermutation<edge_permutation_type>( j,
                                                                                                                   permutation ) );
                                }
                        }
                    else
                        {
                            _oriented_edges.insert( std::make_pair( e.first, _current ) );
                        }

                }
        }
    Debug( 4015 ) << "[Mesh3D::updateEdges] updating edges orientation : " << ti.elapsed() << "\n";
    ti.restart();
}

#if 0
template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::check() const
{
    Debug( 4015 ) << "[Mesh3D::check] numLocalFaces = " << numLocalFaces() << "\n";
    element_iterator iv = this->beginElementWithProcessId( Application::processId() );
    element_iterator en = this->endElementWithProcessId( Application::processId() );
    //boost::tie( iv, en ) = this->elementsRange();
    for ( ;iv != en; ++iv )
        {
            element_type const& __element = *iv;
            for ( size_type j = 0; j < numLocalFaces(); j++ )
                {
                    LIFE_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).error( "invalid element face check" );
                    Debug( 4015 ) << "------------------------------------------------------------\n";
                    Debug( 4015 ) << "Element : " << iv->id() << " face lid: " << j << " face gid:  " << iv->face( j ).id() << "\n";

                }

            size_type counter = 0;
            for (uint16_type ms=0; ms < __element.nNeighbors(); ms++)
                {
                    if ( __element.neighbor(ms).first != invalid_size_type_value )
                        ++counter;

                }
            Debug( 4015 ) << "[Mesh3D::check] element " << __element.id() << " number of neighbors: " << counter << "\n";
            LIFE_ASSERT( counter >= 1 )( __element.id() )( __element.nNeighbors() )( counter ).error( "invalid neighboring data" );

            for ( size_type j = 0; j < (size_type)element_type::numEdges; ++j )
                {
                    LIFE_ASSERT( iv->edgePtr( j ) )( j )( iv->id() ).error( "invalid element edge check" );
                    Debug( 4015 ) << "------------------------------------------------------------\n";
                    Debug( 4015 ) << "Element : " << iv->id() << " edge lid: " << j << " edge gid:  " << iv->edge( j ).id() << "\n";

                }

        }

    // faces check
    location_face_const_iterator itf = this->beginFaceOnBoundary();
    location_face_const_iterator ite = this->endFaceOnBoundary();
    for( ; itf != ite; ++ itf )
        {
            face_type const& __face = *itf;
            LIFE_ASSERT( __face.isConnectedTo0() )
                ( __face.ad_first() )( __face.pos_first() ).error( "invalid face" );

            LIFE_ASSERT( __face.element(0).facePtr( __face.pos_first() ) )
                ( __face.ad_first() )( __face.pos_first() )( __face.element(0).id() ).error( "invalid face in element" );

            LIFE_ASSERT( !__face.isConnectedTo1() )
                ( __face.ad_first() )( __face.pos_first() ).error( "invalid boundary face" );

        }
}
#endif // 0
//
// specialisations
//
template class Mesh3D<GeoEntity<Simplex<3, 1, 3> > >;
template class Mesh3D<GeoEntity<Simplex<3, 2, 3> > >;
template class Mesh3D<GeoEntity<SimplexProduct<3, 1, 3> > >;
template class Mesh3D<GeoEntity<SimplexProduct<3, 2, 3> > >;

} // Life

