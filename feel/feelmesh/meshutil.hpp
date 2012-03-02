/*
 This file is part of the Feel library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#ifndef __MESH_UTIL_BASE__
#define __MESH_UTIL_BASE__

#include <vector>
#include <algorithm>
#include <set>

#include <feel/feelcore/feel.hpp>

#include <feel/feelmesh/entities.hpp>
#include <feel/feelmesh/bareitems.hpp>
#include <feel/feelmesh/marker.hpp>


/// \cond disabled
namespace Feel
{
/*!
  \brief Base utilities operating on meshes

  This file contains a set of base utilities used to test mesh entities or
  operate on them

*/

//! A locally used structure, not meant for general use
typedef std::map<BareFace, std::pair<size_type, size_type >, cmpBareItem<BareFace> > TempFaceContainer;

//! A locally used structure, not meant for general use
typedef std::map<BareEdge, std::pair<size_type, size_type>, cmpBareItem<BareEdge> > TempEdgeContainer;

template<typename MeshType>
struct TempEntityContainer
{
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<3> >,
                              mpl::identity<BareFace>,
                              mpl::identity<BareEdge> >::type::type entity_type;
    typedef std::map<entity_type, std::pair<size_type, size_type >, cmpBareItem<entity_type> > type;
};
template<typename Ele, int Dim>
struct MakeBareEntity
{};
template<typename Ele>
struct MakeBareEntity<Ele, 3>
{
    typedef BareFace entity_type;
    static const int numVertices =  Ele::GeoBShape::numVertices;

    MakeBareEntity( Ele const& ele )
        :
        _M_element( ele )
    {}
    entity_type
    operator()() const
    {
        entity_type bface;

        size_type i1 = ( _M_element.point( 0 ) ).id();
        size_type i2 = ( _M_element.point( 1 ) ).id();
        size_type i3 = ( _M_element.point( 2 ) ).id();
        if ( Ele::face_type::numVertices == 4 )
            {
                size_type i4 = ( _M_element.point( 3 ) ).id();
                bface = makeBareFace( i1, i2, i3, i4 ).first;
            }
        else
            {
                bface = makeBareFace( i1, i2, i3 ).first;
            }
        return bface;
    }
    entity_type
    operator()( uint16_type j ) const
    {
        entity_type bface;
        size_type i1 = _M_element.fToP( j, 0 );
        size_type i2 = _M_element.fToP( j, 1 );
        size_type i3 = _M_element.fToP( j, 2 );
        // go to global
        i1 = ( _M_element.point( i1 ) ).id();
        i2 = ( _M_element.point( i2 ) ).id();
        i3 = ( _M_element.point( i3 ) ).id();
        if ( numVertices == 4 )
        {
            size_type i4 = _M_element.fToP( j, 3 );
            i4 = ( _M_element.point( i4 ) ).id();
            bface = ( makeBareItem( i1, i2, i3, i4 ) ).first;
        }
        else
        {
            bface = ( makeBareItem( i1, i2, i3 ) ).first;
        }
        return bface;
    }
    Ele const& _M_element;
};

template<typename Ele>
struct MakeBareEntity<Ele, 2>
{
    typedef BareEdge entity_type;

    MakeBareEntity( Ele const& ele )
        :
        _M_element( ele )
    {}

    entity_type
    operator()() const
    {
        size_type i1 = ( _M_element.point( 0 ) ).id();
        size_type i2 = ( _M_element.point( 1 ) ).id();
        entity_type bface;
        bface = makeBareEdge( i1, i2 ).first;
        return bface;
    }

    entity_type
    operator()( uint16_type j ) const
    {
        entity_type bface;
        size_type i1 = _M_element.fToP( j, 0 );
        size_type i2 = _M_element.fToP( j, 1 );
        // go to global
        i1 = ( _M_element.point( i1 ) ).id();
        i2 = ( _M_element.point( i2 ) ).id();
        bface = ( makeBareItem( i1, i2 ) ).first;
        return bface;
    }
    Ele const& _M_element;
};

template<typename Ele>
struct MakeBareEntity<Ele, 1>
{
    typedef BarePoint entity_type;

    MakeBareEntity( Ele const& ele )
        :
        _M_element( ele )
    {}

    entity_type
    operator()() const
    {
        size_type i1 = ( _M_element.point( 0 ) ).id();
        Debug( 4015 )  << "[mesh1d::updateFaces] point index in face " << i1 << "\n";
        entity_type bface;
        bface = makeBarePoint( i1 ).first;
        return bface;
    }
    entity_type
    operator()( uint16_type j ) const
    {
        entity_type bface;
        size_type i1 = _M_element.point( j ).id();
        bface = makeBarePoint( i1 ).first;
        return bface;
    }
    Ele const& _M_element;
};


template<typename Ele, int Dim>
struct MakeBareEntityFromFace
{};

template<typename Ele>
struct MakeBareEntityFromFace<Ele, 3>
{
    typedef BareFace entity_type;
    static const int numVertices =  Ele::numVertices;

    MakeBareEntityFromFace( Ele const& ele )
        :
        _M_element( ele )
    {}

    entity_type
    operator()() const
    {
        entity_type bface;
        // go to global
        size_type i1 = ( _M_element.point( 0 ) ).id();
        size_type i2 = ( _M_element.point( 1 ) ).id();
        size_type i3 = ( _M_element.point( 2 ) ).id();
        if ( numVertices == 4 )
        {
            size_type i4 = ( _M_element.point( 3 ) ).id();
            bface = ( makeBareItem( i1, i2, i3, i4 ) ).first;
        }
        else
        {
            bface = ( makeBareItem( i1, i2, i3 ) ).first;
        }
        return bface;
    }
    Ele const& _M_element;
};

template<typename Ele>
struct MakeBareEntityFromFace<Ele, 2>
{
    typedef BareEdge entity_type;
    static const int numVertices =  Ele::numVertices;

    MakeBareEntityFromFace( Ele const& ele )
        :
        _M_element( ele )
    {}

    entity_type
    operator()() const
    {
        entity_type bface;
        // go to global
        size_type i1 = ( _M_element.point( 0 ) ).id();
        size_type i2 = ( _M_element.point( 1 ) ).id();
        bface = ( makeBareItem( i1, i2 ) ).first;
        return bface;
    }
    Ele const& _M_element;
};

/*
*******************************************************************************
FUNCTORS
*******************************************************************************/
//! \defgroup Test_Functors Some useful functors to be used for test mesh entities

/*! \ingroup Test_Functors
  \briefFunctor to check if a Point, Face or Edge is on the boundary.

  \precond It assumes that boundary points in RegionMesh are correctly set.
  \precond   the RegionMesh must export the typenames
  PointType, face_type and EdgeType.
*/
template <typename RegionMesh>
class EnquireBEntity
{
public:
    EnquireBEntity( RegionMesh const & mesh ) : pmesh( &mesh )
    {}
    ;
    typedef typename RegionMesh::face_type face_type;
    typedef typename RegionMesh::EdgeType EdgeType;
    typedef typename RegionMesh::PointType PointType;

    bool operator() ( face_type & face )
    {
        bool isboundary = true;
        for ( size_type k = 0;k < face_type::numVertices;++k )
        {
            isboundary &= face.point( k ).isOnBoundary();
        }
        return isboundary;
    }

    bool operator() ( EdgeType & edge )
    {
        bool isboundary = true;
        for ( size_type k = 0;k < EdgeType::numVertices;++k )
        {
            isboundary &= edge.point( k ).isOnBoundary();
        }
        return isboundary;
    }

    INLINE bool operator() ( PointType & point )
    {
        return point.isOnBoundary();
    }

private:
    EnquireBEntity()
    {}
    RegionMesh const * pmesh;
};

/**
 * \ingroup Test_Functors
 *
 *  Functor to check if a Face is on the boundary, by using the
 *  information contained in a TempFaceContainer produced by
 *  \c findBoundaryFaces(). It does not use the information contained in the
 *  mesh PointList, so it differs from EnquireBEntity.
 *
 *  \warning bfaces have been previously set by a call to FindBoundaryFaces
 */
template <typename RegionMesh>
class EnquireBFace
{
public:
    typedef typename RegionMesh::face_type face_type;
    typedef typename RegionMesh::face_type::GeoShape FaceShape;
    typedef typename TempEntityContainer<RegionMesh>::type container_type;

    typedef MakeBareEntityFromFace<face_type, RegionMesh::nDim> bare_entity_maker_type;
    typedef typename bare_entity_maker_type::entity_type bare_entity_type;

    EnquireBFace( RegionMesh const & mesh, container_type const & bfaces )
        :
        pmesh( &mesh ),
        pbfaces( &bfaces )
    {}

    bool operator() ( face_type const& f )
    {
        size_type i1, i2, i3, i4;
        bare_entity_type bface = bare_entity_maker_type( f )();
        return pbfaces->find( bface ) != pbfaces->end();
    }
private:
    EnquireBFace() {}

private:
    RegionMesh const * pmesh;
    container_type const * pbfaces;
};

//!\ingroup Test_Functors
/*! Functor to check if an edge is on the boundary, by using the information
  contained in a TempFaceContainer produced by findBoundaryEdges(). It does
  not use the information contained in the mesh PointList, so it differs
  from EnquireBEntity.

  \precond bedges have been previously set by a call to FindBoundaryEdges()

*/
template <typename RegionMesh>
class EnquireBEdge
{
public:
    typedef typename RegionMesh::EdgeType EdgeType;
    typedef typename RegionMesh::EdgeShape EdgeShape;

    EnquireBEdge( RegionMesh const & mesh, TempEdgeContainer const & bedges ) :
        pmesh( &mesh ), pbedges( &bedges )
    {}

    bool operator() ( EdgeType & f )
    {
        size_type i1, i2;
        BareEdge bedge;

        i1 = f.point( 0 ).id();
        i2 = f.point( 1 ).id();
        bedge = ( makeBareEdge( i1, i2 ) ).first;
        return pbedges->find( bedge ) != pbedges->end();
    }

private:
    EnquireBEdge()
    {}
    ;
    RegionMesh const * pmesh;
    TempEdgeContainer const * pbedges;
};

/**
 * @class EnquireBPoint
 *
 * Functor to check if a mesh entity with boundary indicator (for instance
 * a GeoPoint) is on the boundary, by enquiring its boundary flag.
 * \warning It assumes that boundary points are correctly set.
 *
 * @ingroup Test_Functors
 */
template <typename RegionMesh>
class EnquireBPoint
{
public:

    EnquireBPoint( RegionMesh & mesh ) : pmesh( &mesh )
    {}

    template<typename ME>
    bool operator() ( ME const& e )
    {
        return e.isOnBoundary();
    }
private:
    EnquireBPoint();

private:
    RegionMesh * pmesh;
};


/*
*******************************************************************************
EDGES/FACES FINDERS
*******************************************************************************
*/
/**
 * @brief Finds boundary faces.
 *
 * A low level routine, not meant to be called directly. It creates a
 * container with all the information needed to set up properly the
 * boundary faces connectivities.
 *
 * \param mesh A 3D mesh.
 *
 * \param NumInternalFaces. A reference to an integer returning the
 * number of internal faces found.
 *
 * \param bfaces This container will eventually contain a map whose
 * key are the BareFace corresponding to the boundary faces and the
 * data a pair of size_types: the size_type of the adjacent element
 * and the relative position of the face in the element.
 *
 * \param allFaces When this bool is set true the function will also
 * construct the set of internale faces, stored in intfaces.
 *
 * \param intfaces A container that will possibly contain a map whose
 * keys are the BareFace corresponding to an internal faces and the
 * data a pair of IDs: them ID of the two elements adjacent to the
 * face.
 *
 *\return Number of boundary faces found
 */
template <typename RegionMesh3D>
size_type
findFaces( const RegionMesh3D & mesh,
           typename TempEntityContainer<RegionMesh3D>::type & bfaces,
           size_type & numInternalFaces,
           typename TempEntityContainer<RegionMesh3D>::type & intfaces,
           bool allFaces = false )
{
    typename RegionMesh3D::element_type::GeoShape ele;
    typedef typename RegionMesh3D::elements_type elements_type;
    typedef typename RegionMesh3D::element_type element_type;

    typedef MakeBareEntity<element_type, RegionMesh3D::nDim> bare_entity_maker_type;
    typedef typename bare_entity_maker_type::entity_type bare_entity_type;

    typename TempEntityContainer<RegionMesh3D>::type::iterator fi;

    // clean first in case it has been alredy used

    bfaces.clear();
    if ( allFaces )
        intfaces.clear();
    numInternalFaces = 0;

    for ( typename elements_type::const_iterator iv = mesh.beginElement();
          iv != mesh.endElement(); ++iv )
    {
        bare_entity_maker_type const maker(  *iv );
        for ( uint16_type j = 0;j < mesh.numLocalFaces();++j )
        {
            bare_entity_type bface = maker( j );
            if ( ( fi = bfaces.find( bface ) ) == bfaces.end() )
            {
                bfaces.insert( std::make_pair( bface, std::make_pair( iv->id(), j ) ) );
            }
            else
            {
                if ( allFaces && bface.first > bface.second )
                {
                    intfaces.insert( ( std::make_pair( bface,
                                                       std::make_pair( iv->id(), j ) ) ) );
                }
                // counted twice: internal face
                bfaces.erase( fi );
                ++numInternalFaces;
            }
        }
    }
    return bfaces.size();
}

template <typename RegionMesh3D>
size_type findBoundaryFaces( const RegionMesh3D & mesh,
                             typename TempEntityContainer<RegionMesh3D>::type & bfaces,
                             size_type & numInternalFaces )
{
    typename TempEntityContainer<RegionMesh3D>::type dummy;
    return findFaces( mesh, bfaces, numInternalFaces, dummy, false );
}



//! Finds boundary edges.
/*!  A low level routine, not meant to be called directly. It creates a
  container with all the information needed to set up properly the boundary
  edges connectivities.

  \param mesh A 3D mesh.

  \param bedges This container contains a set with the BareEdge of the
  boundary edges.

  \return Number of boundary edges found.

  \pre The list of boundary faces must be correctly set.
*/
template <typename RegionMesh3D>
size_type findBoundaryEdges( const RegionMesh3D & mesh, TempEdgeContainer & bedges )
{
    size_type i1, i2;
    BareEdge bedge;
    typedef typename RegionMesh3D::face_type::GeoShape FaceShape;
    typedef typename RegionMesh3D::faces_type faces_type;


    FEELPP_ASSERT( mesh.isUpdatedForUse() ).error( "mesh has not been updated for use" );

    // clean first in case it has been alredy used
    bedges.clear();

    for ( typename faces_type::const_iterator ifa = mesh.beginFace();
          ifa != mesh.endFace(); ++ifa )
    {
        if ( ifa->isOnBoundary() )
        {
            for ( uint16_type j = 0;j < mesh.numLocalEdgesOfFace();++j )
            {
                i1 = FaceShape::eToP( j, 0 );
                i2 = FaceShape::eToP( j, 1 );
                // go to global
                i1 = ( ifa->point( i1 ) ).id();
                i2 = ( ifa->point( i2 ) ).id();
                bedge = ( makeBareEdge( i1, i2 ) ).first;
                bedges.insert( std::make_pair( bedge, std::make_pair( ifa->id(), j ) ) );
            }
        }
    }
    return bedges.size();
}

//! Finds all  edges.
/*!  A low level routine, not meant to be called directly. It creates a
  container with all the information needed to set up properly the edge
  connectivities.
  \param mesh A 3D mesh.

  \param bedges This container contains a set of  BareEdges for all mesh edges.

  \return Number of edges found.

*/

template <typename RegionMesh3D>
size_type findInternalEdges( const RegionMesh3D & mesh,
                             const TempEdgeContainer & boundary_edges,
                             TempEdgeContainer & internal_edges )
{
    size_type i1, i2;
    BareEdge bedge;
    typedef typename RegionMesh3D::element_type::GeoShape VolumeShape;
    typedef typename RegionMesh3D::Volumes Volumes;


    internal_edges.clear();


    for ( typename Volumes::const_iterator ifa = mesh.volumes().begin();
          ifa != mesh.volumes().end(); ++ifa )
    {
        for ( uint16_type j = 0;j < mesh.numLocalEdges();++j )
        {
            i1 = VolumeShape::eToP( j, 0 );
            i2 = VolumeShape::eToP( j, 1 );
            // go to global
            i1 = ( ifa->point( i1 ) ).id();
            i2 = ( ifa->point( i2 ) ).id();
            bedge = ( makeBareEdge( i1, i2 ) ).first;
            if ( boundary_edges.find( bedge ) == boundary_edges.end() )
                internal_edges.insert
                    ( std::make_pair( bedge, std::make_pair( ifa->id(), j ) ) );
        }
    }
    return internal_edges.size();
}
/*
*******************************************************************************
MARKERS HANDLERS
*******************************************************************************/
//! \defgroup marker_handlers Used to manage missing handlers

/**
   \ingroup marker_handlers

   \brief Sets the marker flag of a GeoElement of dimension greater one

   It gets the stronger marker of the GeoElement points. The marker
   hierarchy is defined in the marker.h file.  It returns a bool
   indicating if the flag has changed. If any of the vertices has an
   unset marker the result is an unset flag for the GeoElement.

   \warning It overrides the original marker flag.
*/
template <typename GeoElement>
flag_type
inheritStrongerMarker( GeoElement & fp )
{
    BOOST_STATIC_ASSERT( GeoElement::nDim > 0 );

    fp.setMarker( fp.point( 0 ).marker().value() );
#if 0
    for ( uint16_type j = 1;j < GeoElement::numVertices;++j )
        fp.setStrongerMarker( fp.point( j ).marker() );
#endif
    return fp.marker();

}


/**
   \ingroup marker_handlers
   \brief Sets the marker flag of a GeoElement of dimension greater one

   It gets the weaker marker of the GeoElement points. The marker
   hierarchy is defined in the marker.h file.  It returns a bool
   indicating if the flag has changed. If any of the vertices has an
   unset marker the result is an unset flag for the GeoElement.

   \warning It overrides the original marker flag.
*/
template <typename GeoElement>
flag_type
inheritWeakerMarker( GeoElement & fp )
{
    BOOST_STATIC_ASSERT( GeoElement::nDim > 0 );

    fp.setMarker( fp.point( 0 ).marker().value() );
#if 0
    for ( uint16_type j = 1;j < GeoElement::numVertices;++j )
        fp.setWeakerMarker( fp.point( j ).marker() );
#endif
    return fp.marker().value();

}
/**
   This routine tests if the topological descrption of boundary face
   is sane.  In particular all boundary edges must be adjacent to only
   2 surface elements and the orientation must be correct.

   \param mesh a mesh

   \param numBedges The function also returns the number of boundary
   edges in numBedges.

   \return It it returns 0 the test has been passed. If not it returns
   the number of of wrong boundary edges.

   \warning numBEdges is properly set only if the test has been passed.
*/
template <typename RegionMesh3D>
size_type testClosedDomain_Top( RegionMesh3D const & mesh, size_type & numBEdges )
{

    typedef std::set<BareEdge, cmpBareItem<BareEdge> > TempEdgeContainer2;
    TempEdgeContainer2 bedges;
    size_type i1, i2;
    BareEdge bedge;
    typename RegionMesh3D::BElementShape ele;
    typedef typename RegionMesh3D::faces_type faces_type;
    typedef typename RegionMesh3D::face_type face_type;
    TempEdgeContainer2::iterator ed;


    // clean first in case it has been alredy used

    typename faces_type::const_iterator iv = mesh.faces().begin();

    //for ( size_type k = 0;k < mesh.numBFaces();++k )
    for ( size_type k = 0;k < 0;++k )
    {
        FEELPP_ASSERT( iv != mesh.faces().end() )( k )( mesh.numBFaces() ).error(" Trying to get not existing face" );
        for ( uint16_type j = 0;j < face_type::numEdges;++j )
        {
            i1 = ele.eToP( j, 0 );
            i2 = ele.eToP( j, 1 );
            // go to global
            i1 = ( iv->point( i1 ) ).id();
            i2 = ( iv->point( i2 ) ).id();
            bedge = ( makeBareEdge( i1, i2 ) ).first;

            if ( ( ed = bedges.find( bedge ) ) == bedges.end() )
            {
                bedges.insert( bedge );
                ++numBEdges;
            }
            else
            {
                bedges.erase( ed );
            }
        }
        ++iv;
    }
    return bedges.size();
}
/*
*******************************************************************************
MARKERS FIXING
*******************************************************************************
*/

//! Check wether all markers of a the goemetry entities stored in a list are set
template <typename MeshEntityList>
bool checkMarkerSet( const MeshEntityList & list )
{
    typedef typename MeshEntityList::const_iterator C_Iter;
    bool ok( true );
    for ( C_Iter l = list.begin();l != list.end();++l )
        ok = ( ok & l->isMarkerSet() );
    return ok;
}

//! Sets the marker flag for all boundary edges by inheriting them from boundary points.
/*! The paradigm is that an edge <B>WHOSE MARKER HAS NOT ALREADY BEEN
  SET</B> will get the WEAKER marker flag among its VERTICES. For instance
  is a vertex is assigned to an Essential B.C and the other to a Natural
  B.C. the edge will get the flag related to the Natural B.C.

  /param mesh A mesh
  /param clog ostream to which the logging of the map of the newly assigned marked will be output
  /param err ostream to which error messages will be sent

  /todo better handling of flags: all function handling flags should be
  wrapped into a class
*/

template <typename RegionMesh>
void
setBEdgesMarker( RegionMesh & mesh, std::ostream & clog = std::cout,
                 std::ostream & /*err*/ = std::cerr, bool verbose = true )
{
    typename RegionMesh::EdgeType * fp = 0;
    unsigned int count( 0 );

    if ( verbose )
        Debug( 4100 ) << "NEW EDGE MARKER MAP" << "\n"
                      << " ID->New Marker" << "\n";

    for ( size_type k = 0; k < mesh.numBEdges(); ++k )
    {
        fp = &( mesh.edge( k ) );
        if ( fp->isMarkerUnset() )
        {
            inheritWeakerMarker( *fp );
            if ( verbose )
            {
                Debug( 4100 ) <<  fp->id() << " -> ";
                fp->printFlag( clog );
                Debug( 4100 ) <<  " ";
                if ( ++count % 3 == 0 )
                    Debug( 4100 ) <<  "\n";
            }
        }
    }
    if ( verbose )
        Debug( 4100 ) <<  "\n";
}


//! Sets the marker flag for all boundary faces by inheriting them from boundary points.
/*! The paradigm is that a face WHOSE MARKER HAS NOT ALREADY BEEN SET will
  get the WEAKER marker flag among its VERTICES. For instance if a vertex
  is assigned to a Natural B.C and the others to a Natural B.C. the face
  will get the flag related to the Natural B.C.

  /todo better handling of flags: all function handling flags should be
  wrapped into a class
*/
template <typename RegionMesh>
void
setBFacesMarker( RegionMesh & mesh, std::ostream & /*clog = std::cout*/,
                 std::ostream & /*err*/ = std::cerr, bool verbose = true )
{
    if ( verbose )
        Debug( 4100 ) <<  "NEW FACE MARKER MAP" << "\n"
                      << " ID->New Marker" << "\n";

    typename RegionMesh::location_face_iterator __fit = mesh.beginFaceOnBoundary();
    typename RegionMesh::location_face_iterator __fen = mesh.endFaceOnBoundary();
    for ( ; __fit != __fen; ++__fit )
    {
        if ( __fit->isMarkerUnset() )
        {
            typename RegionMesh::face_type __f ( *__fit );
            inheritWeakerMarker( __f );
            mesh.facesByLocation().replace( __fit, __f );
        }
    }
}

//! It sets the marker flag of boundary points, by inheriting it from boundary elements.
/*! The paradigm is that a point whose marker flag is unset will inherhit
  the strongest marker flag of the surrounding Boundary elements_type, with the
  convention that if the marker flag of one of the surrounding boundary
  elements is null is ignored.
*/
template <typename RegionMesh>
void
setBPointsMarker( RegionMesh & mesh,
                  std::ostream & clog = std::cout,
                  std::ostream& /*err*/ = std::cerr,
                  bool verbose = false )
{
    // First looks at points whose marker has already been set
    std::vector<bool> markset( mesh.storedPoints(), false );

    typedef typename RegionMesh::Points::iterator PointIterator;
    typedef typename RegionMesh::BElementShape BElementShape;

    std::vector<bool>::iterator pm = markset.begin();

    for ( PointIterator p = mesh.points().begin();
          p != mesh.points().end(); ++p )
        *( pm++ ) = p->isMarkerSet();

    typename RegionMesh::location_face_iterator __fit = mesh.beginFaceOnBoundary();
    typename RegionMesh::location_face_iterator __fen = mesh.endFaceOnBoundary();
    for ( ; __fit != __fen; ++__fit )
    {

        if ( __fit->isMarkerSet() )
        {
            typename RegionMesh::face_type __f ( *__fit );

            for ( size_type j = 0;j < BElementShape::numPoints;++j )
            {
                if ( !markset[ ( __f.point( j ).id() ) ] )
                    __f.point( j ).setStrongerMarker( __f.marker() );
            }
            mesh.facesByLocation().replace( __fit, __f );
        }
    }

    unsigned int count( 0 );
    if ( verbose )
    {
        Debug( 4100 ) <<  "**** NEW POINTS MARKERS **************" << "\n";
        Debug( 4100 ) <<  "id->marker    id->marker     id->marker" << "\n";
        pm = markset.begin();
        for ( PointIterator p = mesh.points().begin();
              p != mesh.points().end(); ++p )
        {
            if ( *pm++ )
            {
                Debug( 4100 ) <<  p->id() << " -> ";
                p->printFlag( clog );
                Debug( 4100 ) <<  " ";
                if ( ++count % 3 )
                    Debug( 4100 ) <<  "\n";
            }
        }
        Debug( 4100 ) <<  "\n";
    }
}
/*
*******************************************************************************
FIXING uint16_type AND COUNTERS
*******************************************************************************
*/
//! \brief Verifies if a list of mesh entities have the ID properly set.
/* More precisely, the id() must correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.

   \pre The template argument MeshEntityList must be a stl
   compliant container and its elements must have the method id().
*/
template <typename MeshEntityList>
bool checkIdnumber( const MeshEntityList & list )
{
    typedef typename MeshEntityList::const_iterator C_Iter;
    bool ok( true );
    size_type count = 0;
    for ( C_Iter l = list.begin();l != list.end();++l, ++count )
    {
        FEELPP_ASSERT( l->id() == count )( l->id() )( count ).error( "wrong id" );
        ok &= ( l->id() == count );
    }
    return ok;
}

//! \brief Fixes a a list of mesh entities so that the uint16_type is properly set.
/* \post  The id will correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.

   \pre The template argument MeshEntityList must be a stl
   compliant container and its elements must have the method size_type &id().
*/
template<typename MeshEntityList>
void
fixIdnumber( MeshEntityList & list )
{
    unsigned int count( 0 );
    typedef typename MeshEntityList::iterator Iter;
    for ( Iter l = list.begin() ;l != list.end(); ++l,++count )
    {
        l->setId( count );
    }
}

template<typename RegionMesh>
void
fixIdnumber( typename RegionMesh::Faces & list )
{
    unsigned int count( 0 );
    //typedef typename MeshEntityList::iterator Iter;
    typename RegionMesh::face_iterator __fit = list.begin();
    typename RegionMesh::face_iterator __fen = list.end();

    //for ( Iter l = list.begin() ;l != list.end(); ++l,++count )
    for ( ;__fit != __fen; ++__fit )
    {
        //l->setId( count );
        list.modify_key( __fit, lambda::_1 = count );
    }
}

/*! \brief Fixes boundary points counter
  It fix the boundary points counter by counting
  how many points have te boundary flag set.
  It also reset the Bpoints list.

  \pre It assumes that the points have the boundary flag corretly set
*/

template <typename RegionMesh>
void
setBPointsCounters( RegionMesh & mesh )
{

    unsigned int countBP( 0 );
    unsigned int countBV( 0 );

    mesh.boundaryPoints().clear();

    for ( size_type k = 0;k < mesh.numVertices();++k )
    {
        if ( mesh.isBoundaryPoint( k ) )
        {
            ++countBP;
            ++countBV;
        }
    }

    for ( size_type k = mesh.numVertices();k < mesh.numPoints();++k )
    {
        if ( mesh.isBoundaryPoint( k ) )
        {
            ++countBP;
        }
    }

    mesh.numBVertices() = countBV;
    mesh.setNumBPoints( countBP );
    mesh.boundaryPoints().reserve( countBP );

    for ( size_type k = 0;k < mesh.numPoints();++k )
    {
        if ( mesh.isBoundaryPoint( k ) )
            mesh.boundaryPoints().push_back( &mesh.point( k ) );
    }
}

/*
*******************************************************************************
BOUNDARY INDICATOR FIXING
*******************************************************************************
*/
//! It fixes boundary flag on points laying on boundary faces.
/*!
  \param mesh a mesh
  \param clog logging stream
  \param err error stream
  \param verbose If true you have a verbose output

  \pre mesh point list must exists and boundary face lsist  must have been set properly.
*/
template <typename RegionMesh>
void
fixBPoints( RegionMesh & mesh,
            std::ostream & /*clog*/ = std::cout,
            std::ostream & /*err*/ = std::cerr,
            bool verbose = true )
{
    ASSERT_PRE( mesh.numPoints() > 0, "The point list should not be empty" );
    ASSERT_PRE( mesh.numBelements_type() > 0,
                "The Belements_type list should not be empty" );

    typedef typename RegionMesh::Belements_type Belements_type;
    typedef typename RegionMesh::BElementShape BElementShape;

    if ( verbose )
        Debug( 4100 ) <<  "New BPoints Found " << "\n";

    typename RegionMesh::location_face_iterator __fit = mesh.beginFaceOnBoundary();
    typename RegionMesh::location_face_iterator __fen = mesh.endFaceOnBoundary();
    for ( ; __fit != __fen; ++__fit )
    {
        typename RegionMesh::face_type __f = *__fit;
        for ( size_type j = 0;j < BElementShape::numPoints;++j )
        {
            __f.point( j ).setOnBoundary( true );
        }
        mesh.facesByLocation().replace( __fit, __f );
    }
    // Fix now the number of vertices/points
    setBPointsCounters( mesh );
}

//!It makes sure that boundary edges are stored first
/*!
  \pre It assumes that boundary points are properly stored in the mesh
*/
template <typename RegionMesh>
void
setBoundaryEdgesFirst( RegionMesh & mesh )
{

    typedef typename RegionMesh::Edges Edges;
    // set the functor
    EnquireBEntity<RegionMesh > enquireBEdge( mesh );

    std::partition( mesh.edges().begin(), mesh.edges().end(), enquireBEdge );
    fixIdnumber( mesh.edges() );
}

//!It makes sure that boundary faces are stored first
/*!
  \pre It assumes that boundary points are properly stored in the mesh
*/
template <typename RegionMesh>
void
setBoundaryFacesFirst( RegionMesh & mesh )
{

    typedef typename RegionMesh::Faces Faces;
    // set the functor
    EnquireBEntity<RegionMesh> enquireBFace( mesh );

    std::partition( mesh.faces().begin(), mesh.faces().end(), enquireBFace );
    fixIdnumber( mesh.faces() );
}

//! Tests if boundary faces are stored first
/*! \return true if boundary faces are indeed stored first
  \pre It assumes that boundary points are set */
template <typename RegionMesh>
bool checkBoundaryFacesFirst( const RegionMesh & mesh )
{

    typedef typename RegionMesh::Faces Faces;

    // set the functor
    EnquireBEntity<RegionMesh> enquireBFace( mesh );
    typename RegionMesh::face_type * fp;
    bool ok( true );

    for ( size_type k = 0;k < mesh.numBelements_type();++k )
        ok = ok && enquireBFace( mesh.boundaryFace( k ) );
    for ( size_type k = mesh.numBelements_type();k < mesh.storedFaces();++k )
        ok = ok && ! enquireBFace( mesh.face( k ) );

    return ok;
}

//! Tests if boundary edges are stored first
/*! \return true if boundary edges are indeed stored first
  \pre It assumes that boundary points are set */
template <typename RegionMesh>
bool checkBoundaryEdgesFirst( const RegionMesh & mesh )
{

    typedef typename RegionMesh::Edges Edges;

    // set the functor
    EnquireBEntity<RegionMesh> enquireBEdge( mesh );
    typename RegionMesh::EdgeType * fp;
    bool ok( true );

    for ( size_type k = 0;k < mesh.numBEdges();++k )
        ok = ok && enquireBEdge( mesh.boundaryEdge( k ) );
    for ( size_type k = mesh.numBEdges();k < mesh.storedEdges();++k )
        ok = ok && ! enquireBEdge( mesh.edge( k ) );
    return ok;
}

/****************************************************************************
 *
 * UTILITIES TO VERIFY/CREATE FACES/EDGES
 *
 ****************************************************************************/

/**
   @brief It fixes boundary faces so that they are consistently
   numbered with volumes.

   An important step for building degrees of freedom on faces.  It also
   fixes other face related data.
   \param mesh a mesh
   \param err  ostream for error messages

   \param fixMarker If set to the true value all faces without a markerFlag set will inherit it from the points.

   \param clog ostream that will contain all information on what has been done
   Possible values are
   <ol>
   <li>NUM_FACES_MISMATCH</li>
   <li>FIXED_FACE_COUNTER</li>
   <li>BFACE_MISSING</li>
   <li>BFACE_STORED_MISMATCH</li>
   <li>BELEMENT_COUNTER_UNSET</li>
   <li>BFACE_STORED_MISMATCH</li>
   <li>FIXED_MAX_NUM_FACES</li>
   </ol>

   \param verbose if false nothing is written to clog

   \param numFaces It returns the number of faces found by the function

   \param bfaces_found It returns the number of boundary faces found by the function

   \param ext_container. If not NULL it is a pointer to an external map of bondary faces, already
   produced by a call to findBoundaryFaces(). This parameter may be used to save al lot of computational work, since
   findBoundaryFaces() is rather expensive.

   \pre Boundary faces list must be properly set.
*/
template <class RegionMesh3D>
bool
fixBoundaryFaces( RegionMesh3D & mesh,
                  std::ostream & clog,
                  std::ostream & /*err*/,
                  size_type & numFaces,
                  size_type & bfaces_found,
                  bool /*fixMarker */ ,
                  bool /*verbose*/ ,
                  typename TempEntityContainer<RegionMesh3D>::type * ext_container )
{
    typename RegionMesh3D::element_type::GeoShape ele;
    typedef typename RegionMesh3D::elements_type elements_type;
    typedef typename RegionMesh3D::element_type element_type;
    typedef typename RegionMesh3D::face_type face_type;

    typedef typename TempEntityContainer<RegionMesh3D>::type container_type;
    typedef typename container_type::iterator  container_iterator;

    typedef MakeBareEntityFromFace<face_type, RegionMesh3D::nDim> bare_entity_maker_type;
    typedef typename bare_entity_maker_type::entity_type bare_entity_type;


    container_type * bfaces;
    std::pair<size_type, size_type>info;
    size_type numInternalFaces;
    bool notfound( false );
    bool extcont( false );

    if ( extcont = ( ext_container != 0 ) )
    {
        bfaces = ext_container;
        bfaces_found = bfaces->size();
    }
    else
    {
        bfaces = new container_type;
        bfaces_found = findBoundaryFaces( mesh, *bfaces, numInternalFaces );
        numFaces = bfaces_found + numInternalFaces;
    }


    /*bool notEnough = mesh.storedFaces() < bfaces_found;*/

#if 0

    if ( notEnough )
    {
        Warning() <<  " number of B. Faces stored smaller" << "\n";
        Warning() <<  "         than the number of bfaces found  and build is not set"
                  << "\n";
        Warning() <<  "POSSIBLE ERROR" << "\n";
        clog << "BFACE_STORED_MISMATCH" << std::endl;
    }

    if ( mesh.numBelements_type() == 0 )
    {
        Warning() <<  " Boundary Element counter was not set" << "\n";
        Warning() <<  "I cannot proceed because the situation is ambiguous"
                  << "\n";
        Warning() <<  "Please check and eventually either: (a) call buildBoundaryFaces()" << "\n";
        Warning() <<  "or (b) set the correct number of bfaces in the mesh using mesh.numBelements_type()" << "\n";
        Warning() <<  "ABORT" << "\n";
        clog << "BELEMENT_COUNTER_UNSET" << std::endl;
    }

    if ( mesh.numBFaces() != bfaces_found )
    {
        Warning() <<  " B Face counter in mesh is set to "
                  << mesh.numBFaces();
        Warning() <<  " While I have found " << bfaces_found
                  << " B. elements_type in mesh" << "\n";
        Warning() <<  "Plese check... I continue anyway" << "\n";
        clog << "BFACE_COUNTER_MISMATCH" << std::endl;
    }

    if ( verbose )
    {
        Debug( 4100 ) <<  "**** Marker Flags for Fixed Boundary Faces ***" << "\n";
        Debug( 4100 ) <<  " (it only contains those that were fixed because unset !"
                      << "\n";
        Debug( 4100 ) <<  "id->marker   id->marker  id->marker" << "\n";
    }
#endif
    Debug( 4100 ) << "      total number of faces : " << mesh.faces().size() << "\n";
    Debug( 4100 ) << "number of faces on boundary : " << mesh.facesByLocation().count( ON_BOUNDARY ) << "\n";
    Debug( 4100 ) << "   number of internal faces : " << mesh.facesByLocation().count( INTERNAL ) << "\n";

    typename RegionMesh3D::location_face_const_iterator fit = mesh.beginFaceOnBoundary();
    typename RegionMesh3D::location_face_const_iterator fen = mesh.endFaceOnBoundary();
    //typename RegionMesh3D::face_const_iterator fit = mesh.beginFace();
    //typename RegionMesh3D::face_const_iterator fen = mesh.endFace();

    //std::cWarning() <<  "face ids:  ";
    //std::copy( fit, fen, std::cWarning() <<   lambda::bind( &face_type::id, lambda::_1 ) << " " );
    //std::cWarning() <<  "\n";
    for ( ; fit != fen; ++fit )
    {
        Debug( 4100 ) << "Inspecting face with " << fit->id() << "\n";

        bare_entity_maker_type const maker( *fit );
        bare_entity_type bface = maker();

        container_iterator fi = bfaces->find( bface );
        if ( fi == bfaces->end() )
        {
            notfound = true;
        }
        else
        {
            info = fi->second;
            // Element uint16_type
            size_type vol = info.first;
            // Element
            element_type const&  pv = mesh.element( vol );
            // The local uint16_type of face on element
            size_type j = info.second;


            face_type __f = *fit;
            // Reset face point definition to be consistent with face.
            for ( size_type k = 0;k < face_type::numPoints;++k )
            {
                __f.setPoint( k, pv.point( pv.fToP( j, k ) ) );
            }
            // Correct extra info
            __f.setConnection0( boost::make_tuple( boost::addressof( pv ), vol, j ) );
            Debug( 4100 ) << fit->id() << " -> ("
                          << "T: " << __f.ad_first() << ", "
                          << "L: " << __f.pos_first() << ")\n";
            if ( __f.isMarkerUnset() )
            {
                inheritWeakerMarker( __f );
                Debug( 4100 ) << fit->id() << " -> "
                              << "m: " << __f.marker() << "\n";
            }
            mesh.facesByLocation().replace( fit, __f );
            //mesh.facesById().replace( fit, __f );

            // Take out face from temporary container
            bfaces->erase( fi );
        }
    }

    if ( !extcont )
        delete bfaces;

    if ( notfound )
    {
        Warning() <<  " At least one boundary face has not been found on the list stored in RegionMesh3D\n";
        clog << "BFACE_MISSING" << std::endl;
    }

    Debug( 4100 ) << "  *****  END OF LIST ****" << "\n";

#if 0
    if ( mesh.numFaces() != numFaces )
    {
        Warning() <<  " faces counter in mesh should be " << numFaces
                  << "\n";
        Warning() <<  "         (bfaces->size()+numInternalFaces)" << "\n";
        Warning() <<  "         it is instead " << mesh.numFaces() << "\n";
        clog << "NUM_FACES_MISMATCH" << std::endl;
    }
#endif
    return true;
}

/*!
  @brief Builds faces

  This function may alternatively be used to build the compulsory boundary
  faces, all the mesh faces, or just add to an existing list of just boundary
  faces the internal ones.

  \param mesh A mesh

  \param clog Log file for information on the newly created markers

  \param err  Error stream

  \param buildbounary if true the function builds boundary faces

  \param buildinternal if true the function builds internal faces

  \param verbose. If true markerFrlags info is written on clog.

  \param numInternalFaces It returns the number of internal faces (only if ext_container is not provided!)

  \param bfaces_found It returns the number of boundary faces

  \param ext_container. If not NULL it is a pointer to an external map
  of bondary faces, already produced by a call to
  \c findBoundaryFaces(). This parameter may be used to save al lot of
  computational work, since findBoundaryFaces() is rather expensive.

  \pre If \c buildinternal=true and \c buildboundary=false the mesh must
  contain a proper list of boundary faces

  \note By setting buildinternal=true and buildboundary=true the
  function just fixes the counters with the number of faces in the
  mesh
*/
template <class RegionMesh3D>
bool buildFaces( RegionMesh3D & mesh,
                 std::ostream & clog,
                 std::ostream & /*err*/,
                 size_type & bfaces_found,
                 size_type & numInternalFaces,
                 bool buildboundary = true,
                 bool buildinternal = false,
                 bool verbose = false,
                 TempFaceContainer * ext_container = 0 )
{
    typename RegionMesh3D::VolumeShape ele;
    typedef typename RegionMesh3D::elements_type elements_type;
    typedef typename RegionMesh3D::element_type element_type;
    typedef element_type element_type;
    typedef typename RegionMesh3D::Faces Faces;
    typedef typename RegionMesh3D::face_type face_type;

    element_type * pv;
    TempFaceContainer* bfaces;
    TempFaceContainer::iterator fi;
    bool extcont( false );

    std::pair<size_type, size_type>info;
    size_type j, id;
    size_type vol;

    if ( extcont = ( ext_container != 0 ) )
    {
        bfaces = ext_container;
        bfaces_found = bfaces->size();
    }
    else
    {
        bfaces = new TempFaceContainer;
        bfaces_found = findBoundaryFaces( mesh, *bfaces, numInternalFaces );
    }

    if ( buildboundary )
        mesh.faces().clear();
    mesh.setNumBFaces( bfaces_found );
    if ( !buildinternal )
    {
        mesh.setMaxNumFaces( bfaces_found, false );
        mesh.numFaces() = numInternalFaces + bfaces_found;
    }
    else
    {
        mesh.setMaxNumFaces( numInternalFaces + bfaces_found, true );
    }

    face_type face;

    if ( buildboundary )
    {

        if ( verbose )
        {
            Debug( 4100 ) <<  "**** Marker Flags for Newly Created Boundary Faces ***"
                          << "\n";
            Debug( 4100 ) <<  "id->marker   id->marker  id->marker" << "\n";
        }

        for ( fi = bfaces->begin();fi != bfaces->end();++fi )
        {
            info = fi->second;
            vol = info.first; // Element size_type
            pv = &mesh.volume( vol ); // Element
            j = info.second;       // The local size_type of face on element

            for ( size_type k = 0;k < face_type::numPoints;++k )
                face.setPoint( k, pv->point( ele.fToP( j, k ) ) );
            // Add extra info
            face.setConnection0( boost::make_tuple( pv, vol, j ) );
            // Get marker value
            inheritWeakerMarker( face );
            id = mesh.addFace( face, true ).id();
            if ( verbose )
            {
                if ( id % 3 == 0 )
                    Debug( 4100 ) <<  "\n";
                Debug( 4100 ) <<  id << " -> ";
                face.printFlag( clog );
                Debug( 4100 ) <<  " ";
            }
        }
    }

    if ( !extcont )
        delete bfaces;
#if 0
    if ( ! buildinternal )
        return true;


    if ( !buildboundary )
    {
        if ( mesh.storedFaces() < mesh.numBFaces() )
        {
            Warning() <<  " mesh has not boundary faces, cannot just create internal ones!!!" << "\n";
            Warning() <<  "ABORT CONDITION" << "\n";
            return false;
        }
        else if ( mesh.storedFaces() > mesh.numBFaces() )
        {
            //mesh.faces().resize( mesh.numBFaces() );
        }
    }
#endif

    Debug( 4100 ) << "building internal faces\n";
    //
    // internal faces
    //
    typedef MakeBareEntityFromFace<face_type, RegionMesh3D::nDim> bare_entity_from_face_maker_type;
    typedef typename bare_entity_from_face_maker_type::entity_type bare_entity_from_face_type;

    BareItemsHandler<bare_entity_from_face_type> _be;
    std::pair<size_type, bool> e;
    std::pair<bare_entity_from_face_type, bool> _face;

    typename RegionMesh3D::face_const_iterator __fit = mesh.beginFace();
    typename RegionMesh3D::face_const_iterator __fen = mesh.endFace();
    for ( ;__fit != __fen; ++__fit )
    {

        bare_entity_from_face_maker_type const maker( *__fit );
        bare_entity_from_face_type bface = maker();
        _be.addIfNotThere( bface );
    }

    flag_type mm( mesh.marker() );
    size_type vid;

    typedef MakeBareEntity<element_type, RegionMesh3D::nDim> bare_entity_maker_type;
    typedef typename bare_entity_maker_type::entity_type bare_entity_type;

    for ( typename elements_type::iterator iv = mesh.beginElement();
          iv != mesh.endElement(); ++iv )
    {
        bare_entity_maker_type const maker( *iv );

        vid = iv->id();
        for ( size_type j = 0;j < mesh.numLocalFaces();j++ )
        {
            bare_entity_type bface = maker(j);
            e = _be.addIfNotThere( bface );
            if ( e.second )
            {
                // a new face It must be internal.
                for ( size_type k = 0;k < face_type::numPoints;++k )
                    face.setPoint( k, iv->point( ele.fToP( j, k ) ) );
                face.setConnection0( boost::make_tuple( boost::addressof( *iv ), vid, j ) );
                // gets the marker from the RegionMesh
                face.setMarker( mm );
                mesh.addFace( face, false ); //The id should be correct
            }
            else
            {
                if ( e.first > bfaces_found )  // internal
                {
                    typedef typename RegionMesh3D::face_type face_type;
                    typename RegionMesh3D::face_iterator __fit = mesh.facesById().find( face_type( e.first ) );
                    face_type __f = *__fit;
                    face.setConnection1( boost::make_tuple( boost::addressof( *iv ), vid, j ) );
                    mesh.faces().replace( __fit, __f );
                }
            }
        }
    }
    return true;
}

//! It builds edges.

/*! This function may alternatively be used to build the boundary edges,
  all the mesh faces, or just add the internal edges to an existing list of
  just boundary edges.

  \param mesh A mesh

  \param clog Log file for information on the newly created markers for boundary edges

  \param err  Error stream

  \param bedges_found Returns the number of boundary edges

  \param iedges_found Returns the number of internal edges

  \param buildbounary if true the function builds boundary edges

  \param buildinternal if true the function builds internal edges

  \param verbose. If true markerFlags info is written on clog.

  \param ext_container. If not NULL it is a pointer to an external map of bondary edges, already
  produced by a call to findBoundaryEdges(). This parameter may be used to save al lot of computational work, since
  findBoundaryEdges() is rather expensive.

  \pre If buildinternal=true and buildboundary=false the mesh must contain a proper list
  of boundary edges
  \pre The mesh must copntain a proper list of boundary faces

  \note By setting buildinternal=true and buildboundary=true the function just fixes the counters
  with the number of edges in the mesh
*/

template <typename RegionMesh3D>
bool buildEdges( RegionMesh3D & mesh,
                 std::ostream & /*clog*/,
                 std::ostream & /*err*/,
                 size_type & bedges_found,
                 size_type & iedges_found,
                 bool buildboundary = true,
                 bool buildinternal = false,
                 bool verbose = false,
                 TempEdgeContainer * ext_container = 0 )
{
    typedef typename RegionMesh3D::Volumes Volumes;
    typedef typename RegionMesh3D::Faces Faces;
    typedef typename RegionMesh3D::VolumeType VolumeType;
    typedef typename RegionMesh3D::VolumeShape VolumeShape;
    typedef typename RegionMesh3D::Edges Edges;
    typedef typename RegionMesh3D::EdgeType EdgeType;
    typedef typename RegionMesh3D::face_type face_type;
    typedef typename RegionMesh3D::FaceShape FaceShape;
    typename RegionMesh3D::VolumeType * pv;

    TempEdgeContainer * bedges;
    TempEdgeContainer iedges;
    std::pair<size_type, size_type>info;
    size_type j, id;
    size_type facuint16_type;


    bool extcont( false );


    if ( extcont = ( ext_container != 0 ) )
    {
        bedges = ext_container;
        bedges_found = bedges->size();
    }
    else
    {
        bedges = new TempEdgeContainer;
        bedges_found = findBoundaryEdges( mesh, *bedges );
    }

    iedges_found = findInternalEdges( mesh, *bedges, iedges );
    // free some memory if not needed!
    if ( !buildinternal )
        iedges.clear();
    if ( !buildboundary && buildinternal )
    {
        if ( mesh.storedEdges() < bedges_found )
        {
            Error() << "ERROR in buildedges(): mesh does not contain boundary edges" << "\n";
            Error() << "I need to set buildboundary=true" << "\n";
            Error() << "ABORT CONDITION" << "\n";
            return false;
        }
        else if ( mesh.storedEdges() > bedges_found )
        {
            mesh.edges().resize( bedges_found );
        }
    }
    mesh.setNumBEdges( bedges_found );
    mesh.numEdges() = ( bedges_found + iedges_found );
    if ( buildboundary )
        mesh.edges().clear();
    if ( buildboundary && ! buildinternal )
        mesh.setMaxNumEdges( bedges_found, false );
    if ( buildinternal )
        mesh.setMaxNumEdges( bedges_found + iedges_found, true );
    Debug( 4100 ) << "Building edges from scratch" << "\n";

    EdgeType edge;

    if ( buildboundary )
    {

        if ( verbose )
        {
            Debug( 4100 ) <<  "**** Marker Flags for Newly Created Boundary Edges ***"
                          << "\n";
            Debug( 4100 ) <<  "id->marker   id->marker   id->marker" << "\n";
        }

        // First boundary.
        for ( TempEdgeContainer::iterator ei = bedges->begin();
              ei != bedges->end();++ei )
        {
            info = ei->second;

            // Face uint16_type
            facuint16_type = info.first;

            // Face
            typedef typename RegionMesh3D::face_type face_type;
            typename RegionMesh3D::face_iterator __fit = mesh.faces().find( face_type( facuint16_type ) );
            FEELPP_ASSERT( __fit != mesh.faces().end() )( facuint16_type ).error( "invalid id for the faces" );

            // The local uint16_type of edge on face
            j = info.second;

            for ( size_type k = 0;k < EdgeType::numPoints;++k )
            {
                edge.setPoint( k, __fit->point( FaceShape::eToP( j, k ) ) );
            }

            // Get marker value inheriting from points
            inheritWeakerMarker( edge );

            id = mesh.addEdge( edge, true ).id();
#if 0
            if ( verbose )
            {
                if ( id % 3 == 0 )
                    Debug( 4100 ) <<  "\n";
                Debug( 4100 ) <<  id << " -> ";
                edge.printFlag( clog );
                Debug( 4100 ) <<  " ";
            }
#endif
        }

        if ( verbose )
            Debug( 4100 ) <<  "\n" << "  *****  END OF LIST OF BOUNDARY EDGES ****"
                          << "\n";

    }

    if ( !extcont )
        delete bedges;

    if ( !buildinternal )
    {
        return true;
    }



    // Now internal edges
    // free some memory

    for ( TempEdgeContainer::iterator ei = iedges.begin();
          ei != iedges.end();++ei )
    {
        info = ei->second;
        facuint16_type = info.first; // Volume uint16_type
        pv = &mesh.volume( facuint16_type ); // Volume that generated the edge
        j = info.second;       // The local uint16_type of edge on volume
        for ( size_type k = 0;k < EdgeType::numPoints;++k )
            edge.setPoint( k, pv->point( VolumeShape::eToP( j, k ) ) );
        edge.setMarker( mesh.marker() ); // Get marker value: that of the mesh
        mesh.addEdge( edge, false );
    }

    return true;
}


/*
*******************************************************************************
UTILITIES TO TRANSFORM A MESH
*******************************************************************************
*/
//! It builds a P2 mesh from P1 data.
/*!
  \author L.Formaggia.
  \version Version 1.0
  \pre All compulsory structures in mesh must have been already set: volumes and boundary faces.
  \pre Points list MUST have been dimensioned correctly!!!
  \note the function takes advantage of the fact that
*/
template <typename RegionMesh>
void
p1top2( RegionMesh & mesh, std::ostream & /* out */ = std::cout )
{
    typedef typename RegionMesh::node_type node_type;

    typedef typename RegionMesh::element_type::GeoShape GeoShape;
    typedef typename RegionMesh::BElementShape GeoBShape;
    ASSERT_PRE( GeoShape::numPoints > 4, "p1top2 ERROR: we need a P2 mesh" );

    Debug( 4100 ) <<  "Building P2 mesh points and connectivities from P1 data"
                  << "\n";


    typename RegionMesh::PointType * pp = 0;
    typename RegionMesh::EdgeType * pe = 0;
    typename RegionMesh::element_type * pv = 0;
//     typename RegionMesh::Belement_type * pbe = 0;
    typedef typename RegionMesh::elements_type elements_type;
    typedef typename RegionMesh::Belements_type Belements_type;

    BareItemsHandler<BareEdge> _be;
    std::pair<size_type, bool> _edgeid;
    size_type i1, i2, e_id;
    std::pair<BareEdge, bool> _edge;
    typename RegionMesh::element_type::GeoShape ele;
    Debug( 4100 ) <<  "Processing " << mesh.storedEdges() << " P1 Edges" << "\n";
    size_type nbe = mesh.numBEdges();
    for ( size_type j = 0; j < mesh.storedEdges();++j )
    {
        pe = & mesh.edge( j );
        i1 = ( pe->point( 0 ) ).id();
        i2 = ( pe->point( 1 ) ).id();
        // true for boundary points
        pp = & mesh.addPoint( j < nbe );
        node_type __n( 3 );
        __n = 0.5*( pe->point( 0 ).node()+pe->point( 1 ).node() );
        pp->setNode( __n );

        /*
          If we have set a marker for the boundary edge, that marker is
          inherited by the new created point. Otherwise the edge (and the new
          created point) gets the WORST marker among the two end Vertices
        */
        /*
          if the mesh file do not contain the edges (inria files), they are
          built in fixBoundaryEdges(...), but I suspect that this function
          does not attribute the right marker to the edges (maybe a problem in
          setWorseMarkerOfEntity, or something like that...)
          If you do #undef JFG : no change, if you do #define JFG, we do not
          consider the (wrong) marker of the edge to define the marker of the
          added node (I arbitrarily take the marker of the first node)
        */
        //#define JFG
        //#ifndef JFG
        // original version: DOES NOT work when the edges are not give in the mesh file
        if ( pe->isMarkerUnset() )
            inheritWeakerMarker( *pe );
        pp->setMarker( pe->marker() );
        //#else
        // temporary version that works when the edges are not given in the mesh file
        // pe->setMarker(mesh.point(i1).marker());
        // pp->setMarker(pe->marker());
        //#endif
        // pp->id()=++i; // Indexing from 1
        // if I storealso non B. edges:
        // points()[i].boundary()=
        // (edges()[j].point(1)).boundary() &&(edges()[j].point(2)).boundary();
        pe->setPoint( 2, pp ); //use overloaded version that takes a pointer
        _edge = makeBareEdge( i1, i2 );
        _edgeid = _be.addIfNotThere( _edge.first, pp->id() );
    }
    throw std::logic_error("p1top2 is disabled");
#if 0
    // Now the other edges, of which I do NOT build the global stuff
    if ( GeoShape::nDim == 3 )
    {
        size_type nbf = mesh.numBFaces();
        size_type nbv = GeoBShape::numVertices;
        Debug( 4100 ) <<  "Processing " << mesh.storedFaces() << " Face Edges"
                      << "\n";
        typename RegionMesh::face_iterator __fit = mesh.beginFace();
        typename RegionMesh::face_iterator __fen = mesh.endFace();
        //for ( size_type k = 0; k < mesh.storedFaces(); ++k )
        for ( ; __fit != __fen; ++__fit )
        {
            typename RegionMesh::face_type const* pbe = boost::addressof( *__fit );
            for ( size_type j = 0;j < mesh.numLocalEdgesOfFace();j++ )
            {
                i1 = GeoBShape::eToP( j, 0 );
                i2 = GeoBShape::eToP( j, 1 );
                i1 = ( pbe->point( i1 ) ).id();
                i2 = ( pbe->point( i2 ) ).id();
                _edge = makeBareEdge( i1, i2 );
                e_id = _be.id( _edge.first );
                if ( e_id != 0 )
                {
                    pp = &mesh.point( e_id );
                }
                else
                {
                    // new edge -> new Point
                    pp = &mesh.addPoint( k <= nbf );// true for boundary points
                    _edgeid = _be.addIfNotThere( _edge.first, pp->id() );
                    node_type __n( 3 );
                    __n = 0.5*( mesh.point( i1 ).node()+mesh.point( i2 ).node() );
                    pp->setNode( __n );
                    // If we have set a marker for the face, that marker is
                    // inherited by the new created point
                    pp->setMarker( pbe->marker() );
                }
                pbe->setPoint( nbv + j, pp );
            }
        }
    }
#endif

    Debug( 4100 ) <<  "Processing " << mesh.numElements() << " Mesh elements_type"
                  << "\n";
    size_type nev = GeoShape::numVertices;
    for ( size_type k = 0; k < mesh.numElements(); ++k )
    {
        pv = &mesh.element( k );
        for ( size_type j = 0;j < mesh.numLocalEdges();j++ )
        {
            i1 = ele.eToP( j, 0 );
            i2 = ele.eToP( j, 1 );
            i1 = ( pv->point( i1 ) ).id();
            i2 = ( pv->point( i2 ) ).id();
            _edge = makeBareEdge( i1, i2 );
            e_id = _be.id( _edge.first );
            if ( e_id != 0 )
            {
                pp = &mesh.point( e_id );
            }
            else
            {
                // cannot be on boundary if the mesh is proper!
                pp = &mesh.addPoint( false );
                _edgeid = _be.addIfNotThere( _edge.first, pp->id() );

                node_type __n( 3 );
                __n = 0.5*( mesh.point( i1 ).node()+mesh.point( i2 ).node() );
                pp->setNode( __n );

                pp->setMarker( pe->marker() );
            }
            pv->setPoint( nev + j, pp );
        }
    }
    /*=============================*/
    Debug( 4100 ) <<  " ******* Done Construction of P2 Mmesh *******" << "\n";
}



} // Feel

#include <feel/feelmesh/sphere.hpp>

namespace Feel
{
/**
 * Defines a Cartesian bounding box by the two
 * corner extremum.
 */
typedef std::pair<Point, Point> MeshBoundingBox;

/**
 * @returns two points defining a cartesian box that bounds the
 * mesh.  The first entry in the pair is the mininum, the second
 * is the maximim.
 */
template<typename MeshType>
inline
MeshBoundingBox
boundingBox (const MeshType& mesh)
{
    // processor bounding box with no arguments
    // computes the global bounding box
    return processorBoundingBox(mesh);
}


/**
 * Same, but returns a sphere instead of a box.
 */
template<typename MeshType>
inline
Sphere
boundingSphere (const MeshType& mesh)
{
    MeshBoundingBox bbox = boundingBox(mesh);

    const double diag = Feel::distance( bbox.second, bbox.first);
    const Point cent = Feel::middle( bbox.second, bbox.first);

    return Sphere (cent, .5*diag);

}

/**
 * @returns two points defining a cartesian box that bounds the
 * elements belonging to processor pid.  If no processor id is specified
 * the bounding box for the whole mesh is returned.
 */
template<typename MeshType>
inline
MeshBoundingBox
processorBoundingBox (const MeshType& mesh,
                      const size_type pid = invalid_size_type_value )
{
    FEELPP_ASSERT (mesh.numPoints() != 0).error( "mesh has no points" );

    Point min(1.e30,   1.e30,  1.e30);
    Point max(-1.e30, -1.e30, -1.e30);

    // By default no processor is specified and we compute
    // the bounding box for the whole domain.
    if (pid == invalid_size_type_value)
    {
        Debug( 4100 ) << "[processorBoundingBox] np pid given\n";
        for (unsigned int n=0; n<mesh.numPoints(); n++)
            for (unsigned int i=0; i<mesh.dimension(); i++)
            {
                min(i) = std::min(min(i), mesh.point(n)(i));
                max(i) = std::max(max(i), mesh.point(n)(i));
            }
    }
    // if a specific processor id is specified then we need
    // to only consider those elements living on that processor
    else
    {
        Debug( 4100 ) << "[processorBoundingBox] process bounding box on pid " << pid << "\n";
        typename MeshType::element_iterator it = mesh.beginElementWithProcessId( pid );
        typename MeshType::element_iterator en = mesh.endElementWithProcessId( pid );

        for (; it != en; ++it)
            for (unsigned int n=0; n< MeshType::element_type::numPoints; n++)
                for (unsigned int i=0; i<mesh.dimension(); i++)
                {
                    min(i) = std::min(min(i), mesh.point(n)(i));
                    max(i) = std::max(max(i), mesh.point(n)(i));
                }
    }
    for (unsigned int i=mesh.dimension(); i< min.node().size(); i++)
    {
        min( i ) = 0;
        max( i ) = 0;
    }
    Debug( 4100 ) << "[processorBoundingBox] min= " << min << "\n";
    Debug( 4100 ) << "[processorBoundingBox] max= " << max << "\n";
    const MeshBoundingBox ret_val(min, max);

    return ret_val;
}

/**
 * Same, but returns a sphere instead of a box.
 */
template<typename MeshType>
inline
Sphere
processorBoundingSphere (const MeshType& mesh,
                         const size_type pid = invalid_size_type_value)
{
  MeshBoundingBox bbox = processorBoundingBox(mesh,pid);

  const Real  diag = Feel::distance( bbox.second, bbox.first);
  const Point cent = Feel::middle(bbox.second, bbox.first);

  Debug( 4100 ) << "[processorBoundingSphere] processor " << mesh.comm().rank() << "\n";
  Debug( 4100 ) << "[processorBoundingSphere] center " << cent << "\n";
  Debug( 4100 ) << "[processorBoundingSphere] radius " << 0.5*diag << "\n";
  return Sphere (cent, .5*diag);
}


} // Feel

/// \endcond
#endif
