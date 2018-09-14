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
        M_element( ele )
    {}
    entity_type
    operator()() const
    {
        entity_type bface;

        size_type i1 = ( M_element.point( 0 ) ).id();
        size_type i2 = ( M_element.point( 1 ) ).id();
        size_type i3 = ( M_element.point( 2 ) ).id();

        if ( Ele::face_type::numVertices == 4 )
        {
            size_type i4 = ( M_element.point( 3 ) ).id();
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
        size_type i1 = M_element.fToP( j, 0 );
        size_type i2 = M_element.fToP( j, 1 );
        size_type i3 = M_element.fToP( j, 2 );
        // go to global
        i1 = ( M_element.point( i1 ) ).id();
        i2 = ( M_element.point( i2 ) ).id();
        i3 = ( M_element.point( i3 ) ).id();

        if ( numVertices == 4 )
        {
            size_type i4 = M_element.fToP( j, 3 );
            i4 = ( M_element.point( i4 ) ).id();
            bface = ( makeBareItem( i1, i2, i3, i4 ) ).first;
        }

        else
        {
            bface = ( makeBareItem( i1, i2, i3 ) ).first;
        }

        return bface;
    }
    Ele const& M_element;
};

template<typename Ele>
struct MakeBareEntity<Ele, 2>
{
    typedef BareEdge entity_type;

    MakeBareEntity( Ele const& ele )
        :
        M_element( ele )
    {}

    entity_type
    operator()() const
    {
        size_type i1 = ( M_element.point( 0 ) ).id();
        size_type i2 = ( M_element.point( 1 ) ).id();
        entity_type bface;
        bface = makeBareEdge( i1, i2 ).first;
        return bface;
    }

    entity_type
    operator()( uint16_type j ) const
    {
        entity_type bface;
        size_type i1 = M_element.fToP( j, 0 );
        size_type i2 = M_element.fToP( j, 1 );
        // go to global
        i1 = ( M_element.point( i1 ) ).id();
        i2 = ( M_element.point( i2 ) ).id();
        bface = ( makeBareItem( i1, i2 ) ).first;
        return bface;
    }
    Ele const& M_element;
};

template<typename Ele>
struct MakeBareEntity<Ele, 1>
{
    typedef BarePoint entity_type;

    MakeBareEntity( Ele const& ele )
        :
        M_element( ele )
    {}

    entity_type
    operator()() const
    {
        size_type i1 = ( M_element.point( 0 ) ).id();
        DVLOG(2)  << "[mesh1d::updateFaces] point index in face " << i1 << "\n";
        entity_type bface;
        bface = makeBarePoint( i1 ).first;
        return bface;
    }
    entity_type
    operator()( uint16_type j ) const
    {
        entity_type bface;
        size_type i1 = M_element.point( j ).id();
        bface = makeBarePoint( i1 ).first;
        return bface;
    }
    Ele const& M_element;
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
        M_element( ele )
    {}

    entity_type
    operator()() const
    {
        entity_type bface;
        // go to global
        size_type i1 = ( M_element.point( 0 ) ).id();
        size_type i2 = ( M_element.point( 1 ) ).id();
        size_type i3 = ( M_element.point( 2 ) ).id();

        if ( numVertices == 4 )
        {
            size_type i4 = ( M_element.point( 3 ) ).id();
            bface = ( makeBareItem( i1, i2, i3, i4 ) ).first;
        }

        else
        {
            bface = ( makeBareItem( i1, i2, i3 ) ).first;
        }

        return bface;
    }
    Ele const& M_element;
};

template<typename Ele>
struct MakeBareEntityFromFace<Ele, 2>
{
    typedef BareEdge entity_type;
    static const int numVertices =  Ele::numVertices;

    MakeBareEntityFromFace( Ele const& ele )
        :
        M_element( ele )
    {}

    entity_type
    operator()() const
    {
        entity_type bface;
        // go to global
        size_type i1 = ( M_element.point( 0 ) ).id();
        size_type i2 = ( M_element.point( 1 ) ).id();
        bface = ( makeBareItem( i1, i2 ) ).first;
        return bface;
    }
    Ele const& M_element;
};

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
boundingBox ( const MeshType& mesh )
{
    // processor bounding box with no arguments
    // computes the global bounding box
    return processorBoundingBox( mesh );
}


/**
 * Same, but returns a sphere instead of a box.
 */
template<typename MeshType>
inline
Sphere
boundingSphere ( const MeshType& mesh )
{
    MeshBoundingBox bbox = boundingBox( mesh );

    const double diag = Feel::distance( bbox.second, bbox.first );
    const Point cent = Feel::middle( bbox.second, bbox.first );

    return Sphere ( cent, .5*diag );

}

/**
 * @returns two points defining a cartesian box that bounds the
 * elements belonging to processor pid.  If no processor id is specified
 * the bounding box for the whole mesh is returned.
 */
template<typename MeshType>
inline
MeshBoundingBox
processorBoundingBox ( const MeshType& mesh,
                       const size_type pid = invalid_size_type_value )
{
    FEELPP_ASSERT ( mesh.numPoints() != 0 ).error( "mesh has no points" );

    Point min( 1.e30,   1.e30,  1.e30 );
    Point max( -1.e30, -1.e30, -1.e30 );

    // By default no processor is specified and we compute
    // the bounding box for the whole domain.
    if ( pid == invalid_size_type_value )
    {
        DVLOG(2) << "[processorBoundingBox] np pid given\n";

        for ( unsigned int n=0; n<mesh.numPoints(); n++ )
            for ( unsigned int i=0; i<mesh.dimension(); i++ )
            {
                min( i ) = std::min( min( i ), mesh.point( n )( i ) );
                max( i ) = std::max( max( i ), mesh.point( n )( i ) );
            }
    }

    // if a specific processor id is specified then we need
    // to only consider those elements living on that processor
    else
    {
        DVLOG(2) << "[processorBoundingBox] process bounding box on pid " << pid << "\n";
        typename MeshType::element_iterator it = mesh.beginElementWithProcessId( pid );
        typename MeshType::element_iterator en = mesh.endElementWithProcessId( pid );

        for ( ; it != en; ++it )
            for ( unsigned int n=0; n< MeshType::element_type::numPoints; n++ )
                for ( unsigned int i=0; i<mesh.dimension(); i++ )
                {
                    min( i ) = std::min( min( i ), mesh.point( n )( i ) );
                    max( i ) = std::max( max( i ), mesh.point( n )( i ) );
                }
    }

    for ( unsigned int i=mesh.dimension(); i< min.node().size(); i++ )
    {
        min( i ) = 0;
        max( i ) = 0;
    }

    DVLOG(2) << "[processorBoundingBox] min= " << min << "\n";
    DVLOG(2) << "[processorBoundingBox] max= " << max << "\n";
    const MeshBoundingBox ret_val( min, max );

    return ret_val;
}

/**
 * Same, but returns a sphere instead of a box.
 */
template<typename MeshType>
inline
Sphere
processorBoundingSphere ( const MeshType& mesh,
                          const size_type pid = invalid_size_type_value )
{
    MeshBoundingBox bbox = processorBoundingBox( mesh,pid );

    const Real  diag = Feel::distance( bbox.second, bbox.first );
    const Point cent = Feel::middle( bbox.second, bbox.first );

    DVLOG(2) << "[processorBoundingSphere] processor " << mesh.comm().rank() << "\n";
    DVLOG(2) << "[processorBoundingSphere] center " << cent << "\n";
    DVLOG(2) << "[processorBoundingSphere] radius " << 0.5*diag << "\n";
    return Sphere ( cent, .5*diag );
}


} // Feel

/// \endcond
#endif
