/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-01

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file functors.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-06-01
 */
#ifndef __Functor_H
#define __Functor_H 1

#include <feel/feelalg/glas.hpp>
#include <feel/feelmesh/meshbase.hpp>

namespace Feel
{
class MeshBase;
/// \cond detail
namespace detail
{
struct UpdateMesh
{
    UpdateMesh( MeshBase const* m )
        :
        M_mesh( m )
    {}
    template<typename ElementType>
    void operator()( ElementType& element  ) const
    {
        element.setMesh( M_mesh );
    }
private:
    MeshBase const* M_mesh;

};

struct OnBoundary
{
    OnBoundary( bool is_on_bdy = false )
        :
        M_bdy( is_on_bdy )
    {}
    template<typename ElementType>
    void operator()( ElementType& element  ) const
    {
        element.setOnBoundary( M_bdy );
    }
private:
    bool M_bdy;

};

struct ApplyDisplacement
{
    ApplyDisplacement( uint16_type l, ublas::vector<double> const& displ )
        :
        M_l( l ),
        M_displ( displ )
    {}
    template<typename ElementType>
    void operator()( ElementType& elt )
    {
        elt.applyDisplacement( M_l, M_displ );
    }
    uint16_type M_l;
    ublas::vector<double> const& M_displ;

};

template<typename FaceType>
struct UpdateFace
{
    UpdateFace( FaceType const& face )
        :
        M_face( face )
    {}
    template<typename ElementType>
    void operator()( ElementType& element  ) const
    {
        if ( M_face.ad_first() == element.id() )
        {
            element.setFace( M_face.pos_first(), M_face );

            if ( M_face.isConnectedTo1() )
                element.setNeighbor( M_face.pos_first(), M_face.ad_second(), M_face.proc_second() );

            FEELPP_ASSERT( element.facePtr( M_face.pos_first() ) )
            ( M_face.pos_first() )
            ( M_face.ad_first() ).error( "invalid face" );
        }

        else if ( M_face.ad_second() == element.id() )
        {
            element.setFace( M_face.pos_second(), M_face );
            element.setNeighbor( M_face.pos_second(), M_face.ad_first(), M_face.proc_first() );
            FEELPP_ASSERT( element.facePtr( M_face.pos_second() ) )
            ( M_face.pos_second() )
            ( M_face.ad_second() ).error( "invalid face" );
        }

        else
        {
            FEELPP_ASSERT( 0 )
            ( M_face.ad_first() )( M_face.pos_first() )
            ( M_face.ad_second() )( M_face.pos_second() )
            ( element.id() ).error( "invalid face " );
        }
    }
private:
    FaceType const& M_face;

};

template<typename PermutationType>
struct UpdateFacePermutation
{
    UpdateFacePermutation( uint16_type j, PermutationType const& perm )
        :
        M_localfaceid( j ),
        M_perm( perm )
    {}
    template<typename ElementType>
    void operator()( ElementType& element  ) const
    {
        element.setFacePermutation( M_localfaceid, M_perm );
    }
private:
    uint16_type M_localfaceid;
    PermutationType const& M_perm;

};

template<typename ConnectionType>
struct UpdateFaceConnection0
{
    UpdateFaceConnection0( ConnectionType const& conn )
        :
        M_conn( conn )
    {}
    template<typename FaceType>
    void operator()( FaceType& element  ) const
    {
        element.setConnection0( M_conn );
    }
private:
    ConnectionType const& M_conn;

};
template<typename ConnectionType>
struct UpdateFaceConnection1
{
    UpdateFaceConnection1( ConnectionType const& conn )
        :
        M_conn( conn )
    {}
    template<typename FaceType>
    void operator()( FaceType& element  ) const
    {
        //element.setOnBoundary( false );
        element.setConnection1( M_conn );

        if ( ( element.element0().marker() == element.element1().marker() ) &&
                element.marker().value() == 0 )
            element.setMarker( element.element0().marker().value() );
    }
private:
    ConnectionType const& M_conn;

};

struct UpdateFaceOnBoundary
{
    UpdateFaceOnBoundary( bool onbdy )
        :
        M_onbdy( onbdy )
    {}
    template<typename FaceType>
    void operator()( FaceType& element  ) const
    {
        element.setOnBoundary( M_onbdy );
    }
private:
    bool M_onbdy;

};

template<typename EdgeType>
struct UpdateEdge
{
    UpdateEdge( uint16_type j, EdgeType const& edge )
        :
        M_localedgeid( j ),
        M_edge( edge )
    {}
    template<typename ElementType>
    void operator()( ElementType& element  ) const
    {
        element.setEdge( M_localedgeid, M_edge );
    }
private:
    uint16_type M_localedgeid;
    EdgeType const& M_edge;

};

template<typename PermutationType>
struct UpdateEdgePermutation
{
    UpdateEdgePermutation( uint16_type j, PermutationType const& perm )
        :
        M_localedgeid( j ),
        M_perm( perm )
    {}
    template<typename ElementType>
    void operator()( ElementType& element  ) const
    {
        element.setEdgePermutation( M_localedgeid, M_perm );
    }
private:
    uint16_type M_localedgeid;
    PermutationType const& M_perm;

};

struct updateIdInOthersPartitions
{
    updateIdInOthersPartitions( rank_type pid, size_type id )
        :
        M_pid( pid ),
        M_id( id )
    {}
    template<typename ElementType>
    void operator()( ElementType& element )
    {
        element.setIdInOthersPartitions( M_pid, M_id );
    }
private:
    rank_type M_pid;
    size_type M_id;
};

struct UpdateProcessId
{
    UpdateProcessId( rank_type pid )
        :
        M_pid( pid )
    {}
    template<typename ElementType>
    void operator()( ElementType& element )
    {
        element.setProcessId( M_pid );
    }
private:
    rank_type M_pid;
};

struct UpdateNeighborPartition
{
    UpdateNeighborPartition( rank_type pid )
        :
        M_pid( pid )
    {}
    template<typename ElementType>
    void operator()( ElementType& element )
    {
        element.addNeighborPartitionId( M_pid );
    }
private:
    rank_type M_pid;
};



struct UpdateMarker
{
    UpdateMarker( flag_type v )
        :
        M_v( v )
    {}

    template<typename ElementType>
    void operator()( ElementType& element )
    {
        element.setMarker( M_v );
    }
private:
    flag_type M_v;
};


template<typename PointType>
struct UpdateEltPoint
{
    UpdateEltPoint( uint16_type k, PointType const& pt )
        :
        M_idx( k ),
        M_pt( pt )
    {}
    template<typename EltType>
    void operator()( EltType& element  ) const
    {
        element.setPoint( M_idx, M_pt );
    }
private:
    uint16_type M_idx;
    PointType const& M_pt;
};


} // detail
/// \endcond detail
}
#endif /* __Functor_H */
