//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Author(s) : 
//!     Thibaut Metivet <thibaut.metivet@inria.fr>
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file distancetorange.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 01 Sept 2021
//! @copyright 2021 Feel++ Consortium
//! @copyright 2021 INRIA
//!

#ifndef _DISTANCE_TO_RANGE_HPP
#define _DISTANCE_TO_RANGE_HPP 1

#include <algorithm>

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feells/fastmarching_impl.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/syncdofs.hpp>
#include <feel/feelcore/traits.hpp>

#include "geometryconceptwrappers.hpp"
#include "distancepointtoface.hpp"

namespace Feel {

template< typename FunctionSpaceType >
class DistanceToRange
{
    public:
        typedef DistanceToRange< FunctionSpaceType > self_type;
        typedef std::shared_ptr< self_type > self_ptrtype;

        //--------------------------------------------------------------------//
        // Distance functionspace and mesh
        typedef FunctionSpaceType functionspace_distance_type;
        typedef std::shared_ptr< functionspace_distance_type > functionspace_distance_ptrtype;
        typedef typename functionspace_distance_type::element_type element_distance_type;
        typedef typename functionspace_distance_type::element_ptrtype element_distance_ptrtype;

        static const uint16_type nDofPerEltDistance = functionspace_distance_type::fe_type::nDof;

        typedef typename functionspace_distance_type::mesh_type mesh_distance_type;
        typedef typename functionspace_distance_type::mesh_ptrtype mesh_distance_ptrtype;

        typedef typename MeshTraits<mesh_distance_type>::elements_reference_wrapper_type elements_reference_wrapper_distance_type;
        typedef typename MeshTraits<mesh_distance_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_distance_ptrtype;
        typedef elements_reference_wrapper_t<mesh_distance_type> range_elements_distance_type;

        typedef faces_reference_wrapper_t<mesh_distance_type> range_faces_type;

        //--------------------------------------------------------------------//
        static constexpr uint16_type nRealDim = functionspace_distance_type::nRealDim;
        using size_type = typename functionspace_distance_type::size_type;
        typedef typename functionspace_distance_type::value_type value_type;
        typedef typename node<value_type>::type node_type;
        typedef typename matrix_node<value_type>::type matrix_node_type;

        //--------------------------------------------------------------------//
        // Fast-marching
        typedef FastMarching<functionspace_distance_type> fastmarching_type;
        typedef std::shared_ptr< fastmarching_type > fastmarching_ptrtype;

        //--------------------------------------------------------------------//
        typedef std::unordered_map< size_type, std::vector< size_type > > map_face_dofs_type;

    public:
        //--------------------------------------------------------------------//
        // Constructor
        DistanceToRange( functionspace_distance_ptrtype const& spaceDistance );

        //--------------------------------------------------------------------//
        // Accessors
        functionspace_distance_ptrtype const& functionSpaceDistance() const { return M_spaceDistance; }
        mesh_distance_ptrtype const& meshDistance() const { return this->functionSpaceDistance()->mesh(); }

        //--------------------------------------------------------------------//
        // Fast-marching
        fastmarching_ptrtype const& fastMarching() const;

        //--------------------------------------------------------------------//
        // Geometry
        value_type distanceDofToFace( size_type dofId, size_type faceId ) const;

        //--------------------------------------------------------------------//
        // Result
        template< typename RangeType >
        element_distance_type unsignedDistance( RangeType const& rangeFaces ) const;
        template< typename RangeType >
        element_distance_type signedDistance( RangeType const& rangeFaces ) const;

    private:
        element_distance_type unsignedDistanceToFaces( range_faces_type const& rangeFaces ) const;
        element_distance_type signedDistanceToFaces( range_faces_type const& rangeFaces ) const;

        elements_reference_wrapper_distance_ptrtype eltsTouchingFaces( range_faces_type const& rangeFaces ) const;
        map_face_dofs_type dofsNeighbouringFaces( range_faces_type const& rangeFaces ) const;

    private:
        functionspace_distance_ptrtype M_spaceDistance;

        fastmarching_ptrtype M_fastMarching;

        mutable elements_reference_wrapper_distance_ptrtype M_eltsTouchingFaces;
        mutable std::unordered_map< size_type, std::vector< size_type > > M_dofsNeighbouringFaces;

        mutable std::unordered_map< size_type, Eigen::Matrix<value_type,3,3> > M_faceTriangleTransformationMatrices;
};

template< typename FunctionSpaceType >
DistanceToRange< FunctionSpaceType >::DistanceToRange(
        functionspace_distance_ptrtype const& spaceDistance ) :
    M_spaceDistance( spaceDistance ),
    M_eltsTouchingFaces( new elements_reference_wrapper_distance_type{} )
{}

template< typename FunctionSpaceType >
typename DistanceToRange< FunctionSpaceType >::fastmarching_ptrtype const&
DistanceToRange< FunctionSpaceType >::fastMarching() const
{
    if( !M_fastMarching )
        const_cast<self_type*>(this)->M_fastMarching.reset( new fastmarching_type( this->functionSpaceDistance() ) );

    return M_fastMarching;
}

template< typename FunctionSpaceType >
typename DistanceToRange< FunctionSpaceType >::value_type 
DistanceToRange< FunctionSpaceType >::distanceDofToFace( size_type dofId, size_type faceId ) const
{
    auto const& face = this->meshDistance()->face( faceId );
    auto const& pt = boost::get<0>( this->functionSpaceDistance()->dof()->dofPoint( dofId ) );
    auto const& faceVertices = face.vertices();

    auto P = eigenMap<nRealDim>( pt );
    auto F = eigenMap<nRealDim>( faceVertices );

    if constexpr ( nRealDim == 2 )
    {
        return Feel::detail::geometry::distancePointToSegment( P, F.col(0), F.col(1) );
    }
    else if constexpr ( nRealDim == 3 )
    {
        auto const& P1 = F.col(0);
        auto const& P2 = F.col(1);
        auto const& P3 = F.col(2);
        auto surfaceTriangleTransformationMatrixIt = M_faceTriangleTransformationMatrices.find( faceId );
        if( surfaceTriangleTransformationMatrixIt == M_faceTriangleTransformationMatrices.end() ) 
        {
            surfaceTriangleTransformationMatrixIt = M_faceTriangleTransformationMatrices.insert( { faceId, Feel::detail::geometry::triangleFrameTransformationMatrix( P1, P2, P3 ) } ).first;
        }
        return Feel::detail::geometry::distancePointToTriangleWithTransformationMatrix( P, P1, P2, P3, surfaceTriangleTransformationMatrixIt->second );
    }
}

template< typename FunctionSpaceType >
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::unsignedDistance( RangeType const& rangeFaces ) const
{
    return this->unsignedDistanceToFaces( rangeFaces );
}

template< typename FunctionSpaceType >
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::signedDistance( RangeType const& rangeFaces ) const
{
    return this->signedDistanceToFaces( rangeFaces );
}

template< typename FunctionSpaceType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::unsignedDistanceToFaces( range_faces_type const& rangeFaces ) const
{
    auto unsignedDistance = this->functionSpaceDistance()->element( "unsignedDistance" );
    // Initialise distance with arbitrarily large value
    unsignedDistance.setConstant( std::numeric_limits<value_type>::max() );
    // Set distance=0 on rangeFaces
    unsignedDistance.on( _range=rangeFaces, _expr=cst(0.) );
    // Compute exact distance on neighbouring dofs
    auto const dofsNeighbouringFaces = this->dofsNeighbouringFaces( rangeFaces );
    for( auto const [faceId, dofsIds]: dofsNeighbouringFaces )
    {
        for( size_type dofId: dofsIds )
        {
            auto dist = this->distanceDofToFace( dofId, faceId );
            if( dist < unsignedDistance(dofId) )
            {
                unsignedDistance(dofId) = dist;
            }
        }
    }

    auto eltsTouchingFaces = this->eltsTouchingFaces( rangeFaces );
    range_elements_distance_type rangeEltsTouchingFaces( mpl::size_t<MESH_ELEMENTS>(), eltsTouchingFaces->begin(), eltsTouchingFaces->end(), eltsTouchingFaces );

    // Sync initial distance
    syncDofs( unsignedDistance, rangeEltsTouchingFaces,
            []( value_type valCurrent, std::set<value_type> ghostVals )
            {
                auto const minGhostValueIt = std::min_element( ghostVals.begin(), ghostVals.end() );
                return std::min( valCurrent, *minGhostValueIt );
            }
        );

    // Then perform fast-marching
    unsignedDistance = this->fastMarching()->run( unsignedDistance, rangeEltsTouchingFaces );
    // Return
    return unsignedDistance;
}

template< typename FunctionSpaceType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::signedDistanceToFaces( range_faces_type const& rangeFaces ) const
{
    // Compute unsigned distance
    // Set sign: set negative everywhere, then propagate + sign from the domain boundary
    // note: the element container (VectorUblas) does not provide unary operators at the moment -> TODO: need to improve
    // note2: elements touching the faces are cached by unsignedDistanceToFaces in M_eltsTouchingFaces
    auto signedDistance = this->unsignedDistanceToFaces( rangeFaces );
    signedDistance->scale( -1. );

    // Communication
    rank_type const pidMeshDistance = this->meshDistance()->worldCommPtr()->localRank();
    std::set<rank_type> const& neighborSubdomainsMeshDistance = this->meshDistance()->neighborSubdomains();
    int nRequests = 2*neighborSubdomainsMeshDistance.size();
    mpi::request * mpiRequests = new mpi::request[nRequests];

    std::unordered_set< size_type > eltsToVisit, eltsVisited;
    auto const& intersectingElements = this->intersectingElements();
    auto const rangeMeshBoundaryElements = boundaryelements( this->meshDistance() );
    auto it_boundaryelt = rangeMeshBoundaryElements.template get<1>();
    auto en_boundaryelt = rangeMeshBoundaryElements.template get<2>();
    for( ; it_boundaryelt != en_boundaryelt; it_boundaryelt++ )
    {
        auto const& elt = boost::unwrap_ref( *it_boundaryelt );
        size_type const eltId = elt.id();
        eltsToVisit.insert( eltId );
    }

    bool eltsToVisitIsEmptyOnAllProc = false;
    while( !eltsToVisitIsEmptyOnAllProc )
    {
        /* Contain the ghost elts which require communication */
        std::map<rank_type, std::vector< size_type > > dataToRecv;
        std::map<rank_type, std::vector< size_type > > dataToSend;

        while( !eltsToVisit.empty() )
        {
            auto eltIt = eltsToVisit.begin();
            size_type eltId = *eltIt;
            auto const& elt = this->meshDistance()->element( eltId );
            // visit elt
            // set + sign
            for( uint16_type j = 0; j < nDofPerEltDistance; j++ )
            {
                auto curPhi = signedDistance.localToGlobal( eltId, j, 0 );
                if( curPhi < 0. )
                    signedDistance.assign( eltId, j, 0, -curPhi );
            }
            // add potential neighbors to visit
            for( uint16_type faceId = 0; faceId < elt.nTopologicalFaces(); faceId++ )
            {
                size_type neighId = elt.neighbor( faceId );
                if( neighId == invalid_size_type_value )
                    continue;
                auto const& neigh = this->meshDistance()->element( neighId );
                rank_type const neighPid = neigh.processId();
                // the neighbor elt is a ghost -> send to owner
                if( neighPid != pidMeshDistance )
                {
                    size_type neighIdOnOwner = neigh.idInOthersPartitions( neighPid );
                    dataToSend[neighPid].push_back( neighIdOnOwner );
                    continue;
                }
                if( eltsVisited.find( neighId ) != eltsVisited.end()
                        // stop when reaching an intersecting elt
                        || intersectingElements.find( neighId ) != intersectingElements.end() )
                    continue;
                // need to visit neighbor
                eltsToVisit.insert( neighId );
            }
            eltsVisited.insert( eltId );
            eltsToVisit.erase( eltIt );
        }

        // Perform communications
        int cntRequests = 0;
        for( rank_type neighborRank: neighborSubdomainsMeshDistance )
        {
            mpiRequests[cntRequests++] = this->meshDistance()->worldCommPtr()->localComm().isend( neighborRank, 0, dataToSend[neighborRank] );
            mpiRequests[cntRequests++] = this->meshDistance()->worldCommPtr()->localComm().irecv( neighborRank, 0, dataToRecv[neighborRank] );
        }
        mpi::wait_all( mpiRequests, mpiRequests + cntRequests );

        /* Process received ghosts */
        for( auto const& data: dataToRecv )
        {
            for( size_type const eltId: data.second )
            {
                if( eltsVisited.find( eltId ) != eltsVisited.end()
                        // stop when reaching an intersecting elt
                        || M_eltsTouchingFaces->find( eltId ) != intersectingElements.end() )
                    continue;
                // need to visit neighbor
                eltsToVisit.insert( eltId );
            }
        }

        bool eltsToVisitIsEmpty = eltsToVisit.empty();
        eltsToVisitIsEmptyOnAllProc = mpi::all_reduce(
                this->meshDistance()->worldComm(),
                eltsToVisitIsEmpty,
                std::logical_and<bool>() );
    }

    delete [] mpiRequests;

    sync( signedDistance, "max" );

    // Return
    return signedDistance;
}

template< typename FunctionSpaceType >
typename DistanceToRange< FunctionSpaceType >::elements_reference_wrapper_distance_ptrtype 
DistanceToRange< FunctionSpaceType >::eltsTouchingFaces( range_faces_type const& rangeFaces ) const
{
    // Get ids of elements touching a face in rangeFaces
    // and cache in M_eltsTouchingFaces for reuse
    M_eltsTouchingFaces->clear();

    auto face_it = rangeFaces.template get<1>();
    auto face_en = rangeFaces.template get<2>();
    for( ; face_it != face_en ; ++face_it )
    {
        auto const& face = boost::unwrap_ref( *face_it );
        size_type const faceId = face.id();
        if ( face.isConnectedTo0() )
        {
            M_eltsTouchingFaces->push_back( boost::cref( face.element0() ) );
        }
        if ( face.isConnectedTo1() )
        {
            M_eltsTouchingFaces->push_back( boost::cref( face.element1() ) );
        }
    }
    M_eltsTouchingFaces->shrink_to_fit();

    return M_eltsTouchingFaces;
}

template< typename FunctionSpaceType >
typename DistanceToRange< FunctionSpaceType >::map_face_dofs_type
DistanceToRange< FunctionSpaceType >::dofsNeighbouringFaces( range_faces_type const& rangeFaces ) const
{
    // Get ids of dofs touching a face in rangeFaces via an element, without the face dofs
    // and cache in M_dofsNeighbouringFaces
    M_dofsNeighbouringFaces.clear();

    auto const elementDofsNeighbouringFaces = [&]( size_type const faceId, auto const faceDofsBe, auto const faceDofsEn, size_type const eltId )
    {
        auto const [eltDofsBegin, eltDofsEnd] = this->functionSpaceDistance()->dof()->localDof( eltId );
        for( auto const& [lDof, gDof]: this->functionSpaceDistance()->dof()->localDof( eltId ) )
        {
            size_type const dofId = gDof.index();
            auto const findIt = std::find_if( faceDofsBe, faceDofsEn, 
                    [dofId]( auto const& faceDof ) { return faceDof.index() == dofId; }
                    );
            if( findIt == faceDofsEn )
                M_dofsNeighbouringFaces[faceId].push_back( dofId );
        }
    };

    auto face_it = rangeFaces.template get<1>();
    auto face_en = rangeFaces.template get<2>();
    for( ; face_it != face_en ; ++face_it )
    {
        auto const& face = boost::unwrap_ref( *face_it );
        size_type const faceId = face.id();
        auto const [faceDofsBegin, faceDofsEnd] = this->functionSpaceDistance()->dof()->faceLocalDof( faceId );
        if ( face.isConnectedTo0() )
        {
            elementDofsNeighbouringFaces( faceId, faceDofsBegin, faceDofsEnd, face.element0().id() );
        }
        if ( face.isConnectedTo1() )
        {
            elementDofsNeighbouringFaces( faceId, faceDofsBegin, faceDofsEnd, face.element1().id() );
        }
    }

    return M_dofsNeighbouringFaces;
}

template< typename SpaceType, typename RangeType >
typename Feel::decay_type<SpaceType>::element_type
distanceToRange( SpaceType const& space, RangeType const& range )
{
    DistanceToRange distToRange( space );
    return distToRange.unsignedDistance( range );
}

} // namespace Feel

#endif //_DISTANCE_TO_RANGE_HPP
