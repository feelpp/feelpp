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
#include <feel/feelmesh/ranges.hpp>
#include <feel/feeldiscr/syncdofs.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelvf/cst.hpp>

#include <feel/feells/geometryconceptwrappers.hpp>
#include <feel/feells/distancepointtoface.hpp>

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

        static inline const uint16_type nDofPerEltDistance = functionspace_distance_type::fe_type::nDof;

        typedef typename functionspace_distance_type::mesh_type mesh_distance_type;
        typedef typename functionspace_distance_type::mesh_ptrtype mesh_distance_ptrtype;

        using elements_reference_wrapper_distance_type = Range<mesh_distance_type,MESH_ELEMENTS>;
        using elements_reference_wrapper_distance_ptrtype = std::shared_ptr<Range<mesh_distance_type,MESH_ELEMENTS>>;
        using range_elements_distance_type = Range<mesh_distance_type,MESH_ELEMENTS>;
        using range_faces_type = Range<mesh_distance_type,MESH_FACES> ;

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
        // Options
        /*
         * Maximal computed distance (narrow band)
         * Negative (eg -1) computes full distance
         */
        value_type maxDistance() const { return M_maxDistance; }
        void setMaxDistance( const value_type & maxDist ) { M_maxDistance = maxDist; }
        /*
         * Fast-marching stride: local fast-marching stride before parallel update
         * Negative (eg -1) runs full local marching between each sync
         */
        value_type fastMarchingStride() const { return M_fastMarchingStride; }
        void setFastMarchingStride( const value_type & stride ) { M_fastMarchingStride = stride; }

        //--------------------------------------------------------------------//
        // Geometry
        value_type distanceDofToFace( size_type dofId, size_type faceId ) const;

        //--------------------------------------------------------------------//
        // Result
        template< typename RangeType >
        element_distance_type unsignedDistance( RangeType && rangeFaces ) const;
        template< typename RangeType >
        element_distance_type signedDistance( RangeType && rangeFaces ) const;

    private:
        template< typename RangeType >
        element_distance_type unsignedDistanceToFaces( RangeType && rangeFaces ) const;

        template< typename RangeType >
        element_distance_type signedDistanceToFaces( RangeType && rangeFaces ) const;

        template< typename RangeType >
        range_elements_distance_type const& eltsTouchingFaces( RangeType && rangeFaces ) const;

        template< typename RangeType >
        map_face_dofs_type dofsNeighbouringFaces( RangeType && rangeFaces ) const;

    private:
        functionspace_distance_ptrtype M_spaceDistance;

        fastmarching_ptrtype M_fastMarching;

        value_type M_maxDistance = -1.;
        value_type M_fastMarchingStride = -1.;

        mutable range_elements_distance_type M_eltsTouchingFaces;
        mutable std::unordered_map< size_type, std::vector< size_type > > M_dofsNeighbouringFaces;

        mutable std::unordered_map< size_type, Eigen::Matrix<value_type,3,3> > M_faceTriangleTransformationMatrices;
};

template< typename FunctionSpaceType >
DistanceToRange< FunctionSpaceType >::DistanceToRange(
        functionspace_distance_ptrtype const& spaceDistance ) :
    M_spaceDistance( spaceDistance ),
    M_eltsTouchingFaces( spaceDistance->mesh() )
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
DistanceToRange< FunctionSpaceType >::unsignedDistance( RangeType && rangeFaces ) const
{
    return this->unsignedDistanceToFaces( std::forward<RangeType>(rangeFaces) );
}

template< typename FunctionSpaceType >
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::signedDistance( RangeType && rangeFaces ) const
{
    return this->signedDistanceToFaces( std::forward<RangeType>(rangeFaces) );
}

template< typename FunctionSpaceType >
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::unsignedDistanceToFaces( RangeType && rangeFaces ) const
{
    static_assert(decay_type<RangeType>::entities() == MESH_FACES, "RangeType must have entities() == MESH_FACES");

    auto unsignedDistance = this->functionSpaceDistance()->element( "unsignedDistance" );
    // Initialise distance with arbitrarily large value
    unsignedDistance.setConstant( std::numeric_limits<value_type>::max() );
    // Set distance=0 on rangeFaces
    unsignedDistance.on( _range=std::forward<RangeType>(rangeFaces), _expr=cst(0.) );
    // Compute exact distance on neighbouring dofs
    auto const dofsNeighbouringFaces = this->dofsNeighbouringFaces( std::forward<RangeType>(rangeFaces) );
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

    auto const& rangeEltsTouchingFaces = this->eltsTouchingFaces( std::forward<RangeType>(rangeFaces) );

    // Sync initial distance
    syncDofs( unsignedDistance, rangeEltsTouchingFaces,
            []( value_type valCurrent, std::set<value_type> ghostVals )
            {
                auto const minGhostValueIt = std::min_element( ghostVals.begin(), ghostVals.end() );
                return std::min( valCurrent, *minGhostValueIt );
            }
        );

    // Then perform fast-marching
    this->fastMarching()->setNarrowBandWidth( this->maxDistance() );
    this->fastMarching()->setStride( this->fastMarchingStride() );
    unsignedDistance = this->fastMarching()->run( unsignedDistance, rangeEltsTouchingFaces );
    // Return
    return unsignedDistance;
}

template< typename FunctionSpaceType >
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::element_distance_type
DistanceToRange< FunctionSpaceType >::signedDistanceToFaces( RangeType && rangeFaces ) const
{
    static_assert(decay_type<RangeType>::entities() == MESH_FACES, "RangeType must have entities() == MESH_FACES");

    // Compute unsigned distance
    // Set sign: set negative everywhere, then propagate + sign from the domain boundary
    // note: the element container (VectorUblas) does not provide unary operators at the moment -> TODO: need to improve
    // note2: elements touching the faces are cached by unsignedDistanceToFaces in M_eltsTouchingFaces
    auto signedDistance = this->unsignedDistanceToFaces( std::forward<RangeType>(rangeFaces) );
    signedDistance->scale( -1. );

    // Communication
    rank_type const pidMeshDistance = this->meshDistance()->worldCommPtr()->localRank();
    std::set<rank_type> const& neighborSubdomainsMeshDistance = this->meshDistance()->neighborSubdomains();
    int nRequests = 2*neighborSubdomainsMeshDistance.size();
    mpi::request * mpiRequests = new mpi::request[nRequests];

    std::unordered_set< size_type > eltsToVisit, eltsVisited;
    auto const& intersectingElements = this->intersectingElements();
    auto const rangeMeshBoundaryElements = boundaryelements( this->meshDistance() );
    std::for_each( rangeMeshBoundaryElements.begin(), rangeMeshBoundaryElements.end(),
            [&eltsToVisit]( auto const& elt ) { eltsToVisit.insert( elt.id() ); } );

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
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::range_elements_distance_type const&
DistanceToRange< FunctionSpaceType >::eltsTouchingFaces( RangeType && rangeFaces ) const
{
    // Get ids of elements touching a face in rangeFaces
    // and cache in M_eltsTouchingFaces for reuse
    M_eltsTouchingFaces.clear();

    for( auto const& facew: rangeFaces )
    {
        auto const& face = boost::unwrap_ref( facew );
        if ( face.isConnectedTo0() && !face.element0().isGhostCell() )
        {
            M_eltsTouchingFaces.push_back( face.element0() );
        }
        if ( face.isConnectedTo1() && !face.element1().isGhostCell() )
        {
            M_eltsTouchingFaces.push_back( face.element1() );
        }
    }
    M_eltsTouchingFaces.shrink_to_fit();

    return M_eltsTouchingFaces;
}

template< typename FunctionSpaceType >
template< typename RangeType >
typename DistanceToRange< FunctionSpaceType >::map_face_dofs_type
DistanceToRange< FunctionSpaceType >::dofsNeighbouringFaces( RangeType && rangeFaces ) const
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

    auto face_it = std::forward<RangeType>(rangeFaces).begin();
    auto face_en = std::forward<RangeType>(rangeFaces).end();
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

namespace na::distancetorange {
    using max_distance = NA::named_argument_t<struct max_distance_tag>;
    using fm_stride = NA::named_argument_t<struct fm_stride_tag>;
}
inline constexpr auto& _max_distance = NA::identifier<na::distancetorange::max_distance>;
inline constexpr auto& _fm_stride = NA::identifier<na::distancetorange::fm_stride>;

template< typename ... Args >
auto distanceToRange( Args && ... nargs )
{
    auto args = NA::make_arguments( std::forward<Args>(nargs)... );
    auto && space = args.get( _space );
    auto && range = args.get( _range );
    double maxDistance = args.get_else( _max_distance, -1. );
    double fastMarchingStride = args.get_else( _fm_stride, -1. );
    DistanceToRange distToRange( std::forward<decltype(space)>(space) );
    distToRange.setMaxDistance( maxDistance );
    distToRange.setFastMarchingStride( fastMarchingStride );
    return distToRange.unsignedDistance( std::forward<decltype(range)>(range) );
}


} // namespace Feel

#endif //_DISTANCE_TO_RANGE_HPP
