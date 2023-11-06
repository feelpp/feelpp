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
//! @file distancetomesh.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _DISTANCE_TO_MESH_HPP
#define _DISTANCE_TO_MESH_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feells/reinit_fms_impl.hpp>

#include "geometryconceptwrappers.hpp"
#include "trianglesintersect.hpp"
#include "distancepointtoface.hpp"

namespace Feel {

template< typename MeshType, typename FunctionSpaceType >
class DistanceToMesh
{
    public:
        typedef DistanceToMesh< MeshType, FunctionSpaceType > self_type;
        typedef std::shared_ptr< self_type > self_ptrtype;

        //--------------------------------------------------------------------//
        // Surface mesh
        typedef MeshType mesh_surface_type;
        typedef std::shared_ptr< mesh_surface_type > mesh_surface_ptrtype;

        //--------------------------------------------------------------------//
        // Distance functionspace and mesh
        typedef FunctionSpaceType functionspace_distance_type;
        typedef std::shared_ptr< functionspace_distance_type > functionspace_distance_ptrtype;
        typedef typename functionspace_distance_type::element_type element_distance_type;
        typedef typename functionspace_distance_type::element_ptrtype element_distance_ptrtype;

        static inline const uint16_type nDofPerEltDistance = functionspace_distance_type::fe_type::nDof;

        typedef typename functionspace_distance_type::mesh_type mesh_distance_type;
        typedef typename functionspace_distance_type::mesh_ptrtype mesh_distance_ptrtype;

        typedef typename MeshTraits<mesh_distance_type>::elements_reference_wrapper_type elements_reference_wrapper_distance_type;
        typedef typename MeshTraits<mesh_distance_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_distance_ptrtype;
        typedef elements_reference_wrapper_t<mesh_distance_type> range_elements_distance_type;

        //--------------------------------------------------------------------//
        static constexpr uint16_type nRealDim = functionspace_distance_type::nRealDim;
        using size_type = typename functionspace_type::size_type;
        typedef typename functionspace_distance_type::value_type value_type;
        typedef typename node<value_type>::type node_type;
        typedef typename matrix_node<value_type>::type matrix_node_type;

        //--------------------------------------------------------------------//
        // Fast-marching
        typedef ReinitializerFMS< functionspace_distance_type > fastmarching_type;
        typedef std::shared_ptr<fastmarching_type> fastmarching_ptrtype;

    public:
        //--------------------------------------------------------------------//
        // Constructor
        DistanceToMesh( mesh_surface_ptrtype const& meshSurface, functionspace_distance_ptrtype const& spaceDistance );

        //--------------------------------------------------------------------//
        // Accessors
        mesh_surface_ptrtype const& meshSurface() const { return M_meshSurface; }
        functionspace_distance_ptrtype const& functionSpaceDistance() const { return M_spaceDistance; }
        mesh_distance_ptrtype const& meshDistance() const { return this->functionSpaceDistance()->mesh(); }

        //--------------------------------------------------------------------//
        // Fast-marching
        fastmarching_ptrtype const& fastMarching() const;

        //--------------------------------------------------------------------//
        // Geometry
        static bool segmentsIntersect( matrix_node_type const& seg1, matrix_node_type const& seg2 );
        static bool trianglesIntersect( matrix_node_type const& tri1, matrix_node_type const& tri2 );
        static bool facesIntersect( matrix_node_type const& face1, matrix_node_type const& face2 );

        value_type distanceDofToSurfaceElt( size_type dofId, size_type surfEltId ) const;

        //--------------------------------------------------------------------//
        // Result
        element_distance_ptrtype const& unsignedDistance() const;
        element_distance_ptrtype const& signedDistance() const;

        std::unordered_set< size_type > const& intersectingElements() const;
        range_elements_distance_type rangeIntersectingElements() const;

    private:
        void updateIntersectingElements();
        void updateUnsignedDistance();
        void updateSignedDistance();


    private:
        mesh_surface_ptrtype M_meshSurface;
        functionspace_distance_ptrtype M_spaceDistance;

        fastmarching_ptrtype M_fastMarching;

        mutable std::unordered_map< size_type, Eigen::Matrix<value_type,3,3> > M_meshSurfaceTriangleTransformationMatrices;

        std::unordered_set< size_type > M_intersectingElements;
        std::unordered_map< size_type, std::unordered_set< size_type > > M_eltsIntersectedBySurfaceElt;
        bool M_doUpdateIntersectingElements;

        element_distance_ptrtype M_unsignedDistance;
        bool M_doUpdateUnsignedDistance;
        element_distance_ptrtype M_signedDistance;
        bool M_doUpdateSignedDistance;
};

template< typename MeshType, typename FunctionSpaceType >
DistanceToMesh< MeshType, FunctionSpaceType >::DistanceToMesh(
        mesh_surface_ptrtype const& meshSurface, functionspace_distance_ptrtype const& spaceDistance ) :
    M_meshSurface( meshSurface ),
    M_spaceDistance( spaceDistance ),
    M_doUpdateIntersectingElements( true ),
    M_doUpdateUnsignedDistance( true ),
    M_doUpdateSignedDistance( true )
{}

template< typename MeshType, typename FunctionSpaceType >
typename DistanceToMesh< MeshType, FunctionSpaceType >::fastmarching_ptrtype const&
DistanceToMesh< MeshType, FunctionSpaceType >::fastMarching() const
{
    if( !M_fastMarching )
        const_cast<self_type*>(this)->M_fastMarching.reset( new fastmarching_type( this->functionSpaceDistance() ) );

    return M_fastMarching;
}

template< typename MeshType, typename FunctionSpaceType >
bool
DistanceToMesh< MeshType, FunctionSpaceType >::segmentsIntersect( matrix_node_type const& seg1, matrix_node_type const& seg2 )
{
    return boost::geometry::intersects( Feel::detail::geometry::segmentWrap<nRealDim>( seg1 ), Feel::detail::geometry::segmentWrap<nRealDim>( seg2 ) );
}

template< typename MeshType, typename FunctionSpaceType >
bool
DistanceToMesh< MeshType, FunctionSpaceType >::trianglesIntersect( matrix_node_type const& poly1, matrix_node_type const& poly2 )
{
    typedef Eigen::Matrix<value_type, nRealDim, Eigen::Dynamic, Eigen::ColMajor> PointsMatrix;
    typedef Eigen::Map<PointsMatrix> PointsMatrixMap;
    PointsMatrixMap tri1( const_cast<value_type*>( poly1.data().begin() ), nRealDim, poly1.size2() );
    PointsMatrixMap tri2( const_cast<value_type*>( poly2.data().begin() ), nRealDim, poly2.size2() );
    return Feel::detail::geometry::triangles3DIntersect( 
            tri1.col(0), tri1.col(1), tri1.col(2),
            tri2.col(0), tri2.col(1), tri2.col(2)
            );
}

template< typename MeshType, typename FunctionSpaceType >
bool
DistanceToMesh< MeshType, FunctionSpaceType >::facesIntersect( matrix_node_type const& face1, matrix_node_type const& face2 )
{
    static_assert( nRealDim == 2 || nRealDim == 3, "nRealDim must be 2 or 3" );
    if constexpr ( nRealDim == 2 )
        return segmentsIntersect( face1, face2 );
    else if constexpr ( nRealDim == 3 )
        return trianglesIntersect( face1, face2 );
}

template< typename MeshType, typename FunctionSpaceType >
typename DistanceToMesh< MeshType, FunctionSpaceType >::value_type 
DistanceToMesh< MeshType, FunctionSpaceType >::distanceDofToSurfaceElt( size_type dofId, size_type surfEltId ) const
{
    auto const& surfElt = this->meshSurface()->element( surfEltId );
    auto const& pt = boost::get<0>( this->functionSpaceDistance()->dof()->dofPoint( dofId ) );
    auto const& face = surfElt.vertices();

    auto P = eigenMap<nRealDim>( pt );
    auto F = eigenMap<nRealDim>( face );

    if constexpr ( nRealDim == 2 )
    {
        return Feel::detail::geometry::distancePointToSegment( P, F.col(0), F.col(1) );
    }
    else if constexpr ( nRealDim == 3 )
    {
        auto const& P1 = F.col(0);
        auto const& P2 = F.col(1);
        auto const& P3 = F.col(2);
        auto surfaceTriangleTransformationMatrixIt = M_meshSurfaceTriangleTransformationMatrices.find( surfEltId );
        if( surfaceTriangleTransformationMatrixIt == M_meshSurfaceTriangleTransformationMatrices.end() ) 
        {
            surfaceTriangleTransformationMatrixIt = M_meshSurfaceTriangleTransformationMatrices.insert( { surfEltId, Feel::detail::geometry::triangleFrameTransformationMatrix( P1, P2, P3 ) } ).first;
        }
        return Feel::detail::geometry::distancePointToTriangleWithTransformationMatrix( P, P1, P2, P3, surfaceTriangleTransformationMatrixIt->second );
    }
}

template< typename MeshType, typename FunctionSpaceType >
typename DistanceToMesh< MeshType, FunctionSpaceType >::element_distance_ptrtype const&
DistanceToMesh< MeshType, FunctionSpaceType >::unsignedDistance() const
{
    if( M_doUpdateUnsignedDistance )
        const_cast<self_type*>(this)->updateUnsignedDistance();

    return M_unsignedDistance;
}

template< typename MeshType, typename FunctionSpaceType >
typename DistanceToMesh< MeshType, FunctionSpaceType >::element_distance_ptrtype const&
DistanceToMesh< MeshType, FunctionSpaceType >::signedDistance() const
{
    if( M_doUpdateSignedDistance )
        const_cast<self_type*>(this)->updateSignedDistance();

    return M_signedDistance;
}

template< typename MeshType, typename FunctionSpaceType >
std::unordered_set< size_type > const&
DistanceToMesh< MeshType, FunctionSpaceType >::intersectingElements() const
{
    if( M_doUpdateIntersectingElements )
        const_cast<self_type*>(this)->updateIntersectingElements();
    return M_intersectingElements;
}

template< typename MeshType, typename FunctionSpaceType >
typename DistanceToMesh< MeshType, FunctionSpaceType >::range_elements_distance_type
DistanceToMesh< MeshType, FunctionSpaceType >::rangeIntersectingElements() const
{
    auto const& intersectingElements = this->intersectingElements();
    elements_reference_wrapper_distance_ptrtype intersectingElementsRefWrapper( new elements_reference_wrapper_distance_type() );
    intersectingElementsRefWrapper->reserve( intersectingElements.size() );
    std::transform( 
            intersectingElements.begin(), intersectingElements.end(), 
            std::back_inserter( *intersectingElementsRefWrapper ),
            [this]( size_type id ) { return boost::cref( this->meshDistance()->element( id ) ); }
            );
    return range(_range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
            intersectingElementsRefWrapper->begin(),
            intersectingElementsRefWrapper->end(),
            intersectingElementsRefWrapper
            ), _mesh=this->meshDistance() );
}

template< typename MeshType, typename FunctionSpaceType >
    void
DistanceToMesh< MeshType, FunctionSpaceType >::updateIntersectingElements()
{
    M_intersectingElements.clear();
    M_eltsIntersectedBySurfaceElt.clear();
    std::unordered_map< size_type, std::unordered_set< size_type > > eltsToVisit;

    // Meshes
    auto const& meshSurface = this->meshSurface();
    auto const& meshDistance = this->meshDistance();
    // Localisation tool
    auto locTool = meshDistance->tool_localization();
    locTool->setExtrapolation( false );
    locTool->updateForUse();
    // Communication
    rank_type const pidMeshDistance = meshDistance->worldCommPtr()->localRank();
    std::set<rank_type> const& neighborSubdomainsMeshDistance = meshDistance->neighborSubdomains();
    int nRequests = 2*neighborSubdomainsMeshDistance.size();
    mpi::request * mpiRequests = new mpi::request[nRequests];

    auto it_elt_surf = meshSurface->beginOrderedElement();
    auto en_elt_surf = meshSurface->endOrderedElement();
    for( ; it_elt_surf != en_elt_surf; it_elt_surf++ )
    {
        auto const& surfElt = boost::unwrap_ref( *it_elt_surf );
        size_type const surfEltId = surfElt.id();
        // Localise only the first dof
        node_type ptReal = ublas::column( surfElt.vertices(), 0 );
        //matrix_node_type ptReal( nRealDim, 1 );
        //ublas::column( ptReal, 0 ) = ublas::column( surfElt.vertices(), 0 );
        // Localise
        auto resLocalisation = locTool->searchElement( ptReal );
        //auto resLocalisation = locTool->run_analysis( ptReal, invalid_size_type_value );

        if( resLocalisation.template get<0>() )
            //if( resLocalisation.template get<0>()[0] )
        {
            size_type localisedEltId = resLocalisation.template get<1>();
            M_intersectingElements.insert( localisedEltId );
            eltsToVisit[surfEltId].insert( localisedEltId );
            M_eltsIntersectedBySurfaceElt[surfEltId].insert( localisedEltId );
        }
    }

    bool eltsToVisitIsEmptyOnAllProc = false;
    while( !eltsToVisitIsEmptyOnAllProc )
    {
        /* Contain the ghost elts which require communication */
        std::map<rank_type, std::vector< std::pair<size_type, size_type> > > dataToSend;
        std::map<rank_type, std::vector< std::pair<size_type, size_type> > > dataToRecv;

        for( auto & eltsToVisitPerSurfElt: eltsToVisit )
        {
            size_type const surfEltId = eltsToVisitPerSurfElt.first;
            auto const& surfElt = meshSurface->element( surfEltId );
            auto & eltsToVisitInSurfElt = eltsToVisitPerSurfElt.second;

            while( !eltsToVisitInSurfElt.empty() )
            {
                size_type const eltId = *(eltsToVisitInSurfElt.begin());
                auto const& elt = meshDistance->element( eltId );
                // Add elements connected via the "faces" which intersect meshSurface element
                for( uint16_type faceLocId = 0; faceLocId < elt.nTopologicalFaces(); faceLocId++ )
                {
                    size_type const neighId = elt.neighbor( faceLocId );

                    if( neighId == invalid_size_type_value )
                        continue;

                    auto const& neigh = meshDistance->element( neighId );
                    rank_type const neighPid = neigh.processId();

                    auto const& face = elt.face( faceLocId );
                    bool intersect = self_type::facesIntersect( face.vertices(), surfElt.vertices() );
                    if( intersect )
                    {
                        // the neighbor elt is a ghost -> send to owner
                        if( neighPid != pidMeshDistance )
                        {
                            size_type neighIdOnOwner = neigh.idInOthersPartitions( neighPid );
                            dataToSend[neighPid].push_back( { neighIdOnOwner, surfEltId } );
                        }
                        // else process
                        else if( M_eltsIntersectedBySurfaceElt[surfEltId].find( neighId ) == M_eltsIntersectedBySurfaceElt[surfEltId].end() )
                        {
                            M_eltsIntersectedBySurfaceElt[surfEltId].insert( neighId );
                            eltsToVisitInSurfElt.insert( neighId );
                        }
                    }
                }
                eltsToVisitInSurfElt.erase( eltId );
            }
        }

        /* Perform ghost elts communications */
        // Send and recv
        int cntRequests = 0;
        for( rank_type neighborRank: neighborSubdomainsMeshDistance )
        {
            mpiRequests[cntRequests++] = meshDistance->worldCommPtr()->localComm().isend( neighborRank, 0, dataToSend[neighborRank] );
            mpiRequests[cntRequests++] = meshDistance->worldCommPtr()->localComm().irecv( neighborRank, 0, dataToRecv[neighborRank] );
        }
        mpi::wait_all( mpiRequests, mpiRequests + cntRequests );

        /* Process received ghosts */
        for( auto const& data: dataToRecv )
        {
            for( auto const& eltRecv: data.second )
            {
                size_type const eltId = eltRecv.first;
                size_type const surfEltId = eltRecv.second;
                if( M_eltsIntersectedBySurfaceElt[surfEltId].find( eltId ) == M_eltsIntersectedBySurfaceElt[surfEltId].end() )
                {
                    M_eltsIntersectedBySurfaceElt[surfEltId].insert( eltId );
                    eltsToVisit[surfEltId].insert( eltId );
                }
            }
        }

        bool eltsToVisitIsEmpty = std::reduce( 
                eltsToVisit.begin(), eltsToVisit.end(),
                true,
                [](bool empty, auto const& set) { return empty && set.second.empty(); }
                );
        eltsToVisitIsEmptyOnAllProc = mpi::all_reduce(
                meshDistance->worldComm(),
                eltsToVisitIsEmpty,
                std::logical_and<bool>()
                );
    }

    delete [] mpiRequests;

    // Merge intersecting elements
    for( auto const& eltsIntersected: M_eltsIntersectedBySurfaceElt )
    {
        M_intersectingElements.insert( 
                eltsIntersected.second.begin(),
                eltsIntersected.second.end()
                );
    }

    M_doUpdateIntersectingElements = false;
}

template< typename MeshType, typename FunctionSpaceType >
    void
DistanceToMesh< MeshType, FunctionSpaceType >::updateUnsignedDistance()
{
    if( !M_unsignedDistance )
        M_unsignedDistance.reset( new element_distance_type( this->functionSpaceDistance(), "unsignedDistance" ) );
    // Find intersecting elements
    if( M_doUpdateIntersectingElements )
        this->updateIntersectingElements();
    // Initialise distance with arbitrarily large value
    //M_unsignedDistance->setConstant( std::numeric_limits<value_type>::max() );
    M_unsignedDistance->setConstant(1e8);
    // Compute exact distance on intersecting elements
    auto const& dofTable = this->functionSpaceDistance()->dof();
    for( auto const& eltsIntersectedBySurfaceEltPair: M_eltsIntersectedBySurfaceElt )
    {
        size_type const surfEltId = eltsIntersectedBySurfaceEltPair.first;
        auto const& surfElt = this->meshSurface()->element( surfEltId );
        auto const& eltsIntersected = eltsIntersectedBySurfaceEltPair.second;

        for( size_type const eltId: eltsIntersected )
        {
            auto const elt = this->meshDistance()->element( eltId );
            for( uint16_type j = 0; j < nDofPerEltDistance; j++ )
            {
                size_type const eltDofGlobalId = dofTable->localToGlobal( elt, j, 0 ).index();

                auto dist = this->distanceDofToSurfaceElt( eltDofGlobalId, surfEltId );
                if( dist < M_unsignedDistance->localToGlobal( eltId, j, 0 ) )
                    M_unsignedDistance->assign( eltId, j, 0, dist );
            }
        }
    }
    // Then perform fast-marching
    *M_unsignedDistance = this->fastMarching()->march( *M_unsignedDistance, this->rangeIntersectingElements() );
    sync( *M_unsignedDistance, "min" );

    M_doUpdateUnsignedDistance = false;
}

template< typename MeshType, typename FunctionSpaceType >
    void
DistanceToMesh< MeshType, FunctionSpaceType >::updateSignedDistance()
{
    if( !M_signedDistance )
        M_signedDistance.reset( new element_distance_type( this->functionSpaceDistance(), "signedDistance" ) );

    // Compute unsigned distance
    if( M_doUpdateUnsignedDistance )
        this->updateUnsignedDistance();

    // Communication
    rank_type const pidMeshDistance = this->meshDistance()->worldCommPtr()->localRank();
    std::set<rank_type> const& neighborSubdomainsMeshDistance = this->meshDistance()->neighborSubdomains();
    int nRequests = 2*neighborSubdomainsMeshDistance.size();
    mpi::request * mpiRequests = new mpi::request[nRequests];

    // Set sign: set negative everywhere, then propagate + sign from the domain boundary
    // note: the element container (VectorUblas) does not provide unary operators at the moment -> TODO: need to improve
    *M_signedDistance = *M_unsignedDistance;
    M_signedDistance->scale( -1. );

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
                auto curPhi = M_signedDistance->localToGlobal( eltId, j, 0 );
                if( curPhi < 0. )
                    M_signedDistance->assign( eltId, j, 0, -curPhi );
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
                        || intersectingElements.find( eltId ) != intersectingElements.end() )
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

    sync( *M_signedDistance, "max" );

    M_doUpdateSignedDistance = false;
}

} // namespace Feel


#endif // _DISTANCE_TO_MESH_HPP
