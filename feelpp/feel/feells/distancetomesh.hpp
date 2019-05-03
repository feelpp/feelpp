#ifndef _DISTANCE_TO_MESH_HPP
#define _DISTANCE_TO_MESH_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feells/reinit_fms_impl.hpp>

#include "geometryconceptwrappers.hpp"

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

    static const uint16_type nDofPerEltDistance = functionspace_distance_type::fe_type::nDof;

    typedef typename functionspace_distance_type::mesh_type mesh_distance_type;
    typedef typename functionspace_distance_type::mesh_ptrtype mesh_distance_ptrtype;

    typedef typename MeshTraits<mesh_distance_type>::elements_reference_wrapper_type elements_reference_wrapper_distance_type;
    typedef typename MeshTraits<mesh_distance_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_distance_ptrtype;
    typedef elements_reference_wrapper_t<mesh_distance_type> range_elements_distance_type;

    //--------------------------------------------------------------------//
    static const uint16_type nRealDim = functionspace_distance_type::nRealDim;
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

    static value_type distancePointToSegment( node_type const& pt, matrix_node_type const& seg );

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
    //TODO
    return false;
}

template< typename MeshType, typename FunctionSpaceType >
bool
DistanceToMesh< MeshType, FunctionSpaceType >::facesIntersect( matrix_node_type const& face1, matrix_node_type const& face2 )
{
    if constexpr ( nRealDim == 2 )
        return segmentsIntersect( face1, face2 );
    else if constexpr ( nRealDim == 3 )
        return trianglesIntersect( face1, face2 );
    else
        CHECK( false ) << "nRealDim must be 2 or 3" << std::endl;
}

template< typename MeshType, typename FunctionSpaceType >
typename DistanceToMesh< MeshType, FunctionSpaceType >::value_type
DistanceToMesh< MeshType, FunctionSpaceType >::distancePointToSegment( node_type const& pt, matrix_node_type const& seg )
{
    return boost::geometry::distance( Feel::detail::geometry::pointWrap<nRealDim>( pt ), Feel::detail::geometry::segmentWrap<nRealDim>( seg ) );
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
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
            intersectingElementsRefWrapper->begin(),
            intersectingElementsRefWrapper->end(),
            intersectingElementsRefWrapper
            );
}

template< typename MeshType, typename FunctionSpaceType >
void
DistanceToMesh< MeshType, FunctionSpaceType >::updateIntersectingElements()
{
    M_intersectingElements.clear();
    M_eltsIntersectedBySurfaceElt.clear();

    auto const& meshSurface = this->meshSurface();
    auto const& meshDistance = this->meshDistance();
    auto locTool = meshDistance->tool_localization();
    locTool->updateForUse();

    auto it_elt_surf = meshSurface->beginOrderedElement();
    auto en_elt_surf = meshSurface->endOrderedElement();
    for( ; it_elt_surf != en_elt_surf; it_elt_surf++ )
    {
        auto const& surfElt = boost::unwrap_ref( *it_elt_surf );
        size_type const surfEltId = surfElt.id();
        // Localise only the first dof
        node_type ptReal = ublas::column( surfElt.vertices(), 0 );
        // Localise
        //auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,theImageElt.vertices()[>theImageElt.G()<],mpl::int_<interpolation_type::isConforming()>());
        auto resLocalisation = locTool->searchElement( ptReal );

        std::unordered_set< size_type > eltsToVisit;
        if( !resLocalisation.template get<0>() )
            Feel::cout << "point localisation failed !" << std::endl;
        else
        {
            size_type localisedEltId = resLocalisation.template get<1>();
            M_intersectingElements.insert( localisedEltId );
            eltsToVisit.insert( localisedEltId );
            M_eltsIntersectedBySurfaceElt[surfEltId].insert( localisedEltId );
        }

        while( !eltsToVisit.empty() )
        {
            size_type const eltId = *(eltsToVisit.begin());
            auto const& elt = meshDistance->element( eltId );
            // Add elements connected via the "faces" which intersect meshSurface element
            for( uint16_type faceLocId = 0; faceLocId < elt.nTopologicalFaces(); faceLocId++ )
            {
                size_type const neighId = elt.neighbor( faceLocId );

                if( neighId == invalid_size_type_value )
                    continue;

                auto const& face = elt.face( faceLocId );
                size_type const faceId = face.id();
                bool intersect = self_type::facesIntersect( face.vertices(), surfElt.vertices() );

                if( intersect && M_eltsIntersectedBySurfaceElt[surfEltId].find( neighId ) == M_eltsIntersectedBySurfaceElt[surfEltId].end() )
                {
                    M_eltsIntersectedBySurfaceElt[surfEltId].insert( neighId );
                    auto const& neigh = meshDistance->element( neighId );
                    eltsToVisit.insert( neighId );
                }
            }
            eltsToVisit.erase( eltId );
        }
        M_intersectingElements.insert( 
                M_eltsIntersectedBySurfaceElt[surfEltId].begin(),
                M_eltsIntersectedBySurfaceElt[surfEltId].end()
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
                auto const pt = boost::get<0>( dofTable->dofPoint( eltDofGlobalId ) );

                auto dist = distancePointToSegment( pt, surfElt.vertices() );
                if( dist < M_unsignedDistance->localToGlobal( eltId, j, 0 ) )
                    M_unsignedDistance->assign( eltId, j, 0, dist );
            }
        }
    }
    // Then perform fast-marching
    *M_unsignedDistance = this->fastMarching()->march( *M_unsignedDistance, this->rangeIntersectingElements() );

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
            if( neighId == invalid_size_type_value 
                    || eltsVisited.find( neighId ) != eltsVisited.end()
                    // stop when reaching an intersecting elt
                    || intersectingElements.find( neighId ) != intersectingElements.end() )
                continue;
            // need to visit neighbor
            eltsToVisit.insert( neighId );
        }
        eltsVisited.insert( eltId );
        eltsToVisit.erase( eltIt );
    }

    M_doUpdateSignedDistance = false;
}

} // namespace Feel


#endif // _DISTANCE_TO_MESH_HPP
