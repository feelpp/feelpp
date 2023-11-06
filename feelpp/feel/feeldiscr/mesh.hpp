//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
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
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 21 May 2011
//! @copyright 2011-2017 Feel++ Consortium
//!
#ifndef FEELPP_MESH_HPP
#define FEELPP_MESH_HPP 1

#include <bitset>

#include <boost/unordered_map.hpp>
#include <boost/version.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/signals2/signal.hpp>

#if defined( __clang__ )
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdivision-by-zero"
#endif
#include <boost/archive/text_oarchive.hpp>
#if defined( __clang__ )
#pragma clang diagnostic pop
#endif
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/mpi/operations.hpp>
#include <feel/feelcore/disablewarnings.hpp>
#include <feel/feelcore/reenablewarnings.hpp>

#include <feel/feelcore/context.hpp>

#include <feel/feelcore/functors.hpp>

#include <feel/feeldiscr/traits.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/mesh0d.hpp>
#include <feel/feelmesh/mesh1d.hpp>
#include <feel/feelmesh/mesh2d.hpp>
#include <feel/feelmesh/mesh3d.hpp>
#include <feel/feeldiscr/mesh_fwd.hpp>

#include <feel/feelalg/boundingbox.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/geomapinv.hpp>

#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/facilities/identity.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/list/cat.hpp>
#include <boost/preprocessor/list/for_each_product.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>

#include <boost/enable_shared_from_this.hpp>

#if defined( FEELPP_HAS_VTK )
#include <feel/feelcore/disablewarnings.hpp>
#include <feel/feelcore/reenablewarnings.hpp>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#endif


namespace Feel
{
enum class IOStatus
{
    isLoading,
    isSaving
};

const size_type EXTRACTION_KEEP_POINTS_IDS = ( 1 << 0 );
const size_type EXTRACTION_KEEP_EDGES_IDS = ( 1 << 1 );
const size_type EXTRACTION_KEEP_FACES_IDS = ( 1 << 2 );
const size_type EXTRACTION_KEEP_VOLUMES_IDS = ( 1 << 3 );
const size_type EXTRACTION_KEEP_ALL_IDS = ( EXTRACTION_KEEP_POINTS_IDS |
                                            EXTRACTION_KEEP_EDGES_IDS |
                                            EXTRACTION_KEEP_FACES_IDS |
                                            EXTRACTION_KEEP_VOLUMES_IDS );
const size_type EXTRACTION_KEEP_MESH_RELATION = ( 1 << 4 );
const size_type EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT = ( 1 << 5 );

} // namespace Feel
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/localization.hpp>

namespace Feel
{
struct PeriodicEntity
{
    PeriodicEntity( int d, int s, int m )
        : dim( d ),
          slave( s ),
          master( m )
    {
    }
    int dim;
    int slave;
    int master;
    std::map<int, int> correspondingVertices;
};
struct MeshMarkerName
{
    std::string name;
    std::vector<int> ids;
};

std::vector<MeshMarkerName> markerMap( int Dim );

//!  partitioner class
template <typename Mesh>
class Partitioner;

class DummySharedFromThis {};
//!
//!  @brief unifying mesh class
//!  @details This structure is an aggregation of elements, faces,
//!  edges(3D) and points data structures and provides a unified
//!  interface with respect to the dimension.
//!
//!  @tparam GeoShape Geometric entities type
//!  @tparam T numerical type for coordinates
//!  @tparam Tag = 0 to discriminate between meshes of the same type
//!
template <typename GeoShape, typename T, int Tag, typename IndexT, bool _EnableSharedFromThis>
class Mesh
    : public mp11::mp_if_c< is_3d_v<GeoShape>,
                            Mesh3D<GeoShape, T, IndexT>,
                            mp11::mp_if_c<
                                is_2d_v<GeoShape>,
                                Mesh2D<GeoShape, T, IndexT>,
                                mp11::mp_if_c<
                                    is_1d_v<GeoShape>,
                                    Mesh1D<GeoShape, T, IndexT>,
                                    Mesh0D<GeoShape, T, IndexT>
                                >
                            >
                          >,
      public boost::addable<Mesh<GeoShape, T, Tag, IndexT>>,
      public std::conditional_t<_EnableSharedFromThis, 
                                std::enable_shared_from_this<Mesh<GeoShape, T, Tag, IndexT, _EnableSharedFromThis>>, 
                                DummySharedFromThis>
{
    using super = mp11::mp_if_c< is_3d_v<GeoShape>,
                            Mesh3D<GeoShape, T, IndexT>,
                            mp11::mp_if_c<
                                is_2d_v<GeoShape>,
                                Mesh2D<GeoShape, T, IndexT>,
                                mp11::mp_if_c<
                                    is_1d_v<GeoShape>,
                                    Mesh1D<GeoShape, T, IndexT>,
                                    Mesh0D<GeoShape, T, IndexT>
                                >
                            >
                        >;

  public:
    //!  @name Constants
    //!
    //! @{

    static inline const uint16_type nDim = GeoShape::nDim;
    static inline const uint16_type nRealDim = GeoShape::nRealDim;
    static inline const uint16_type Shape = GeoShape::Shape;
    static inline const uint16_type nOrder = GeoShape::nOrder;
    static inline const uint16_type tag = Tag;

    //! @}
    //!  @name Typedefs
    //!
    //! @{
    typedef Mesh<GeoShape, T, Tag, IndexT> type;
    typedef std::shared_ptr<type> ptrtype;

    typedef T value_type;
    using index_type = typename super::index_type;
    using size_type = typename super::size_type;

    typedef GeoShape shape_type;
    typedef typename super::return_type return_type;
    typedef typename node<double>::type node_type;

    typedef typename super::super_elements super_elements;
    typedef typename super::elements_type elements_type;
    using element_ptrtype = typename super::element_ptrtype;
    typedef typename super::element_type element_type;
    typedef typename super::element_iterator element_iterator;
    typedef typename super::element_const_iterator element_const_iterator;

    typedef typename super::super_faces super_faces;
    typedef typename super::faces_type faces_type;
    typedef typename super::face_type face_type;
    typedef typename super::face_iterator face_iterator;
    typedef typename super::face_const_iterator face_const_iterator;

    typedef typename super::edge_type edge_type;

    typedef typename super::points_type points_type;
    typedef typename super::point_type point_type;
    typedef typename super::point_iterator point_iterator;
    typedef typename super::point_const_iterator point_const_iterator;

    typedef typename element_type::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;

    typedef typename element_type::gm1_type gm1_type;
    typedef std::shared_ptr<gm1_type> gm1_ptrtype;

    // using gmc_type = typename gm_type::template Context<element_type>;
    // using gmc_ptrtype = std::shared_ptr<gmc_type>;

    typedef Mesh<shape_type, T, Tag, IndexT, EnableSharedFromThis> self_type;
    typedef self_type mesh_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    typedef self_ptrtype mesh_ptrtype;

    typedef typename element_type::template reference_convex<T>::type reference_convex_type;

    //! typedef Partitioner<self_type> partitioner_type;
    //! typedef std::shared_ptr<partitioner_type> partitioner_ptrtype;

    typedef typename super::face_processor_type face_processor_type;
    typedef typename super::face_processor_type element_edge_type;

    using P1_mesh_type = boost::mp11::mp_if_c<is_simplex_v<GeoShape>,
                                       Mesh<Simplex<GeoShape::nDim, 1, GeoShape::nRealDim>, value_type, Tag, IndexT, _EnableSharedFromThis>,
                                       Mesh<Hypercube<GeoShape::nDim, 1, GeoShape::nRealDim>, value_type, Tag, IndexT, _EnableSharedFromThis>>;
    typedef std::shared_ptr<P1_mesh_type> P1_mesh_ptrtype;

    template<int TheTag = 0>
    using parent_mesh_type = parent_mesh_t<mesh_type, TheTag>;
    template<int TheTag = 0>
    using parent_mesh_ptrtype = parent_mesh_ptr_t<mesh_type, TheTag>;

    template<int TheTag = 0>
    using trace_mesh_type = trace_mesh_t<mesh_type, TheTag>;
    template<int TheTag = 0>
    using trace_mesh_ptrtype = trace_mesh_ptr_t<mesh_type, TheTag>;

    template<int TheTag = 0>
    using trace_trace_mesh_type = trace_trace_mesh_t<mesh_type, TheTag>;
    template<int TheTag = 0>
    using trace_trace_mesh_ptrtype = trace_trace_mesh_ptr_t<mesh_type, TheTag>;


    //@}

    //!
    //!  Default mesh constructor
    //
    explicit Mesh( std::string const& name,
                   worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                   std::string const& props = "00001" );

    explicit Mesh( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(), std::string const& props = "00001"  )
        : Mesh( "", worldComm, props ) {}

    ~Mesh() override {}

    void clear() override
    {
        VLOG( 1 ) << "Mesh clear()";
        M_gm.reset();
        M_gm1.reset();
        M_tool_localization.reset();
        super::clear();
    }

    //! return current shared_ptr of type MeshBase
    std::shared_ptr<MeshBase<IndexT>> shared_from_this_meshbase() override 
    { 
        if constexpr ( _EnableSharedFromThis )
            return std::dynamic_pointer_cast<MeshBase<IndexT>>( this->shared_from_this() );  
        else
            return std::shared_ptr<MeshBase<IndexT>>{};
    }

    //! return current shared_ptr of type MeshBase
    std::shared_ptr<const MeshBase<IndexT>> shared_from_this_meshbase() const override 
    { 
        if constexpr ( _EnableSharedFromThis )
            return std::dynamic_pointer_cast<const MeshBase<IndexT>>( this->shared_from_this() );  
        else
            return std::shared_ptr<MeshBase<IndexT>>{};
    }

    //!
    //! @brief allocate a new Mesh
    //! @param worldcomm communicator defaulting to Environment::worldComm()
    //! @return the Mesh shared pointer
    //!
    static mesh_ptrtype New( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
    {
        return std::make_shared<mesh_type>( worldComm );
    }

    self_type& operator+=( self_type const& m );

    //!  @name Accessors
    //!
    //! @{

    //!
    //!  @brief get the number of active elements,faces,edges,points for each partition
    //!  @return the number of active elements
    //!
    size_type statNumElementsActive( rank_type p = invalid_rank_type_value ) const { return std::get<0>( M_statElements[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumElementsAll( rank_type p = invalid_rank_type_value ) const { return std::get<1>( M_statElements[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumFacesActive( rank_type p = invalid_rank_type_value ) const { return std::get<0>( M_statFaces[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumFacesMarkedAll( rank_type p = invalid_rank_type_value ) const { return std::get<1>( M_statFaces[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumEdgesActive( rank_type p = invalid_rank_type_value ) const { return std::get<0>( M_statEdges[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumEdgesMarkedAll( rank_type p = invalid_rank_type_value ) const { return std::get<1>( M_statEdges[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumPointsActive( rank_type p = invalid_rank_type_value ) const { return std::get<0>( M_statPoints[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumPointsAll( rank_type p = invalid_rank_type_value ) const { return std::get<1>( M_statPoints[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumPointsMarkedAll( rank_type p = invalid_rank_type_value ) const { return std::get<2>( M_statPoints[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )] ); }
    size_type statNumVerticesAll( rank_type p = invalid_rank_type_value ) const { return M_statVertices[( ( p != invalid_rank_type_value ) ? p : MeshBase<>::worldComm().localRank() )]; }

    //!
    //!  @brief get the global number of elements
    //!  @details it requires communication in parallel to
    //!  retrieve and sum the number of elements in each subdomain.
    //!  @return the global number of elements
    //!
    size_type numGlobalElements() const override { return M_numGlobalElements; }

    //!
    //!  @return the maximum number of elements over all subdomains
    //!
    size_type maxNumElements() const { return M_maxNumElements; }

    //!
    //!  @brief get the global number of faces
    //!  @details it requires communication in parallel to
    //!  retrieve and sum the number of faces in each subdomain.
    //!  @return the global number of faces
    //!
    size_type numGlobalFaces() const override { return M_numGlobalFaces; }

    //!
    //!  @return the maximum number of faces over all subdomains
    //!
    size_type maxNumFaces() const { return M_maxNumFaces; }

    //!
    //!  @brief get the global number of edges
    //!  @details it requires communication in parallel to
    //!  retrieve and sum the number of edges in each subdomain.
    //!  @return the global number of edges
    //!
    size_type numGlobalEdges() const override { return M_numGlobalEdges; }

    //!
    //!  @return the maximum number of edges over all subdomains
    //!
    size_type maxNumEdges() const { return M_maxNumEdges; }

    //!
    //!  @brief get the global number of points
    //!  @details it requires communication in parallel to
    //!  retrieve and sum the number of points in each subdomain.
    //!  @return the global number of points
    //!
    size_type numGlobalPoints() const override { return M_numGlobalPoints; }

    //!
    //!  @return the maximum number of edges over all subdomains
    //!
    size_type maxNumPoints() const { return M_maxNumPoints; }

    //!
    //!  @brief get the global number of vertices
    //!  @details it requires communication in parallel to
    //!  retrieve and sum the number of vertices in each subdomain.
    //!  @return the global number of vertices
    //!
    size_type numGlobalVertices() const override { return M_numGlobalVertices; }

    //!
    //!  @return the maximum number of vertices over all subdomains
    //!
    size_type maxNumVertices() const { return M_maxNumVertices; }

    //!
    //!  set the global number of points, edges, faces and elements in the mesh
    //!
    void setNumGlobalElements( std::vector<size_type> I )
    {
        CHECK( I.size() == 4 ) << "Invalid information data on elements: num points, num edges, num faces, num elements ";
        M_numGlobalPoints = I[0];
        M_numGlobalEdges = I[1];
        M_numGlobalFaces = I[2];
        M_numGlobalElements = I[3];
    }

    struct UpdateNumGlobalEntitiesForAllReduce : public std::function<boost::tuple<std::vector<size_type>, std::vector<size_type>>(
                                                                        boost::tuple<std::vector<size_type>, std::vector<size_type>> const&,
                                                                        boost::tuple<std::vector<size_type>, std::vector<size_type>> const& )>
    {
        typedef boost::tuple<std::vector<size_type>, std::vector<size_type>> cont_type;
        cont_type operator()( cont_type const& x, cont_type const& y ) const
        {
            auto const& numEltVector_x = boost::get<0>( x );
            auto const& numEltVector_y = boost::get<0>( y );
            auto const& maxVector_x = boost::get<1>( x );
            auto const& maxVector_y = boost::get<1>( y );
            CHECK( numEltVector_x.size() == numEltVector_y.size() ) << "invalid size";
            CHECK( maxVector_x.size() == maxVector_y.size() ) << "invalid size";

            std::vector<size_type> numEltVector_res;
            for ( int k = 0; k < numEltVector_x.size(); ++k )
                numEltVector_res.push_back( numEltVector_x[k] + numEltVector_y[k] );
            std::vector<size_type> maxVector_res;
            for ( int k = 0; k < maxVector_x.size(); ++k )
                maxVector_res.push_back( std::max( maxVector_x[k], maxVector_y[k] ) );

            return boost::make_tuple( numEltVector_res, maxVector_res );
        }
    };

    //!
    //!  @brief compute the global number of elements,faces,points and vertices
    //!  @details it requires communications in parallel to
    //!  retrieve and sum the contribution of each subdomain.
    //!
    template <typename MT>
    void updateNumGlobalElements( typename std::enable_if<is_3d<MT>::value>::type* = nullptr )
    {
        rank_type nProc = MeshBase<>::worldComm().localSize();
        rank_type currentRank = MeshBase<>::worldComm().localRank();

        auto rangeElements = this->elementsWithProcessId( currentRank );
        size_type ne = std::distance( std::get<0>( rangeElements ), std::get<1>( rangeElements ) );
        auto rangeFaces = this->facesWithProcessId( currentRank );
        size_type nf = std::distance( std::get<0>( rangeFaces ), std::get<1>( rangeFaces ) );
        auto rangeEdges = this->edgesWithProcessId( currentRank );
        size_type ned = std::distance( std::get<0>( rangeEdges ), std::get<1>( rangeEdges ) );
        auto rangePoints = this->pointsWithProcessId( currentRank );
        size_type np = std::distance( std::get<0>( rangePoints ), std::get<1>( rangePoints ) );
        size_type nv = this->numVertices();

        size_type neall = this->numElements();
        size_type nfall = this->numFaces();
        size_type nedall = this->numEdges();
        size_type npall = this->numPoints();
        size_type nvall = this->numVertices();

        size_type nfmarkedall = std::count_if( this->beginFace(), this->endFace(), []( auto const& theface ) { return theface.second.hasMarker(); } );
        size_type nedmarkedall = std::count_if( this->beginEdge(), this->endEdge(), []( auto const& theedge ) { return theedge.second.hasMarker(); } );
        size_type npmarkedall = std::count_if( this->beginPoint(), this->endPoint(), []( auto const& thepoint ) { return thepoint.second.hasMarker(); } );

        M_statElements.resize( nProc, std::make_tuple( 0, 0 ) );
        M_statFaces.resize( nProc, std::make_tuple( 0, 0 ) );
        M_statEdges.resize( nProc, std::make_tuple( 0, 0 ) );
        M_statPoints.resize( nProc, std::make_tuple( 0, 0, 0 ) );
        M_statVertices.resize( nProc, 0 );

        if ( MeshBase<>::worldComm().localSize() > 1 )
        {
            std::vector<boost::tuple<boost::tuple<size_type, size_type>, boost::tuple<size_type, size_type>, boost::tuple<size_type, size_type>,
                                     boost::tuple<size_type, size_type, size_type>, size_type>>
                dataRecvFromAllGather;
            auto dataSendToAllGather = boost::make_tuple( boost::make_tuple( ne, neall ), boost::make_tuple( nf, nfmarkedall ), boost::make_tuple( ned, nedmarkedall ),
                                                          boost::make_tuple( np, npall, npmarkedall ), nv );
            mpi::all_gather( MeshBase<>::worldComm(),
                             dataSendToAllGather,
                             dataRecvFromAllGather );
            for ( rank_type p = 0; p < nProc; ++p )
            {
                auto const& dataOnProc = dataRecvFromAllGather[p];
                M_statElements[p] = std::make_tuple( boost::get<0>( boost::get<0>( dataOnProc ) ), boost::get<1>( boost::get<0>( dataOnProc ) ) );
                M_statFaces[p] = std::make_tuple( boost::get<0>( boost::get<1>( dataOnProc ) ), boost::get<1>( boost::get<1>( dataOnProc ) ) );
                M_statEdges[p] = std::make_tuple( boost::get<0>( boost::get<2>( dataOnProc ) ), boost::get<1>( boost::get<2>( dataOnProc ) ) );
                M_statPoints[p] = std::make_tuple( boost::get<0>( boost::get<3>( dataOnProc ) ), boost::get<1>( boost::get<3>( dataOnProc ) ), boost::get<2>( boost::get<3>( dataOnProc ) ) );
                M_statVertices[p] = boost::get<4>( dataOnProc );
            }

            size_type numFaceGlobalCounter = nf, numEdgeGlobalCounter = ned, numPointGlobalCounter = np, numVerticeGlobalCounter = 0;

            for ( auto it = std::get<0>( rangeFaces ), en = std::get<1>( rangeFaces ); it != en; ++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( !face.isInterProcessDomain() )
                    continue;
                if ( face.partition1() < face.partition2() )
                    continue;
                --numFaceGlobalCounter;
            }
            for ( auto it = std::get<0>( rangeEdges ), en = std::get<1>( rangeEdges ); it != en; ++it )
            {
                auto const& edge = unwrap_ref( *it );
                bool countThisEntity = true;
                for ( auto const& ghostData : edge.elementsGhost() )
                {
                    if ( ghostData.first < currentRank )
                    {
                        countThisEntity = false;
                        break;
                    }
                }
                if ( countThisEntity )
                    continue;
                --numEdgeGlobalCounter;
            }
            for ( auto it = std::get<0>( rangePoints ), en = std::get<1>( rangePoints ); it != en; ++it )
            {
                auto const& point = unwrap_ref( *it );
                bool countThisEntity = true;
                for ( auto const& ghostData : point.elementsGhost() )
                {
                    if ( ghostData.first < currentRank )
                    {
                        countThisEntity = false;
                        break;
                    }
                }
                if ( countThisEntity )
                    continue;
                --numPointGlobalCounter;
            }

            std::vector<size_type> numEntitiesGlobalCounter = {numFaceGlobalCounter, numEdgeGlobalCounter, numPointGlobalCounter};
            if ( nOrder > 1 )
                numEntitiesGlobalCounter.push_back( numVerticeGlobalCounter );
            std::vector<size_type> maxNumEntities = {neall, nfall, nedall, npall};
            if ( nOrder > 1 )
                maxNumEntities.push_back( nvall );

            auto dataAllReduce = boost::make_tuple( numEntitiesGlobalCounter, maxNumEntities );
            mpi::all_reduce( MeshBase<>::worldComm(), mpi::inplace( dataAllReduce ), UpdateNumGlobalEntitiesForAllReduce() );
            auto const& numEntitiesGlobalCounterGlobal = boost::get<0>( dataAllReduce );
            auto const& maxNumEntitiesGlobal = boost::get<1>( dataAllReduce );

            auto opBinaryPlusTuple2 = []( std::tuple<size_type, size_type> const& cur, std::tuple<size_type, size_type> const& res ) { return std::make_tuple( std::get<0>( cur ) + std::get<0>( res ), std::get<1>( cur ) + std::get<1>( res ) ); };
            M_numGlobalElements = std::get<0>( std::accumulate( M_statElements.begin(), M_statElements.end(), std::make_tuple( 0, 0 ), opBinaryPlusTuple2 ) );
            M_numGlobalFaces = numEntitiesGlobalCounterGlobal[0];
            M_numGlobalEdges = numEntitiesGlobalCounterGlobal[1];
            M_numGlobalPoints = numEntitiesGlobalCounterGlobal[2];
            M_numGlobalVertices = ( nOrder > 1 ) ? numEntitiesGlobalCounterGlobal[3] : M_numGlobalPoints;

            M_maxNumElements = maxNumEntitiesGlobal[0];
            M_maxNumFaces = maxNumEntitiesGlobal[1];
            M_maxNumEdges = maxNumEntitiesGlobal[2];
            M_maxNumPoints = maxNumEntitiesGlobal[3];
            M_maxNumVertices = ( nOrder > 1 ) ? maxNumEntitiesGlobal[4] : M_maxNumPoints;
        }
        else
        {
            M_statElements[0] = std::make_tuple( ne, neall );
            M_statFaces[0] = std::make_tuple( nf, nfmarkedall );
            M_statEdges[0] = std::make_tuple( ned, nedmarkedall );
            M_statPoints[0] = std::make_tuple( np, npall, npmarkedall );
            M_statVertices[0] = nv;

            M_maxNumElements = ne;
            M_maxNumFaces = nf;
            M_maxNumEdges = ned;
            M_maxNumPoints = np;
            M_maxNumVertices = nv;

            M_numGlobalElements = ne;
            M_numGlobalFaces = nf;
            M_numGlobalEdges = ned;
            M_numGlobalPoints = np;
            M_numGlobalVertices = nv;
        }
    }

    template <typename MT>
    void updateNumGlobalElements( typename std::enable_if<mpl::not_<is_3d<MT>>::value>::type* = nullptr )
    {
        rank_type nProc = MeshBase<>::worldComm().localSize();
        rank_type currentRank = MeshBase<>::worldComm().localRank();
        auto rangeElements = this->elementsWithProcessId( currentRank );
        size_type ne = std::distance( std::get<0>( rangeElements ), std::get<1>( rangeElements ) );
        auto rangeFaces = this->facesWithProcessId( currentRank );
        size_type nf = std::distance( std::get<0>( rangeFaces ), std::get<1>( rangeFaces ) );
        size_type ned = 0;

        auto rangePoints = this->pointsWithProcessId( currentRank );
        size_type np = std::distance( std::get<0>( rangePoints ), std::get<1>( rangePoints ) );
        size_type nv = this->numVertices();

        size_type neall = this->numElements();
        size_type nfall = this->numFaces();
        size_type nedall = 0;
        size_type npall = this->numPoints();
        size_type nvall = this->numVertices();

        size_type nfmarkedall = std::count_if( this->beginFace(), this->endFace(), []( auto const& theface ) { return theface.second.hasMarker(); } );
        //size_type nfmarkedall = std::count_if( this->beginFace(),this->endFace(),
        //                                       [this]( face_type const& theface ) { return theface.marker().isOn() && this->hasFaceMarker( this->markerName( theface.marker().value() ) ) ; } );

        size_type nedmarkedall = 0;
        size_type npmarkedall = std::count_if( this->beginPoint(), this->endPoint(), []( auto const& thepoint ) { return thepoint.second.hasMarker(); } );

        M_statElements.resize( nProc, std::make_tuple( 0, 0 ) );
        M_statFaces.resize( nProc, std::make_tuple( 0, 0 ) );
        M_statEdges.resize( nProc, std::make_tuple( 0, 0 ) );
        M_statPoints.resize( nProc, std::make_tuple( 0, 0, 0 ) );
        M_statVertices.resize( nProc, 0 );

        if ( nProc > 1 )
        {
            std::vector<boost::tuple<boost::tuple<size_type, size_type>, boost::tuple<size_type, size_type>,
                                     boost::tuple<size_type, size_type, size_type>, size_type>>
                dataRecvFromAllGather;
            auto dataSendToAllGather = boost::make_tuple( boost::make_tuple( ne, neall ), boost::make_tuple( nf, nfmarkedall ),
                                                          boost::make_tuple( np, npall, npmarkedall ), nv );
            mpi::all_gather( MeshBase<>::worldComm().localComm(),
                             dataSendToAllGather,
                             dataRecvFromAllGather );

            for ( rank_type p = 0; p < nProc; ++p )
            {
                auto const& dataOnProc = dataRecvFromAllGather[p];
                M_statElements[p] = std::make_tuple( boost::get<0>( boost::get<0>( dataOnProc ) ), boost::get<1>( boost::get<0>( dataOnProc ) ) );
                M_statFaces[p] = std::make_tuple( boost::get<0>( boost::get<1>( dataOnProc ) ), boost::get<1>( boost::get<1>( dataOnProc ) ) );
                M_statPoints[p] = std::make_tuple( boost::get<0>( boost::get<2>( dataOnProc ) ), boost::get<1>( boost::get<2>( dataOnProc ) ), boost::get<2>( boost::get<2>( dataOnProc ) ) );
                M_statVertices[p] = boost::get<3>( dataOnProc );
            }

            size_type numFaceGlobalCounter = nf, numPointGlobalCounter = np, numVerticeGlobalCounter = 0;

            for ( auto it = std::get<0>( rangeFaces ), en = std::get<1>( rangeFaces ); it != en; ++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( !face.isInterProcessDomain() )
                    continue;
                if ( face.partition1() < face.partition2() )
                    continue;
                --numFaceGlobalCounter;
            }
            for ( auto it = std::get<0>( rangePoints ), en = std::get<1>( rangePoints ); it != en; ++it )
            {
                auto const& point = unwrap_ref( *it );
                bool countThisEntity = true;
                for ( auto const& ghostData : point.elementsGhost() )
                {
                    if ( ghostData.first < currentRank )
                    {
                        countThisEntity = false;
                        break;
                    }
                }
                if ( countThisEntity )
                    continue;
                --numPointGlobalCounter;
            }

            std::vector<size_type> numEntitiesGlobalCounter = {numFaceGlobalCounter, numPointGlobalCounter};
            if ( nOrder > 1 )
                numEntitiesGlobalCounter.push_back( numVerticeGlobalCounter );
            std::vector<size_type> maxNumEntities = {neall, nfall, npall};
            if ( nOrder > 1 )
                maxNumEntities.push_back( nvall );

            auto dataAllReduce = boost::make_tuple( numEntitiesGlobalCounter, maxNumEntities );
            mpi::all_reduce( MeshBase<>::worldComm().localComm(), mpi::inplace( dataAllReduce ), UpdateNumGlobalEntitiesForAllReduce() );
            auto const& numEntitiesGlobalCounterGlobal = boost::get<0>( dataAllReduce );
            auto const& maxNumEntitiesGlobal = boost::get<1>( dataAllReduce );

            auto opBinaryPlusTuple2 = []( std::tuple<size_type, size_type> const& cur, std::tuple<size_type, size_type> const& res ) { return std::make_tuple( std::get<0>( cur ) + std::get<0>( res ), std::get<1>( cur ) + std::get<1>( res ) ); };
            M_numGlobalElements = std::get<0>( std::accumulate( M_statElements.begin(), M_statElements.end(), std::make_tuple( 0, 0 ), opBinaryPlusTuple2 ) );
            M_numGlobalFaces = numEntitiesGlobalCounterGlobal[0];
            M_numGlobalEdges = 0;
            M_numGlobalPoints = numEntitiesGlobalCounterGlobal[1];
            M_numGlobalVertices = ( nOrder > 1 ) ? numEntitiesGlobalCounterGlobal[2] : M_numGlobalPoints;

            M_maxNumElements = maxNumEntitiesGlobal[0];
            M_maxNumFaces = maxNumEntitiesGlobal[1];
            M_maxNumEdges = 0;
            M_maxNumPoints = maxNumEntitiesGlobal[2];
            M_maxNumVertices = ( nOrder > 1 ) ? maxNumEntitiesGlobal[3] : M_maxNumPoints;
        }
        else
        {
            M_statElements[0] = std::make_tuple( ne, neall );
            M_statFaces[0] = std::make_tuple( nf, nfmarkedall );
            M_statEdges[0] = std::make_tuple( ned, nedmarkedall );
            M_statPoints[0] = std::make_tuple( np, npall, npmarkedall );
            M_statVertices[0] = nv;

            M_maxNumElements = ne;
            M_maxNumFaces = nf;
            M_maxNumEdges = ned;
            M_maxNumPoints = np;
            M_maxNumVertices = nv;

            M_numGlobalElements = ne;
            M_numGlobalFaces = nf;
            M_numGlobalEdges = ned;
            M_numGlobalPoints = np;
            M_numGlobalVertices = nv;
        }
    }


    template <typename ContType>
    struct UpdateSetForAllReduce : public std::function<ContType(ContType const& ,ContType const&)>
    {
        using cont_type = ContType;
        cont_type operator()( cont_type const& x, cont_type const& y ) const
            {
                cont_type ret = x;
                for ( auto const& yVal : y )
                    ret.insert( yVal );
                return ret;
            }
    };


    void updateMeshFragmentation()
        {
            using _marker_element_type = typename element_type::marker_type;

            std::array<std::set<_marker_element_type>,nDim+1> collectMarkerIds;
            std::vector<ElementsType> ets(collectMarkerIds.size());
            _marker_element_type emptyMarker;

            if constexpr ( nDim >= 1 )
            {
                auto & collectEltMarkerIds = collectMarkerIds[0];
                ets[0] = ElementsType::MESH_ELEMENTS;
                auto it = this->beginOrderedElement();
                auto en = this->endOrderedElement();
                for ( ; it != en; ++it )
                {
                    auto const& elt = unwrap_ref( *it );
#if 0
                    if ( elt.isGhostCell() )
                        continue;
#endif
                    _marker_element_type const& eltMarkers = elt.hasMarker()? elt.marker() : emptyMarker;
                    collectMarkerIds[0].insert( eltMarkers );
                }
            }

            if constexpr ( nDim >= 2 )
            {
                auto & collectFaceMarkerIds = collectMarkerIds[1];
                ets[1] = ElementsType::MESH_FACES;
                auto itf = this->beginOrderedFace();
                auto enf = this->endOrderedFace();
                for ( ; itf != enf; ++itf )
                {
                    auto const& face = unwrap_ref( *itf );
#if 0
                    if ( face.isGhostCell() )
                        continue;
#endif
                    _marker_element_type const& faceMarkers = face.hasMarker()? face.marker() : emptyMarker;
                    collectFaceMarkerIds.insert( faceMarkers );
                }
            }

            if constexpr ( nDim >= 3 )
            {
                auto & collectEdgeMarkerIds = collectMarkerIds[2];
                ets[2] = ElementsType::MESH_EDGES;
                auto ited = this->beginOrderedEdge();
                auto ened = this->endOrderedEdge();
                for ( ; ited != ened; ++ited )
                {
                    auto const& edge = unwrap_ref( *ited );
#if 0
                    if ( edge.isGhostCell() )
                        continue;
#endif
                    _marker_element_type const& edgeMarkers = edge.hasMarker()? edge.marker() : emptyMarker;
                    collectEdgeMarkerIds.insert( edgeMarkers );
                }
            }

            auto & collectPointMarkerIds = collectMarkerIds[nDim];
            ets[nDim] = ElementsType::MESH_POINTS;
            auto itp = this->beginOrderedPoint();
            auto enp = this->endOrderedPoint();
            for ( ; itp != enp; ++itp )
            {
                auto const& point = unwrap_ref( *itp );
#if 0
                if ( point.isGhostCell() )
                    continue;
#endif
                _marker_element_type const& pointMarkers = point.hasMarker()? point.marker() : emptyMarker;
                collectPointMarkerIds.insert( pointMarkers );
            }

            mpi::all_reduce( MeshBase<>::worldComm().localComm(), mpi::inplace( collectMarkerIds.data() ), collectMarkerIds.size(), UpdateSetForAllReduce<std::set<_marker_element_type>>() );

            for (int cd=0;cd<ets.size();++cd )
            {
                ElementsType et = ets[cd];
                M_meshFragmentationByMarker[et].clear();
                auto & mfbym = M_meshFragmentationByMarker[et];
                int fragmentId = 0;
                for ( _marker_element_type const& mIds : collectMarkerIds[cd] )
                    mfbym.emplace( fragmentId++, mIds );
            }
        }

    //! return mesh fragmentation : mapping fragment id to elements marker ids
    std::map<int,typename element_type::marker_type> const& meshFragmentationByMarker( ElementsType et = ElementsType::MESH_ELEMENTS ) const { return M_meshFragmentationByMarker.find( et )->second; }

    //! return mesh fragmentation : mapping fragment id to elements marker ids for all entities
    std::map<ElementsType, std::map<int,typename element_type::marker_type>> const& meshFragmentationByMarkerByEntity() const { return M_meshFragmentationByMarker; }

    //! !
    //! ! @return the topological dimension
    //! ! @code
    //! ! auto mesh = new Mesh<Simplex<2>>;
    //! ! CHECK( mesh->dimension() == 2 ); // true
    //! ! auto mesh2 = new Mesh<Simplex<2,1,3>>;
    //! ! CHECK( mesh2->dimension() == 2 ); // true
    //! ! @endcode
    //! !
    constexpr uint16_type dimension() const noexcept
    {
        return nDim;
    }

    //! !
    //! ! @return the real space dimension
    //! ! @code
    //! ! auto mesh = new Mesh<Simplex<2>>;
    //! ! CHECK( mesh->realDimension() == 2 ); // true
    //! ! auto mesh2 = new Mesh<Simplex<2,1,3>>;
    //! ! CHECK( mesh2->realDimension() == 3 ); // true
    //! ! @endcode
    //! !
    constexpr uint16_type realDimension() const noexcept
    {
        return nRealDim;
    }

    //!
    //!  @return geometric mapping
    //!
    gm_ptrtype const& gm() const
    {
        return M_gm;
    }

    //!
    //!  @return geometric mapping of order 1
    //!
    gm1_ptrtype const& gm1() const
    {
        return M_gm1;
    }

    //!
    //!  @return the geometric mapping
    //!
    gm_ptrtype& gm()
    {
        return M_gm;
    }

    //!
    //!  @return the geometric mapping of order 1
    //!
    gm1_ptrtype& gm1()
    {
        return M_gm1;
    }

    //!
    //!  @return the reference convex associated with the element of the
    //!  mesh
    //!
    reference_convex_type referenceConvex() const
    {
        return reference_convex_type();
    }

    //!
    //!  @return the face index of the face \p n in the element \p e
    //!
    FEELPP_DEPRECATED
    face_processor_type const& localFaceId( element_type const& e,
                                            size_type const n ) const
    {
        return M_e2f.find( std::make_pair( e.id(), n ) )->second;
    }

    //!
    //!  @return the face index of the face \p n in the element \p e
    //!
    FEELPP_DEPRECATED
    face_processor_type const& localFaceId( size_type const e,
                                            size_type const n ) const
    {
        return M_e2f.find( std::make_pair( e, n ) )->second;
    }

    //!
    //!  \return true if the mesh is substructured, false otherwise
    //!
    bool subStructuring() const
    {
        return M_substructuring;
    }

    //!
    //! return true if mesh is unstructured, false otherwise
    //!
    bool isUnstructured() const
    {
        return M_structure_property.test( 0 );
    }
    //!
    //! return true if mesh is structured, false otherwise
    //!
    bool isStructured() const
        {
            return M_structure_property.test( 1 ) || isCartesian();
        }
    //!
    //! return true if mesh is cartesian, false otherwise
    //!
    bool isCartesian() const
        {
            return M_structure_property.test( 2 );
        }
    //! @}

    //!  @name  Mutators
    //!
    //! @{

    void setSubStructuring( bool s )
    {
        M_substructuring = s;
    }
protected:
    void setStructureProperty( std::string const& prop )
        {
            M_structure_property = std::bitset<5>( prop );
        }
public:
    //!
    //!  set the partitioner to \p partitioner
    //!
    void setPartitioner( std::string partitioner )
    {
        //! M_part = partitioner_ptrtype( partitioner_type::New( partitioner ) );
    }

    //!
    //!  @return the face index of the face \p n in the element \p e
    //!
    FEELPP_DEPRECATED
    face_processor_type& localFaceId( element_type const& e,
                                      size_type const n )
    {
        return M_e2f[std::make_pair( e.id(), n )];
    }

    //!
    //!  @return the face index of the face \p n in the element \p e
    //!
    FEELPP_DEPRECATED
    face_processor_type& localFaceId( size_type const e,
                                      size_type const n )
    {
        return M_e2f[std::make_pair( e, n )];
    }

    //!
    //!  @return a localization tool
    //!
    std::shared_ptr<Localization<self_type>> tool_localization()
    {
        return M_tool_localization;
    }

    //!
    //!  @brief get the average h
    //!  @return the average h
    //!
    value_type hAverage() const { return M_h_avg; }

    //!
    //!  @brief get the min h
    //!  @return the min h
    //!
    value_type hMin() const { return M_h_min; }

    //!
    //!  @brief get the max h
    //!  @return the max h
    //!
    value_type hMax() const { return M_h_max; }

    //!
    //!  @return the measure of the mesh (sum of the measure of the elements)
    //!
    value_type measure( bool parallel = true ) const
    {
        if ( parallel )
            return M_meas;
        return M_local_meas;
    }

    //!
    //!  @return the measure of the mesh (sum of the measure of the elements)
    //!
    value_type measureBoundary() const
    {
        return M_measbdy;
    }

    //! @}

    //!
    //!  set the periodic entities
    //!
    void setPeriodicEntities( std::vector<PeriodicEntity> const& e ) { M_periodic_entities = e; }

    //!
    //!  @return true if the mesh has periodic entities
    //!
    bool isPeriodic() const { return M_periodic_entities.empty() == false; }

    //!  @name  Methods
    //!
    //! @{

    //!
    //!  erase element at position \p position
    //!
    //!  @param position \p position is a valid dereferenceable iterator of the index.
    //!  @param modify if true, update mesh data structure and in particular faces
    //!
    //!  @return An iterator pointing to the element immediately
    //!  following the one that was deleted, or \c end() if no such element
    //!  exists.
    //!
    element_iterator eraseElement( element_iterator position, bool modify = true );

    //!
    //! add a new element in the mesh
    //! @param f a new point
    //! @return the new point from the list
    //!
    std::pair<element_iterator,bool> addElement( element_type& f, bool setid = true )
    {
        auto ret = super::addElement( f, setid );
        if ( ret.second )
        {
            if ( !M_geondEltCommon )
                M_geondEltCommon = std::make_shared<GeoNDCommon<typename element_type::super>>( this, this->gm(), this->gm1() );
            auto & eltInserted = ret.first->second;
            eltInserted.setCommonData( M_geondEltCommon.get() );
        }
        return ret;
    }

    //!
    //! move an element into the mesh
    //! @param f a new point
    //! @return the new point from the list
    //!
    std::pair<element_iterator,bool> addElement( element_type&& f )
    {
        auto ret = super::addElement( f );
        if ( ret.second )
        {
            if ( !M_geondEltCommon )
                M_geondEltCommon = std::make_shared<GeoNDCommon<typename element_type::super>>( this, this->gm(), this->gm1() );
            auto & eltInserted = ret.first->second;
            eltInserted.setCommonData( M_geondEltCommon.get() );
        }
        return ret;
    }

    //!
    //!  @brief compute the trace mesh
    //!  @details compute the trace mesh associated to the
    //!  range \p range and tag \p TheTag. The \p range provides
    //!  the iterators over the boundary faces of the mesh.
    //!  \@code
    //!  auto trace_mesh1 = mesh->trace(boundaryfaces(mesh));
    //!  auto trace_mesh2 = mesh->trace(markedfaces(mesh,"marker"));
    //!  @endcode
    //!
    //!  @param range iterator range
    //!  @return the computed trace mesh
    //!
    template <typename RangeT, int TheTag = 0>
    trace_mesh_ptrtype<TheTag>
    trace( RangeT && range, mpl::int_<TheTag> ) const
    {
        DVLOG( 2 ) << fmt::format("[trace] extracting range: {}", range ); 
        return Feel::createSubmesh( _mesh=this->shared_from_this(), _range=std::forward<RangeT>(range) );
    }

    //!
    //!  creates a mesh by iterating over the elements between
    //!  \p begin_elt and \p end_elt and adding them to the
    //!  \p mesh
    //!
    //!  \param mesh new mesh to construct
    //!  \param begin_elt begin iterator
    //!  \param begin_elt end iterator
    //!  \param extraction_policies not in use yet
    //!
    //!  \todo make use of \c extraction_policies
    //!
    template <int TheTag = 0>
    trace_mesh_ptrtype<TheTag>
    trace() const
    {
        return trace( boundaryfaces( this->shared_from_this() ), mpl::int_<TheTag>() );
    }

    template <typename RangeT>
    trace_mesh_ptrtype<Tag>
    trace( RangeT && range ) const
    {
        DVLOG( 2 ) << fmt::format("[trace] extracting range: {}", range ); 
        return Feel::createSubmesh( _mesh=this->shared_from_this(), _range=std::forward<RangeT>(range) );
    }

    template <int TheTag = 0>
    trace_trace_mesh_ptrtype<TheTag>
    wireBasket() const
    {
        return wireBasket( boundaryfaces( this->shared_from_this() ), mpl::int_<TheTag>() );
    }

    template <typename RangeT, int TheTag = 0>
    trace_trace_mesh_ptrtype<TheTag>
    wireBasket( RangeT && range, mpl::int_<TheTag> ) const
    {
        DVLOG( 2 ) << fmt::format("[trace] extracting range: {}", range ); 
        return Feel::createSubmesh( _mesh=this->shared_from_this(), _range=std::forward<RangeT>(range) );
    }

    template <typename RangeT>
    trace_trace_mesh_ptrtype<Tag>
    wireBasket( RangeT && range ) const
    {
        DVLOG( 2 ) << fmt::format("[wirebasked] extracting range: {}", range ); 
        return Feel::createSubmesh( _mesh=this->shared_from_this(), _range=std::forward<RangeT>(range) );
    }

    template <typename Iterator>
    void createSubmesh( self_type& mesh,
                        Iterator const& begin_elt,
                        Iterator const& end_elt,
                        size_type extraction_policies = EXTRACTION_KEEP_ALL_IDS ) const;

    //!
    //!  A special case of \p createSubmesh() : it creates a mesh by
    //!  iterating over all elements having the process id \p pid
    //!
    //!  \param mesh new mesh to construct
    //!  \param pid process id that will select the elements to add
    //!
    FEELPP_DEPRECATED void createSubmeshByProcessId( self_type& mesh, uint16_type pid ) const
    {
#if 0
        this->createSubmesh( mesh,
                             this->beginElementWithProcessId( pid ),
                             this->endElementWithProcessId( pid ) );
#else
        CHECK( 0 ) << "createSubmeshByProcessId is not yet implemented";
#endif
    }

    //!
    //!  Create a P1 mesh from the HO mesh
    //!
    P1_mesh_ptrtype createP1mesh( size_type ctxExtraction = EXTRACTION_KEEP_MESH_RELATION, size_type ctxMeshUpdate = MESH_UPDATE_EDGES | MESH_UPDATE_FACES ) const
    {
        return this->createP1mesh( elements(this->shared_from_this()), ctxExtraction, ctxMeshUpdate );
    }

    template <typename RangeType>
        P1_mesh_ptrtype createP1mesh( RangeType const& range, size_type ctxExtraction = EXTRACTION_KEEP_MESH_RELATION|EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT, size_type ctxMeshUpdate = MESH_UPDATE_EDGES | MESH_UPDATE_FACES ) const;

#if defined( FEELPP_HAS_VTK )
    //!
    //! exporter to VTK data structure
    //!
    typename MeshBase<>::vtk_export_type exportVTK( bool exportMarkers, std::string const& vtkFieldNameMarkers ) const override;
#endif // FEELPP_HAS_VTK

    //!
    //!  update the Marker2 with a range of elements or faces
    //!  if elements -> update marker2 for this elements
    //!  if faces -> update marker2 for this faces
    //!
    template <typename IteratorRange>
    void updateMarker2WithRange( IteratorRange const& range, flag_type flag )
    {
        const size_type iDim = boost::tuples::template element<0, IteratorRange>::type::value;
        this->updateMarker2WithRange( range, flag, mpl::int_<iDim>() );
    }

    //!
    //!  sub method of updateMarker2WithRange : MESH_ELEMENTS
    //!
    template <typename IteratorRange>
    void updateMarker2WithRange( IteratorRange const& range, flag_type flag, mpl::int_<MESH_ELEMENTS> /**/ )
    {
        this->updateMarker2WithRangeElements( range, flag );
    }

    //!
    //!  sub method of updateMarker2WithRange : MESH_FACES
    //!
    template <typename IteratorRange>
    void updateMarker2WithRange( IteratorRange const& range, flag_type flag, mpl::int_<MESH_FACES> /**/ )
    {
        this->updateMarker2WithRangeFaces( range, flag );
    }

    //!
    //!  update the Marker3 with a range of elements or faces
    //!  if elements -> update marker3 for this elements
    //!  if faces -> update marker3 for this faces
    //!
    template <typename IteratorRange>
    void updateMarker3WithRange( IteratorRange const& range, flag_type flag )
    {
        this->updateMarker3WithRange( range, flag, mpl::int_<IteratorRange::entities()>() );
    }

    //!
    //!  sub method of updateMarker3WithRange : MESH_ELEMENTS
    //!
    template <typename IteratorRange>
    void updateMarker3WithRange( IteratorRange const& range, flag_type flag, mpl::int_<MESH_ELEMENTS> /**/ )
    {
        this->updateMarker3WithRangeElements( range, flag );
    }

    //!
    //!  sub method of updateMarker3WithRange : MESH_FACES
    //!
    template <typename IteratorRange>
    void updateMarker3WithRange( IteratorRange const& range, flag_type flag, mpl::int_<MESH_FACES> /**/ )
    {
        this->updateMarker3WithRangeFaces( range, flag );
    }

    //!
    //!  Call the default partitioner (currently \p metis_partition()).
    //!
    void partition( const uint16_type n_parts = 1 ) override;

    //!
    //!  After loading/defining a mesh, we want to have as much locality
    //!  as possible (elements/faces/nodes to be contiguous). In order
    //!  to do that the mesh elements/faces/nodes are renumbered. That
    //!  will be then most helpful when generating the \p Dof table.
    //!  This procedure should work also with
    //!  \p comm().size() == 1
    //!
    void renumber() override
    {
        renumber( mpl::bool_<( nDim > 1 )>() );
    }
    void renumber( std::vector<size_type> const& node_map, mpl::int_<1> );
    void renumber( std::vector<size_type> const& node_map, mpl::int_<2> );
    void renumber( std::vector<size_type> const& node_map, mpl::int_<3> );

    //!
    //!  This function only take sense in the 3D modal case with a simplex mesh.
    //!  In order to construct a global C^0 expansion, we need to assure that
    //!  two contiguous elements share the same top vertex !
    //!  This can be achieve by a local renumbering of the vertex.
    //!

    void localrenumber();

    //!
    //!  check elements permutation and fix it if needed
    //!
    //void checkAndFixPermutation();

    //!
    //!  send the mesh data structure to processor \p p with  \p tag
    //!
    void send( int p, int tag );

    //!
    //!  receive the mesh data structure to processor \p p with  \p tag
    //!
    void recv( int p, int tag );

    FEELPP_DEPRECATED
    void saveMD( std::ostream& out )
    {
        out << "| Shape              |" << Shape << "|\n";
        out << "|---|---|\n";
        out << "| DIM              |" << dimension() << "|\n";
        out << "| Order              |" << nOrder << "|\n";
        out << "| hMin              |" << hMin() << "|\n";
        out << "| hMax              |" << hMax() << "|\n";
        out << "| hAverage              |" << hAverage() << "|\n";
        out << "| nPoints              |" << this->numPoints() << "|\n";
        out << "| nEdges              |" << this->numEdges() << "|\n";
        out << "| nFaces              |" << this->numFaces() << "|\n";
        out << "| nVertices              |" << this->numVertices() << "|\n\n";
    }
#if defined( FEELPP_HAS_HDF5 )
    //!
    //!  load mesh in hdf5
    //!
    void loadHDF5( std::string const& filename, size_type ctxMeshUpdate = MESH_UPDATE_EDGES | MESH_UPDATE_FACES, double scale = 1 )
    {
        ioHDF5( IOStatus::isLoading, filename, ctxMeshUpdate, scale );
    }

    //!
    //!  save mesh in hdf5
    //!
    void saveHDF5( std::string const& filename, double scale = 1 )
    {
        ioHDF5( IOStatus::isSaving, filename, 0, scale );
    }
#endif

  private:
    //! Private Methods
    //! @{

    //! update information for the current object
    void updateInformationObject( nl::json& p ) const override;

#if defined( FEELPP_HAS_HDF5 )
    //!
    //!  save mesh in hdf5
    //!
    void ioHDF5( IOStatus status, std::string const& filename, size_type ctxMeshUpdate = MESH_UPDATE_EDGES | MESH_UPDATE_FACES, double scale = 1 );
#endif

    //! @}
  public:
    //!
    //!  encode the mesh data structure into a tighter data structure and to avoid
    //!  pointers in order to serialize it for saving/loading and
    //!  sending/receiving the mesh
    //!
    void encode();

    //!
    //!  decode the mesh data structure from a tighter data structure and to avoid
    //!  pointers in order to serialize it for saving/loading and
    //!  sending/receiving the mesh
    //!
    void decode();

    template <typename ... Ts>
    void save( Ts && ... v ) const
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        std::string const& name = args.get(_name);
        auto && path = args.get(_path);
        std::string const& type = args.get_else(_type,"binary");
        std::string const& suffix = args.get_else(_suffix,"");
        std::string const& sep = args.get_else(_sep,"");


        if ( !fs::exists( fs::path( path ) ) )
        {
            fs::create_directories( fs::path( path ) );
        }

        std::ostringstream os1;
        os1 << name << sep << suffix << "-" << MeshBase<>::worldComm().globalSize() << "." << MeshBase<>::worldComm().globalRank() << ".fdb";
        fs::path p = fs::path( path ) / os1.str();
        fs::ofstream ofs( p );

        if ( type == "binary" )
        {
            boost::archive::binary_oarchive oa( ofs );
            oa << *this;
        }

        else if ( type == "text" )
        {
            boost::archive::text_oarchive oa( ofs );
            oa << *this;
        }

        else if ( type == "xml" )
        {
            //! boost::archive::xml_oarchive oa(ofs);
            //! oa << *this;
        }
    }
    template <typename ... Ts>
    bool load( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        std::string const& name = args.get(_name );
        auto && path = args.get(_path);
        size_type update = args.get_else(_update,MESH_CHECK | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        std::string const& type = args.get_else(_type,"binary");
        std::string const& suffix = args.get_else(_suffix,"");
        std::string const& sep = args.get_else(_sep,"");

        std::ostringstream os1;
        os1 << name << sep << suffix << "-" << MeshBase<>::worldComm().globalSize() << "." << MeshBase<>::worldComm().globalRank() << ".fdb";
        fs::path p = fs::path( path ) / os1.str();
        if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            std::cout << "try loading " << p.native() << "\n";
        if ( !fs::exists( p ) )
        {
            LOG( INFO ) << "[mesh::load] failed loading " << p.native() << "\n";
            std::ostringstream os2;
            os2 << name << sep << suffix << "-" << MeshBase<>::worldComm().globalSize() << "." << MeshBase<>::worldComm().globalRank();
            p = fs::path( path ) / os2.str();
            if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                std::cout << " now try loading " << p.native() << "\n";

            if ( !fs::exists( p ) )
            {
                LOG( INFO ) << "[mesh::load] failed loading " << p.native() << "\n";
                return false;
            }
        }

        if ( !fs::is_regular_file( p ) )
        {
            LOG( INFO ) << "[mesh::load] failed loading " << p.native() << "\n";
            return false;
        }

        fs::ifstream ifs( p );

        if ( type == "binary" )
        {
            boost::archive::binary_iarchive ia( ifs );
            ia >> *this;
        }

        else if ( type == "text" )
        {
            boost::archive::text_iarchive ia( ifs );
            ia >> *this;
        }

        else if ( type == "xml" )
        {
            //! boost::archive::xml_iarchive ia(ifs);
            //! ia >> *this;
        }
        if ( update )
        {
            this->components().reset();
            this->components().set( update );
            this->updateForUse();
        }

        else
        {
            this->components().reset();
        }
        return true;
    }

    FEELPP_DEFINE_VISITABLE();
    //! @}

    //! private:

    //!
    //!  Finds all the processors that may contain
    //!  elements that neighbor my elements.  This list
    //!  is guaranteed to include all processors that border
    //!  any of my elements, but may include additional ones as
    //!  well.  This method computes bounding boxes for the
    //!  elements on each processor and checks for overlaps.
    //!
    //! void findNeighboringProcessors();

    //!
    //!  This function checks if the local numbering of the mesh elements
    //!  is anticlockwise oriented. For the time being, this function
    //!  only applies to tetrahedra meshes
    //!
    void checkLocalPermutation( mpl::bool_<false> ) const {}
    void checkLocalPermutation( mpl::bool_<true> ) const;

    //!
    //!  update the mesh data structure before using it
    //!
    //!  by default if the number of processors if > 1, we partition
    //!  the mesh. A different behaviour is controlled by setting
    //!  properly \p setComponents(), \p components()
    //!
    void updateForUse() override;

    //! update the mesh when nodes have moved
    void updateForUseAfterMovingNodes( bool upMeasures = true ) override;

    //!
    //!  update hAverage, hMin, hMax, measure of the mesh and measure of the boundary mesh
    //!
    void updateMeasures();

    void meshModified() override
    {
        for ( auto fit = this->beginFace(), fen = this->endFace(); fit != fen; ++fit )
        {
            auto& faceModified = fit->second;
            if ( faceModified.isOnBoundary() && !faceModified.isConnectedTo0() )
            {
                std::cout << "erase boundary face...\n";
                this->eraseFace( fit );
            }
            if ( !faceModified.isOnBoundary() && faceModified.isConnectedTo0() && !faceModified.isConnectedTo1() )
            {
                std::cout << "found boundary face...\n";
                faceModified.setOnBoundary( true );
            }
        }
        this->setUpdatedForUse( false );
        this->updateForUse();
    }

    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void removeMarkerNameWithoutEntity( std::enable_if_t<TheShape::nDim == 0>* = nullptr );
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void removeMarkerNameWithoutEntity( std::enable_if_t<TheShape::nDim == 1>* = nullptr );
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void removeMarkerNameWithoutEntity( std::enable_if_t<TheShape::nDim == 2>* = nullptr );
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void removeMarkerNameWithoutEntity( std::enable_if_t<TheShape::nDim == 3>* = nullptr );

  private:
    FEELPP_NO_EXPORT void propagateMarkers( mpl::int_<1> ) {}
    FEELPP_NO_EXPORT void propagateMarkers( mpl::int_<2> );
    FEELPP_NO_EXPORT void propagateMarkers( mpl::int_<3> );

    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 0>* = nullptr );
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 1>* = nullptr );
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 2>* = nullptr );
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 3>* = nullptr );

    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        if ( Archive::is_saving::value )
        {
            DVLOG( 2 ) << "Serializing mesh(saving) ...\n";
            DVLOG( 2 ) << "encoding...\n";
            encode();
            DVLOG( 2 ) << "loading markers...\n";
            ar & this->M_markername;
            DVLOG( 2 ) << "loading pts...\n";
            ar& M_enc_pts;
            DVLOG( 2 ) << "loading faces...\n";
            ar& M_enc_faces;
            DVLOG( 2 ) << "loading elts...\n";
            ar& M_enc_elts;
        }
        if ( Archive::is_loading::value )
        {
            DVLOG( 2 ) << "Serializing mesh(loading) ...\n";
            DVLOG( 2 ) << "loading markers...\n";
            ar & this->M_markername;
            DVLOG( 2 ) << "loading pts...\n";
            ar& M_enc_pts;
            DVLOG( 2 ) << "loading faces...\n";
            ar& M_enc_faces;
            DVLOG( 2 ) << "loading elts...\n";
            ar& M_enc_elts;
            decode();
        }
    }

  public:
    

    //!  @name  Signals
    //!
    //! @{

#if !defined( __INTEL_COMPILER )
    //!
    //!  mesh changed its connectivity
    //!
    boost::signals2::signal<void( MESH_CHANGES )> meshChanged;

    template <typename meshObserver>
    void addObserver( meshObserver& obs )
    {
        LOG( INFO ) << "Observer attached ! \n";
        meshChanged.connect( obs );
    }
#endif // __INTEL_COMPILER

    void removeFacesFromBoundary( std::initializer_list<uint16_type> markers );

    typename std::set<rank_type>::const_iterator beginNeighborSubdomains() const { return M_neighbor_processors.begin(); }
    typename std::set<rank_type>::const_iterator endNeighborSubdomains() const { return M_neighbor_processors.end(); }
    std::set<rank_type> const& neighborSubdomains() const { return M_neighbor_processors; }
    void addNeighborSubdomain( rank_type p ) { M_neighbor_processors.insert( p ); }

    typename std::set<rank_type>::const_iterator beginFaceNeighborSubdomains() const { return M_face_neighbor_processors.begin(); }
    typename std::set<rank_type>::const_iterator endFaceNeighborSubdomains() const { return M_face_neighbor_processors.end(); }
    std::set<rank_type> const& faceNeighborSubdomains() const { return M_face_neighbor_processors; }
    void addFaceNeighborSubdomain( rank_type p ) { M_face_neighbor_processors.insert( p ); }

    //! @}

  protected:
    //!
    //!  update the adjacency graph elements (useful if coDimensionOne not updated)
    //!
    void updateAdjacencyElements();

    //!
    //!  Update connectivity of entities of codimension 1
    //!
    void updateEntitiesCoDimensionOne() override;
    void updateEntitiesCoDimensionOne( mpl::bool_<true> );
    void updateEntitiesCoDimensionOne( mpl::bool_<false> );

    //!
    //!  Update in ghost cells of entities of codimension 1. Done only for marked entities and ghost faces
    //!
    void updateEntitiesCoDimensionOneMinimal();
    /**
     * Update in ghost cells of entities of codimension 1
     */
    void updateEntitiesCoDimensionGhostCellByUsingBlockingComm();
    void updateEntitiesCoDimensionGhostCellByUsingNonBlockingComm();

    //!
    //!  check mesh connectivity
    //!
    void check() const override;

  private:
    //!
    //!  \sa renumber()
    //!
    void renumber( mpl::bool_<false> ) {}

    //!
    //!  \sa renumber()
    //!
    void renumber( mpl::bool_<true> );

    /**
     * modify edges on boundary in 3D
     */
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void modifyEdgesOnBoundary( face_type& face, std::enable_if_t<TheShape::nDim == 3>* = nullptr );

    /**
     * modify edges on boundary in 2D or 1D
     */
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT void modifyEdgesOnBoundary( face_type& face, std::enable_if_t<TheShape::nDim != 3>* = nullptr );

    /**
     * modify element that may touch the boundary through one of its edge in 1D or 2D
     */
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT bool modifyElementOnBoundaryFromEdge( element_type& elt, std::enable_if_t<TheShape::nDim != 3>* = nullptr );

    /**
     * modify element that may touch the boundary through one of its edge in 3D
     */
    template <typename TheShape = GeoShape>
    FEELPP_NO_EXPORT bool modifyElementOnBoundaryFromEdge( element_type& elt, std::enable_if_t<TheShape::nDim == 3>* = nullptr );

    //!
    //!  update entities on boundary (point, edge, face and element)
    //!
    FEELPP_NO_EXPORT void updateOnBoundary();

    /**
     * fix duplication of point in connection1 with 3d mesh at order 3 and 4
     */
    FEELPP_NO_EXPORT void fixPointDuplicationInHOMesh( element_type& elt, face_type const& face, mpl::true_ );
    FEELPP_NO_EXPORT void fixPointDuplicationInHOMesh( element_type& elt, face_type const& face, mpl::false_ );

  private:

    // entity type -> ( fragment id to elements marker ids )
    std::map<ElementsType, std::map<int,typename element_type::marker_type>> M_meshFragmentationByMarker;

    //! ! communicator
    size_type M_numGlobalElements, M_numGlobalFaces, M_numGlobalEdges, M_numGlobalPoints, M_numGlobalVertices;
    size_type M_maxNumElements, M_maxNumFaces, M_maxNumEdges, M_maxNumPoints, M_maxNumVertices;
    std::vector<std::tuple<size_type, size_type>> M_statElements;          // ( actives,all )
    std::vector<std::tuple<size_type, size_type, size_type>> M_statPoints; // ( actives,all,marked_all )
    std::vector<std::tuple<size_type, size_type>> M_statFaces;             // (actives,marked_all)
    std::vector<std::tuple<size_type, size_type>> M_statEdges;             // (actives,marked_all)
    std::vector<size_type> M_statVertices;

    bool M_is_gm_cached = false;
    gm_ptrtype M_gm;
    gm1_ptrtype M_gm1;

    value_type M_h_avg;
    value_type M_h_min;
    value_type M_h_max;

    //! ! measure of the mesh
    value_type M_meas, M_local_meas;

    //! ! measure of the boundary of the mesh
    value_type M_measbdy, M_local_measbdy;

    //! !sub structuring
    bool M_substructuring;

    //!
    //! 0: unstructured
    //! 1: structured
    //! 2: cartesian
    //! 3: semistructured
    //! 4: boundary layer
    //!
    std::bitset<5> M_structure_property;
    
    //!
    //!  The processors who neighbor the current
    //!  processor
    //!
    //! std::vector<uint16_type> M_neighboring_processors;
    std::set<rank_type> M_neighbor_processors;
    std::set<rank_type> M_face_neighbor_processors;

    //! partitioner_ptrtype M_part;

    //!
    //!  Arrays containing the global ids of Faces of each element
    //!
    boost::unordered_map<std::pair<int, int>, face_processor_type> M_e2f;

    //!
    //!  Arrays containing the global ids of edges of each element
    //!
    boost::multi_array<element_edge_type, 2> M_e2e;

    //!
    //!  periodic entities
    //!
    std::vector<PeriodicEntity> M_periodic_entities;

    //!
    //!  to encode points coordinates
    //!
    std::map<uint64_type, boost::tuple<bool, std::vector<int>, std::vector<double>>> M_enc_pts;

    //!
    //!  to encode elements
    //!
    std::map<uint64_type, std::vector<int>> M_enc_faces;

    //!
    //!  to encode elements
    //!
    std::map<uint64_type, std::vector<int>> M_enc_elts;

    //!
    //!  tool for localize point in the mesh
    //!
    std::shared_ptr<Localization<self_type>> M_tool_localization;

    //! data accessibles in each elements
    std::shared_ptr<GeoNDCommon<typename element_type::super>> M_geondEltCommon;
    //! data accessibles in each faces
    std::shared_ptr<GeoNDCommon<typename face_type::super>> M_geondFaceCommon;
    //! data accessibles in each edges
    std::shared_ptr<GeoNDCommon<typename edge_type::super>> M_geondEdgeCommon;
};


template <int TheTag>
struct Tag : public mpl::int_<TheTag>
{
};

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename Iterator>
void Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::createSubmesh( self_type& new_mesh,
                                                                       Iterator const& begin_elt,
                                                                       Iterator const& end_elt,
                                                                       size_type extraction_policies ) const
{
#if 0
    new_mesh = Feel::createSubmesh( this->shared_from_this(), boost::make_tuple(mpl::int_<MESH_ELEMENTS>(),begin_elt, end_elt), extraction_policies );
#else
    Context policies( extraction_policies );

    DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] start\n";
    //!  make sure it is all empty
    new_mesh.clear();
    //!  inherit all the markers
    new_mesh.M_markername = this->markerNames();
    for ( auto marker : new_mesh.M_markername )
    {
        LOG( INFO ) << "marker name " << marker.first
                    << " id: " << marker.second[0]
                    << " geoe: " << marker.second[1] << "\n";
    }
    //!  How the nodes on this mesh will be renumbered to nodes
    //!  on the new_mesh.
    std::vector<size_type> new_node_numbers( this->numPoints(), invalid_v<size_type> );
    std::vector<size_type> new_vertex( this->numPoints(), 0 );

    //!  the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem = 0;
    size_type n_new_faces = 0;

    for ( Iterator it = begin_elt; it != end_elt; ++it )
    {

        element_type const& old_elem = *it;

        //!  copy element so that we can modify it
        element_type new_elem = old_elem;

        //!  Loop over the nodes on this element.
        for ( unsigned int n = 0; n < old_elem.nPoints(); n++ )
        {
            FEELPP_ASSERT( old_elem.point( n ).id() < new_node_numbers.size() ).error( "invalid point id()" );

            if ( new_node_numbers[old_elem.point( n ).id()] == invalid_v<size_type> )
            {
                new_node_numbers[old_elem.point( n ).id()] = n_new_nodes;

                DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] insert point " << old_elem.point( n ) << "\n";

                point_type pt( old_elem.point( n ) );
                pt.setId( n_new_nodes );

                //!  Add this node to the new mesh
                new_mesh.addPoint( pt );

                DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] number of  points " << new_mesh.numPoints() << "\n";

                //!  Increment the new node counter
                n_new_nodes++;

                if ( n < element_type::numVertices )
                {
                    FEELPP_ASSERT( new_vertex[old_elem.point( n ).id()] == 0 ).error( "already seen this point?" );
                    new_vertex[old_elem.point( n ).id()] = 1;
                }
            }

            //!  Define this element's connectivity on the new mesh
            FEELPP_ASSERT( new_node_numbers[old_elem.point( n ).id()] < new_mesh.numPoints() ).error( "invalid connectivity" );

            DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] adding point old(" << old_elem.point( n ).id()
                       << ") as point new(" << new_node_numbers[old_elem.point( n ).id()]
                       << ") in element " << new_elem.id() << "\n";

            new_elem.setPoint( n, new_mesh.point( new_node_numbers[old_elem.point( n ).id()] ) );
        }

        //!  set id of element
        new_elem.setId( n_new_elem );

        //!  increment the new element counter
        n_new_elem++;

        //!  Add an equivalent element type to the new_mesh
        new_mesh.addElement( new_elem );

        //!  Maybe add faces for this element
        for ( unsigned int s = 0; s < old_elem.numTopologicalFaces; s++ )
        {
            if ( !old_elem.facePtr( s ) ) continue;

#if 0
            std::cout << "local face id: " << s
                      << " global face id: " << old_elem.face( s ).id() << "\n";
#endif
            //!  only add face on the boundary: they have some data
            //!  (boundary ids) which cannot be retrieved otherwise
            //! if ( old_elem.neighbor(s) == invalid_v<size_type> )
            size_type global_face_id = old_elem.face( s ).id();

            if ( this->hasFace( global_face_id ) )
            {

                //! std::cout << "found face " << global_face_id << "\n";
                //!  get the corresponding face
                face_type const& old_face = old_elem.face( s );
                face_type new_face = old_face;

                //!  disconnect from elements of old mesh,
                //!  the connection will be redone in
                //!  \c updateForUse()
                new_face.disconnect();

                //! std::cout << "disconnect face\n";
                //!  update points info
                for ( uint16_type p = 0; p < new_face.nPoints(); ++p )
                {
#if 0
                    std::cout << "add new point " << new_face.point( p ).id() << " to face \n";
                    std::cout << "add old point " << old_face.point( p ).id() << " to face \n";
                    std::cout << "new point id " << new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] << "\n";
#endif
                    new_face.setPoint( p, new_mesh.point( new_node_numbers[old_elem.point( old_elem.fToP( s, p ) ).id()] ) );
                }

                new_face.setId( n_new_faces++ );
#if 0
                std::cout << "face id" << new_face.id()
                          << " marker1 : " << new_face.marker()
                          << " old marker1 : " << old_face.marker()
                          << "\n";
#endif
                //!  add it to the list of faces
                new_mesh.addFace( new_face );
            }
        }
    }

    new_mesh.setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0 ) );

    DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] update face/edge info if necessary\n";
    //!  Prepare the new_mesh for use
    new_mesh.components().set( MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES | MESH_CHECK );
    new_mesh.updateForUse();

    DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] stop\n";
#endif
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename RangeType>
typename Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::P1_mesh_ptrtype
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::createP1mesh( RangeType const& range, size_type ctxExtraction, size_type ctxMeshUpdate ) const
{
    if constexpr (false /*nOrder == 1*/ )
         return Feel::createSubmesh( _mesh=this->shared_from_this(), _range=elements(this->shared_from_this()), _context=ctxExtraction, _update=ctxMeshUpdate );
    else
    {
    
    std::shared_ptr<SubMeshData<>> smd;
    Context c( ctxExtraction );
    bool keepMeshRelation = c.test( EXTRACTION_KEEP_MESH_RELATION );
    if ( keepMeshRelation )
        smd.reset( new SubMeshData<>( this->shared_from_this() ) );

    P1_mesh_ptrtype new_mesh{std::make_shared<P1_mesh_type>( this->worldCommPtr() )};

    //!  How the nodes on this mesh will be renumbered to nodes on the new_mesh.
    std::unordered_map<size_type, size_type> new_node_numbers;
    std::unordered_map<size_type, int> new_vertex;
    std::unordered_map<size_type, size_type> new_element_numbers; 

    const int nProc = new_mesh->worldComm().localSize();

    //!  the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem = 0;
    size_type n_new_faces = 0;

    //!  inherit the table of markersName
    for( auto const& itMark : this->markerNames() )
    {
        new_mesh->addMarkerName( itMark.first, itMark.second[0], itMark.second[1] );
    }

    //!  data useful for parallism
    std::map<int, std::set<boost::tuple<size_type, size_type>>> memoryGhostId;
    //! std::set< int > setOfRecvProc;
    std::vector<int> nbMsgToRecv( nProc, 0 );

    auto it = this->beginElement();
    auto const en = this->endElement();
    for ( ; it != en; ++it )
    {
        element_type const& old_elem = it->second;

        //!  create a new element
        typename P1_mesh_type::element_type new_elem;
        //!  set id of element
        new_elem.setId( n_new_elem );
        new_element_numbers[old_elem.id()] = n_new_elem;
        // set element markers
        new_elem.setMarkers( old_elem.markers() );
        // partitioning update
        new_elem.setProcessIdInPartition( old_elem.pidInPartition() );
        new_elem.setProcessId( old_elem.processId() );
        new_elem.setNeighborPartitionIds( old_elem.neighborPartitionIds() );

        //!  Loop over the P1 nodes on this element.
        for ( uint16_type n = 0; n < element_type::numVertices; n++ )
        {
            auto const& old_point = old_elem.point( n );

            //! if ( !new_node_numbers[old_point.id()] )
            if ( new_node_numbers.find( old_point.id() ) == new_node_numbers.end() )
            {
                new_node_numbers[old_point.id()] = n_new_nodes;
                DVLOG( 2 ) << "[Mesh<Shape,T>::createP1mesh] insert point " << old_point << "\n";
                //typename P1_mesh_type::point_type pt( old_point );
                //pt.setId( n_new_nodes );
                typename P1_mesh_type::point_type pt( n_new_nodes, old_point, false, false );
                pt.setProcessId( old_point.processId() );
                pt.clearElementsGhost();

                //!  Add this node to the new mesh
                new_mesh->addPoint( pt );
                DVLOG( 2 ) << "[Mesh<Shape,T>::createSubmesh] number of  points " << new_mesh->numPoints() << "\n";
                //!  Increment the new node counter
                n_new_nodes++;
                FEELPP_ASSERT( !new_vertex[old_point.id()] ).error( "already seen this point?" );
                new_vertex[old_point.id()] = 1;
            }
            //!  Define this element's connectivity on the new mesh
            //! FEELPP_ASSERT ( new_node_numbers[old_elem.point( n ).id()] < new_mesh->numPoints() ).error( "invalid connectivity" );
            DVLOG( 2 ) << "[Mesh<Shape,T>::createP1mesh] adding point old(" << old_point.id()
                       << ") as point new(" << new_node_numbers[old_point.id()]
                       << ") in element " << new_elem.id() << "\n";
            //!  add point in element
            new_elem.setPoint( n, new_mesh->point( new_node_numbers[old_point.id()] ) );
        } //for ( uint16_type n=0; n < element_type::numVertices; n++ )

        //!  Add an equivalent element type to the new_mesh
        auto eit = new_mesh->addElement( new_elem );
        auto const& e = eit.first->second;
        if ( keepMeshRelation )
            smd->bm.insert( typename SubMeshData<>::bm_type::value_type( e.id(), old_elem.id() ) );

        //!  increment the new element counter
        n_new_elem++;

#if 0
        //!  Maybe add faces for this element
        for ( unsigned int s=0; s<old_elem.numTopologicalFaces; s++ )
        {
            if ( !old_elem.facePtr( s ) ) continue;
            //!  only add face on the boundary: they have some data
            //!  (boundary ids) which cannot be retrieved otherwise
            const size_type global_face_id = old_elem.face( s ).id();
            if ( this->hasFace( global_face_id ) )
            {
                //!  get the corresponding face
                face_type const& old_face = old_elem.face( s );
                //! if ( old_face.marker().isOff() ) continue;
                typename P1_mesh_type::face_type new_face;
                //!  disconnect from elements of old mesh,
                //!  the connection will be redone in updateForUse()
                new_face.disconnect();
                //!  is on boundary
                new_face.setOnBoundary( old_face.isOnBoundary() );
                //!  set id of face
                new_face.setId( n_new_faces );
                // set face markers
                new_face.setMarkers( old_face.markers() );
                // partitioning update
                new_face.setProcessIdInPartition( old_face.pidInPartition() );
                new_face.setProcessId(old_face.processId());
                new_face.clearIdInOthersPartitions();
                new_face.setNeighborPartitionIds(old_face.neighborPartitionIds());
                //!  update P1 points info
                for ( uint16_type p = 0; p < face_type::numVertices; ++p )
                {
                    //! new_face.setPoint( p, new_mesh->point( new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] ) );
                    new_face.setPoint( p, new_mesh->point( new_node_numbers[ old_face.point(p).id()] ) );
                }
                //!  add it to the list of faces
                new_mesh->addFace( new_face );
                //!  increment the new face counter
                ++n_new_faces;
            } // if ( this->hasFace( global_face_id ) )
        } // for ( unsigned int s=0; s<old_elem.numTopologicalFaces; s++ )

#endif
        if ( old_elem.isGhostCell() )
        {
            DVLOG( 2 ) << "element " << old_elem.id() << " is a ghost cell\n";
            for ( auto it_pid = old_elem.idInOthersPartitions().begin(), en_pid = old_elem.idInOthersPartitions().end(); it_pid != en_pid; ++it_pid )
            {
                DVLOG( 2 ) << " " << it_pid->first << "-" << it_pid->second << "-" << old_elem.pidInPartition() << "-" << new_mesh->worldComm().localRank();
                const int procToSend = it_pid->first;
                DCHECK( procToSend != old_elem.pidInPartition() ) << "invalid\n";
                memoryGhostId[procToSend].insert( boost::make_tuple( new_elem.id(), it_pid->second ) );
            }
        }
        else if ( old_elem.numberOfNeighborPartitions() /*old_elem.numberOfPartitions()*/ > 0 )
        {
#if 0
            setOfRecvProc.insert( old_elem.neighborPartitionIds().begin(), old_elem.neighborPartitionIds().end() );
#else
            auto itneighbor = old_elem.neighborPartitionIds().begin();
            auto const enneighbor = old_elem.neighborPartitionIds().end();
            for ( ; itneighbor != enneighbor; ++itneighbor )
                nbMsgToRecv[*itneighbor]++;
#endif
        }
    } // end for it

    //!  add marked faces in P1 mesh
    auto face_it = this->beginFace();
    auto face_en = this->endFace();
    for ( ; face_it != face_en; ++face_it )
    {
        auto const& old_face = face_it->second;
        if ( !old_face.hasMarker() ) continue;

        typename P1_mesh_type::face_type new_face;
        //!  is on boundary
        new_face.setOnBoundary( old_face.isOnBoundary() );
        //!  set id of face
        new_face.setId( n_new_faces );
        // set face markers
        new_face.setMarkers( old_face.markers() );
        // partitioning update
        new_face.setProcessIdInPartition( old_face.pidInPartition() );
        new_face.setProcessId( old_face.processId() );
        new_face.clearIdInOthersPartitions();
        new_face.setNeighborPartitionIds( old_face.neighborPartitionIds() );
        //!  update P1 points info
        for ( uint16_type p = 0; p < face_type::numVertices; ++p )
        {
            //! new_face.setPoint( p, new_mesh->point( new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] ) );
            new_face.setPoint( p, new_mesh->point( new_node_numbers[old_face.point( p ).id()] ) );
        }
        //!  add it to the list of faces
        new_mesh->addFace( new_face );
        //!  increment the new face counter
        ++n_new_faces;
    }

#if 0
    if ( nProc > 1 )
    {
        std::map< int, std::vector<size_type> > memoryMpiMsg;

        auto itghostproc = memoryGhostId.begin();
        auto const enghostproc = memoryGhostId.end();
        for ( ; itghostproc!=enghostproc ; ++itghostproc )
        {
            const int procToSend = itghostproc->first;
            const int sizeMsgToSend = itghostproc->second.size();
            std::vector<size_type> dataToSend( sizeMsgToSend );
            memoryMpiMsg[procToSend].resize( sizeMsgToSend );

            auto itghostelt = itghostproc->second.begin();
            auto const enghostelt = itghostproc->second.end();
            for ( int k=0 ; itghostelt!=enghostelt ; ++itghostelt,++k )
            {
                dataToSend[k] = itghostelt->template get<1>();
                memoryMpiMsg[procToSend][k] = itghostelt->template get<0>();
            }
            new_mesh->worldComm().localComm().send(procToSend, 0, dataToSend);
        }

        auto itrecvproc = setOfRecvProc.begin();
        auto const enrecvproc = setOfRecvProc.end();
        for ( ; itrecvproc!=enrecvproc ; ++itrecvproc )
        {
            const int procToRecv = *itrecvproc;
            std::vector<size_type> dataToRecv;
            new_mesh->worldComm().localComm().recv(procToRecv, 0, dataToRecv);

            const int nbDataToTreat = dataToRecv.size();
            std::vector<size_type> dataToSend( nbDataToTreat );
            for ( int k=0;k<nbDataToTreat;++k )
            {
                CHECK( dataToRecv[k]!=invalid_v<size_type> ) << "invalid id recv \n";
                dataToSend[k] = new_element_numbers[ dataToRecv[k] ];
            }
            new_mesh->worldComm().localComm().send(procToRecv, 1, dataToSend);
        }

        itghostproc = memoryGhostId.begin();
        for ( ; itghostproc!=enghostproc ; ++itghostproc )
        {
            const int procToRecv = itghostproc->first;
            std::vector<size_type> dataToRecv;
            new_mesh->worldComm().localComm().recv(procToRecv, 1, dataToRecv);
            const int nbDataToTreat = dataToRecv.size();
            for ( int k=0;k<nbDataToTreat;++k )
            {
                auto eltToUpdate = new_mesh->elementIterator( memoryMpiMsg[procToRecv][k]/*e.id()*/ );
                new_mesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( procToRecv, dataToRecv[k]/*idEltAsked*/ ) );
            }
        }

    } // if ( nProc > 1 )
#else
    if ( nProc > 1 )
    {
        std::map<int, std::vector<size_type>> memoryMpiMsg;
        std::vector<int> nbMsgToSend( nProc, 0 );

        auto itghostproc = memoryGhostId.begin();
        auto const enghostproc = memoryGhostId.end();
        for ( ; itghostproc != enghostproc; ++itghostproc )
        {
            const int procToSend = itghostproc->first;
            const int sizeMsgToSend = itghostproc->second.size();
            memoryMpiMsg[procToSend].resize( sizeMsgToSend );

            auto itghostelt = itghostproc->second.begin();
            auto const enghostelt = itghostproc->second.end();
            for ( int k = 0; itghostelt != enghostelt; ++itghostelt, ++k )
            {
                const size_type idInOtherPart = itghostelt->template get<1>();
                new_mesh->worldComm().localComm().send( procToSend, nbMsgToSend[procToSend], idInOtherPart );
                ++nbMsgToSend[procToSend];
                memoryMpiMsg[procToSend][k] = itghostelt->template get<0>();
            }
            //! CHECK( nbMsgToSend[procToSend] == sizeMsgToSend ) << "invalid data to send\n";
        }

#if !defined( NDEBUG )
        //!  check nbMsgToRecv computation
        std::vector<int> nbMsgToRecv2( nProc, 0 );
        mpi::all_to_all( new_mesh->worldComm().localComm(),
                         nbMsgToSend,
                         nbMsgToRecv2 );
        for ( int proc = 0; proc < nProc; ++proc )
            CHECK( nbMsgToRecv[proc] == nbMsgToRecv2[proc] ) << "partitioning data incorect "
                                                             << "myrank " << MeshBase<>::worldComm().localRank() << " proc " << proc
                                                             << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
                                                             << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc] << "\n";
#endif

        //!  recv dof asked and re-send dof in this proc
        for ( int procToRecv = 0; procToRecv < nProc; ++procToRecv )
        {
            for ( int cpt = 0; cpt < nbMsgToRecv[procToRecv]; ++cpt )
            {
                //! recv
                size_type idEltRecv;
                new_mesh->worldComm().localComm().recv( procToRecv, cpt, idEltRecv );

                const size_type idEltAsked = new_element_numbers[idEltRecv];
                DCHECK( idEltAsked != invalid_v<size_type> ) << "invalid elt id\n";

                new_mesh->worldComm().localComm().send( procToRecv, cpt, idEltAsked );
            }
        }

        itghostproc = memoryGhostId.begin();
        for ( ; itghostproc != enghostproc; ++itghostproc )
        {
            const int procToRecv = itghostproc->first;

            auto itghostelt = itghostproc->second.begin();
            auto const enghostelt = itghostproc->second.end();
            for ( int k = 0; itghostelt != enghostelt; ++itghostelt, ++k )
            {
                size_type idEltAsked;
                new_mesh->worldComm().localComm().recv( procToRecv, k, idEltAsked );

                auto& eltModified = new_mesh->elementIterator( memoryMpiMsg[procToRecv][k] /*e.id()*/ )->second;
                eltModified.setIdInOtherPartitions( procToRecv, idEltAsked );
            }
        }

    } // if ( nProc > 1 )
#endif

#if 0
    new_mesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0,
                                               []( size_type lhs, std::pair<size_type,int> const& rhs )
                                               {
                                                   return lhs+rhs.second;
                                               } ) );
#endif

    //!  Prepare the new_mesh for use
    new_mesh->components().reset();
    new_mesh->components().set( ctxMeshUpdate ); //MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    //!  run intensive job
    new_mesh->updateForUse();

    if ( keepMeshRelation )
        new_mesh->setSubMeshData( smd );

    return new_mesh;
    }
}

#if defined( FEELPP_HAS_VTK )
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
typename MeshBase<>::vtk_export_type
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::exportVTK( bool exportMarkers, std::string const& vtkFieldNameMarkers ) const
{
    /* Compute element type from the parameters */
    typedef typename
        /* if (Mdim == 1) */
        mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<1>>,
                 mpl::identity<vtkLine>,
                 /* if (Mdim == 2) */
                 typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2>>,
                                   /* if(MShape == SHAPE_TRIANGLE) */
                                   typename mpl::if_<mpl::equal_to<mpl::int_<Shape>, mpl::size_t<SHAPE_TRIANGLE>>,
                                                     mpl::identity<vtkTriangle>,
                                                     mpl::identity<vtkQuad>>::type,
                                   /* if (Mdim == 3) */
                                   typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3>>,
                                                     /* if(MShape == SHAPE_TETRA) */
                                                     typename mpl::if_<mpl::equal_to<mpl::int_<Shape>, mpl::size_t<SHAPE_TETRA>>,
                                                                       mpl::identity<vtkTetra>,
                                                                       mpl::identity<vtkHexahedron>>::type,
                                                     /* We should normally not reach this case */
                                                     /* anyway we set a default vtkTetra for face type */
                                                     mpl::identity<vtkTetra>>::type>::type>::type::type vtkelement_type;

    vtkSmartPointer<vtkUnstructuredGrid> out = vtkSmartPointer<vtkUnstructuredGrid>::New();

    rank_type currentRank = MeshBase<>::worldComm().localRank();

    auto rangePoints = this->pointsWithProcessId( currentRank );
    auto itPoint = std::get<0>( rangePoints );
    auto enPoint = std::get<1>( rangePoints );

    size_type nPoints = std::distance( itPoint, enPoint );
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToFloat();
    points->SetNumberOfPoints( nPoints );

    auto mappingWithVTK = std::make_shared<typename MeshBase<>::MappingDataWithVTK>();
    auto& mapPointsFeelIdToVTKId = mappingWithVTK->mapPointsFeelIdToVTKId;
    size_type cptPoint = 0;
    std::vector<double> node( 3, 0. );
    for ( ; itPoint != enPoint; ++itPoint, ++cptPoint )
    {
        auto const& pt = unwrap_ref( *itPoint );
        if ( nDim >= 1 )
            node[0] = pt.node()[0];
        if ( nDim >= 2 )
            node[1] = pt.node()[1];
        if ( nDim >= 3 )
            node[2] = pt.node()[2];
        points->SetPoint( ( vtkIdType )( cptPoint ), node.data() );
        mapPointsFeelIdToVTKId[pt.id()] = cptPoint;
    }
    out->SetPoints( points );

    auto rangeElements = this->elementsWithProcessId( currentRank );
    auto itElement = std::get<0>( rangeElements );
    auto enElement = std::get<1>( rangeElements );

    auto& mapElementsFeelIdToVTKId = mappingWithVTK->mapElementsFeelIdToVTKId;
    vtkSmartPointer<vtkelement_type> cell = vtkSmartPointer<vtkelement_type>::New();
    size_type nCell = std::distance( itElement, enElement );
    for ( ; itElement != enElement; ++itElement )
    {
        auto const& elt = unwrap_ref( *itElement );
        for ( uint16_type p = 0; p < element_type::numVertices /*numPoints*/; ++p )
        {
            size_type ptIdFeel = elt.point( p ).id();
            size_type ptIdVTK = mapPointsFeelIdToVTKId.find( ptIdFeel )->second;
            cell->GetPointIds()->SetId( p, ptIdVTK );
        }
        vtkIdType newCellId = out->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
        mapElementsFeelIdToVTKId[elt.id()] = ( size_type )( newCellId );
    }

    if ( exportMarkers )
    {
        vtkSmartPointer<vtkFloatArray> da = vtkSmartPointer<vtkFloatArray>::New();
        da->SetName( vtkFieldNameMarkers.c_str() );

        da->SetNumberOfComponents( 1 );
        da->SetNumberOfTuples( nCell );

        float* arrayValue = new float[1];

        for ( itElement = std::get<0>( rangeElements ); itElement != enElement; ++itElement )
        {
            auto const& elt = unwrap_ref( *itElement );
            if ( elt.hasMarker() )
                arrayValue[0] = (float)( elt.marker().value() );
            else
                arrayValue[0] = 0;
            size_type vtkEltId = mapElementsFeelIdToVTKId.find( elt.id() )->second;
            da->SetTuple( vtkEltId, arrayValue );
        }

        delete[] arrayValue;

        out->GetCellData()->AddArray( da );
    }

    return std::make_pair( out, mappingWithVTK );
}
#endif // FEELPP_HAS_VTK

//! Fill mesh properties in journal publication
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateInformationObject( nl::json& p ) const
{
    if ( p.contains( "shape" ) )
        return;
    p.emplace( "shape", Shape::name() );
    p.emplace( "dim", this->dimension() );
    p.emplace( "order", this->nOrder );
    p.emplace( "real_dim", this->realDimension() );
    p.emplace( "h_min", this->hMin() );
    p.emplace( "h_max", this->hMax() );
    p.emplace( "h_average", this->hAverage() );
    p.emplace( "n_points", this->numGlobalPoints() );
    if ( nOrder > 1 )
        p.emplace( "n_vertices", this->numGlobalVertices() );
    if ( this->dimension() == 3 )
        p.emplace( "n_edges", this->numGlobalEdges() );
    p.emplace( "n_faces", this->numGlobalFaces() );
    p.emplace( "n_elements", this->numGlobalElements() );

    rank_type nProc = MeshBase<>::worldComm().localSize();
    p.emplace( "n_partition", nProc );

    if ( nProc > 1 )
    {
        nl::json::array_t ptEltActive, ptEltAll, ptFacesActives, ptEdgesActives, ptPointsActives;
        for ( rank_type p = 0; p < nProc; ++p )
        {
            ptEltActive.push_back( this->statNumElementsActive( p ) );
            ptEltAll.push_back( this->statNumElementsAll( p ) );;
            if ( this->dimension() > 1 )
                ptFacesActives.push_back( this->statNumFacesActive( p ) );
            if ( this->dimension() > 2 )
                ptEdgesActives.push_back( this->statNumEdgesActive( p ) );
            if ( this->dimension() > 0 )
                ptPointsActives.push_back( this->statNumPointsActive( p ) );
        }
        p["/partitioning/n_elements"_json_pointer] = ptEltActive;
        p["/partitioning/n_elements_with_ghost"_json_pointer] = ptEltAll;
        if ( this->dimension() > 1 )
            p["/partitioning/n_faces"_json_pointer] = ptFacesActives;
        if ( this->dimension() > 2 )
            p["/partitioning/n_edges"_json_pointer] = ptEdgesActives;
        if ( this->dimension() > 0 )
            p["/partitioning/n_points"_json_pointer] = ptPointsActives;
    }
}

template<typename MeshType>
class MeshInverse
        : public boost::mp11::mp_if_c<
              MeshType::shape_type::is_simplex,
              GeoMapInverse<MeshType::nDim, MeshType::nOrder, MeshType::nRealDim, typename MeshType::value_type, Simplex>,
              GeoMapInverse<MeshType::nDim, MeshType::nOrder, MeshType::nRealDim, typename MeshType::value_type, Hypercube>
            >
{
public:
    using mesh_t = decay_type<MeshType>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    using value_type = typename mesh_t::value_type;
    using size_type = typename mesh_t::size_type;
    using super = boost::mp11::mp_if_c<
        MeshType::shape_type::is_simplex,
        GeoMapInverse<MeshType::nDim, MeshType::nOrder, MeshType::nRealDim, typename MeshType::value_type, Simplex>,
        GeoMapInverse<MeshType::nDim, MeshType::nOrder, MeshType::nRealDim, typename MeshType::value_type, Hypercube>
    >;
    using gic_type = typename super::gic_type;
    using gm_type = typename super::gm_type;

    using super1 = boost::mp11::mp_if_c<
        MeshType::shape_type::is_simplex,
        GeoMapInverse<MeshType::nDim, 1, MeshType::nRealDim, typename MeshType::value_type, Simplex>,
        GeoMapInverse<MeshType::nDim, 1, MeshType::nRealDim, typename MeshType::value_type, Hypercube>
    >;

    using gic1_type = typename super1::gic_type;
    using gm1_type = typename super1::gm_type;
    explicit MeshInverse( std::shared_ptr<mesh_t> const& m )
        : super(),
            M_mesh( m )
    {
    }

    ~MeshInverse() override
    {
    }

    size_type nPointsInConvex( size_type i ) const
    {
        auto itFind = M_pts_cvx.find( i );
        if ( itFind != M_pts_cvx.end() )
            return M_pts_cvx.find( i )->second.size();
        else
            return 0;
    }
    void pointsInConvex( size_type i, std::vector<boost::tuple<size_type, uint16_type>>& itab ) const
    {
        const size_type nPts = this->nPointsInConvex( i );
        itab.resize( nPts );
        if ( nPts == 0 ) return;

        auto it = M_pts_cvx.find( i )->second.begin();
        auto const en = M_pts_cvx.find( i )->second.end();
        for ( size_type j = 0; it != en; ++it )
            itab[j++] = boost::make_tuple( it->first, it->second );
    }

    const boost::unordered_map<size_type, node_type>& referenceCoords( void )
    {
        return M_ref_coords;
    }

    //!
    //!  distribute the points of the mesh in a kdtree
    //!
    //!  extrapolation = false : Only the points inside the mesh are distributed.
    //!  extrapolation = true  : Try to project the exterior points.
    //!  \todo : for extrapolation, verify that all the points have been taken
    //!         into account, else test them on the frontiere convexes.
    //!
    void distribute( bool extrapolation = false );

    private:
    mesh_ptr_t M_mesh;
    boost::unordered_map<size_type, boost::unordered_map<size_type, uint16_type>> M_pts_cvx;
    typedef typename std::map<size_type, uint16_type>::const_iterator map_iterator;
    //! typedef typename node<value_type>::type node_type;
    boost::unordered_map<size_type, node_type> M_ref_coords;
    boost::unordered_map<size_type, double> M_dist;
    boost::unordered_map<size_type, size_type> M_cvx_pts;

}; // MeshInverse


template <typename MeshType>
void  
MeshInverse<MeshType>::distribute( bool extrapolation )
{
    auto rangeElements = elements(M_mesh);
    auto el_it = rangeElements.begin();
    auto el_en = rangeElements.end();
    const size_type nActifElt = std::distance( el_it, el_en );
    if ( nActifElt == 0 ) return;

    using element_type = element_t<mesh_t>;
    BoundingBox<> bb;

    typename gm_type::reference_convex_type refelem;
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( M_mesh->gm(),
                                                                                         refelem.points() ) );
    boost::unordered_map<size_type, bool> npt;

    M_dist.clear();
    M_ref_coords.clear();
    M_cvx_pts.clear();
    M_pts_cvx.clear();
    M_pts_cvx.clear();

    KDTree::points_type boxpts;

    auto __c = M_mesh->gm()->template context<vm::JACOBIAN | vm::KB | vm::POINT>( unwrap_ref( *el_it ),__geopc );
    VLOG( 2 ) << "[Mesh::Inverse] distribute mesh points ion kdtree\n";

    for ( ; el_it != el_en; ++el_it )
    {
        auto const& elt = boost::unwrap_ref( *el_it );
        // get geometric transformation
        __c->template update<vm::JACOBIAN | vm::KB | vm::POINT>( elt );
        gic_type gic( M_mesh->gm(), elt );

        // create bounding box
        //bb.make( el_it->points() );
        bb.make( elt.G() );

        for ( size_type k = 0; k < bb.min.size(); ++k )
        {
            bb.min[k] -= 1e-10;
            bb.max[k] += 1e-10;
        }

        DVLOG( 2 ) << "G = " << elt.G() << " min = " << bb.min << ", max = " << bb.max << "\n";

        // check if the points
        this->pointsInBox( boxpts, bb.min, bb.max );

        DVLOG( 2 ) << "boxpts size = " << boxpts.size() << "\n";

        for ( size_type i = 0; i < boxpts.size(); ++i )
        {
            size_type index = boost::get<1>( boxpts[i] );

            if ( ( !npt[index] ) || M_dist[index] < 0 )
            {
                // check if we are in
                gic.setXReal( boost::get<0>( boxpts[i] ) );
                bool isin;
                value_type dmin;
                boost::tie( isin, dmin ) = refelem.isIn( gic.xRef() );
                bool tobeadded = extrapolation || isin;

                DVLOG( 2 ) << "i = " << i << " index = " << index << " isin = " << ( isin >= -1e-10 )
                           << " xref = " << gic.xRef() << " xreal = " << boost::get<0>( boxpts[i] )
                           << " tobeadded= " << tobeadded << " dist=" << dmin << "\n";

                if ( tobeadded && npt[index] )
                {
                    if ( dmin > M_dist[index] )
                        M_pts_cvx[M_cvx_pts[index]].erase( index );

                    else
                        tobeadded = false;
                }

                if ( tobeadded )
                {
                    M_ref_coords[index] = gic.xRef();
                    M_dist[index] = dmin;
                    M_cvx_pts[index] = elt.id();
                    M_pts_cvx[elt.id()][index] = boost::get<2>( boxpts[i] );
                    npt[index] = true;
                }
            }
        }
    }

    VLOG( 2 ) << "[Mesh::Inverse] distribute mesh points in kdtree done\n";
}


} // namespace Feel

//#if !defined(FEELPP_INSTANTIATION_MODE)
#include <feel/feeldiscr/meshimpl.hpp>
//#endif

#endif /* FEELPP_MESH_HPP */
