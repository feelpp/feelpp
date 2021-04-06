/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
  \file geond.hpp
*/
#ifndef _GEOND_HH_
#define _GEOND_HH_

#include <boost/numeric/ublas/storage.hpp>

#include <feel/feelmesh/geo0d.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/im.hpp>

namespace Feel
{
template <typename IndexT>
class MeshBase;

template <int Dim, int Order, int RealDim, template <uint16_type, uint16_type, uint16_type> class Entity, typename T>
struct GT_Lagrange;

template <class Convex, uint16_type O, typename T2>
class Gauss;
template <int IMORDER,
          int DIM,
          template <uint16_type, uint16_type, uint16_type> class Entity,
          template <class Convex, uint16_type O, typename T2> class QPS,
          typename T>
struct IMGeneric;

/// \cond detail
namespace detail
{
/**
 * \class ReversePoint
 * An utility to  invert point numbering on a GeoShape
 */
template <typename GeoShape>
class ReversePoint
{
  public:
    uint16_type operate( uint16_type const& point )
    {
        return ( point < GeoShape::numVertices ) ? GeoShape::numVertices - point : GeoShape::numPoints - point + GeoShape::numVertices;
    }
};
} // namespace detail
/// \endcond

/**
 * @class GeoNDCommon
 * @ingroup Mesh
 * @brief Common data shared in a collection of multi-dimensional geometrical entity.
 */
template <typename GeoNDType>
class GeoNDCommon
{
  public:
    typedef GeoNDType geond_type;
    typedef typename geond_type::gm_ptrtype gm_ptrtype;
    typedef typename geond_type::gm1_ptrtype gm1_ptrtype;
    typedef typename geond_type::index_type index_type;

    GeoNDCommon()
        : M_mesh( nullptr )
    {
    }
    explicit GeoNDCommon( MeshBase<index_type> const* mesh )
        : M_mesh( mesh )
    {
    }
    GeoNDCommon( MeshBase<index_type> const* mesh, gm_ptrtype const& gm, gm1_ptrtype const& gm1 )
        : M_mesh( mesh ),
          M_gm( gm ),
          M_gm1( gm1 )
    {
    }
    GeoNDCommon( GeoNDCommon const& e ) = default;
    GeoNDCommon( GeoNDCommon&& e ) = default;

    MeshBase<index_type> const* mesh() const { return M_mesh; }
    gm_ptrtype const& gm() const { return M_gm; }
    gm1_ptrtype const& gm1() const { return M_gm1; }

    void setMesh( MeshBase<index_type> const* mesh ) { M_mesh = mesh; }
    void setGm( gm_ptrtype const& gm, gm1_ptrtype const& gm1 )
    {
        M_gm = gm;
        M_gm1 = gm1;
    }

  private:
    // mesh to which the geond element belongs to
    MeshBase<index_type> const* M_mesh;
    // geometric mapping
    gm_ptrtype M_gm;
    gm1_ptrtype M_gm1;
};

template< typename T>
class GeoNDMeasuresStorage
{
    typedef T value_type;
public :
    explicit GeoNDMeasuresStorage( uint16_type nF ) : M_face_measures( nF ) {}

    GeoNDMeasuresStorage( GeoNDMeasuresStorage const& ) = default;
    GeoNDMeasuresStorage( GeoNDMeasuresStorage && ) = default;

    //! \return the measure of the element
    value_type measureBIS() const
    {
        return M_measure.value_or( 0. );
    }

    //! \return the measure of the element face \p f
    value_type faceMeasureBIS( uint16_type f ) const
    {
        return M_face_measures[f].value_or( 0. );
    }

    //! set the measure of point element neighbors
    void setMeasurePointElementNeighborsBIS( value_type meas )
    {
        M_measurePointElementNeighbors = meas;
    }
    //! \return the measure of point element neighbors
    value_type measurePointElementNeighbors() const
    {
        return M_measurePointElementNeighbors.value_or( 0. );
    }

    void setMeasure( value_type meas ) const
        {
            M_measure = meas;
        }

    void setFaceMeasure( uint16_type f, value_type meas ) const
        {
            M_face_measures[f] = meas;
        }
private :

    mutable std::optional<value_type> M_measure;
    mutable std::vector<std::optional<value_type>> M_face_measures;
    mutable std::optional<value_type> M_measurePointElementNeighbors;
};

template< typename T>
class GeoNDNoMeasuresStorage
{
public :
    explicit GeoNDNoMeasuresStorage( uint16_type nF ) {}
    GeoNDNoMeasuresStorage( GeoNDNoMeasuresStorage const& ) = default;
    GeoNDNoMeasuresStorage( GeoNDNoMeasuresStorage && ) = default;
private :
};

/**
 * @class GeoND
 * @ingroup Mesh
 * @brief Base class for Multi-dimensional basis Geometrical Entities.
 *
 */
template <uint16_type Dim,
          typename GEOSHAPE,
          typename T = double,
          typename IndexT = uint32_type,
          typename POINTTYPE = Geo0D<Dim, T, IndexT>,
          bool UseMeasuresStorage = false >
class GeoND
    : public GeoEntity<GEOSHAPE, T>,
      public mpl::if_<mpl::bool_<UseMeasuresStorage>,  mpl::identity<GeoNDMeasuresStorage<T>>,  mpl::identity<GeoNDNoMeasuresStorage<T>>>::type::type
{
    typedef GeoEntity<GEOSHAPE, T> super;
    typedef typename mpl::if_<mpl::bool_<UseMeasuresStorage>,  mpl::identity<GeoNDMeasuresStorage<T>>,  mpl::identity<GeoNDNoMeasuresStorage<T>>>::type::type super2;

  public:
    typedef T value_type;
    using index_type = IndexT;
    using size_type = index_type;
    typedef GeoND<Dim, GEOSHAPE, value_type, index_type, POINTTYPE, UseMeasuresStorage> self_type;
    typedef self_type element_type;

    typedef GEOSHAPE GeoShape;
    typedef POINTTYPE PointType;

    //typedef PointType point_type;
    typedef typename super::face_type face_type;
    using point_type = PointType;

    static inline const size_type Shape = super::Shape;
    static inline const uint16_type numPoints = super::numPoints;
    static inline const uint16_type numVertices = super::numVertices;
    static inline const uint16_type numLocalPoints = super::numPoints;
    static inline const uint16_type numLocalEdges = super::numEdges;
    static inline const uint16_type numLocalVertices = super::numVertices;
    static inline const int numFaces = super::numFaces;
    static inline const int numEdges = super::numEdges;
    static inline const int numTopologicalFaces = super::numTopologicalFaces;
    static inline const uint16_type numNeighbors = super::numTopologicalFaces;

    typedef typename std::vector<point_type*>::iterator point_iterator;
    typedef typename std::vector<point_type*>::const_iterator point_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename node<value_type>::type node_type;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nRealDim = super::nRealDim;

    template <int GmOrder>
    struct GetGm
    {

        typedef typename mpl::if_<mpl::bool_<GeoShape::is_hypercube>,
                                  mpl::identity<GT_Lagrange<nDim, GmOrder, nRealDim, Hypercube, T>>,
                                  mpl::identity<GT_Lagrange<nDim, GmOrder, nRealDim, Simplex, T>>>::type::type type;
        typedef std::shared_ptr<type> ptrtype;
    };
    typedef typename GetGm<nOrder>::type gm_type;
    typedef typename GetGm<nOrder>::ptrtype gm_ptrtype;

    typedef typename GetGm<1>::type gm1_type;
    typedef typename GetGm<1>::ptrtype gm1_ptrtype;

    typedef typename gm_type::super::reference_convex_type reference_convex_type;
    typedef typename gm1_type::super::reference_convex_type reference_convex1_type;

    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>,
                                            mpl::int_<1>>,
                              mpl::identity<vertex_permutation_type>,
                              typename mpl::if_<mpl::equal_to<mpl::int_<nDim>,
                                                              mpl::int_<2>>,
                                                mpl::identity<edge_permutation_type>,
                                                mpl::identity<face_permutation_type>>::type>::type::type permutation_type;

    static inline constexpr uint16_type meas_quad_order = ( nOrder - 1 ) * nDim;

    template <int GeoOrder>
    struct GetImMeasure
    {
        // quadrature formula used in entity measure (for ho geo, need to check)
        static inline const uint16_type quad_order = ( nOrder - 1 ) * nDim;
        typedef typename mpl::if_<mpl::bool_<GeoShape::is_hypercube>,
                                  mpl::identity<typename IMGeneric<quad_order, Dim, Hypercube, Gauss, value_type /*double*/>::type>,
                                  mpl::identity<typename IMGeneric<quad_order, Dim, Simplex, Gauss, value_type /*double*/>::type>>::type::type type;
    };
    typedef typename GetImMeasure<nOrder>::type quad_meas_type;
    typedef typename GetImMeasure<1>::type quad_meas1_type;

    /**
     * default constructor
     */
    GeoND()
        : super( 0 ),
          super2( numTopologicalFaces ),
          M_points( numPoints, nullptr ),
        M_neighbors( 0 ),
        M_commonData( nullptr )
    {
    }

    /**
     * constructor from an id
     *
     * @param id identifier for the element to store
     *
     */
    explicit GeoND( size_type id )
        : super( id ),
          super2( numTopologicalFaces ),
          M_points( numPoints, nullptr ),
        M_neighbors( 0 ),
        M_commonData( nullptr )
    {
    }

    GeoND( GeoND const& e ) = default;
    GeoND( GeoND&& e )
        : super( std::move( e ) ),
          super2( std::move( e ) ),
          M_points( std::move( e.M_points ) ),
          M_neighbors( std::move( e.M_neighbors ) ),
          M_markers( std::move( e.M_markers ) ),
          M_commonData( std::move( e.M_commonData ) )
    {
        //std::cout << "GeoND move ctor\n";
    }

    GeoND& operator=( GeoND const& ) = default;
    GeoND& operator=( GeoND&& e )
    {
        super::operator=( std::move( e ) );
        super2::operator=( std::move( e ) );
        M_points = std::move( e.M_points );
        M_neighbors = std::move( e.M_neighbors );
        M_markers = std::move( e.M_markers );
        M_commonData = std::move( e.M_commonData );
        //M_face_measures = std::move( e.M_face_measures );
        //std::cout << "GeoND move assign\n";
        return *this;
    }

    /**
     * destructor, make it virtual for derived classes
     */
    ~GeoND() override
    {
    }

    /**
     * set common data
     */
    void setCommonData( GeoNDCommon<self_type> * commonData )
    {
        M_commonData = commonData;
    }

    //! return the geometric mapping if a mesh was set
    gm_ptrtype gm() const
    {
        if ( M_commonData != nullptr )
            return M_commonData->gm();
        return gm_ptrtype();
    }

    //! return the geometric mapping if a mesh was set
    gm1_ptrtype gm1() const
    {
        if ( M_commonData != nullptr )
            return M_commonData->gm1();
        return gm1_ptrtype();
    }

    /**
     * \return the mesh to which this geometric entity belongs to
     */
    MeshBase<> const* mesh() const
    {
        if ( M_commonData != nullptr )
            return M_commonData->mesh();
        return nullptr;
    }

    /**
     * \return true if points have been inserted in elements, false
     * otherwise
     */
    bool hasPoints() const
    {
        auto it = std::find( M_points.begin(), M_points.end(), nullptr );
        return it == M_points.end();
    }
    //! @return true if the element has at least a point with marker1 active
    bool hasPointWithMarker() const
    {
        auto it = std::find_if( M_points.begin(), M_points.end(),
                                []( auto const* p ) { return p->hasMarker(); } );
        return it != M_points.end();
    }

    /**
     * \return the number of points in convex
     */
    uint16_type nPoints() const
    {
        return numPoints;
    }

    /**
     * the number of neighbors is equal to the number of
     * faces. Neighbors are stored as pointers and some of them can be
     * null pointers if the corresponding face is on the boundary.
     *
     *\return the number of neighbors
     */
    uint16_type nNeighbors() const
    {
        return M_neighbors.size(); //numNeighbors;
    }

    /**
     * Neighbors are stored as pointers and some of them can be null
     * pointers if the corresponding face is on the boundary.
     *
     * \return the pair neighbor \p n index and process \p id it belongs to
     */
    size_type neighbor( uint16_type n ) const
    {
        return M_neighbors[n];
    }

    /**
     *
     * \return true if a neighbor has been registered
     */
    size_type hasNeighbor( uint16_type n ) const
    {
        if ( M_neighbors.empty() )
            return false;
        if ( n >= M_neighbors.size() )
            return false;
        return M_neighbors[n] != invalid_v<size_type>;
    }

    /**
     * set the \p n -th neighbor with \p neigh
     */
    void setNeighbor( uint16_type n, size_type neigh_id )
    {
        if ( M_neighbors.empty() )
        {
            M_neighbors.reserve( numNeighbors );
            M_neighbors.resize( numNeighbors, invalid_v<size_type> );
        }
        M_neighbors[n] = neigh_id;
    }

    bool isNeighbor( self_type const& G ) const
    {
        for ( uint16_type i = 0; i < this->nNeighbors(); ++i )
            if ( this->neighbor( i ) == G.id() ) return true;

        return false;
    }

    /**
     * \return the barycenter of the element
     */
    node_type barycenter() const
    {
        auto M = glas::average( this->G() );
        return ublas::column( M, 0 );
    }

    /**
     * \return the barycenter at the faces of the element
     */
    node_type faceBarycenter( uint16_type f ) const
    {
        constexpr int nPtsInFace = GEOSHAPE::topological_face_type::numPoints;

        node_type n( nRealDim );
        em_node_type<value_type> en( n.data().begin(), n.size() );
        en.setZero();
        for ( uint16_type p = 0; p < nPtsInFace; ++p )
        {
            // get pt id in element  from local pt id  in face
            int ptid = this->fToP( f, p );
            auto const& node = this->point( ptid ).node();
            em_node_type<value_type> emnode( const_cast<typename point_type::value_type /*double*/*>( node.data().begin() ), node.size() );
            en += emnode;
        }
        en /= nPtsInFace;
        return n;
    }

    /**
     * \return the barycenters at the faces of the element
     */
    FEELPP_DEPRECATED matrix_node_type faceBarycenters() const
    {
        matrix_node_type n;
        return n;
    }

    /**
     * \return permutation
     */
    permutation_type permutation( uint16_type /*f*/ ) const
    {
        return permutation_type();
    }

    //! \return identity permutation of entity co-dim 1 (topological faces)
    virtual typename super::template PermutationSubEntity<1> permutation( uint16_type /*i*/, mpl::int_<1> ) const
    {
        return typename super::template PermutationSubEntity<1>();
    }

    //!
    //! @return points set with constness
    //!
    std::vector<point_type*> const& points() const { return M_points; }

    //!
    //! @return points set
    //!
    std::vector<point_type*>& points() { return M_points; }

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    point_type& point( uint16_type i )
    {
        return *( static_cast<point_type*>( M_points[i] ) );
    }

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    point_type const& point( uint16_type i ) const
    {
        return *( static_cast<point_type*>( M_points[i] ) );
    }
    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    point_type* pointPtr( uint16_type i )
    {
        return M_points[i];
    }

    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    point_type const* pointPtr( uint16_type i ) const
    {
        return M_points[i];
    }

    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    point_type& reversepoint( uint16_type const i )
    {
        return *( static_cast<point_type*>( M_points[Feel::detail::ReversePoint<GEOSHAPE>::operate( i )] ) );
    }

    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    point_type const& reversepoint( uint16_type const i ) const
    {
        return *( static_cast<point_type*>( M_points[Feel::detail::ReversePoint<GEOSHAPE>::operate( i )] ) );
    }

    /**
     * Inserts a point.  Uses point references
     * put point
     */
    void setPoint( uint16_type i, point_type& p );

    /**
     * show information about the geoND
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c output stream
     *
     * @return the output stream
     */
    std::ostream& showMe( bool verbose = false, std::ostream& c = std::cout ) const;

    /**
     * Swap Points
     *
     *  This is a member function to be used ONLY by routines for checking or
     *  amending meshes. You must give the local id
     *
     * @param pt1 1st point to swap with 2nd point
     * @param pt2 2nd point to swap with 1st point
     */
    void swapPoints( const uint16_type& pt1, const uint16_type& pt2 );

    /**
     *  Exchange Points
     *
     *  Exchanges points according to a list of old2new local id numbering !
     *  old2new[i] is the new local id of a point whose old local id was ! i+1
     *  (remeber the numbering from 1 of the id's!. This is a member function
     *  to be used ONLY by routines for checking or amending meshes. You must
     *  give uint16_type (which start
     */
    void exchangePoints( const uint16_type otn[numPoints] );

    /**
     * matrix of geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of geometric nodes
     */
    matrix_node_type G() const
    {
        matrix_node_type G( nRealDim, numPoints );
        for ( uint16_type i = 0; i < numPoints; ++i )
            ublas::column( G, i ) = M_points[i]->node();

        return G;
    }

    /**
     * matrix of vertices geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of vertices geometric nodes
     */
    //matrix_node_type  vertices() const { return ublas::subrange( M_G, 0, nRealDim, 0, numVertices ); }
    matrix_node_type vertices() const
    {
        matrix_node_type G( nRealDim, numVertices );
        for ( uint16_type i = 0; i < numVertices; ++i )
            ublas::column( G, i ) = M_points[i]->node();
        return G;
        // return ublas::subrange( M_G, 0, nRealDim, 0, numVertices );
    }

    //! update G with all points of the element
    template <typename TheNodesType>
    void updateG( TheNodesType & G,
                  std::enable_if_t<is_eigen_matrix_v<TheNodesType> >* = nullptr ) const
    {
        DCHECK( nRealDim == G.rows() && G.cols() >= numPoints ) << "G is not compatible";
        for ( uint16_type i = 0; i < numPoints; ++i )
            G.col( i ) = em_fixed_size_cmatrix_t<nRealDim,1,value_type>( M_points[i]->node().data().begin() );
    }

    //! update G with all vertices of the element
    template <typename TheNodesType>
    void updateVertices( TheNodesType & G,
                         std::enable_if_t<is_eigen_matrix_v<TheNodesType> >* = nullptr ) const
    {
        DCHECK( nRealDim == G.rows() && G.cols() >= numVertices ) << "G is not compatible";
        for ( uint16_type i = 0; i < numVertices; ++i )
            G.col( i ) = em_fixed_size_cmatrix_t<nRealDim,1,value_type>( M_points[i]->node().data().begin() );
    }

    point_iterator beginPoint()
    {
        return M_points.begin();
    }
    point_const_iterator beginPoint() const
    {
        return M_points.begin();
    }
    point_iterator endPoint()
    {
        return M_points.end();
    }
    point_const_iterator endPoint() const
    {
        return M_points.end();
    }

    /**
     * get the max length of the edges of the element
     *
     *
     * @return the max length of the edges of the element
     */
    value_type h() const
    {
        matrix_node_type nodesG = this->G();
        em_matrix_col_type<value_type> G( const_cast<value_type*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );
        value_type res = 0.;
        for ( uint16_type __e = 0; __e < numLocalEdges; ++__e )
        {
            int col1 = this->eToP( __e, 0 );
            int col2 = this->eToP( __e, 1 );
            value_type r = ( G.col( col1 ) - G.col( col2 ) ).norm();

            res = ( res > r ) ? res : r;
        }
        return res;
    }
    /**
     * @brief get the minimum edge length in the element
     * @return the minimum edge length in the element
     */
    value_type hMin() const
    {
        matrix_node_type nodesG = this->G();
        em_matrix_col_type<value_type> G( const_cast<value_type*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );

        value_type res = 1e10;
        for ( uint16_type __e = 0; __e < numLocalEdges; ++__e )
        {
            int col1 = this->eToP( __e, 0 );
            int col2 = this->eToP( __e, 1 );
            value_type r = ( G.col( col1 ) - G.col( col2 ) ).norm();
            res = ( res > r ) ? r : res;
        }
        return res;
    }
    /**
     * get the max length of the edge in the local face \c f
     *
     * @param f local id of the face
     *
     * @return the max length of the edges of the local face
     */
    value_type hFace( uint16_type f ) const
    {
        if ( nRealDim == 1 )
            return 1;

        constexpr int nEdges = GEOSHAPE::topological_face_type::numEdges;
        matrix_node_type nodesG = this->G();
        em_matrix_col_type<value_type> G( const_cast<value_type*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );

        value_type res = 0.;
        for ( uint16_type e = 0; e < nEdges; ++e )
        {
            // get edge id in face from local edge id
            int edg = this->f2e( f, ( nDim == 2 ) ? f : e );
            int col1 = this->eToP( edg, 0 );
            int col2 = this->eToP( edg, 1 );
            value_type r = ( G.col( col1 ) - G.col( col2 ) ).norm();
            res = ( res > r ) ? res : r;
        }
        return res;
    }

    value_type hEdge( uint16_type f ) const
    {
        int col1 = this->eToP( f, 0 );
        int col2 = this->eToP( f, 1 );
        auto const& node1 = this->point( col1 ).node();
        em_node_type<value_type> emnode1( const_cast<value_type*>( node1.data().begin() ), node1.size() );
        auto const& node2 = this->point( col2 ).node();
        em_node_type<value_type> emnode2( const_cast<value_type*>( node2.data().begin() ), node2.size() );
        return ( emnode1 - emnode2 ).norm();
    }

    struct tt
    {
        static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
        {
            return super::eToP( _localFace, _point );
        }
    };


    /**
     * \return the normals at the barycenter of the faces
     */
    matrix_node_type normals() const
    {
        matrix_node_type _normals( nRealDim, numTopologicalFaces );
        if ( nDim != nRealDim )
        {
            CHECK( false ) << "normal when nDim != nRealDim is not implemented";
            return _normals;
        }
        gm1_ptrtype gm1 = this->gm1();
        if ( !gm1.use_count() )
            gm1 = gm1_ptrtype( new gm1_type );

        auto const& baryOnRefFaces = gm1->referenceConvex().barycenterFaces();
        matrix_node_type baryOnRefFace( nRealDim, 1 );

        std::vector<std::map<uint16_type, matrix_node_type>> ctxPtsOnRefFaces( numTopologicalFaces );
        for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
        {
            ublas::column( baryOnRefFace, 0 ) = ublas::column( baryOnRefFaces, f );
            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                ctxPtsOnRefFaces[f][__p.value()] = baryOnRefFace;
        }
        auto pcf = gm1->preComputeOnFaces( gm1, ctxPtsOnRefFaces );
        static const size_type gmc_context_v = /*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN;
        auto ctx = gm1->template context<gmc_context_v>( *this, pcf, 0 );
        em_matrix_col_type<value_type> ns( _normals.data().begin(), _normals.size1(), _normals.size2() );
        for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
        {
            ctx->template update<gmc_context_v>( *this, f );
            ns.col( f ) = ctx->unitNormal( 0 );
            //ublas::column( _normals, f ) = ctx->unitNormal( 0 );
        }
        return _normals;
    }

    /**
     * \return the normal at the barycenter of the face \p f
     */
    node_type normal( uint16_type f ) const
    {
        if ( nDim != nRealDim )
        {
            CHECK( false ) << "normal when nDim != nRealDim is not implemented";
            return node_type( 0 );
        }

        gm1_ptrtype gm1 = this->gm1();
        if ( !gm1.use_count() )
            gm1 = gm1_ptrtype( new gm1_type );

        matrix_node_type baryOnFace( nRealDim, 1 );
        ublas::column( baryOnFace, 0 ) = gm1->referenceConvex().faceBarycenter( f );
        auto pcf = gm1->preComputeOnFaces( gm1, baryOnFace );
        static const size_type gmc_context_v =/*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN;
        auto ctx = gm1->template context<gmc_context_v>( *this, pcf, f );
        ctx->template update<gmc_context_v>( *this, f );

        node_type n( nRealDim );
        em_node_type<value_type> en( n.data().begin(), n.size() );
        en = ctx->unitNormal( 0 );
        return n;
    }

    /**
     * Get the local id of the point in the element
     *
     * @param _localFace local id of a face in the element
     * @param _point local id of a point in the face
     *
     * @return the local id of the point in the element
     */
    static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
    {
#if 1
        typedef typename mpl::if_<mpl::not_equal_to<mpl::int_<super::nDim>, mpl::int_<2>>,
                                  mpl::identity<super>,
                                  mpl::identity<tt>>::type the_type;
        return the_type::type::fToP( _localFace, _point );
#else
        return super::fToP( _localFace, _point );
#endif
    }

    /**
     * get the number of opposite points per face
     *
     * @return the number of opposite points per face
     */
    uint16_type nOppositePointsPerFace() const
    {
        return super::nbOppositePointsPerFace;
    }

    /**
     * faceToOppositePoint(i,j) = localId of jth opposite point to ith local face
     *
     * @return the localId of _point-th opposite point to _localFace-th local face
     */
    uint16_type faceToOppositePoint( uint16_type const _localFace, uint16_type const _point ) const
    {
        return super::faceToOppositePoint( _localFace, _point );
    }

    /**
     * Determines if the local numbering of a 2D element
     * is oriented anticlockwise
     */
    bool isAnticlockwiseOriented() const
    {

        // Calculate vectors originating from vertex zero
        ublas::matrix<T> orientation_matrix( nRealDim, nRealDim );

        for ( int i = 0; i < nRealDim; ++i )
        {
            matrix_node_type nodesG = this->G();
            ublas::row( orientation_matrix, i ) = ( ublas::column( nodesG, i + 1 ) -
                                                    ublas::column( nodesG, 0 ) );
        }

        LU<ublas::matrix<T>> lu( orientation_matrix );
        T sgn = lu.det();

        return ( sgn > 0 ) ? 1 : 0;
    }

    FEELPP_DEPRECATED
    void setPointCoordG( int i, ublas::vector<double> const& u )
    {
        // ublas::column( M_G, i ) = u;
    }
    void applyDisplacement( int i, ublas::vector<double> const& u )
    {
        // ublas::column( M_G, i ) += u;
        ( *M_points[i] ) += u;
    }
    FEELPP_DEPRECATED
    void applyDisplacementG( int i, ublas::vector<double> const& u )
    {
        // ublas::column( M_G, i ) += u;
    }
    /**
     * set the tags associated to the points
     * - tags[0] physical region
     * - tags[1] elementary region
     * - tags[2] particular region
     */
    void setTags( std::vector<int> const& tags )
    {
        M_markers[1].assign( tags[0] );
        if ( tags.size() > 1 )
            M_markers[2].assign( tags[1] );

        if ( tags.size() > 2 )
        {
            this->setProcessId( tags[3] );

            if ( tags[2] > 1 )
            {
                // ghosts
                std::vector<rank_type> p( tags[2] - 1 );

                for ( size_type i = 0; i < p.size(); ++i )
                {
                    p[i] = tags[4 + i];
                }

                this->setNeighborPartitionIds( p );
            }
        }
    }

    std::map<uint16_type, Marker1> const&
    markers() const
    {
        return M_markers;
    }
    void setMarkers( std::map<uint16_type, Marker1> const& markers )
    {
        M_markers = markers;
    }
    bool hasMarker( uint16_type k ) const
    {
        auto itFindMarker = M_markers.find( k );
        if ( itFindMarker == M_markers.end() )
            return false;
        if ( itFindMarker->second.isOff() )
            return false;
        return true;
    }
    Marker1 const& marker( uint16_type k ) const
    {
        DCHECK( this->hasMarker( k ) ) << "no marker type " << k;
        return M_markers.find( k )->second;
    }
    Marker1& marker( uint16_type k )
    {
        return M_markers[k];
    }
    void setMarker( uint16_type k, flag_type v )
    {
        M_markers[k].assign( v );
    }

    bool hasMarker() const
    {
        return this->hasMarker( 1 );
    }
    Marker1 const& marker() const
    {
        DCHECK( this->hasMarker( 1 ) ) << "no marker type 1";
        return M_markers.find( 1 )->second;
    }
    Marker1& marker()
    {
        return M_markers[1];
    }
    Marker1 markerOr( uint16_type v = 0 ) const
    {
        if ( hasMarker() )
            return M_markers.find( 1 )->second;
        else
            return Marker1{v};
    }
    Marker1 markerOr( uint16_type v = 0 )
    {
        if ( hasMarker() )
            return M_markers[1];
        else
            return Marker1{v};
    }
    void setMarker( flag_type v )
    {
        M_markers[1].assign( v );
    }

    bool hasMarker2() const
    {
        return this->hasMarker( 2 );
    }
    Marker1 const& marker2() const
    {
        DCHECK( this->hasMarker( 2 ) ) << "no marker type 2";
        return M_markers.find( 2 )->second;
    }
    Marker1& marker2()
    {
        return M_markers[2];
    }
    void setMarker2( flag_type v )
    {
        M_markers[2].assign( v );
    }

    bool hasMarker3() const
    {
        return this->hasMarker( 3 );
    }
    Marker1 const& marker3() const
    {
        DCHECK( this->hasMarker( 3 ) ) << "no marker type 3";
        return M_markers.find( 3 )->second;
    }
    Marker1& marker3()
    {
        return M_markers[3];
    }
    void setMarker3( flag_type v )
    {
        M_markers[3].assign( v );
    }

    //! \return the number of point element neighbors
    size_type numberOfPointElementNeighbors() const
    {
        return pointElementNeighborIds().size();
    }
    //! \return the set of ids of point element neighbors
    std::set<size_type> pointElementNeighborIds() const
    {
        std::set<size_type> n;
        for ( uint16_type __p = 0; __p < numPoints; ++__p )
        {
            std::for_each( M_points[__p]->elements().begin(),
                           M_points[__p]->elements().end(),
                           [&n]( auto const& e ) { n.insert( e.first ); } );
        }
        return n;
    }

    value_type measure() const override
        {
            if constexpr ( UseMeasuresStorage )
                return this->measureBIS();
            else
            {
                gm_ptrtype gm = this->gm();
                if ( !gm.use_count() )
                    gm = gm_ptrtype( new gm_type );

                quad_meas_type thequad( meas_quad_order );
                auto pc = gm->preCompute( gm, thequad.points() );
                auto ctx = gm->template context<vm::JACOBIAN>( *this, pc );
                return this->computeMeasureImpl( thequad, ctx );
            }
        }
    value_type faceMeasure( uint16_type f ) const
        {
            if constexpr ( UseMeasuresStorage )
                return this->faceMeasureBIS( f );
            else
            {
                gm_ptrtype gm = this->gm();
                if ( !gm.use_count() )
                    gm = gm_ptrtype( new gm_type );

                quad_meas_type thequad( meas_quad_order );
                auto pcf = gm->preComputeOnFaces( gm, thequad.allfpoints() );
                auto ctxf = gm->template context<vm::NORMAL | vm::KB | vm::JACOBIAN>( *this, pcf, f );
                return this->computeFaceMeasureImpl( thequad, ctxf, f );
            }
        }
    Eigen::Matrix<value_type,numTopologicalFaces,1> faceMeasures() const
        {
            Eigen::Matrix<value_type,numTopologicalFaces,1> m;
            for( int f = 0; f < numTopologicalFaces; ++f )
                m( f ) = faceMeasure( f );
            return m;
        }

    void setMeasurePointElementNeighbors( value_type meas )
    {
        if constexpr ( UseMeasuresStorage )
            this->setMeasurePointElementNeighborsBIS( meas );
    }

    /**
     * @return a quadrature to compute the measure of the element
     */
    static quad_meas_type imMeasure() { return quad_meas_type( meas_quad_order ); }

    /**
     * @return a quadrature to compute the measure of the order 1 approximation
     * of the element
     */
    static quad_meas1_type imMeasureOrder1() { return quad_meas1_type( 0 ); }

    //!
    //! check if an element is valid: check points are != nullptr
    //!
    bool isValid() const;

    void update();

    template <typename QuadType, typename CtxType, typename CtxFaceType>
    void updateWithCtx( QuadType const& thequad,
                        std::shared_ptr<CtxType>& ctx,
                        std::shared_ptr<CtxFaceType>& ctxf ) const;
    template <typename CtxType, typename CtxFaceType>
    void updateWithCtx1( quad_meas1_type const& thequad,
                         std::shared_ptr<CtxType>& ctx,
                         std::shared_ptr<CtxFaceType>& ctxf ) const;

  private:

    template <typename QuadType, typename CtxType>
        value_type computeMeasureImpl( QuadType const& thequad,
                                       std::shared_ptr<CtxType>& ctx ) const;
    template <typename QuadType, typename CtxType>
        value_type computeFaceMeasureImpl( QuadType const& thequad,
                                           std::shared_ptr<CtxType>& ctxf,
                                           uint16_type f ) const;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        DVLOG( 2 ) << "Serializing GeoND...\n";
        DVLOG( 2 ) << "  - base class...\n";
        ar& boost::serialization::base_object<super>( *this );
        DVLOG( 2 ) << "  - points...\n";
        //ar & M_points;
        for ( uint16_type i = 0; i < numPoints; ++i )
        {
            if ( Archive::is_loading::value )
            {
                // WARNING ! a memory leak is introduce here
                // we need to register in M_commonData for example : in the mesh if exist else in this struct
                if ( !M_points[i] )
                    M_points[i] = new point_type;
            }
            ar&*( M_points[i] );
        }
        // DVLOG(2) << "  - G...\n";
        // ar & M_G;
        DVLOG( 2 ) << "  - markers...\n";
        ar& M_markers;
    }

  private:
    /** geometric nodes of the element */
    std::vector<point_type*> M_points;

    //! store neighbor element id
    std::vector<size_type> M_neighbors;

    //! mapping from marker index to marker flag
    std::map<uint16_type, Marker1> M_markers;

    //! common data shared in a collection of multi-dimensional geometrical entity
    mutable GeoNDCommon<self_type> * M_commonData;
};

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline void
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::setPoint( uint16_type i, point_type& p )
{
    M_points[i] = std::addressof( p );
    DCHECK( M_points[i] != nullptr ) << "invalid point added to element at index " << i;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
std::ostream&
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::showMe( bool verbose, std::ostream& out ) const
{
    out << "----- BEGIN OF GeoND data ---" << std::endl
        << std::endl;
    out << " GeoND object of shape " << Shape << std::endl;
    out << " Number of Vertices = " << numVertices << std::endl;
    out << "   Number of Points = " << numPoints << std::endl;
    out << "                 id = " << this->id() << std::endl;
    out << "                  G = " << this->G() /*M_G*/ << "\n";

    for ( int i = 0; i < numVertices; i++ )
    {
        out << "POINT id = " << i << std::endl;
        point( i ).showMe( verbose, out );
    }

    out << "----- END OF GeoND data ---" << std::endl
        << std::endl;
    return out;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::swapPoints( const uint16_type& pt1, const uint16_type& pt2 )
{
    point_type* tmp( M_points[pt1] );
    M_points[pt1] = M_points[pt2];
    M_points[pt2] = tmp;

    // swap also the entries in G
    // ublas::column( M_G, pt1 ).swap( ublas::column( M_G, pt2 ) );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::exchangePoints( const uint16_type otn[numPoints] )
{
    point_type* tmp[numPoints];

    for ( unsigned int i = 0; i < numPoints; ++i )
    {
        tmp[i] = M_points[i];
    }

    for ( unsigned int i = 0; i < numPoints; ++i )
    {
        M_points[i] = tmp[otn[i]];
        // ublas::column( M_G, i ) = M_points[i]->node();
    }
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
bool GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::isValid() const
{
    auto it = std::find( M_points.begin(), M_points.end(), nullptr );
    DCHECK( it == M_points.end() ) << "element with invalid point : " << this->id();
    return it == M_points.end();
}
template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::update()
{
    if constexpr ( !UseMeasuresStorage )
       return;

    gm_ptrtype gm = this->gm();
    if ( !gm.use_count() )
        gm = gm_ptrtype( new gm_type );

    quad_meas_type thequad( meas_quad_order );
    auto pc = gm->preCompute( gm, thequad.points() );
    auto pcf = gm->preComputeOnFaces( gm, thequad.allfpoints() );

    auto ctx = gm->template context<vm::JACOBIAN>( *this, pc );
    auto ctxf = gm->template context</*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN>( *this, pcf, 0 );
    this->updateWithCtx( thequad, ctx, ctxf );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename QuadType, typename CtxType>
typename GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::value_type
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::computeMeasureImpl( QuadType const& thequad,
                                                                                    std::shared_ptr<CtxType>& ctx ) const
{
    value_type meas = 0.;
    if constexpr ( CtxType::is_linear )
    {
        meas = thequad.weightsSum() * ctx->J( 0 );
    }
    else
    {
        for ( int q = 0; q < thequad.nPoints(); ++q )
            meas += thequad.weight( q ) * ctx->J( q );
    }
    return meas;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename QuadType, typename CtxFaceType>
typename GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::value_type
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::computeFaceMeasureImpl( QuadType const& thequad,
                                                                                        std::shared_ptr<CtxFaceType>& ctxf,
                                                                                        uint16_type f ) const
{
    value_type measFace = 0.;
    if constexpr ( nDim == nRealDim )
    {
        for ( int q = 0; q < thequad.nPointsOnFace( f ); ++q )
            measFace += thequad.weight( f, q ) * ctxf->J( q ) * ctxf->normalNorm( q );
    }
    return measFace;
}


template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename QuadType, typename CtxType, typename CtxFaceType>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::updateWithCtx( QuadType const& thequad,
                                                        std::shared_ptr<CtxType>& ctx,
                                                        std::shared_ptr<CtxFaceType>& ctxf ) const
{
    if constexpr ( UseMeasuresStorage )
    {
        this->setMeasure( this->computeMeasureImpl( thequad,ctx ) );

        if constexpr ( nDim == nRealDim )
        {
            for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
            {
                ctxf->template update</*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN>( dynamic_cast<typename CtxFaceType::element_type const&>( *this ), f );
                this->setFaceMeasure( f, computeFaceMeasureImpl( thequad,ctxf, f ) );
            }
        }
    }
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename CtxType, typename CtxFaceType>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::updateWithCtx1( quad_meas1_type const& thequad,
                                                         std::shared_ptr<CtxType>& ctx,
                                                         std::shared_ptr<CtxFaceType>& ctxf ) const
{
    updateWithCtx( thequad, ctx, ctxf );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline FEELPP_DEPRECATED
    DebugStream&
operator<<( DebugStream& __os, GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE,UseMeasuresStorage> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.showMe( true, __str );

        __os << __str.str() << "\n";
    }

    return __os;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline NdebugStream& FEELPP_DEPRECATED
operator<<( NdebugStream& __os, GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage> const& __n )
{
    return __os;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline std::ostream&
operator<<( std::ostream& __os, GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE,UseMeasuresStorage> const& __n )
{
    return __n.showMe( true, __os );
}

} // namespace Feel
#endif
