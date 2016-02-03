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

#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/geo0d.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelpoly/im.hpp>

namespace Feel
{
class MeshBase;

template<int Dim, int Order, int RealDim, template<uint16_type,uint16_type,uint16_type> class Entity, typename T> struct GT_Lagrange;

template<class Convex, uint16_type O, typename T2> class Gauss;
template<int IMORDER,
         int DIM,
         template<uint16_type, uint16_type, uint16_type> class Entity,
         template<class Convex, uint16_type O, typename T2> class QPS,
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
    uint16_type operate( uint16_type const & point )
    {
        return ( point < GeoShape::numVertices ) ?
               GeoShape::numVertices - point :
               GeoShape::numPoints - point + GeoShape::numVertices;
    }

};
}
/// \endcond


/**
 * @class GeoND
 * @brief Base class for Multi-dimensional basis Geometrical Entities.
 *
 */
template <uint16_type Dim,
         typename GEOSHAPE,
         typename T = double,
         typename POINTTYPE = Geo0D<Dim, T> >
class GeoND
    :
public GeoEntity<GEOSHAPE>
{
    typedef GeoEntity<GEOSHAPE> super;
public:

    typedef T value_type;
    typedef GeoND<Dim, GEOSHAPE, value_type, POINTTYPE> self_type;
    typedef self_type element_type;

    typedef GEOSHAPE GeoShape;
    typedef POINTTYPE PointType;

    typedef PointType point_type;
    typedef typename super::face_type face_type;

    static const size_type Shape = super::Shape;
    static const uint16_type numPoints = super::numPoints;
    static const uint16_type numVertices = super::numVertices;
    static const uint16_type numLocalPoints = super::numPoints;
    static const uint16_type numLocalEdges = super::numEdges;
    static const uint16_type numLocalVertices = super::numVertices;
    static const int numFaces = super::numFaces;
    static const int numEdges = super::numEdges;
    static const int numTopologicalFaces = super::numTopologicalFaces;
    static const uint16_type numNeighbors = super::numTopologicalFaces;


    typedef typename ublas::bounded_array<point_type*, numPoints>::iterator point_iterator;
    typedef typename ublas::bounded_array<point_type*, numPoints>::const_iterator point_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename node<value_type>::type node_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    template<int GmOrder>
    struct GetGm
    {

        typedef typename mpl::if_<mpl::bool_<GeoShape::is_hypercube>,
                mpl::identity<GT_Lagrange<nDim, GmOrder, nRealDim, Hypercube, T> >,
                mpl::identity<GT_Lagrange<nDim, GmOrder, nRealDim, Simplex, T> > >::type::type type;
        typedef boost::shared_ptr<type> ptrtype;
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
            mpl::int_<1> >,
            mpl::identity<vertex_permutation_type>,
            typename mpl::if_<mpl::equal_to<mpl::int_<nDim>,
            mpl::int_<2> >,
            mpl::identity<edge_permutation_type>,
            mpl::identity<face_permutation_type> >::type>::type::type permutation_type;

    template<int GeoOrder>
    struct GetImMeasure
    {
        // quadrature formula used in entity measure (for ho geo, need to check)
        static const uint16_type quad_order = (nOrder-1)*nDim;
        typedef typename mpl::if_<mpl::bool_<GeoShape::is_hypercube>,
                                  mpl::identity<typename IMGeneric<quad_order,Dim,Hypercube,Gauss,double>::type >,
                                  mpl::identity<typename IMGeneric<quad_order,Dim,Simplex,Gauss,double>::type > >::type::type type;
    };
    typedef typename GetImMeasure<nOrder>::type quad_meas_type;
    typedef typename GetImMeasure<1>::type quad_meas1_type;

    /**
     * default constructor
     */
    GeoND()
        :
        super( 0 ),
        M_points( numPoints ),
        M_G( nRealDim, numPoints ),
        M_neighbors( 0 ),
        M_marker1(),
        M_marker2(),
        M_marker3(),
        M_gm(),
        M_gm1()
    {
    }

    /**
     * constructor from an id
     *
     * @param id identifier for the element to store
     *
     */
    explicit GeoND( size_type id )
        :
        super( id ),
        M_points( numPoints ),
        M_G( nRealDim, numPoints ),
        M_neighbors( 0 ),
        M_marker1(),
        M_marker2(),
        M_marker3(),
        M_gm(),
        M_gm1()
    {
    }

    GeoND( GeoND const& e ) = default;
    GeoND( GeoND && e ) = default;

    GeoND& operator=( GeoND const& ) = default;
    GeoND& operator=( GeoND && ) = default;

    /**
     * destructor, make it virtual for derived classes
     */
    ~GeoND()
    {
    }

    /**
     * set the mesh to which this geometric entity belongs to
     */
    void setMeshAndGm( MeshBase const* m, gm_ptrtype const& gm, gm1_ptrtype const& gm1 ) const
    {
        M_mesh = m;
        M_gm = gm;
        M_gm1 = gm1;
    }

    void setMesh( MeshBase const* m ) const
    {
        M_mesh = m;
    }

    //! return the geometric mapping if a mesh was set
    gm_ptrtype gm() const
    {
        return M_gm;
    }

    //! return the geometric mapping if a mesh was set
    gm1_ptrtype gm1() const
    {
        return M_gm1;
    }

    /**
     * \return the mesh to which this geometric entity belongs to
     */
    MeshBase const* mesh() const
    {
        return M_mesh;
    }

    /**
     * \return true if points have been inserted in elements, false
     * otherwise
     */
    bool hasPoints() const
    {
        for ( int i = 0; i < numPoints; ++i )
            if ( M_points[ i ] == nullptr )
                return false;
        return numPoints > 0;
    }
    //! @return true if the element has at least a point with marker1 active
    bool hasPointWithMarker() const
        {
            bool pt_with_marker = false;
            for ( int i = 0; i < numPoints; ++i )
            {
                pt_with_marker = this->point( i ).marker().isOn();
                if ( pt_with_marker )
                    break;
            }
            return pt_with_marker;
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
        return M_neighbors.size();//numNeighbors;
    }

    /**
     * Neighbors are stored as pointers and some of them can be null
     * pointers if the corresponding face is on the boundary.
     *
     * \return the pair neighbor \p n index and process \p id it belongs to
     */
    std::pair<size_type,rank_type> const& neighbor( uint16_type n ) const
    {
        return M_neighbors[n];
    }

    /**
     * set the \p n -th neighbor with \p neigh
     */
    void setNeighbor( uint16_type n, size_type neigh_id, rank_type proc_id )
    {
        if ( M_neighbors.empty() )
        {
            M_neighbors.reserve( numNeighbors );
            M_neighbors.resize( numNeighbors, std::make_pair( invalid_size_type_value, invalid_rank_type_value ) );
        }
        M_neighbors[n] = std::make_pair( neigh_id, proc_id );
    }

    bool isNeighbor( self_type const& G ) const
    {
        for ( uint16_type i = 0; i< this->nNeighbors() ; ++i )
            if ( this->neighbor( i ).first==G.id() ) return true;

        return false;
    }

    /**
     * \return the barycenter of the element
     */
    node_type barycenter() const
    {
        auto M = glas::average( M_G );
        return ublas::column( M, 0 );
    }

    /**
     * \return the barycenter at the faces of the element
     */
    node_type faceBarycenter( uint16_type f ) const
    {
        constexpr int nPtsInFace = GEOSHAPE::topological_face_type::numPoints;
        em_matrix_col_type G( const_cast<double*>(M_G.data().begin()), M_G.size1(), M_G.size2() );

        node_type n( nRealDim );
        em_node_type en( n.data().begin(), n.size() );
        en.setZero();
        for ( uint16_type p =  0;  p < nPtsInFace; ++p )
        {
            // get pt id in element  from local pt id  in face
            int ptid = this->fToP( f, p );
            en += G.col(ptid);
        }
        return en / nPtsInFace;
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

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    PointType & point( uint16_type i )
    {
        return *( static_cast<POINTTYPE *>( M_points[ i ] ) );
    }


    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    PointType const & point ( uint16_type i ) const
    {
        return *( static_cast<POINTTYPE *>( M_points[ i ] ) );
    }
    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    PointType* pointPtr( uint16_type i )
    {
        return M_points[ i ];
    }

    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    PointType const* pointPtr ( uint16_type i ) const
    {
        return M_points[ i ];
    }

    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    PointType & reversepoint( uint16_type const i )
    {
        return *( static_cast<POINTTYPE *>( M_points[ Feel::detail::ReversePoint<GEOSHAPE>::operate( i ) ] ) );
    }


    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    PointType const & reversepoint ( uint16_type const i ) const
    {
        return *( static_cast<POINTTYPE *>( M_points[ Feel::detail::ReversePoint<GEOSHAPE>::operate( i ) ] ) );
    }


    /**
     * Inserts a point.  Uses point references
     * put point
     */
    void setPoint( uint16_type const i, point_type const & p );

    /**
     * show information about the geoND
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c output stream
     *
     * @return the output stream
     */
    std::ostream & showMe( bool verbose = false, std::ostream & c = std::cout ) const;

    /**
     * Swap Points
     *
     *  This is a member function to be used ONLY by routines for checking or
     *  amending meshes. You must give the local id
     *
     * @param pt1 1st point to swap with 2nd point
     * @param pt2 2nd point to swap with 1st point
     */
    void swapPoints( const uint16_type & pt1, const uint16_type & pt2 );

    /**
     *  Exchange Points
     *
     *  Exchanges points according to a list of old2new local id numbering !
     *  old2new[i] is the new local id of a point whose old local id was ! i+1
     *  (remeber the numbering from 1 of the id's!. This is a member function
     *  to be used ONLY by routines for checking or amending meshes. You must
     *  give uint16_type (which start
     */
    void exchangePoints( const uint16_type otn[ numPoints ] );

    /**
     * matrix of geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of geometric nodes
     */
    matrix_node_type const& G() const
    {
        return M_G;
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
    matrix_node_type  vertices() const
    {
        return ublas::subrange( M_G, 0, nRealDim, 0, numVertices );
    }

    /**
     * matrix of geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of geometric nodes
     */
    matrix_node_type & G()
    {
        return M_G;
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
    double h() const
    {
        em_matrix_col_type G( const_cast<double*>(M_G.data().begin()), M_G.size1(), M_G.size2() );
        double res = 0.;
        for ( uint16_type __e = 0; __e < numLocalEdges; ++__e )
        {
            int col1 = this->eToP( __e, 0 );
            int col2 = this->eToP( __e, 1 );
            double r = (G.col(col1)-G.col(col2)).norm();
            res = ( res > r )?res:r;
        }
        return res;
    }
    /**
     * @brief get the minimum edge length in the element
     * @return the minimum edge length in the element
     */
    double hMin() const
    {
        em_matrix_col_type G( const_cast<double*>(M_G.data().begin()), M_G.size1(), M_G.size2() );

        double res = 1e10;
        for ( uint16_type __e = 0; __e < numLocalEdges; ++__e )
        {
            int col1 = this->eToP( __e, 0 );
            int col2 = this->eToP( __e, 1 );
            double r = (G.col(col1)-G.col(col2)).norm();
            res = ( res > r )?r:res;
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
    double hFace( uint16_type f ) const
    {
        if ( nRealDim==1 )
            return 1;

        constexpr int nEdges = GEOSHAPE::topological_face_type::numEdges;
        em_matrix_col_type G( const_cast<double*>(M_G.data().begin()), M_G.size1(), M_G.size2() );

        double res = 0.;
        for ( uint16_type e =  0;  e < nEdges; ++e )
        {
            // get edge id in face from local edge id
            int edg = this->f2e( f, (nDim==2)?f:e );
            int col1 = this->eToP( edg, 0 );
            int col2 = this->eToP( edg, 1 );
            double r = (G.col(col1)-G.col(col2)).norm();
            res = ( res > r )?res:r;
        }
        return res;
    }

    double hEdge( uint16_type f ) const
    {
        int col1 = this->eToP( f, 0 );
        int col2 = this->eToP( f, 1 );
        em_matrix_col_type G( const_cast<double*>(M_G.data().begin()), M_G.size1(), M_G.size2() );
        return (G.col(col1)-G.col(col2)).norm();
    }

    struct tt
    {
        static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
        {
            return super::eToP( _localFace, _point );
        }
    };

    /**
     * \return the measure of the element
     */
    double measure() const
    {
        //return M_measure;
        auto itFindMeasure = M_measures.find( GEOND_MEASURES::MEAS_ELEMENT );
        if ( itFindMeasure != M_measures.end() )
            return M_measures.find( GEOND_MEASURES::MEAS_ELEMENT )->second[0];
        else
        {
            LOG(WARNING) << "element measure in mesh not updated : return 0 value";
            return 0.;
        }
    }

    /**
     * \return the measure of the element face \p f
     */
    double faceMeasure( uint16_type f ) const
    {
        //return M_measurefaces[f];
        auto itFindMeasure = M_measures.find( GEOND_MEASURES::MEAS_FACES );
        if ( itFindMeasure != M_measures.end() )
            return M_measures.find( GEOND_MEASURES::MEAS_FACES )->second[f];
        else
        {
            LOG(WARNING) << "faces measure in mesh not updated : return 0 value";
            return 0.;
        }
    }

    /**
     * \return the measure of the element faces
     */
    std::vector<double> const& faceMeasures() const
    {
        //return M_measurefaces;
        CHECK( M_measures.find( GEOND_MEASURES::MEAS_FACES ) != M_measures.end() ) << "FACE_MEASURES is not computed";
        return M_measures.find( GEOND_MEASURES::MEAS_FACES )->second;
    }

    /**
     * \return the normals at the barycenter of the faces
     */
    matrix_node_type normals() const
    {
        matrix_node_type _normals( nRealDim, numTopologicalFaces);
        if  ( nDim != nRealDim )
        {
            CHECK( false ) << "normal when nDim != nRealDim is not implemented";
            return _normals;
        }

        if ( !M_gm1.use_count() )
            M_gm1 = gm1_ptrtype( new gm1_type );

        auto const& baryOnRefFaces = M_gm1->referenceConvex().barycenterFaces();
        matrix_node_type baryOnRefFace( nRealDim, 1 );

        std::vector<std::map<uint16_type,matrix_node_type > > ctxPtsOnRefFaces( numTopologicalFaces );
        for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
        {
            ublas::column( baryOnRefFace, 0 ) = ublas::column( baryOnRefFaces,f );
            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                ctxPtsOnRefFaces[f][__p.value()] = baryOnRefFace;
        }
        auto pcf =  M_gm1->preComputeOnFaces( M_gm1, ctxPtsOnRefFaces );
        auto ctx = M_gm1->template context</*vm::POINT|*/vm::NORMAL|vm::KB|vm::JACOBIAN>( *this, pcf, 0 );
        for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
        {
            ctx->update( *this, f );
            ublas::column( _normals, f ) = ctx->unitNormal( 0 );
        }
        return _normals;
    }

    /**
     * \return the normal at the barycenter of the face \p f
     */
    node_type normal( uint16_type f ) const
    {
        if  ( nDim != nRealDim )
        {
            CHECK( false ) << "normal when nDim != nRealDim is not implemented";
            return node_type(0);
        }

        if ( !M_gm1.use_count() )
            M_gm1 = gm1_ptrtype( new gm1_type );

        matrix_node_type baryOnFace( nRealDim, 1 );
        ublas::column( baryOnFace, 0 ) = M_gm1->referenceConvex().faceBarycenter(f);
        auto pcf =  M_gm1->preComputeOnFaces( M_gm1, baryOnFace );
        auto ctx = M_gm1->template context</*vm::POINT|*/vm::NORMAL|vm::KB|vm::JACOBIAN>( *this, pcf, f );
        ctx->update( *this, f );
        return ctx->unitNormal( 0 );
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
        typedef typename mpl::if_<mpl::not_equal_to<mpl::int_<super::nDim>, mpl::int_<2> >,
                mpl::identity<super>,
                mpl::identity<tt> >::type the_type;
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
        ublas::matrix<T> orientation_matrix ( nRealDim,nRealDim );

        for ( int i = 0; i < nRealDim ; ++i )
        {
            ublas::row( orientation_matrix, i ) = ( ublas::column( this->G(), i+1 ) -
                                                    ublas::column( this->G(),   0 ) );

        }

        LU< ublas::matrix<T> > lu( orientation_matrix );
        T sgn=lu.det();

        return ( sgn > 0 ) ? 1 : 0;
    }

    void setPointCoordG( int i, ublas::vector<double> const& u )
    {
        ublas::column( M_G, i ) = u;
    }
    void applyDisplacement( int i, ublas::vector<double> const& u )
    {
        ublas::column( M_G, i ) += u;
        ( *M_points[ i ] ) += u;
    }
    void applyDisplacementG( int i, ublas::vector<double> const& u )
    {
        ublas::column( M_G, i ) += u;
    }
    /**
     * set the tags associated to the points
     * - tags[0] physical region
     * - tags[1] elementary region
     * - tags[2] particular region
     */
    void setTags( std::vector<int> const& tags )
    {
        M_marker1.assign( tags[0] );

        if ( tags.size() > 1 )
            M_marker2.assign( tags[1] );

        if ( tags.size() > 2 )
        {
            this->setProcessId( tags[3] );

            if ( tags[2] > 1 )
            {
                // ghosts
                std::vector<rank_type> p( tags[2]-1 );

                for ( size_type i = 0; i < p.size(); ++i )
                {
                    p[i] = tags[4+i];
                }

                this->setNeighborPartitionIds( p );
            }

        }
    }
    Marker1 const& marker() const
    {
        return M_marker1;
    }
    Marker1& marker()
    {
        return M_marker1;
    }
    void setMarker( flag_type v )
    {
        return M_marker1.assign( v );
    }

    Marker2 const& marker2() const
    {
        return M_marker2;
    }
    Marker2& marker2()
    {
        return M_marker2;
    }
    void setMarker2( flag_type v )
    {
        return M_marker2.assign( v );
    }

    Marker3 const& marker3() const
    {
        return M_marker3;
    }
    Marker3& marker3()
    {
        return M_marker3;
    }
    void setMarker3( flag_type v )
    {
        return M_marker3.assign( v );
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
    //! set the measure of point element neighbors
    void setMeasurePointElementNeighbors( value_type meas )
    {
        if ( M_measures.find( GEOND_MEASURES::MEAS_NEIGHBORS_ELEMENT ) == M_measures.end() )
        {
            M_measures[GEOND_MEASURES::MEAS_NEIGHBORS_ELEMENT].reserve( 1 );
            M_measures[GEOND_MEASURES::MEAS_NEIGHBORS_ELEMENT].resize( 1 );
        }

        //M_meas_pneighbors = meas;
        M_measures[GEOND_MEASURES::MEAS_NEIGHBORS_ELEMENT][0] = meas;
    }
    //! \return the measure of point element neighbors
    value_type measurePointElementNeighbors() const
    {
        //return M_meas_pneighbors;
        CHECK( M_measures.find( GEOND_MEASURES::MEAS_NEIGHBORS_ELEMENT ) != M_measures.end() ) << "MEAS_NEIGHBORS_ELEMENT is not computed";
        return M_measures.find(GEOND_MEASURES::MEAS_NEIGHBORS_ELEMENT)->second[0];
    }

    void update();
    void updateWithPc( typename gm_type::precompute_ptrtype const& pc,
                       typename gm_type::faces_precompute_type & pcf,
                       quad_meas_type const& thequad );

    void updateWithPc1( typename gm1_type::precompute_ptrtype const& pc,
                        typename gm1_type::faces_precompute_type & pcf,
                        quad_meas1_type const& thequad );
private:

    template<typename GmType, typename QuadType>
    void updateMeasureImpl( boost::shared_ptr<GmType> gm, typename GmType::precompute_ptrtype const& pc,
                            QuadType const& thequad );
    template<typename GmType, typename QuadType>
    void updateMeasureFaceImpl( boost::shared_ptr<GmType> gm, typename GmType::faces_precompute_type & pcf,
                                QuadType const& thequad, mpl::bool_<true> );
    template<typename GmType, typename QuadType>
    void updateMeasureFaceImpl( boost::shared_ptr<GmType> gm, typename GmType::faces_precompute_type & pcf,
                                QuadType const& thequad, mpl::bool_<false> );

private:

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            DVLOG(2) << "Serializing GeoND...\n";
            DVLOG(2) << "  - base class...\n";
            ar & boost::serialization::base_object<super>( *this );
            DVLOG(2) << "  - points...\n";
            ar & M_points;
            DVLOG(2) << "  - G...\n";
            ar & M_G;
            DVLOG(2) << "  - marker1...\n";
            ar & M_marker1;
            DVLOG(2) << "  - marker1: " << M_marker1.value() << "...\n";
            DVLOG(2) << "  - marker2...\n";
            ar & M_marker2;
            DVLOG(2) << "  - marker2: " << M_marker2.value() << "...\n";
            DVLOG(2) << "  - marker3...\n";
            ar & M_marker3;
            DVLOG(2) << "  - marker3: " << M_marker3.value() << "...\n";
        }

private:
    /** geometric nodes of the element */
    std::vector<point_type*> M_points;

    /**< matrix of the geometric nodes */
    matrix_node_type M_G;

    enum GEOND_MEASURES
    {
        MEAS_ELEMENT          = 0,
        MEAS_FACES            = 1,
        MEAS_NEIGHBORS_ELEMENT= 2
    };

    std::map<GEOND_MEASURES,std::vector<value_type> > M_measures;

    //double M_measure;
    //std::vector<double> M_measurefaces;

    //! store neighbor element id
    std::vector<std::pair<size_type,rank_type> > M_neighbors;

    //! measure of the set of point element neighbors
    //value_type M_meas_pneighbors;

    Marker1 M_marker1;
    Marker2 M_marker2;
    Marker3 M_marker3;

    // mesh to which the geond element belongs to
    mutable MeshBase const* M_mesh;
    mutable gm_ptrtype M_gm;
    mutable gm1_ptrtype M_gm1;
};

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
const uint16_type GeoND<Dim,GEOSHAPE, T, POINTTYPE>::numLocalPoints;

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
const uint16_type GeoND<Dim,GEOSHAPE, T, POINTTYPE>::numLocalVertices;

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::setPoint( uint16_type const i, point_type const & p )
{
    M_points[ i ] = const_cast<point_type *>( &p );
    //VLOG(1) << "[setPoint] üpdate point index " << i << " with "<< M_points[i]->id() << "\n";
    FEELPP_ASSERT( const_cast<point_type *>( &p ) != 0 ).error( "invalid Geo0D<>" );
    DCHECK( M_G.size1() == M_points[i]->node().size()) << "Invalid dimension " << M_G.size1() << "  vs "  << M_points[i]->node().size()
                                                       << " n=" << M_points[i]->node();
    ublas::column( M_G, i ) = M_points[i]->node();
}


template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
std::ostream &
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::showMe( bool verbose, std::ostream & out ) const
{
    out << "----- BEGIN OF GeoND data ---" << std::endl << std::endl;
    out << " GeoND object of shape " << Shape << std::endl;
    out << " Number of Vertices = " << numVertices << std::endl;
    out << "   Number of Points = " << numPoints << std::endl;
    out << "                 id = " << this->id() << std::endl;
    out << "                  G = " << M_G << "\n";

    for ( int i = 0; i < numVertices; i++ )
    {
        out << "POINT id = " << i << std::endl;
        point( i ).showMe( verbose, out );
    }

    out << "----- END OF GeoND data ---" << std::endl << std::endl;
    return out;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void GeoND<Dim,GEOSHAPE, T, POINTTYPE>::swapPoints( const uint16_type & pt1, const uint16_type & pt2 )
{
    point_type * tmp( M_points[ pt1 ] );
    M_points[ pt1 ] = M_points[ pt2 ];
    M_points[ pt2 ] = tmp;

    // swap also the entries in G
    ublas::column( M_G, pt1 ).swap( ublas::column( M_G, pt2 ) );
}


template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void GeoND<Dim,GEOSHAPE, T, POINTTYPE>::exchangePoints( const uint16_type otn[ numPoints ] )
{
    point_type * tmp[ numPoints ];

    for ( unsigned int i = 0; i < numPoints; ++i )
    {
        tmp[ i ] = M_points[ i ];
    }

    for ( unsigned int i = 0; i < numPoints; ++i )
    {
        M_points[ i ] = tmp[ otn[ i ] ];
        ublas::column( M_G, i ) = M_points[i]->node();
    }
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::update()
{
    if ( !M_gm.use_count() )
        M_gm = gm_ptrtype( new gm_type );

    if ( !M_gm1.use_count() )
        M_gm1 = gm1_ptrtype( new gm1_type );

    quad_meas_type thequad;
    auto pc = M_gm->preCompute( M_gm, thequad.points() );
    auto pcf =  M_gm->preComputeOnFaces( M_gm, thequad.allfpoints() );

    updateWithPc( pc, pcf, thequad );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::updateWithPc( typename gm_type::precompute_ptrtype const& pc,
                                                 typename gm_type::faces_precompute_type& pcf,
                                                 quad_meas_type const& thequad )
{
    updateMeasureImpl( M_gm, pc, thequad );
    updateMeasureFaceImpl( M_gm, pcf, thequad, typename mpl::equal_to<mpl::int_<nDim>, mpl::int_<nRealDim> >::type() );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::updateWithPc1( typename gm1_type::precompute_ptrtype const& pc,
                                                  typename gm1_type::faces_precompute_type & pcf,
                                                  quad_meas1_type const& thequad )
{
    updateMeasureImpl( M_gm1, pc, thequad );
    updateMeasureFaceImpl( M_gm1, pcf, thequad, typename mpl::equal_to<mpl::int_<nDim>, mpl::int_<nRealDim> >::type() );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
template<typename GmType, typename QuadType>
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::updateMeasureImpl( boost::shared_ptr<GmType> gm,
                                                      typename GmType::precompute_ptrtype const& pc,
                                                      QuadType const& thequad )
{
    if ( M_measures.find( GEOND_MEASURES::MEAS_ELEMENT ) == M_measures.end() )
    {
        M_measures[GEOND_MEASURES::MEAS_ELEMENT].reserve( 1 );
        M_measures[GEOND_MEASURES::MEAS_ELEMENT].resize( 1 );
    }

    auto ctx = gm->template context<vm::JACOBIAN>( *this, pc );
    double meas = 0.;
    for ( int q=0 ; q < thequad.nPoints() ; ++q )
        meas += thequad.weight(q)*ctx->J( q );
    M_measures[GEOND_MEASURES::MEAS_ELEMENT][0] = meas;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
template<typename GmType, typename QuadType>
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::updateMeasureFaceImpl( boost::shared_ptr<GmType> gm,
                                                          typename GmType::faces_precompute_type& pcf,
                                                          QuadType const& thequad,
                                                          mpl::bool_<true> )
{
    if ( M_measures.find( GEOND_MEASURES::MEAS_FACES ) == M_measures.end() )
    {
        // must be resize before ctx because used into
        M_measures[GEOND_MEASURES::MEAS_FACES].reserve( numTopologicalFaces );
        M_measures[GEOND_MEASURES::MEAS_FACES].resize( numTopologicalFaces );
    }

    auto ctx = gm->template context</*vm::POINT|*/vm::NORMAL|vm::KB|vm::JACOBIAN>( *this,pcf,0 );
    for ( int f = 0; f < numTopologicalFaces; ++f )
    {
        ctx->update( *this, f );
        double meas = 0.;
        for ( int q=0 ; q < thequad.nPointsOnFace( f ) ; ++q )
            meas += thequad.weight(f,q)*ctx->J( q )*ctx->normalNorm( q );
        M_measures[GEOND_MEASURES::MEAS_FACES][f] = meas;
    }
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
template<typename GmType, typename QuadType>
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::updateMeasureFaceImpl( boost::shared_ptr<GmType> gm,
                                                          typename GmType::faces_precompute_type& pcf,
                                                          QuadType const& thequad,
                                                          mpl::bool_<false> )
{
    // need because M_measurefaces is used in Geomap::Context
    if ( M_measures.find( GEOND_MEASURES::MEAS_FACES ) == M_measures.end() )
    {
        M_measures[GEOND_MEASURES::MEAS_FACES].reserve( numTopologicalFaces );
        M_measures[GEOND_MEASURES::MEAS_FACES].resize( numTopologicalFaces );
    }

}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline FEELPP_DEPRECATED
DebugStream&
operator<<( DebugStream& __os, GeoND<Dim,GEOSHAPE, T, POINTTYPE> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.showMe( true, __str );

        __os << __str.str() << "\n";
    }

    return __os;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
NdebugStream& FEELPP_DEPRECATED
operator<<( NdebugStream& __os, GeoND<Dim,GEOSHAPE, T, POINTTYPE> const& __n )
{
    return __os;
}


template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
std::ostream&
operator<<( std::ostream& __os, GeoND<Dim,GEOSHAPE, T, POINTTYPE> const& __n )
{
    return __n.showMe( true, __os );
}

} // Feel
#endif
