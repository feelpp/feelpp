/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Gilles Steiner <gilles.steiner@epfl.ch>
 * SPDX-FileContributor: Goncalo Pena <goncalo.pena@epfl.ch>
 *
 * SPDX-FileCopyrightText: 2006-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file equispaced.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-09-26
*/
#ifndef __PointSetEquiSpaced_H
#define __PointSetEquiSpaced_H 1

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>

#include <feel/feelcore/visitor.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelpoly/order.hpp>

#include <stdexcept>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelalg/glas.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

/**
 * @brief Equispaced point set on a convex.
 *
 * @tparam Convex The convex type (Simplex or Hypercube).
 * @tparam Order Polynomial order. Use Dynamic (-1) for runtime order.
 * @tparam T Value type for coordinates.
 *
 * Supports both static and dynamic order following the Eigen-style pattern:
 * - Static order: PointSetEquiSpaced<Simplex<2,1,2>, 3, double>
 * - Dynamic order: PointSetEquiSpaced<Simplex<2,1,2>, Dynamic, double>
 */
template<class Convex, int Order, typename T>
class PointSetEquiSpaced :  public PointSet< Convex, T>, public OrderBase<Order>
{

public :

    typedef PointSet<Convex, T> super;
    using order_base_type = OrderBase<Order>;

    typedef T value_type;

    typedef typename super::return_type return_type;
    typedef typename super::self_type self_type;

    typedef typename super::node_type node_type;

    typedef typename super::nodes_type nodes_type;

    typedef typename matrix_node<value_type>::type points_type;

    static inline const uint32_type Dim = Convex::nDim;
    static inline const uint32_type convexOrder = Convex::nOrder;
    static inline const uint32_type topological_dimension = Convex::topological_dimension;
    static inline const uint32_type nRealDim = Convex::nRealDim;

    static const size_type Shape = Convex::Shape;

    static inline const bool is_simplex = Convex::is_simplex;
    static inline const bool is_hypercube = Convex::is_hypercube;
    using cardinality_type = std::conditional_t<
        is_simplex,
        ::Feel::detail::SimplexCardinality,
        ::Feel::detail::HypercubeCardinality>;

    //! @name Static/Dynamic Order Support
    //! @{
    static constexpr int order_spec = order_base_type::order_spec;
    static constexpr int static_order = order_spec;
    static constexpr bool is_order_static = order_base_type::is_order_static;
    static constexpr bool is_order_dynamic = order_base_type::is_order_dynamic;
    //! @}

    // For static order, use compile-time convex type; for dynamic, use order 1 as placeholder
    static constexpr int order_for_convex = is_order_static ? order_spec : 1;
    typedef typename std::conditional_t<is_simplex,
            Simplex<Dim, order_for_convex, Dim>,
            Hypercube<Dim, order_for_convex, Dim>> conv_order_type;

    typedef Reference<Convex, Dim, convexOrder, Dim/*nRealDim*/, value_type> RefElem;

    // Compile-time point counts (only valid for static order)
    static inline const uint32_type numPoints = is_order_static ? conv_order_type::numPoints : 0;
    static inline const uint32_type nbPtsPerVertex = is_order_static ? conv_order_type::nbPtsPerVertex : 0;
    static inline const uint32_type nbPtsPerEdge = is_order_static ? conv_order_type::nbPtsPerEdge : 0;
    static inline const uint32_type nbPtsPerFace = is_order_static ? conv_order_type::nbPtsPerFace : 0;
    static inline const uint32_type nbPtsPerVolume = is_order_static ? conv_order_type::nbPtsPerVolume : 0;

    typedef typename Convex::edge_to_point_t edge_to_point_t;
    typedef typename Convex::face_to_point_t face_to_point_t;
    typedef typename Convex::face_to_edge_t face_to_edge_t;

    typedef std::pair<uint16_type, uint16_type> range_type;
    typedef std::vector< std::vector<size_type> > index_map_type;

    RefElem RefConv;

    using order_base_type::order;
    using super::nPoints;

    /**
     * @brief Return the mathematical number of points attached to one vertex.
     */
    [[nodiscard]] constexpr uint32_type nPointsOnVertex() const noexcept
    {
        return cardinality_type::perVertex( Dim, this->order() );
    }

    /**
     * @brief Return the mathematical number of points attached to one edge.
     */
    [[nodiscard]] constexpr uint32_type nPointsOnEdge() const noexcept
    {
        return cardinality_type::perEdge( Dim, this->order() );
    }

    /**
     * @brief Return the mathematical number of points attached to one face.
     */
    [[nodiscard]] constexpr uint32_type nPointsOnFace() const noexcept
    {
        return cardinality_type::perFace( Dim, this->order() );
    }

    /**
     * @brief Return the mathematical number of points attached to one volume.
     */
    [[nodiscard]] constexpr uint32_type nPointsOnVolume() const noexcept
    {
        return cardinality_type::perVolume( Dim, this->order() );
    }

    //! @}

    /**
     * @brief Constructor for static order (Order >= 0).
     * @param interior If 1, only interior points are generated.
     */
    PointSetEquiSpaced( int interior = 0 )
        requires ( is_order_static )
        :
        super( 0, Dim ),
        order_base_type(),
        M_eid()
    {
        initialize( interior );
    }

    /**
     * @brief Constructor for dynamic order (Order == Dynamic).
     * @param order The runtime order.
     * @param interior If 1, only interior points are generated.
     */
    PointSetEquiSpaced( RuntimeOrder ro, int interior = 0 )
        requires ( is_order_dynamic )
        :
        super( 0, Dim ),  // Will be resized in body
        order_base_type(
            ::Feel::detail::validatedNodalOrder<cardinality_type>(
                Dim, ro ) ),
        M_eid()
    {
        initialize( interior );
    }

    ~PointSetEquiSpaced() override {}

    ublas::matrix_range<nodes_type const> pointsByEntity( uint16_type e ) const
    {
        FEELPP_ASSERT( M_eid[e].size() )( e ).error( "no points defined on this entity" );

        return ublas::project( this->points(),
                               ublas::range( 0,Dim ),
                               ublas::range( *M_eid[e].begin(), *M_eid[e].rbegin() + 1 ) );
    }

    std::pair<uint16_type, uint16_type> interiorRangeById( uint16_type e, uint16_type id ) const
    {
        uint16_type numEntities = 1;

        if ( e == 0 )
            numEntities = Convex::numVertices;

        else if ( e == 1 )
            numEntities = Convex::numEdges;

        else if ( e == 2 )
            numEntities = Convex::numFaces;

        uint16_type N = M_eid[e].size()/numEntities;

        return std::make_pair( *M_eid[e].begin() + id*N, *M_eid[e].begin() + ( id+1 )*N );
    }

    ublas::matrix_range<nodes_type const> interiorPointsById( uint16_type e, uint16_type id ) const
    {
        std::pair<uint16_type, uint16_type> position = interiorRangeById( e, id );

        ublas::matrix_range<nodes_type const> G = ublas::project( this->points(),
                ublas::range( 0, Dim ),
                ublas::range( position.first, position.second ) );

        return G;
    }

    uint32_type entityIds( int i, int j ) const
    {
        return M_eid[i][j];
    }

    uint32_type numEntities( int i ) const
    {
        return M_eid[i].size();
    }

    std::pair<uint16_type, uint16_type> const& pointToEntity( int p ) const
    {
        return M_pt_to_entity[p];
    }

    //Returns the local indices of all the subentities that compose the entity
    index_map_type entityToLocal ( uint16_type top_dim, uint16_type local_id, bool boundary = 0 ) const
    {
        index_map_type indices( top_dim+1 );

        if ( top_dim == 0 && boundary )
        {
            range_type pair = interiorRangeById( top_dim, local_id );

            indices[0].push_back( pair.first );
        }

        else
        {
            if ( !boundary && top_dim > 0 )
            {
                if ( numEntities( top_dim ) != 0 )
                    indices[top_dim].push_back( local_id );
            }

            else
            {
                if ( top_dim > 0 )
                {
                    //For the time being, this definition works, but in a more general framework,
                    //a entity shoulkd be created here with the respective top_dim
                    //in order to retrieve the number of vertices, edges, faces.

                    //number of vertices, number of edges and number of faces in the volume
                    std::map<uint16_type, uint16_type> numPointsInEntity;

                    numPointsInEntity[0] = ( uint16_type ) ( is_simplex )?( top_dim+1 ):( 2 + ( top_dim-1 )*( 2+3*( top_dim-2 ) ) );
                    numPointsInEntity[1] = ( uint16_type ) ( 2 + ( top_dim-1 )*( 1+2*is_simplex + ( 5-4*is_simplex )*( top_dim-1 ) ) )/2;
                    numPointsInEntity[2] = ( uint16_type ) 6-2*is_simplex;

                    for ( uint16_type i=0; i < numPointsInEntity[0] ; i++ )
                    {
                        if ( top_dim == 1 )
                            indices[0].push_back( RefConv.e2p( local_id, i ) );

                        else if ( top_dim == 2 )
                            indices[0].push_back( RefConv.f2p( local_id, i ) );
                    }

                    if ( ( top_dim == 2 ) && ( M_eid[1].size() != 0 ) )
                    {
                        for ( uint16_type i=0; i < numPointsInEntity[1] ; i++ )
                            indices[1].push_back( RefConv.f2e( local_id, i ) );
                    }

                    if ( top_dim == 3 )
                    {
                        for ( uint16_type k=0; k<numPointsInEntity.size(); k++ )
                            for ( uint16_type i=0; i < numPointsInEntity[k]; i++ )
                                indices[k].push_back( i );
                    }

                    if ( M_eid[top_dim].size() )
                        indices[top_dim].push_back( local_id );
                }
            }
        }

        return indices;
    }


    points_type pointsBySubEntity( uint16_type top_dim, uint16_type local_id, bool boundary = 0 ) const
    {
        index_map_type index_list = entityToLocal( top_dim, local_id, boundary );

        // Use runtime point counts for dynamic order support
        const auto pts_per_vertex = nPointsOnVertex();
        const auto pts_per_edge = nPointsOnEdge();
        const auto pts_per_face = nPointsOnFace();
        const auto pts_per_volume = nPointsOnVolume();

        uint16_type matrix_size = 0;

        if ( index_list[0].size() != 0 )
            matrix_size +=index_list[0].size()*pts_per_vertex;

        if ( ( top_dim >= 1 ) && ( index_list[1].size() != 0 ) )
            matrix_size +=index_list[1].size()*pts_per_edge;

        if ( ( top_dim >= 2 ) && ( index_list[2].size() != 0 ) )
            matrix_size +=index_list[2].size()*pts_per_face;

        if ( ( top_dim == 3 ) && ( index_list[3].size() != 0 ) )
            matrix_size +=pts_per_volume;

        points_type G ( Dim, matrix_size );

        for ( uint16_type i=0, p=0; i < top_dim+1; i++ )
        {
            if ( index_list[i].size() )
            {
                for ( uint16_type j=0; j < index_list[i].size(); j++ )
                {
                    points_type aux = interiorPointsById( i, index_list[i][j] );

                    ublas::subrange( G, 0, Dim, p, p+aux.size2() ) = aux;

                    p+=aux.size2();
                }
            }
        }

        return G;
    }

    index_map_type getEid ()
    {
        return M_eid;
    }

    std::vector<range_type> getPtE()
    {
        return M_pt_to_entity;
    }

    void setEid ( index_map_type eid )
    {
        M_eid = eid;
    }

    void setPtE ( std::vector<range_type> pt_ent )
    {
        M_pt_to_entity = pt_ent;
    }

    void addToEid ( uint16_type p, uint16_type q )
    {
        M_eid[p].push_back( q );
    }

    void addToPtE ( uint16_type p, range_type q )
    {
        M_pt_to_entity[p] = q;
    }

    FEELPP_DEFINE_VISITABLE();

private:

    /**
     * @brief Return the complete cell cardinality used while constructing a
     * full lattice.
     *
     * This is deliberately private because `nPoints()` reports the points
     * actually stored by the object, which may be an interior-only subset.
     */
    [[nodiscard]] constexpr uint32_type fullPointCount() const noexcept
    {
        return static_cast<uint32_type>(
            cardinality_type::total( Dim, this->order() ) );
    }

    /** @brief Build complete or interior-only storage for either order mode. */
    void initialize( int interior )
    {
        if ( interior != 0 && interior != 1 )
            throw std::invalid_argument(
                "PointSetEquiSpaced interior selector must be 0 or 1" );

        M_eid.assign( topological_dimension + 1, {} );

        if ( this->order() == 0 )
        {
            this->setPoints( glas::average( RefConv.vertices() ) );
            M_pt_to_entity.resize( this->nPoints() );
            addToEid( Dim, 0 );
            addToPtE( 0, std::make_pair( Dim, 0 ) );
        }
        else if ( interior == 1 )
        {
            this->setPoints( makePoints( Dim, 0 ) );
            M_pt_to_entity.resize( this->nPoints() );
            for ( uint16_type point = 0;
                  point < this->nPoints(); ++point )
            {
                addToEid( Dim, point );
                addToPtE( point, std::make_pair( Dim, 0 ) );
            }
        }
        else
        {
            const uint32_type pointCount = fullPointCount();
            M_pt_to_entity.resize( pointCount );
            nodes_type points( Dim, pointCount );

            for ( uint16_type entityDim = 0, point = 0;
                  entityDim < topological_dimension + 1;
                  ++entityDim )
            {
                for ( int entity = RefConv.entityRange( entityDim ).begin();
                      entity < RefConv.entityRange( entityDim ).end();
                      ++entity )
                {
                    nodes_type entityPoints(
                        makePoints( entityDim, entity ) );
                    if ( entityPoints.size2() == 0 )
                        continue;

                    ublas::subrange(
                        points, 0, Dim, point,
                        point + entityPoints.size2() ) = entityPoints;
                    for ( size_type localPoint = 0;
                          localPoint < entityPoints.size2();
                          ++localPoint )
                    {
                        addToEid(
                            entityDim, point + localPoint );
                        addToPtE(
                            point + localPoint,
                            std::make_pair( entityDim, entity ) );
                    }
                    point += entityPoints.size2();
                }
            }
            this->setPoints( points );
        }

        this->setName( "equispaced", this->order() );
    }

    index_map_type M_eid;
    std::vector<range_type> M_pt_to_entity;

    points_type makePoints( uint16_type topo_dim, uint16_type __id, int interior = 1 )
    {
        // vertices
        if ( topo_dim == 0 )
        {
            points_type G( RefConv.vertices().size1(), 1 );
            ublas::column( G, 0 ) = ublas::column( RefConv.vertices(), __id );
            return G;
        }

        // interior points of the convex
        else if ( topo_dim == topological_dimension )
        {
            if ( __id == 0 )
                return makeLattice<Shape>( interior );

            throw std::logic_error( "cannot make those points" );
            return points_type();
        }

        // all the other points
        else
        {
            points_type G;
            points_type Gret;

            if ( topo_dim == 1 )
            {
                G = makeLattice<SHAPE_LINE>( 1 );
                Gret.resize( nRealDim, G.size2() );

                if ( is_simplex )
                {
                    pt_to_entity_tetrahedron<Shape, 1> p_to_e( __id );

                    for ( size_type i = 0; i < G.size2(); ++i )
                        ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );
                }

                else
                {
                    pt_to_entity_hexahedron<Shape, 1> p_to_e( __id );

                    for ( size_type i = 0; i < G.size2(); ++i )
                        ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );
                }

                return Gret;
            }

            else if ( topo_dim == 2 )
            {
                if ( is_simplex )
                {
                    G = makeLattice<SHAPE_TRIANGLE>( 1 );
                    Gret.resize( nRealDim, G.size2() );
                    pt_to_entity_tetrahedron<Shape, 2> p_to_e( __id );

                    for ( size_type i = 0; i < G.size2(); ++i )
                        ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );
                }

                else
                {
                    G = makeLattice<SHAPE_QUAD>( 1 );
                    Gret.resize( nRealDim, G.size2() );
                    pt_to_entity_hexahedron<Shape, 2> p_to_e( __id );

                    for ( size_type i = 0; i < G.size2(); ++i )
                        ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );
                }

                return Gret;
            }
        }

        return points_type();
    }

    template<size_type shape>
    points_type makeLattice( uint16_type interior = 0 )
    {
        points_type G;
        const auto runtime_order = this->order();

        if ( runtime_order > 0 )
        {
            if ( shape == SHAPE_LINE )
                G = make_line_points( interior );

            else if ( shape == SHAPE_TRIANGLE )
                G = make_triangle_points( interior );

            else if ( shape == SHAPE_TETRA )
                G = make_tetrahedron_points( interior );

            else if ( shape == SHAPE_QUAD )
                return make_quad_points( interior );

            else if ( shape == SHAPE_HEXA )
                return make_hexa_points( interior );
        }

        else if ( runtime_order == 0 )
            G = glas::average( RefConv.vertices() );

        return G;
    }

    //---------------------------------------------------------------------------------------------
    // Point count methods - use semantic order accessor for both static and dynamic cases
    //---------------------------------------------------------------------------------------------
    int n_line_points( int interior = 0 ) const
    {
        const int runtime_order = static_cast<int>( this->order() );
        return std::max( 0, runtime_order+1-2*interior );
    }
    int n_triangle_points( int interior = 0 ) const
    {
        const int runtime_order = static_cast<int>( this->order() );
        if ( interior == 1 )
            return std::max( 0, ( runtime_order+1-2*interior )*( runtime_order-2*interior )/2 );

        return ( runtime_order+1 )*( runtime_order+2 )/2;
    }
    int n_tetrahedron_points( int interior = 0 ) const
    {
        const int runtime_order = static_cast<int>( this->order() );
        if ( interior == 1 )
            return std::max( 0, ( runtime_order+1-2*interior )*( runtime_order-2*interior )*( runtime_order-1-2*interior )/6 );

        return ( runtime_order+1 )*( runtime_order+2 )*( runtime_order+3 )/6;
    }

    int n_quad_points( int interior = 0 ) const
    {
        const int runtime_order = static_cast<int>( this->order() );
        if ( interior == 1 )
            return std::max( 0, ( runtime_order+1-2*interior )*( runtime_order+1-2*interior ) );

        return ( runtime_order+1 )*( runtime_order+1 );
    }

    int n_hexa_points( int interior = 0 ) const
    {
        const int runtime_order = static_cast<int>( this->order() );
        if ( interior == 1 )
            return std::max( 0, ( runtime_order+1-2*interior )*( runtime_order+1-2*interior )*( runtime_order+1-2*interior ) );

        return ( runtime_order+1 )*( runtime_order+1 )*( runtime_order+1 );
    }

    points_type
    make_line_points( int interior = 0 )
    {
        points_type G;
        const int runtime_order = static_cast<int>( this->order() );

        if ( runtime_order > 0 )
        {
            ublas::vector<node_type> h ( 1 );
            h( 0 ) = RefConv.vertex( 1 ) - RefConv.vertex( 0 );

            G.resize( Dim, n_line_points( interior ) );

            for ( int i = interior, indp = 0; i < runtime_order+1-interior; ++i, ++indp )
            {
                ublas::column( G, indp ) = RefConv.vertex( 0 ) + ( h( 0 ) * value_type( i ) )/value_type( runtime_order );
            }
        }

        else
            G = glas::average( RefConv.vertices() );

        return G;

    }

    points_type
    make_triangle_points( int interior = 0 )
    {
        points_type G;
        const int runtime_order = static_cast<int>( this->order() );

        if ( runtime_order > 0 )
        {
            ublas::vector<node_type> h ( 2 );
            h( 0 ) = RefConv.vertex( 1 ) - RefConv.vertex( 0 );
            h( 1 ) = RefConv.vertex( 2 ) - RefConv.vertex( 0 );

            G.resize( Dim, n_triangle_points( interior ) );

            for ( int i = interior, p = 0; i < runtime_order+1-interior; ++i )
            {
                for ( int j = interior; j < runtime_order + 1 - i-interior; ++j, ++p )
                {
                    ublas::column( G, p ) = RefConv.vertex( 0 ) + ( value_type( i ) * h( 1 )  +
                                            value_type( j ) * h( 0 ) )/ value_type( runtime_order );
                }
            }
        }

        else
            G = glas::average( RefConv.vertices() );

        return G;
    }

    points_type
    make_tetrahedron_points( int interior = 0 )
    {
        points_type G;
        const int runtime_order = static_cast<int>( this->order() );

        if ( runtime_order > 0 )
        {
            ublas::vector<node_type> h ( 3 );
            h( 0 ) = RefConv.vertex( 1 ) - RefConv.vertex( 0 );
            h( 1 ) = RefConv.vertex( 2 ) - RefConv.vertex( 0 );
            h( 2 ) = RefConv.vertex( 3 ) - RefConv.vertex( 0 );

            G.resize( Dim, n_tetrahedron_points( interior ) );

            for ( int i = interior, p = 0; i < runtime_order+1-interior; ++i )
            {
                for ( int j = interior; j < runtime_order + 1 - i - interior; ++j )
                {
                    for ( int k = interior; k < runtime_order + 1 - i - j - interior; ++k, ++p )
                    {
                        ublas::column( G, p ) = RefConv.vertex( 0 ) + ( value_type( i ) * h( 2 ) +
                                                value_type( j ) * h( 1 ) +
                                                value_type( k ) * h( 0 ) ) / value_type( runtime_order );

                    }
                }
            }
        }

        else
            G = glas::average( RefConv.vertices() );

        return G;
    }

    points_type
    make_quad_points( int interior = 0 )
    {
        const int runtime_order = static_cast<int>( this->order() );
        if ( runtime_order > 0 )
        {
            ublas::vector<node_type> h ( 2 );
            h( 0 ) = RefConv.vertex( 1 ) - RefConv.vertex( 0 );
            h( 1 ) = RefConv.vertex( 3 ) - RefConv.vertex( 0 );

            DVLOG(2) << "n quad pts = " << n_quad_points( interior ) << "\n";
            points_type G( Dim, n_quad_points( interior ) );

            for ( int i = interior, p = 0; i < runtime_order+1-interior; ++i )
            {
                for ( int j = interior; j < runtime_order + 1 -interior; ++j, ++p )
                {
                    ublas::column( G, p ) = RefConv.vertex( 0 ) + ( value_type( i ) * h( 0 )  +
                                            value_type( j ) * h( 1 ) )/ value_type( runtime_order );
                }
            }

            return G;
        }

        else
            return glas::average( RefConv.vertices() );
    }

    points_type
    make_hexa_points( int interior = 0 )
    {
        const int runtime_order = static_cast<int>( this->order() );
        if ( runtime_order > 0 )
        {
            ublas::vector<node_type> h ( 3 );
            h( 0 ) = RefConv.vertex( 1 ) - RefConv.vertex( 0 );
            h( 1 ) = RefConv.vertex( 3 ) - RefConv.vertex( 0 );
            h( 2 ) = RefConv.vertex( 4 ) - RefConv.vertex( 0 );

            points_type G( Dim, n_hexa_points( interior ) );
            DVLOG(2) << "n hexa pts = " << n_hexa_points( interior ) << "\n";

            for ( int i = interior, p = 0; i < runtime_order+1-interior; ++i )
            {
                for ( int j = interior; j < runtime_order + 1 - interior; ++j )
                {
                    for ( int k = interior; k < runtime_order + 1 - interior; ++k, ++p )
                    {
                        ublas::column( G, p ) = RefConv.vertex( 0 ) + ( value_type( i ) * h( 0 ) +
                                                value_type( j ) * h( 1 ) +
                                                value_type( k ) * h( 2 ) ) / value_type( runtime_order );

                    }
                }
            }

            return G;
        }

        else
            return glas::average( RefConv.vertices() );
    }

    template<size_type shape>
    struct pt_to_edge
    {
        pt_to_edge( std::vector<uint16_type> vert_ids )
            :
            h( 1.0 ),
            a( Entity<SHAPE_LINE, value_type>().vertex( 0 ) ),
            b( Entity<SHAPE_LINE, value_type>().vertex( 1 ) ),
            u( Entity<shape, value_type>().vertex( vert_ids[ 0 ] ) ),
            v( Entity<shape, value_type>().vertex( vert_ids[ 1 ] ) ),
            diff( v-u )
        {
            h = 1.0/( b[0]-a[0] );
        }
        node_type
        operator()( node_type const& x ) const
        {
            return u + h * ( x[ 0 ] - a[ 0 ] ) * diff;
        }
        value_type h;
        node_type a, b;
        node_type u, v, diff;
    };

    //
    // pt_to_face tetrahedron
    //
    template<size_type shape>
    struct pt_to_face_tetrahedron
    {
        pt_to_face_tetrahedron( std::vector<uint16_type> vert_ids )
            :
            u( Entity<shape, value_type>().vertex( vert_ids[ 0 ] ) ),
            v( Entity<shape, value_type>().vertex( vert_ids[ 1 ] ) ),
            w( Entity<shape, value_type>().vertex( vert_ids[ 2 ] ) ),
            diff( 2 )
        {
            diff[0] = v-u;
            diff[1] = w-u;
        }
        node_type
        operator()( node_type const& x ) const
        {
            return u + 0.5*( x[ 0 ]+1.0 ) * diff[ 0 ] + 0.5*( x[ 1 ]+1.0 ) * diff[ 1 ];
        }
        node_type u, v, w;
        ublas::vector<node_type> diff;
    };

    //
    // pt_to_face hexa
    //
    template<size_type shape>
    struct pt_to_face_hexahedron
    {
        pt_to_face_hexahedron( std::vector<uint16_type> vert_ids )
            :
            u( Entity<shape, value_type>().vertex( vert_ids[ 0 ] ) ),
            v( Entity<shape, value_type>().vertex( vert_ids[ 1 ] ) ),
            w( Entity<shape, value_type>().vertex( vert_ids[ 3 ] ) ),
            diff( 2 )
        {
            diff[0] = v-u;
            diff[1] = w-u;
        }
        node_type
        operator()( node_type const& x ) const
        {
            return u + 0.5*( x[ 0 ]+1.0 ) * diff[ 0 ] + 0.5*( x[ 1 ]+1.0 ) * diff[ 1 ];
        }
        node_type u, v, w;
        ublas::vector<node_type> diff;
    };

    template<size_type shape>
    struct pt_to_element
    {
        pt_to_element() {}
        pt_to_element( std::vector<uint16_type> const& ) {}
        node_type operator()( node_type const& x ) const
        {
            return x;
        }
    };

    template<size_type shape, uint16_type topo_dim>
    struct pt_to_entity_tetrahedron
    {
        typedef typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_LINE> >,
                mpl::identity<mpl::vector<boost::none_t,pt_to_edge<shape>,pt_to_edge<shape> > >,
                typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_TRIANGLE> >,
                mpl::identity<mpl::vector<boost::none_t, pt_to_edge<shape>, pt_to_element<shape> > >,
                mpl::identity<mpl::vector<boost::none_t, pt_to_edge<shape>, pt_to_face_tetrahedron<shape>, pt_to_element<shape> > >
                >::type // 2
                >::type::type _type;
        typedef typename mpl::at<_type, mpl::int_<topo_dim> >::type mapping_type;
        typedef mpl::vector<boost::none_t, edge_to_point_t, face_to_point_t> list_v;

        pt_to_entity_tetrahedron( uint16_type entity_id )
            :
            mapping( typename mpl::at<list_v, mpl::int_<topo_dim> >::type().entity( topo_dim, entity_id ) )
        {}

        node_type operator()( node_type const& x ) const
        {
            return mapping( x );
        }
        mapping_type mapping;
    };

    template<size_type shape,uint16_type topo_dim>
    struct pt_to_entity_hexahedron
    {
        typedef typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_LINE> >,
                mpl::identity<mpl::vector<boost::none_t,pt_to_edge<shape>,pt_to_edge<shape> > >,
                typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_QUAD> >,
                mpl::identity<mpl::vector<boost::none_t, pt_to_edge<shape>, pt_to_element<shape> > >,
                mpl::identity<mpl::vector<boost::none_t, pt_to_edge<shape>, pt_to_face_hexahedron<shape>, pt_to_element<shape> > >
                >::type // 2
                >::type::type _type;
        typedef typename mpl::at<_type, mpl::int_<topo_dim> >::type mapping_type;
        typedef mpl::vector<boost::none_t, edge_to_point_t, face_to_point_t> list_v;

        pt_to_entity_hexahedron( uint16_type entity_id )
            :
            mapping( typename mpl::at<list_v, mpl::int_<topo_dim> >::type().entity( topo_dim, entity_id ) )
        {}

        node_type operator()( node_type const& x ) const
        {
            return mapping( x );
        }
        mapping_type mapping;
    };

};

} // Feel
#endif /* __PointSetEquiSpaced_H */
