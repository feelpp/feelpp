/* -*- Mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-02-20

  Copyright (C) 2006 EPFL

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
   \file hypercube.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-02-20
 */
#ifndef __Hypercube_H
#define __Hypercube_H 1

#include <boost/detail/identifier.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelmesh/entities.hpp>
#include <feel/feelmesh/convex.hpp>
#include <feel/feelmesh/hypercubeordering.hpp>
#include <feel/feelpoly/order.hpp>

#include <variant>

namespace Feel
{
class HypercubeBase {};

/**
 * @class Hypercube
 * @brief Hypercube (tensor-product) convex of dimension \c Dim
 *
 * Supports both static order (default, backward compatible) and dynamic order.
 * For dynamic order, use Hypercube<Dim, Dynamic, RDim> and construct with RuntimeOrder(n).
 *
 * @tparam Dim Topological dimension
 * @tparam Order Polynomial order (or Dynamic = -1 for runtime order)
 * @tparam RDim Real/spatial dimension
 */
template<int Dim, int Order=1, int RDim = Dim>
class Hypercube  : public Convex<Dim, (Order >= 0 ? Order : 1), RDim>, HypercubeBase
{
public:
    //! @name Order type detection
    //! @{
    static constexpr bool is_order_static = (Order != Dynamic);
    static constexpr bool is_order_dynamic = !is_order_static;
    //! @}

private:
    typedef mpl::vector_c<size_type, SHAPE_POINT, SHAPE_LINE, SHAPE_QUAD, SHAPE_HEXA, SHAPE_SP4, SHAPE_SP5> shapes_t;
    typedef mpl::vector_c<size_type, GEOMETRY_POINT, GEOMETRY_LINE, GEOMETRY_SURFACE, GEOMETRY_VOLUME, GEOMETRY_4, GEOMETRY_5> geometries_t;

    typedef mpl::vector_c<size_type, 1, 2, 4, 8, 16, 32> vertices_t;
    typedef mpl::vector_c<size_type, 0, 1, 4, 12, 32> edges_t;
    typedef mpl::vector_c<size_type, 0, 2, 4, 6, 24> faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 1, 6, 24> geo_faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 4, 6, 24> normals_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 0, 1, 8> volumes_t;

    // Use effective order for mpl computations (fallback to 1 for dynamic case)
    static constexpr int EffectiveOrder = (Order >= 0) ? Order : 1;

    typedef mpl::vector_c<size_type, 1, EffectiveOrder+1, ( EffectiveOrder+1 )*( EffectiveOrder+1 ), ( EffectiveOrder+1 )*( EffectiveOrder+1 )*( EffectiveOrder+1 ) , ( EffectiveOrder+1 )*( EffectiveOrder+1 )*( EffectiveOrder+1 )*( EffectiveOrder+1 )> points_t;
    typedef mpl::vector_c<size_type, 0, EffectiveOrder+1-2, ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ), ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ) , ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 )> points_interior_t;
    typedef mpl::vector_c<size_type, 0, EffectiveOrder+1-2,                   EffectiveOrder+1-2,                         ( EffectiveOrder+1-2 ), ( EffectiveOrder+1-2 )> points_edge_t;
    typedef mpl::vector_c<size_type, 0,         0,     ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ),             ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ), ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ) > points_face_t;
    typedef mpl::vector_c<size_type, 0,         0,                           0, ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ), ( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 )*( EffectiveOrder+1-2 ) > points_volume_t;

    template<uint16_type rdim>
    struct faces_t
    {
        typedef mpl::vector<Hypercube<0, EffectiveOrder, rdim>,
                            Hypercube<0, EffectiveOrder, rdim>,
                            Hypercube<1, EffectiveOrder, rdim>,
                            Hypercube<2, EffectiveOrder, rdim> > type;
    };
    typedef mpl::vector<boost::none_t,Hypercube<1, EffectiveOrder,1>, Hypercube<1, EffectiveOrder, 2>, Hypercube<1, EffectiveOrder, 3>, boost::none_t > v_edges_t;
    typedef mpl::vector<Hypercube<1, EffectiveOrder>, Hypercube<2, EffectiveOrder>, Hypercube<3, EffectiveOrder>, Hypercube<4, EffectiveOrder>, boost::none_t > elements_t;

    typedef mpl::vector_c<uint16_type, 0, 1, 2, 8> permutations_t;

public:

    static inline const bool is_simplex = false;
    static inline const bool is_hypercube = true;

    static const size_type Shape = mpl::at<shapes_t, mpl::int_<Dim> >::type::value;
    static const size_type Geometry = mpl::at<geometries_t, mpl::int_<Dim> >::type::value;

    static inline const uint16_type nDim = Dim;
    //! Static order value (backward compatibility only).
    //! WARNING: for dynamic order types, this is a placeholder (0). Use order().
    static inline const uint16_type nOrder = is_order_static ? static_cast<uint16_type>(Order) : 0;
    //! Template order parameter value (may be Dynamic = -1)
    static constexpr int nOrder_v = Order;
    static inline const uint16_type nRealDim = RDim;

    static inline const uint16_type topological_dimension = nDim;
    static inline const uint16_type real_dimension = RDim;

    //! Runtime order storage (zero-cost for static order via [[no_unique_address]])
    [[no_unique_address]]
    std::conditional_t<is_order_dynamic, uint16_type, std::monostate> M_runtime_order{};

    typedef typename mpl::at<elements_t, mpl::int_<nDim> >::type element_type;
    typedef typename mpl::at<typename faces_t<real_dimension>::type, mpl::int_<nDim> >::type topological_face_type;
    typedef typename mpl::at<v_edges_t, mpl::int_<real_dimension> >::type edge_type;
    
    static inline const uint16_type numVertices = mpl::at<vertices_t, mpl::int_<Dim> >::type::value;
    static inline const uint16_type numEdges = mpl::at<edges_t, mpl::int_<Dim> >::type::value;
    static inline const uint16_type numFaces = mpl::at<geo_faces_index_t, mpl::int_<Dim> >::type::value;
    static inline const uint16_type numGeometricFaces = mpl::at<geo_faces_index_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type numTopologicalFaces = mpl::at<faces_index_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type numNormals = mpl::at<normals_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type numVolumes = mpl::at<volumes_t, mpl::int_<nDim> >::type::value;

    //! Static point counts (using EffectiveOrder for mpl, only meaningful for static order)
    static inline const uint16_type nbPtsPerVertex_static = ( EffectiveOrder==0 )?( ( nDim==0 )?1:0 ):1;
    static inline const uint16_type nbPtsPerEdge_static = ( EffectiveOrder==0 )?( ( nDim==1 )?1:0 ):(uint16_type)mpl::at<points_edge_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type nbPtsPerFace_static = ( EffectiveOrder==0 )?( ( nDim==2 )?1:0 ):(uint16_type)mpl::at<points_face_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type nbPtsPerVolume_static = ( EffectiveOrder == 0 ) ? ( ( nDim == 3 ) ? 1 : 0 ) : (uint16_type)mpl::at<points_volume_t, mpl::int_<nDim>>::type::value;
    static inline const uint16_type numPoints_static = ( numVertices * nbPtsPerVertex_static +
                                           numEdges * nbPtsPerEdge_static +
                                           numFaces * nbPtsPerFace_static +
                                           numVolumes * nbPtsPerVolume_static );

    //! Backward-compatible static constants.
    //! WARNING: for dynamic order types, these are placeholder P1-layout values.
    //! Use nPointsOn*()/nPointsTotal() for runtime-correct values.
    static inline const uint16_type nbPtsPerVertex = nbPtsPerVertex_static;
    static inline const uint16_type nbPtsPerEdge = nbPtsPerEdge_static;
    static inline const uint16_type nbPtsPerFace = nbPtsPerFace_static;
    static inline const uint16_type nbPtsPerVolume = nbPtsPerVolume_static;
    static inline const uint16_type numPoints = numPoints_static;

    static inline const uint16_type orderSquare = boost::mpl::if_<boost::mpl::greater< boost::mpl::int_<EffectiveOrder>,
                             boost::mpl::int_<5> >,
                             boost::mpl::int_<5>,
                             typename boost::mpl::if_<boost::mpl::less< boost::mpl::int_<EffectiveOrder>,
                             boost::mpl::int_<1> >,
                             boost::mpl::int_<1>,
                             boost::mpl::int_<EffectiveOrder>
                             >::type
                             >::type::value;


    typedef mpl::vector<details::point<orderSquare>, details::line<orderSquare>, details::quad<orderSquare>, details::hexa<orderSquare> > map_entity_to_point_t;
    typedef typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type edge_to_point_t;
    typedef typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type face_to_point_t;
    typedef typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type face_to_edge_t;


    typedef no_permutation vertex_permutation_type;

    typedef typename mpl::if_<mpl::greater_equal<mpl::int_<nDim>, mpl::int_<2> >,
            mpl::identity<line_permutations>,
            mpl::identity<no_permutation> >::type::type edge_permutation_type;


    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
            mpl::identity<quadrangular_faces>,
            mpl::identity<no_permutation> >::type::type face_permutation_type;

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
            mpl::identity<edge_permutation_type>,
            typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
            mpl::identity<face_permutation_type>,
            mpl::identity<no_permutation> >::type>::type::type permutation_type;

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
            mpl::identity<mpl::vector<edge_permutation_type, vertex_permutation_type, vertex_permutation_type> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
            mpl::identity<mpl::vector<face_permutation_type, edge_permutation_type, vertex_permutation_type> >,
            mpl::identity<mpl::vector<vertex_permutation_type, vertex_permutation_type, vertex_permutation_type> > >::type>::type::type permutation_by_subentity_type;

    template<int N>
    using PermutationSubEntity =  typename mpl::at_c<permutation_by_subentity_type,N-1>::type;

    template<int shape_dim, int O = Order,  int R=nDim>
    using shape = Hypercube<shape_dim, O, R>;
    template<int shape_dim, int O = Order,  int R=nDim>
    using shape_t = Hypercube<shape_dim, O, R>;

    //! @name Constructors
    //! @{

    //! Default constructor (only for static order)
    Hypercube() requires( is_order_static ) = default;

    //! Constructor with runtime order (only for dynamic order)
    explicit Hypercube( Feel::RuntimeOrder o ) requires( is_order_dynamic )
        : M_runtime_order( o.value )
    {
    }

    Hypercube( Hypercube const& ) = default;
    Hypercube( Hypercube && ) = default;
    Hypercube& operator=( Hypercube const& ) = default;
    Hypercube& operator=( Hypercube && ) = default;

    //! @}

    //! @name Order Accessors
    //! @{

    /**
     * @brief Returns the geometric order
     * @return The order (static compile-time value or runtime value)
     */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        if constexpr ( is_order_static )
            return nOrder;
        else
            return M_runtime_order;
    }

    //! Static-only order accessor (deleted for dynamic types)
    static constexpr uint16_type staticOrder() requires( is_order_static )
    {
        return static_cast<uint16_type>( Order );
    }
    static constexpr uint16_type staticOrder() requires( is_order_dynamic ) = delete;

    //! @}

    //! @name Dimension Accessors
    //! @{

    /**
     * \return the topological dimension of the hypercube
     */
    uint16_type topologicalDimension() const
    {
        return topological_dimension;
    }

    /**
     * \return the dimension of the space where the hypercube resides
     */
    uint16_type dimension() const
    {
        return real_dimension;
    }

    //! @}

    //! @name Point Count Methods (work for both static and dynamic order)
    //! @{

    /**
     * Returns the number of points per vertex
     */
    [[nodiscard]] uint16_type nPointsOnVertex() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerVertex_static;
        else
            return detail::hypercubePerVertex( nDim, order() );
    }

    /**
     * Returns the number of points per edge
     */
    [[nodiscard]] uint16_type nPointsOnEdge() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerEdge_static;
        else
            return detail::hypercubePerEdge( nDim, order() );
    }

    /**
     * Returns the number of points per face
     */
    [[nodiscard]] uint16_type nPointsOnFace() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerFace_static;
        else
            return detail::hypercubePerFace( nDim, order() );
    }

    /**
     * Returns the number of points per volume
     */
    [[nodiscard]] uint16_type nPointsOnVolume() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerVolume_static;
        else
            return detail::hypercubePerVolume( nDim, order() );
    }

    /**
     * Returns the total number of points (works for both static and dynamic order)
     */
    [[nodiscard]] uint16_type nPointsTotal() const noexcept
    {
        if constexpr ( is_order_static )
            return numPoints_static;
        else
            return static_cast<uint16_type>( detail::hypercubeTotal( nDim, order() ) );
    }

    //! Static-only point-count accessors (deleted for dynamic types)
    static constexpr uint16_type staticPointsPerVertex() requires( is_order_static ) { return nbPtsPerVertex_static; }
    static constexpr uint16_type staticPointsPerVertex() requires( is_order_dynamic ) = delete;
    static constexpr uint16_type staticPointsPerEdge() requires( is_order_static ) { return nbPtsPerEdge_static; }
    static constexpr uint16_type staticPointsPerEdge() requires( is_order_dynamic ) = delete;
    static constexpr uint16_type staticPointsPerFace() requires( is_order_static ) { return nbPtsPerFace_static; }
    static constexpr uint16_type staticPointsPerFace() requires( is_order_dynamic ) = delete;
    static constexpr uint16_type staticPointsPerVolume() requires( is_order_static ) { return nbPtsPerVolume_static; }
    static constexpr uint16_type staticPointsPerVolume() requires( is_order_dynamic ) = delete;
    static constexpr uint16_type staticNumPoints() requires( is_order_static ) { return numPoints_static; }
    static constexpr uint16_type staticNumPoints() requires( is_order_dynamic ) = delete;

    //! @}

    /**
     * \return the number of polynomials of total degree \c n on the
     * shape:
     *
     * -# (n+1)^2 over the quadrangle
     * -# (n+1)^3 over the hexahedron
     */
    template<int N>
    struct PolyDims
    {
        static inline const uint32_type value = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<3> >,
                                 mpl::identity<mpl::int_<( N+1 )*( N+1 )*( N+1 )> >,
                                 typename mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,
                                 mpl::identity<mpl::int_<( N+1 )*( N+1 )> >,
                                 mpl::identity<mpl::int_<( N+1 )> > >::type>::type::value;
    };
    static uint32_type polyDims( int n )
    {
        return uint32_type( math::pow( double( n+1 ), double( nDim ) ) );
    }

    /**
     * Given an edge \p e in the element and the local index \p p (0
     * or 1) of a point in the edge \p e , \return the index in the
     * element of the point.
     */
    static int e2p( int e,  int p )
    {
        return edge_to_point_t::e2p( e, p );
    }

    /**
     * Given a face \p f in the element and the local index \p e of an
     * edge in the face \p f, \return the index in the element of the
     * edge.
     */
    static int f2e( int f,  int e )
    {
        return face_to_edge_t::f2e( f, e );
    }

    /**
     * Given a face \p f in the element and the local index \p p of a
     * point in the face \p f , \return the index in the element of
     * the point.
     */
    static int f2p( int f,  int p )
    {
        return face_to_point_t::f2p( f, p );
    }

    /**
     * \return the name of the hypercube (static version for backward compat)
     */
    static std::string name() requires( is_order_static )
    {
        std::ostringstream ostr;
        ostr << "Hypercube"
             << "_"
             << nDim
             << "_"
             << nOrder
             << "_"
             << nRealDim;
        return ostr.str();
    }

    /**
     * \return the name of the hypercube (instance version for dynamic order)
     */
    [[nodiscard]] std::string name() const requires( is_order_dynamic )
    {
        std::ostringstream ostr;
        ostr << "Hypercube"
             << "_"
             << nDim
             << "_"
             << order()
             << "_"
             << nRealDim;
        return ostr.str();
    }

    static std::string type()
    {
        return "hypercube";
    }
};

}
#endif /* __Hypercube_H */
