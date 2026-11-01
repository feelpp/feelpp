/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2006-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
class Hypercube : public Convex<
                      Dim,
                      ( OrderBase<Order>::order_spec >= 0
                            ? OrderBase<Order>::order_spec
                            : 1 ),
                      RDim>,
                  public OrderBase<Order>, HypercubeBase
{
public:
    using order_base_type = OrderBase<Order>;
    using order_base_type::order;

    //! @name Order type detection
    //! @{
    static constexpr bool is_order_static = order_base_type::is_order_static;
    static constexpr bool is_order_dynamic = order_base_type::is_order_dynamic;
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
    static constexpr int EffectiveOrder =
        is_order_static ? order_base_type::order_spec : 1;

    template<uint16_type rdim>
    struct faces_t
    {
        typedef mpl::vector<Hypercube<0, Order, rdim>,
                            Hypercube<0, Order, rdim>,
                            Hypercube<1, Order, rdim>,
                            Hypercube<2, Order, rdim> > type;
    };
    typedef mpl::vector<boost::none_t,Hypercube<1, Order,1>, Hypercube<1, Order, 2>, Hypercube<1, Order, 3>, boost::none_t > v_edges_t;
    typedef mpl::vector<Hypercube<1, Order>, Hypercube<2, Order>, Hypercube<3, Order>, Hypercube<4, Order>, boost::none_t > elements_t;

    typedef mpl::vector_c<uint16_type, 0, 1, 2, 8> permutations_t;

public:

    static inline const bool is_simplex = false;
    static inline const bool is_hypercube = true;

    static const size_type Shape = mpl::at<shapes_t, mpl::int_<Dim> >::type::value;
    static const size_type Geometry = mpl::at<geometries_t, mpl::int_<Dim> >::type::value;

    static constexpr int order_spec = order_base_type::order_spec;
    static inline const uint16_type nDim = Dim;
    //! Static order value (backward compatibility only).
    //! WARNING: for dynamic order types, this is a placeholder (0). Use order().
    static inline const uint16_type nOrder = is_order_static ? static_cast<uint16_type>( order_spec ) : 0;
    //! Template order parameter value (may be Dynamic = -1)
    static constexpr int nOrder_v = order_spec;
    static inline const uint16_type nRealDim = RDim;

    static inline const uint16_type topological_dimension = nDim;
    static inline const uint16_type real_dimension = RDim;

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

    //! Static point counts derived from the shared cardinality kernel.
    static constexpr uint16_type nbPtsPerVertex_static =
        detail::HypercubeCardinality::perVertex( nDim, EffectiveOrder );
    static constexpr uint16_type nbPtsPerEdge_static =
        detail::HypercubeCardinality::perEdge( nDim, EffectiveOrder );
    static constexpr uint16_type nbPtsPerFace_static =
        detail::HypercubeCardinality::perFace( nDim, EffectiveOrder );
    static constexpr uint16_type nbPtsPerVolume_static =
        detail::HypercubeCardinality::perVolume( nDim, EffectiveOrder );
    static constexpr uint16_type numPoints_static = static_cast<uint16_type>(
        detail::HypercubeCardinality::total( nDim, EffectiveOrder ) );

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


    typedef mpl::vector<Feel::details::point<orderSquare>, Feel::details::line<orderSquare>, Feel::details::quad<orderSquare>, Feel::details::hexa<orderSquare> > map_entity_to_point_t;
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
        : order_base_type(
              Feel::detail::validatedNodalOrder<
                  Feel::detail::HypercubeCardinality>( Dim, o ) )
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
     * @brief Clear the construction-only order after mesh attachment.
     */
    void clearDetachedOrder() noexcept
        requires( is_order_dynamic )
    {
        this->setOrder( RuntimeOrder{ 0 } );
    }

    //! Static-only order accessor (deleted for dynamic types)
    static constexpr uint16_type staticOrder() requires( is_order_static )
    {
        return order_base_type::staticOrder();
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
    [[nodiscard]] constexpr uint16_type nPointsOnVertex() const noexcept
    {
        return detail::HypercubeCardinality::perVertex( nDim, order() );
    }

    /**
     * Returns the number of points per edge
     */
    [[nodiscard]] constexpr uint16_type nPointsOnEdge() const noexcept
    {
        return detail::HypercubeCardinality::perEdge( nDim, order() );
    }

    /**
     * Returns the number of points per face
     */
    [[nodiscard]] constexpr uint16_type nPointsOnFace() const noexcept
    {
        return detail::HypercubeCardinality::perFace( nDim, order() );
    }

    /**
     * Returns the number of points per volume
     */
    [[nodiscard]] constexpr uint16_type nPointsOnVolume() const noexcept
    {
        return detail::HypercubeCardinality::perVolume( nDim, order() );
    }

    /**
     * Returns the total number of points (works for both static and dynamic order)
     */
    [[nodiscard]] constexpr uint16_type nPointsTotal() const noexcept
    {
        return static_cast<uint16_type>(
            detail::HypercubeCardinality::total( nDim, order() ) );
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
