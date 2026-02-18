/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
   \file simplex.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-02-20
 */
#ifndef __Simplex_H
#define __Simplex_H 1

#include <boost/detail/identifier.hpp>
#include <feel/feelmesh/entities.hpp>
#include <feel/feelmesh/convex.hpp>
#include <feel/feelmesh/simplexordering.hpp>
#include <feel/feelpoly/order.hpp>

#include <variant>

namespace Feel
{
namespace details
{
template<int Order>
struct points
{
    typedef mpl::vector_c<size_type, 1, Order+1, ( Order+1 )*( Order+2 )/2, ( Order+1 )*( Order+2 )*( Order+3 )/6 > type;
    typedef mpl::vector_c<size_type, 0, Order+1-2, ( Order+1-2 )*( Order+2-2 )/2, ( Order+1-2 )*( Order+2-2 )*( Order+3-2 )/6 > interior_type;
    typedef mpl::vector_c<size_type, 0, Order+1-2,                 Order+1-2,                             Order+1-2 > edge_type;
    typedef mpl::vector_c<size_type, 0,         0, ( Order-1 )*( Order-2 )/2,             ( Order-1 )*( Order-2 )/2 > face_type;
    typedef mpl::vector_c<size_type, 0,         0,                         0, ( Order-1 )*( Order-2 )*( Order-3 )/6 > volume_type;
};
template<>
struct points<0>
{
    typedef mpl::vector_c<size_type, 1, 1, 1, 1> type;
    typedef mpl::vector_c<size_type, 0, 1, 1, 1> interior_type;
    typedef mpl::vector_c<size_type, 0, 1, 0, 0> edge_type;
    typedef mpl::vector_c<size_type, 0, 0, 1, 0> face_type;
    typedef mpl::vector_c<size_type, 0, 0, 0, 1> volume_type;
};


}
class SimplexBase {};

/**
 * @class Simplex
 *  @brief simplex of dimension \c Dim
 *
 * Supports both static order (default, backward compatible) and dynamic order.
 * For dynamic order, use Simplex<Dim, Dynamic, RDim> and construct with Order(n).
 *
 *  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<int Dim,
         int Order = 1,
         int RDim = Dim>
class Simplex : public Convex<Dim, (Order >= 0 ? Order : 1), RDim>, SimplexBase
{
public:
    //! @name Order type detection
    //! @{
    static constexpr bool is_order_static = (Order != Dynamic);
    static constexpr bool is_order_dynamic = !is_order_static;
    //! @}

private:
    //! Placeholder order used for internal type computations when Order is Dynamic
    static constexpr int OrderPlaceholder = is_order_static ? Order : 1;
    /**
     * for Dim >= 3 : n edges = n(vertices) + n(faces) - 2
     * thanks to Euler formula
     */
    typedef mpl::vector_c<uint16_type, 0, 1, 3, ( 4 ) + ( 4 ) - 2> edges_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 1, 4> geo_faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 3, 4> faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 0, 1> volumes_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 3, 4> normals_t;

    // Use OrderPlaceholder for mpl-based types (works for both static and dynamic)
    typedef typename details::points<OrderPlaceholder>::type points_t;
    typedef typename details::points<OrderPlaceholder>::interior_type points_interior_t;
    typedef typename details::points<OrderPlaceholder>::edge_type points_edge_t;
    typedef typename details::points<OrderPlaceholder>::face_type points_face_t;
    typedef typename details::points<OrderPlaceholder>::volume_type points_volume_t;

    typedef mpl::vector_c<size_type, SHAPE_POINT, SHAPE_LINE, SHAPE_TRIANGLE, SHAPE_TETRA> shapes_t;
    typedef mpl::vector_c<size_type, GEOMETRY_POINT, GEOMETRY_LINE, GEOMETRY_SURFACE, GEOMETRY_VOLUME> geometries_t;

    static constexpr int computeOrderTriangle()
        {
            if constexpr ( OrderPlaceholder > 5 ) return 5;
            if constexpr ( OrderPlaceholder < 1 ) return 1;
            return OrderPlaceholder;
        }
    inline static constexpr int orderTriangle = computeOrderTriangle();

    typedef mpl::vector<details::point<orderTriangle>, details::line<orderTriangle>, details::triangle<orderTriangle>, details::tetra<orderTriangle> > map_entity_to_point_t;

    typedef mpl::vector_c<uint16_type, 0, 1, 2, 6> permutations_t;

    template<uint16_type rdim>
    struct faces_t
    {
        typedef mpl::vector<Simplex<0, OrderPlaceholder, rdim>,
                            Simplex<0, OrderPlaceholder, rdim>,
                            Simplex<1, OrderPlaceholder, rdim>,
                            Simplex<2, OrderPlaceholder, rdim> > type;
    };

    typedef mpl::vector<Simplex<0, OrderPlaceholder, 0>, Simplex<1, OrderPlaceholder, 1>, Simplex<1, OrderPlaceholder, 2>, Simplex<1, OrderPlaceholder, 3>, boost::none_t > v_edges_t;
    typedef mpl::vector<Simplex<1, OrderPlaceholder>, Simplex<2, OrderPlaceholder>, Simplex<3, OrderPlaceholder>, boost::none_t > elements_t;

    //! Runtime order storage (only used when is_order_dynamic)
    using order_storage_type = std::conditional_t<is_order_dynamic, uint16_type, std::monostate>;
    [[no_unique_address]] order_storage_type M_runtime_order{};

public:

    static inline const bool is_simplex = true;
    static inline const bool is_hypercube = false;

    static inline const uint16_type nDim = Dim;
    //! Static order value (backward compatibility only).
    //! WARNING: for dynamic order types, this is a placeholder (0). Use order().
    static inline const uint16_type nOrder = is_order_static ? static_cast<uint16_type>(Order) : 0;
    //! Template order parameter value (may be Dynamic = -1)
    static constexpr int nOrder_v = Order;
    static inline const uint16_type nRealDim = RDim;

    static inline const uint16_type topological_dimension = nDim;
    static inline const uint16_type real_dimension = nRealDim;

    /**
     * handle a shape 
     */
    static constexpr int getShape()
        {
            if constexpr ( nDim == 0 ) return SHAPE_POINT;
            else if constexpr ( nDim == 1 ) return SHAPE_LINE;
            else if constexpr ( nDim == 2 ) return SHAPE_TRIANGLE;
            else if constexpr ( nDim == 3 ) return SHAPE_TETRA;
        }
    static constexpr int Shape = getShape();
    template<int shape_dim = nDim, int O = Order,  int R=nDim>
    using shape = Simplex<shape_dim, O, R>;
    template<int shape_dim, int O = Order,  int R=nDim>
    using shape_t = Simplex<shape_dim, O, R>;
    template<int shape_dim, int  O = Order,  int  R=nDim>
    static constexpr int shape_v = Simplex<shape_dim, O, R>::Shape;

    
    static const size_type Geometry = mpl::at<geometries_t, mpl::int_<nDim> >::type::value;

    typedef typename mpl::at<elements_t, mpl::int_<nDim> >::type element_type;
    typedef typename mpl::at<typename faces_t<real_dimension>::type, mpl::int_<nDim> >::type topological_face_type;
    typedef typename mpl::at<v_edges_t, mpl::int_<real_dimension> >::type edge_type;
    typedef topological_face_type GeoBShape;

    static constexpr int numberOfVertices()
        {
            return nDim+1;
        }
    static inline const int numVertices = numberOfVertices();
    
    static constexpr int numberOfGeometricFaces()
        {
            if constexpr ( nDim == 3 ) return 4;
            else if constexpr ( nDim == 2 ) return 1;
            else return 0;
        }

    static inline const int numFaces = numberOfGeometricFaces();
    static inline const int numGeometricFaces = numberOfGeometricFaces();

    static constexpr int numberOfTopologicalFaces()
        {
            if constexpr ( nDim == 3 ) return 4;
            else if constexpr ( nDim == 2 ) return 3;
            else if constexpr ( nDim == 1 ) return 2;
            else return 0;
        }
    static inline const int numTopologicalFaces = numberOfTopologicalFaces();

    static constexpr int numberOfEdges()
        {
            if constexpr ( nDim == 3 ) return 4+4-2;
            else if constexpr ( nDim == 2 ) return 3;
            else if constexpr ( nDim == 1 ) return 1;
            else return 0;
        }
    static inline const int numEdges = numberOfEdges();

    static constexpr int numberOfVolumes()
        {
            if constexpr ( nDim == 3 ) return 1;
            else return 0;
        }
    static const int numVolumes = numberOfVolumes();

    static inline const uint16_type numNormals = mpl::at<normals_t, mpl::int_<nDim> >::type::value;

    //! Static point counts (only valid for static order; use runtime methods for dynamic)
    static inline const uint16_type nbPtsPerVertex_static = ( OrderPlaceholder == 0 ) ? ( ( nDim == 0 ) ? 1 : 0 ) : 1;
    static inline const uint16_type nbPtsPerEdge_static = mpl::at<points_edge_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type nbPtsPerFace_static = mpl::at<points_face_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type nbPtsPerVolume_static = mpl::at<points_volume_t, mpl::int_<nDim> >::type::value;
    static inline const uint16_type numPoints_static = ( numVertices * nbPtsPerVertex_static +
                                                         numEdges * nbPtsPerEdge_static +
                                                         numFaces * nbPtsPerFace_static +
                                                         numVolumes * nbPtsPerVolume_static );

    //! Backward compatible static constants.
    //! WARNING: for dynamic order types, these are placeholder P1-layout values.
    //! Use nPointsOn*()/nPointsTotal() for runtime-correct values.
    static inline const uint16_type nbPtsPerVertex = nbPtsPerVertex_static;
    static inline const uint16_type nbPtsPerEdge = nbPtsPerEdge_static;
    static inline const uint16_type nbPtsPerFace = nbPtsPerFace_static;
    static inline const uint16_type nbPtsPerVolume = nbPtsPerVolume_static;
    static inline const uint16_type numPoints = numPoints_static;

    typedef typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type edge_to_point_t;
    typedef typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type face_to_point_t;
    typedef typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type face_to_edge_t;


    typedef no_permutation vertex_permutation_type;

    typedef typename mpl::if_<mpl::greater_equal<mpl::int_<nDim>, mpl::int_<2> >,
            mpl::identity<line_permutations>,
            mpl::identity<no_permutation> >::type::type edge_permutation_type;


    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
            mpl::identity<triangular_faces_type>,
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
    using PermutationSubEntity = typename mpl::at_c<permutation_by_subentity_type,N-1>::type;


    //! @name Constructors
    //! @{

    //! Default constructor (only for static order)
    Simplex() requires( is_order_static ) = default;

    //! Constructor with runtime order (only for dynamic order)
    explicit Simplex( Feel::RuntimeOrder o ) requires( is_order_dynamic )
        : M_runtime_order( o.value )
    {
    }

    Simplex( Simplex const& ) = default;
    Simplex( Simplex&& ) = default;
    Simplex& operator=( Simplex const& ) = default;
    Simplex& operator=( Simplex&& ) = default;

    //! @}

    //! @name Order accessors
    //! @{

    /**
     * \return the geometric order of the simplex
     */
    [[nodiscard]] uint16_type order() const noexcept
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

    /**
     * \return the topological dimension of the simplex
     */
    [[nodiscard]] uint16_type topologicalDimension() const noexcept
    {
        return topological_dimension;
    }

    /**
     * \return the dimension of the space where the simplex resides
     */
    [[nodiscard]] uint16_type dimension() const noexcept
    {
        return real_dimension;
    }

    //! @name Point count accessors (work for both static and dynamic order)
    //! @{

    /**
     * Returns the number of points per vertex
     */
    [[nodiscard]] uint16_type nPointsOnVertex() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerVertex_static;
        else
            return detail::simplexPerVertex( nDim, order() );
    }

    /**
     * Returns the number of points per edge
     */
    [[nodiscard]] uint16_type nPointsOnEdge() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerEdge_static;
        else
            return detail::simplexPerEdge( nDim, order() );
    }

    /**
     * Returns the number of points per face
     */
    [[nodiscard]] uint16_type nPointsOnFace() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerFace_static;
        else
            return detail::simplexPerFace( nDim, order() );
    }

    /**
     * Returns the number of points per volume
     */
    [[nodiscard]] uint16_type nPointsOnVolume() const noexcept
    {
        if constexpr ( is_order_static )
            return nbPtsPerVolume_static;
        else
            return detail::simplexPerVolume( nDim, order() );
    }

    /**
     * Returns the total number of points (works for both static and dynamic order)
     */
    [[nodiscard]] uint16_type nPointsTotal() const noexcept
    {
        if constexpr ( is_order_static )
            return numPoints_static;
        else
            return static_cast<uint16_type>( detail::simplexTotal( nDim, order() ) );
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
     * -# (n+1) over the line
     * -# (n+1)(n+2)/2 over the triangle,
     * -# (n+1)(n+2)(n+3)/6 in three dimensions.
     */
    static int polyDims( int n )
    {
        if constexpr ( nDim == 0 )
            return (n>0)?0:1;

        else if constexpr ( nDim == 1 )
            return std::max( 0, n + 1 );

        else if constexpr ( nDim == 2 )
            return std::max( 0, ( n+1 )*( n+2 )/2 );

        else if constexpr ( nDim == 3 )
            return std::max( 0, ( n+1 )*( n+2 )*( n+3 )/6 );

        return -1;
    }

    /**
     * Given an edge \p e in the element and the local index \p p (0
     * or 1) of a point in the edge \p e , \return the index in the
     * element of the point.
     */
    static uint16_type e2p( uint16_type e,  uint16_type p )
    {
        return edge_to_point_t::e2p( e, p );
    }

    /**
     * Given a face \p f in the element and the local index \p e of an
     * edge in the face \p f, \return the index in the element of the
     * edge.
     */
    static uint16_type f2e( uint16_type f,  uint16_type e )
    {
        return face_to_edge_t::f2e( f, e );
    }

    /**
     * Given a face \p f in the element and the global index \p e of an
     * edge in the face \p f, \return the local index in the element of the
     * edge.
     */
    static uint16_type f2eLoc( uint16_type f,  uint16_type e )
    {
        return face_to_edge_t::f2eLoc( f, e );
    }

    /**
     * Given a face \p f in the element and the local index \p p of a
     * point in the face \p f , \return the index in the element of
     * the point.
     */
    static uint16_type f2p( uint16_type f,  uint16_type p )
    {
        return face_to_point_t::f2p( f, p );
    }

    /**
     * \return the name of the simplex (static version for backward compat)
     */
    static std::string name()
        requires( is_order_static )
    {
        std::ostringstream ostr;
        ostr << "Simplex"
             << "_"
             << nDim
             << "_"
             << nOrder
             << "_"
             << nRealDim;
        return ostr.str();
    }

    /**
     * \return the name of the simplex (instance method for dynamic order)
     */
    [[nodiscard]] std::string name() const
        requires( is_order_dynamic )
    {
        std::ostringstream ostr;
        ostr << "Simplex"
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
        return "simplex";
    }
};

template<int Dim> struct Line : public Simplex<1, Dim> {};
template<int Dim> struct Triangle : public Simplex<2, Dim> {};
template<int Dim> struct Tetrahedron : public Simplex<3, Dim> {};


}

#endif /* __Simplex_H */
