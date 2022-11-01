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

#include <feel/feelalg/eigen.hpp>
#include <boost/detail/identifier.hpp>
#include <feel/feelmesh/entities.hpp>
#include <feel/feelmesh/convex.hpp>
#include <feel/feelmesh/simplexordering.hpp>


namespace Feel
{
namespace details
{
template<int Order>
struct points
{
    static constexpr int n2D(int i) { return ( Order+1-i )*( Order+2-i )/2; }
    static constexpr int n3D(int i) { return ( Order+1-i )*( Order+2-i )*( Order+3-i )/6; }
    static constexpr int nPerFacet() { return ( Order-1 )*( Order-2 )/2; }
    static constexpr int nPerVolume() { return ( Order-1 )*( Order-2 )*( Order-3 )/6; }
    typedef mpl::vector_c<int, 1, Order+1  , n2D(0),n3D(0) > type;
    typedef mpl::vector_c<int, 0, Order+1-2, n2D(2),n3D(2) > interior_type;
    typedef mpl::vector_c<int, 0, Order+1-2, Order+1-2,Order+1-2 > edge_type;
    typedef mpl::vector_c<int, 0,         0, nPerFacet(), nPerFacet() > face_type;
    typedef mpl::vector_c<int, 0,         0, 0, nPerVolume() > volume_type;
};
template<>
struct points<0>
{
    typedef mpl::vector_c<int, 1, 1, 1, 1> type;
    typedef mpl::vector_c<int, 0, 1, 1, 1> interior_type;
    typedef mpl::vector_c<int, 0, 1, 0, 0> edge_type;
    typedef mpl::vector_c<int, 0, 0, 1, 0> face_type;
    typedef mpl::vector_c<int, 0, 0, 0, 1> volume_type;
};


}
class SimplexBase {};

/**
 * @class Simplex
 *  @brief simplex of dimension \c Dim
 *
 *  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<uint16_type Dim,
         int Order = 1,
         uint16_type RDim = Dim>
class Simplex : public Convex<Dim,Order,RDim>, SimplexBase
{
public:
    /**
     * for Dim >= 3 : n edges = n(vertices) + n(faces) - 2
     * thanks to Euler formula
     */
    inline static constexpr std::array<std::array<int,4>,4> topology = {{
                { 1, 2, 3, 4  }, // vertices
                { 0, 1, 3, ( 4 ) + ( 4 ) - 2 }, // edges
                { 0, 0, 1, 4 }, // faces
                { 0, 0, 0, 1 } // volumes
    }};
    /**
     * @brief number of facets per topological dimension
     * 
     */
    inline static constexpr std::array<int,4> facets = { 0, 2, 3, 4 };

    static constexpr std::array<std::array<int,4>,5> computeNumberPointsPerEntity(int O)
    {
        return std::array<std::array<int,4>,5>{{ 
            { 1,   O+1,     ( O+1)*( O+2 )/2,       ( O+1 )*( O+2 )*( O+3 )/6 }, // all points
            { 0, O+1-2,                O+1-2,                           O+1-2 }, // edges
            { 0,     0,    ( O-1 )*( O-2 )/2,               ( O-1 )*( O-2 )/2 }, // faces
            { 0,     0,                    0,       ( O-1 )*( O-2 )*( O-3 )/6 }, // volumes
            { 1, O+1-2,  (O+1-2)*( O+2-2 )/2, ( O+1-2 )*( O+2-2 )*( O+3-2 )/6 }, // interior points
        }};
    }

    /**
     * @brief number of points per topological entity as well as total and total interior points
     * 
     */
    template<int O=Order>
    inline static constexpr std::array<std::array<int,4>,5> numPointsPerEntityAtCompileTime = computeNumberPointsPerEntity(O);
    template<>
    inline static constexpr std::array<std::array<int,4>,5> numPointsPerEntityAtCompileTime<0> = {{ 
        { 1, 1, 1, 1}, // all points
        { 0, 1, 0, 0 }, // edges
        { 0, 0, 1, 0 }, // faces
        { 0, 0, 0, 1 }, // volumes
        { 0, 1, 1, 1 }, // interior points
    }};

    typedef mpl::vector_c<int, 0, 1, 3, ( 4 ) + ( 4 ) - 2> edges_t;
    typedef mpl::vector_c<int, 0, 0, 1, 4> geo_faces_index_t;
    typedef mpl::vector_c<int, 0, 2, 3, 4> faces_index_t;
    typedef mpl::vector_c<int, 0, 0, 0, 1> volumes_t;
    typedef mpl::vector_c<int, 0, 2, 3, 4> normals_t;

    typedef typename details::points<Order>::type points_t;
    typedef typename details::points<Order>::interior_type points_interior_t;
    typedef typename details::points<Order>::edge_type points_edge_t;
    typedef typename details::points<Order>::face_type points_face_t;
    typedef typename details::points<Order>::volume_type points_volume_t;

    inline static constexpr std::array<int,4> shapes{ { SHAPE_POINT, SHAPE_LINE, SHAPE_TRIANGLE, SHAPE_TETRA } };
    inline static constexpr std::array<int,4> geometries{ { GEOMETRY_POINT, GEOMETRY_LINE, GEOMETRY_SURFACE, GEOMETRY_VOLUME } };


    static constexpr int computeOrderTriangle() 
        {
            if constexpr ( Order > 5 ) return 5;
            if constexpr ( Order < 1 ) return 1;
            return Order;
        }
    inline static constexpr int orderTriangle = computeOrderTriangle();

    typedef mpl::vector<details::point<orderTriangle>, details::line<orderTriangle>, details::triangle<orderTriangle>, details::tetra<orderTriangle> > map_entity_to_point_t;

    typedef mpl::vector_c<uint16_type, 0, 1, 2, 6> permutations_t;

    template<uint16_type rdim>
    struct faces_t
    {
        typedef mpl::vector<Simplex<0, Order, rdim>,
                            Simplex<0, Order, rdim>,
                            Simplex<1, Order, rdim>,
                            Simplex<2, Order, rdim> > type;
    };

    typedef mpl::vector<Simplex<0, Order,0>, Simplex<1, Order,1>, Simplex<1, Order, 2>, Simplex<1, Order, 3>, boost::none_t > v_edges_t;
    typedef mpl::vector<Simplex<1, Order>, Simplex<2, Order>, Simplex<3, Order>, boost::none_t > elements_t;

public:

    inline static const bool is_simplex = true;
    inline static const bool is_hypercube = false;

    inline static const uint16_type nDim = Dim;
    inline static const int nOrder = Order;
    inline static const uint16_type nRealDim = RDim;
    inline static const bool is_order_dynamic = (Order == Dynamic);
    inline static const uint16_type topological_dimension = nDim;
    inline static const uint16_type real_dimension = nRealDim;

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

    
    inline static const int Geometry = geometries[nDim];

    typedef typename mpl::at<elements_t, mpl::int_<nDim> >::type element_type;
    typedef typename mpl::at<typename faces_t<real_dimension>::type, mpl::int_<nDim> >::type topological_face_type;
    typedef typename mpl::at<v_edges_t, mpl::int_<real_dimension> >::type edge_type;
    typedef topological_face_type GeoBShape;

    static constexpr int numberOfVertices()
        {
            return nDim+1;
        }
    static const int numVertices = numberOfVertices();
    
    static constexpr int numberOfGeometricFaces()
        {
            if constexpr ( nDim == 3 ) return 4;
            else if constexpr ( nDim == 2 ) return 1;
            else return 0;
        }

    static const int numFaces = numberOfGeometricFaces();
    static const int numGeometricFaces = numberOfGeometricFaces();

    static constexpr int numberOfTopologicalFaces()
        {
            if constexpr ( nDim == 3 ) return 4;
            else if constexpr ( nDim == 2 ) return 3;
            else if constexpr ( nDim == 1 ) return 2;
            else return 0;
        }
    static const int numTopologicalFaces = numberOfTopologicalFaces();

    static constexpr int numberOfEdges()
        {
            if constexpr ( nDim == 3 ) return 4+4-2;
            else if constexpr ( nDim == 2 ) return 3;
            else if constexpr ( nDim == 1 ) return 1;
            else return 0;
        }
    static const int numEdges = numberOfEdges();

    static constexpr int numberOfVolumes()
        {
            if constexpr ( nDim == 3 ) return 1;
            else return 0;
        }
    inline static const int numVolumes = numberOfVolumes();
    inline static const int numNormals = mpl::at<normals_t, mpl::int_<nDim> >::type::value;

    
    static constexpr int numberOfPointsPerVertexAtCompileTime()
    {
        if constexpr ( !is_order_dynamic )
        {
            return ( nOrder == 0 ) ? 0 : 1;
        }
        else
            return -1;
    }
    static constexpr int numberOfPointsPerEdgeAtCompileTime()
    {
        if constexpr ( !is_order_dynamic )
        {
            return numPointsPerEntityAtCompileTime<Order>[1][nDim];
        }
        else
            return -1;
    }
    static constexpr int numberOfPointsPerFacetAtCompileTime()
    {
        if constexpr ( !is_order_dynamic )
        {
            return numPointsPerEntityAtCompileTime<Order>[2][nDim];
        }
        else
            return -1;
    }
    static constexpr int numberOfPointsPerVolumeAtCompileTime()
    {
        if constexpr ( !is_order_dynamic )
        {
            return numPointsPerEntityAtCompileTime<Order>[3][nDim];
        }
        else
            return -1;
    }
    static constexpr int numberOfPointsAtCompileTime()
    {
        if constexpr ( !is_order_dynamic )
        {
            //static_assert(numPoints==numPointsPerEntityAtCompileTime<Order>[0][nDim],"invalid number of points at compile time");
            return numPointsPerEntityAtCompileTime<Order>[0][nDim];
        }
        else
            return -1;
    }
    inline static constexpr int nbPtsPerVertex = numberOfPointsPerVertexAtCompileTime();
    inline static constexpr int nbPtsPerEdge = numberOfPointsPerEdgeAtCompileTime();
    inline static constexpr int nbPtsPerFace = numberOfPointsPerFacetAtCompileTime();
    inline static constexpr int nbPtsPerVolume = numberOfPointsPerVolumeAtCompileTime();
    inline static constexpr int numPoints = numberOfPointsAtCompileTime();
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


    Simplex(): Simplex( Order ) {}
    Simplex( int order ): M_order( order ) 
    {
        if constexpr ( is_order_dynamic )
        {
            if ( order == 0 )
                numPointsPerEntity = numPointsPerEntityAtCompileTime<0>;
            else
            {
                numPointsPerEntity = computeNumberPointsPerEntity(order);
            }
        }
    }
    Simplex( Simplex const& ) = default;
    Simplex( Simplex && ) = default;
    Simplex& operator=( Simplex const& ) = default;
    Simplex& operator=( Simplex && ) = default;

    /**
     * \return the topological dimension of the simplex
     */
    uint16_type topologicalDimension() const
    {
        return topological_dimension;
    }

    /**
     * \return the dimension of the space where the simplex resides
     */
    uint16_type dimension() const
    {
        return real_dimension;
    }

    /**
     * Returns the number of points per vertex
     */
    constexpr int numberOfPointsPerVertex() const
    {
        if constexpr( !is_order_dynamic )
            return nbPtsPerVertex;
        else
            return ( M_order == 0 ) ? 0 : 1;
    }

    /**
     * Returns the number of points per edge without the vertices
     */
    constexpr int numberOfPointsPerEdge() const
    {
        if constexpr( !is_order_dynamic )
            return nbPtsPerEdge;
        else
            return  numPointsPerEntity[1][nDim];
    }

    /**
     * Returns the number of points per face
     */
    constexpr int numberOfPointsPerFacet() const
    {
        if constexpr( !is_order_dynamic )
            return nbPtsPerFace;
        else
            return numPointsPerEntity[2][nDim];
    }

    /**
     * Returns the number of points per volume
     */
    constexpr int numberOfPointsPerVolume() const
    {
        if constexpr( !is_order_dynamic )
            return nbPtsPerVolume;
        else
            return numPointsPerEntity[3][nDim];
    }

    /**
     * @brief get the number of points in the Simplex 
     * 
     * @return int the number of points
     */
    constexpr int numberOfPoints() const
    {
        if constexpr ( !is_order_dynamic )
            return numPoints;
        else
            return numPointsPerEntity[0][nDim];
    }
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
        
        static_assert( nDim == 0 || nDim == 1 || nDim == 2 || nDim == 3, "Invalid simplex topological dimension" );
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
     * \return the name of the simplex
     */
    static std::string name()
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
    static std::string type()
    {
        return "simplex";
    }
private:
    int M_order;
    std::array<std::array<int,4>,5> numPointsPerEntity;
};

template<int Dim> struct Line : public Simplex<1, Dim> {};
template<int Dim> struct Triangle : public Simplex<2, Dim> {};
template<int Dim> struct Tetrahedron : public Simplex<3, Dim> {};


}

#endif /* __Simplex_H */
