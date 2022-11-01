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

#include <feel/feelalg/eigen.hpp>
#include <boost/detail/identifier.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelmesh/entities.hpp>
#include <feel/feelmesh/convex.hpp>
#include <feel/feelmesh/hypercubeordering.hpp>


namespace Feel
{
class HypercubeBase {};
template<uint16_type Dim, int Order=1, uint16_type RDim = Dim>
class Hypercube  : public Convex<Dim,Order,RDim>, HypercubeBase
{ 
public:    
    /**
     * for Dim >= 3 : n edges = n(vertices) + n(faces) - 2
     * thanks to Euler formula
     */
    inline static constexpr std::array<std::array<int,6>,4> topology = {{
                { 1, 2, 4, 8, 16, 32  }, // vertices
                { 0, 1, 4, 12, 32 }, // edges
                { 0, 0, 1, 6, 24 }, // faces
                { 0, 0, 0, 1, 8 } // volumes
    }};
    /**
     * @brief number of facets per topological dimension
     * 
     */
    inline static constexpr std::array<int,5> facets = { 0, 2, 4, 6, 24 };
    inline static constexpr std::array<int, 5> normals = { 0, 2, 4, 6, 24 };

    static constexpr std::array<std::array<int,4>,5> computeNumberPointsPerEntity(int O)
    {
        return std::array<std::array<int,4>,5>{{ 
            { 1,   O+1,     ( O+1)*(O+1),       ( O+1 )*( O+1 )*( O+1 ) }, // all points
            { 0, O+1-2,          (O+1-2),                         O+1-2 }, // edges
            { 0,     0,    ( O-1 )*( O-1 ),               ( O-1 )*( O-1 ) }, // faces
            { 0,     0,                    0,       ( O-1 )*( O-1 )*( O-1 ) }, // volumes
            { 1, O+1-2,  (O+1-2)*( O+1-2 ), ( O+1-2 )*( O+1-2 )*( O+1-2 ) }, //( Order+1-2 )*( Order+1-2 )*( Order+1-2 )*( Order+1-2 ) }, // interior points
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
  inline static constexpr std::array<int, 4> shapes{ { SHAPE_POINT, SHAPE_LINE, SHAPE_TRIANGLE, SHAPE_TETRA } };
  inline static constexpr std::array<int, 4> geometries{ { GEOMETRY_POINT, GEOMETRY_LINE, GEOMETRY_SURFACE, GEOMETRY_VOLUME } };

  inline static const bool is_simplex = false;
  inline static const bool is_hypercube = true;
  inline static const size_type Shape = shapes[Dim];
  inline static const size_type Geometry = geometries[Dim];
  inline static const uint16_type nDim = Dim;
  inline static const int nOrder = Order;
  inline static const uint16_type nRealDim = RDim;
  inline static const uint16_type topological_dimension = nDim;
  inline static const uint16_type real_dimension = RDim;
  inline static const bool is_order_dynamic = ( Order == Dynamic );

  typedef typename mpl::at<elements_t, mpl::int_<nDim>>::type element_type;
  typedef typename mpl::at<typename faces_t<real_dimension>::type, mpl::int_<nDim>>::type topological_face_type;
  typedef typename mpl::at<v_edges_t, mpl::int_<real_dimension>>::type edge_type;

  static constexpr int numberOfVertices()
  {
      return computeNumberPointsPerEntity( 1 )[0][Dim];
  }
  static const int numVertices = numberOfVertices();

  static constexpr int numberOfGeometricFaces()
  {
      if constexpr ( nDim == 3 )
          return 4;
      else if constexpr ( nDim == 2 )
          return 1;
      else
          return 0;
  }

  static const int numFaces = numberOfGeometricFaces();
  static const int numGeometricFaces = numberOfGeometricFaces();

  static constexpr int numberOfTopologicalFaces()
  {
    return topology[2][Dim];
  }
  static const int numTopologicalFaces = numberOfTopologicalFaces();

  static constexpr int numberOfEdges()
  {
    return topology[1][Dim];
  }
  static const int numEdges = numberOfEdges();

  static constexpr int numberOfVolumes()
  {
    return topology[3][Dim];
  }
  inline static const int numVolumes = numberOfVolumes();
  inline static const int numNormals = normals[Dim];

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
          // static_assert(numPoints==numPointsPerEntityAtCompileTime<Order>[0][nDim],"invalid number of points at compile time");
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
  
  static constexpr int orderSquare()
  {
      if constexpr ( Order > 5 )
          return 5;
      else if constexpr ( Order < 1 )
          return 1;
      return Order;
    }
    inline static const uint16_type nOrderSquare = orderSquare();


    typedef mpl::vector<details::point<orderSquare()>, details::line<orderSquare()>, details::quad<orderSquare()>, details::hexa<orderSquare()> > map_entity_to_point_t;
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

    template<uint16_type shape_dim, uint16_type O = Order,  uint16_type R=nDim>
    using shape = Hypercube<shape_dim, O, R>;
    template<uint16_type shape_dim, uint16_type O = Order,  uint16_type R=nDim>
    using shape_t = Hypercube<shape_dim, O, R>;

    Hypercube()
        : Hypercube( Order ) {}
    Hypercube( int order )
        : M_order( order )
    {
        if constexpr ( is_order_dynamic )
        {
            if ( order == 0 )
                numPointsPerEntity = numPointsPerEntityAtCompileTime<0>;
            else
            {
                numPointsPerEntity = computeNumberPointsPerEntity( order );
            }
        }
    }
    Hypercube( Hypercube const& ) = default;
    Hypercube( Hypercube && ) = default;
    Hypercube& operator=( Hypercube const& ) = default;
    Hypercube& operator=( Hypercube && ) = default;
    
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
        if constexpr ( !is_order_dynamic )
            return nbPtsPerVertex;
        else
            return ( M_order == 0 ) ? 0 : 1;
    }

    /**
     * Returns the number of points per edge without the vertices
     */
    constexpr int numberOfPointsPerEdge() const
    {
        if constexpr ( !is_order_dynamic )
            return nbPtsPerEdge;
        else
            return numPointsPerEntity[1][nDim];
    }

    /**
     * Returns the number of points per face
     */
    constexpr int numberOfPointsPerFacet() const
    {
        if constexpr ( !is_order_dynamic )
            return nbPtsPerFace;
        else
            return numPointsPerEntity[2][nDim];
    }

    /**
     * Returns the number of points per volume
     */
    constexpr int numberOfPointsPerVolume() const
    {
        if constexpr ( !is_order_dynamic )
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
     * Returns the number of points per vertex
     */
    static uint16_type nPointsOnVertex()
    {
        return nbPtsPerVertex;
    }

    /**
     * Returns the number of points per edge
     */
    static uint16_type nPointsOnEdge()
    {
        return nbPtsPerEdge;
    }

    /**
     * Returns the number of points per face
     */
    static uint16_type nPointsOnFace()
    {
        return nbPtsPerFace;
    }

    /**
     * Returns the number of points per volume
     */
    static uint16_type nPointsOnVolume()
    {
        return nbPtsPerVolume;
    }

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
        static const uint32_type value = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<3> >,
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
     * \return the name of the simplex product
     */
    static std::string name()
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
    static std::string type()
    {
        return "hypercube";
    }
private:
    int M_order;
    std::array<std::array<int,4>,5> numPointsPerEntity;
};


}
#endif /* __Hypercube_H */
