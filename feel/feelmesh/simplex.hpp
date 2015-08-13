/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-02-20

  Copyright (C) 2006 EPFL
  Copyright (C) 2010-2015 Feel++ Consortium

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


namespace Feel
{
namespace details
{
template<int Order>
struct points
{
#if 0
    static constexpr numberOfPoints( uint16_type Dim, uint16_type Order ) 
        {
            uint16_type pts[4] = { 1, Order+1, ( Order+1 )*( Order+2 )/2, ( Order+1 )*( Order+2 )*( Order+3 )/6};
            return pts[Dim];
        }
#endif
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
 *  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<uint16_type Dim,
         uint16_type Order = 1,
         uint16_type RDim = Dim>
class Simplex : public Convex, SimplexBase
{
private:
    using super = Convex;
    using super2 = SimplexBase;
    /**
     * for Dim >= 3 : n edges = n(vertices) + n(faces) - 2
     * thanks to Euler formula
     */
    typedef mpl::vector_c<uint16_type, 0, 1, 3, ( 4 ) + ( 4 ) - 2> edges_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 1, 4> geo_faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 3, 4> faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 0, 1> volumes_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 3, 4> normals_t;

#if 0
    typedef typename details::points<Order>::type points_t;
    typedef typename details::points<Order>::interior_type points_interior_t;
    typedef typename details::points<Order>::edge_type points_edge_t;
    typedef typename details::points<Order>::face_type points_face_t;
    typedef typename details::points<Order>::volume_type points_volume_t;
#endif

    typedef mpl::vector_c<size_type, SHAPE_POINT, SHAPE_LINE, SHAPE_TRIANGLE, SHAPE_TETRA> shapes_t;
    typedef mpl::vector_c<size_type, GEOMETRY_POINT, GEOMETRY_LINE, GEOMETRY_SURFACE, GEOMETRY_VOLUME> geometries_t;

    
    typedef mpl::vector<details::vertex, details::line, details::triangle, details::tetra > map_entity_to_point_t;

    typedef mpl::vector_c<uint16_type, 0, 1, 2, 6> permutations_t;

    template<int rdim>
    using faces_t = mpl::vector<Simplex<0, Order, rdim>,
                                Simplex<0, Order, rdim>,
                                Simplex<1, Order, rdim>,
                                Simplex<2, Order, rdim> >;
    using v_edges_t = mpl::vector<Simplex<0, Order,0>,Simplex<1, Order,1>, Simplex<1, Order, 2>, Simplex<1, Order, 3>, boost::none_t >; 

    
    typedef mpl::vector<Simplex<1, Order>, Simplex<2, Order>, Simplex<3, Order>, boost::none_t > elements_t;

public:

    static constexpr bool is_simplex = true;
    static constexpr bool is_hypercube = false;

    static constexpr uint16_type nDim = Dim;
    //static constexpr uint16_type nOrder = Order;
    static constexpr uint16_type nRealDim = RDim;

    static constexpr uint16_type topological_dimension = nDim;
    static constexpr uint16_type real_dimension = nRealDim;

    static constexpr size_type Shape = mpl::at<shapes_t, mpl::int_<nDim> >::type::value;
    static constexpr size_type Geometry = mpl::at<geometries_t, mpl::int_<nDim> >::type::value;

    using element_type = typename mpl::at<elements_t, mpl::int_<nDim> >::type;
    using topological_face_type = typename mpl::at<faces_t<real_dimension>, mpl::int_<nDim> >::type;
    using edge_type = typename mpl::at<v_edges_t, mpl::int_<real_dimension> >::type;
    typedef topological_face_type GeoBShape;

    static constexpr uint16_type numVertices = nDim+1;
    static constexpr uint16_type numFaces = mpl::at<geo_faces_index_t, mpl::int_<nDim> >::type::value;
    static constexpr uint16_type numGeometricFaces = mpl::at<geo_faces_index_t, mpl::int_<nDim> >::type::value;
    static constexpr uint16_type numTopologicalFaces = mpl::at<faces_index_t, mpl::int_<nDim> >::type::value;
    static constexpr uint16_type numEdges = mpl::at<edges_t, mpl::int_<nDim> >::type::value;
    static constexpr uint16_type numVolumes = mpl::at<volumes_t, mpl::int_<nDim> >::type::value;

    static constexpr uint16_type numNormals = mpl::at<normals_t, mpl::int_<nDim> >::type::value;

    static constexpr uint16_type numberOfPointsPerVertex( uint16_type _Order ) 
        {
            return ( _Order==0 )?0:1;
        }
    static constexpr uint16_type numberOfPointsPerEdge( uint16_type _Order ) 
        {
            int pts[4] = { 0, _Order+1-2,                 _Order+1-2,                             _Order+1-2};
            return (pts[Dim]>=0)?pts[Dim]:0;
        }
    static constexpr uint16_type numberOfPointsPerFace( uint16_type _Order ) 
        {
            int pts[4] = { 0,         0, ( _Order-1 )*( _Order-2 )/2,             ( _Order-1 )*( _Order-2 )/2 };
            return (pts[Dim]>=0)?pts[Dim]:0;
        }
    static constexpr uint16_type numberOfPointsPerVolume( uint16_type _Order ) 
        {
            int pts[4] = { 0,         0,                         0, ( _Order-1 )*( _Order-2 )*( _Order-3 )/6 };
            return (pts[Dim]>=0)?pts[Dim]:0;
        }
#if 0
    typedef mpl::vector_c<size_type, 0, _Order+1-2, ( _Order+1-2 )*( _Order+2-2 )/2, ( _Order+1-2 )*( _Order+2-2 )*( _Order+3-2 )/6 > interior_type;
    typedef mpl::vector_c<size_type, 0, _Order+1-2,                 _Order+1-2,                             _Order+1-2 > edge_type;
    typedef mpl::vector_c<size_type, 0,         0, ( _Order-1 )*( _Order-2 )/2,             ( _Order-1 )*( _Order-2 )/2 > face_type;
    typedef mpl::vector_c<size_type, 0,         0,                         0, ( _Order-1 )*( _Order-2 )*( _Order-3 )/6 > volume_type;
#endif
    static constexpr uint16_type numberOfPoints( uint16_type _Order )
        { 
            return ( numVertices * numberOfPointsPerVertex( _Order ) +
                     numEdges * numberOfPointsPerEdge( _Order ) +
                     numFaces * numberOfPointsPerFace( _Order ) + 
                     numVolumes * numberOfPointsPerVolume( _Order ) ); 
        }
    
    constexpr uint16_type numberOfPointsPerVertex() const { return numberOfPointsPerVertex( order() ); }
    constexpr uint16_type numberOfPointsPerEdge() const { return numberOfPointsPerEdge( order() ); }
    constexpr uint16_type numberOfPointsPerFace() const { return  numberOfPointsPerFace( order() ); }
    constexpr uint16_type numberOfPointsPerVolume() const { return  numberOfPointsPerVolume( order() ); }
    constexpr uint16_type numberOfPoints() const { return numberOfPoints( order() ); }


    using edge_to_point_t = typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type;
    using face_to_point_t = typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type;
    using face_to_edge_t =  typename mpl::at<map_entity_to_point_t, mpl::int_<nDim> >::type;


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

    template<uint16_type shape_dim, uint16_type O = Order,  uint16_type R=nDim>
    struct shape
    {
        typedef Simplex<shape_dim, O, R> type;
    };


    constexpr Simplex() : Simplex( Dim, Order, RDim ) {}
    constexpr Simplex( uint16_type O ) : Convex( Dim, O, RDim ) {}
 
    static constexpr uint16_type numberOfVertices() { return numVertices; }
    static constexpr uint16_type numberOfEdges() { return numEdges; }
    static constexpr uint16_type numberOfFaces() { return numFaces; }
    static constexpr uint16_type numberOfTopologicalFaces() { return numTopologicalFaces; }
    static constexpr uint16_type numberOfVolumes() { return numVolumes; }

    /**
     * \return the number of polynomials of total degree \c n on the
     * shape:
     *
     * -# (n+1) over the line
     * -# (n+1)(n+2)/2 over the triangle,
     * -# (n+1)(n+2)(n+3)/6 in three dimensions.
     */
    static constexpr uint32_type polyDims( int n )
    {
#if 0
        if ( nDim == 1 )
            return std::max( 0, n + 1 );

        if ( nDim == 2 )
            return std::max( 0, ( n+1 )*( n+2 )/2 );

        if ( nDim == 3 )
            return std::max( 0, ( n+1 )*( n+2 )*( n+3 )/6 );

        BOOST_STATIC_ASSERT( nDim == 1 || nDim == 2 || nDim == 3 );
        return uint32_type( -1 );
#else
        return numberOfPoints( n );
#endif
    }

    /**
     * Given an edge \p e in the element and the local index \p p (0
     * or 1) of a point in the edge \p e , \return the index in the
     * element of the point.
     */
    constexpr uint16_type e2p( uint16_type e,  uint16_type p )
    {
        return edge_to_point_t::e2p( order(), e, p );
    }

    /**
     * Given a face \p f in the element and the local index \p e of an
     * edge in the face \p f, \return the index in the element of the
     * edge.
     */
    constexpr uint16_type f2e( uint16_type f,  uint16_type e )
    {
        return face_to_edge_t::f2e( f, e );
    }

    /**
     * Given a face \p f in the element and the global index \p e of an
     * edge in the face \p f, \return the local index in the element of the
     * edge.
     */
    constexpr uint16_type f2eLoc( uint16_type f,  uint16_type e )
    {
        return face_to_edge_t::f2eLoc( f, e );
    }

    /**
     * Given a face \p f in the element and the local index \p p of a
     * point in the face \p f , \return the index in the element of
     * the point.
     */
    constexpr uint16_type f2p( uint16_type f,  uint16_type p )
    {
        return face_to_point_t::f2p( order(), f, p );
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
             << nRealDim;
        return ostr.str();
    }
    static std::string type()
    {
        return "simplex";
    }
};

template<uint16_type Dim, uint16_type Order, uint16_type RDim >
const uint16_type Simplex<Dim, Order, RDim>::topological_dimension;

template<int RDim> using Line = Simplex<1, 1, RDim>;
template<int RDim> using Triangle = Simplex<2, 1, RDim>;
template<int RDim> using Tetrahedron = Simplex<3, 1>;


template<uint16_type D, uint16_type O, uint16_type R>
inline std::ostream& 
operator<<( std::ostream& os, Feel::Simplex<D,O,R> const& s )
{
    os << "Simplez<" << s.topologicalDimension() << "," << s.order() << "," << s.realDimension() << ">\n"
       << " - Topology:\n"
       << "   - number of vertices : " << s.numberOfVertices() << "\n"
       << "   - number of edges : " << s.numberOfEdges() << "\n"
       << "   - number of faces : " << s.numberOfFaces() << "\n"
       << "   - number of topological faces : " << s.numberOfTopologicalFaces() << "\n"
       << "   - number of volumes : " << s.numberOfVolumes() << "\n"
       << " - Points:\n"
       << "   - number of points per vertex : " << s.numberOfPointsPerVertex() << "\n"
       << "   - number of points per edge : " << s.numberOfPointsPerEdge() << "\n"
       << "   - number of points per faec : " << s.numberOfPointsPerFace() << "\n"
       << "   - number of points per volume : " << s.numberOfPointsPerVolume() << "\n"
       << "   - number of points : " << s.numberOfPoints() << "\n"
       << std::endl;
    return os;
}
}

#endif /* __Simplex_H */
