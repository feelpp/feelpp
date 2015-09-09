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


namespace Feel
{
class HypercubeBase {};

template<uint16_type Dim, uint16_type RDim = Dim>
class Hypercube  : public Convex, HypercubeBase
{
    using super = Convex;
    using super2 = HypercubeBase;
    
    typedef mpl::vector_c<size_type, SHAPE_POINT, SHAPE_LINE, SHAPE_QUAD, SHAPE_HEXA, SHAPE_SP4, SHAPE_SP5> shapes_t;
    typedef mpl::vector_c<size_type, GEOMETRY_POINT, GEOMETRY_LINE, GEOMETRY_SURFACE, GEOMETRY_VOLUME, GEOMETRY_4, GEOMETRY_5> geometries_t;

    typedef mpl::vector_c<size_type, 1, 2, 4, 8, 16, 32> vertices_t;
    typedef mpl::vector_c<size_type, 0, 1, 4, 12, 32> edges_t;
    typedef mpl::vector_c<size_type, 0, 2, 4, 6, 24> faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 1, 6, 24> geo_faces_index_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 4, 6, 24> normals_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 0, 1, 8> volumes_t;


    template<uint16_type rdim>
    struct faces_t
    {
        typedef mpl::vector<boost::none_t,
                Hypercube<0, rdim>,
                Hypercube<1, rdim>,
                Hypercube<2, rdim> > type;
    };
    typedef mpl::vector<boost::none_t,Hypercube<1,1>, Hypercube<1, 2>, Hypercube<1, 3>, boost::none_t > v_edges_t;
    typedef mpl::vector<Hypercube<1>, Hypercube<2>, Hypercube<3>, Hypercube<4>, Hypercube<4>  > elements_t;

    typedef mpl::vector_c<uint16_type, 0, 1, 2, 8> permutations_t;

public:

    static const bool is_simplex = false;
    static const bool is_hypercube = true;

    static const size_type Shape = mpl::at<shapes_t, mpl::int_<Dim> >::type::value;
    static const size_type Geometry = mpl::at<geometries_t, mpl::int_<Dim> >::type::value;

    static const uint16_type nDim = Dim;
    //static const uint16_type nOrder = Order;
    static const uint16_type nRealDim = RDim;

    static const uint16_type topological_dimension = nDim;
    static const uint16_type real_dimension = RDim;

    typedef typename mpl::at<elements_t, mpl::int_<nDim> >::type element_type;
    typedef typename mpl::at<typename faces_t<real_dimension>::type, mpl::int_<nDim> >::type topological_face_type;
    typedef typename mpl::at<v_edges_t, mpl::int_<real_dimension> >::type edge_type;
    
    static const uint16_type numVertices = mpl::at<vertices_t, mpl::int_<Dim> >::type::value;
    static const uint16_type numEdges = mpl::at<edges_t, mpl::int_<Dim> >::type::value;
    static const uint16_type numFaces = mpl::at<geo_faces_index_t, mpl::int_<Dim> >::type::value;
    static const uint16_type numGeometricFaces = mpl::at<geo_faces_index_t, mpl::int_<nDim> >::type::value;
    static const uint16_type numTopologicalFaces = mpl::at<faces_index_t, mpl::int_<nDim> >::type::value;
    static const uint16_type numNormals = mpl::at<normals_t, mpl::int_<nDim> >::type::value;
    static const uint16_type numVolumes = mpl::at<volumes_t, mpl::int_<nDim> >::type::value;

#if 0
    static const uint16_type nbPtsPerVertex = ( nOrder==0 )?0:1;
    static const uint16_type nbPtsPerEdge = ( nOrder==0 )?( ( nDim==1 )?1:0 ):mpl::at<points_edge_t, mpl::int_<nDim> >::type::value;
    static const uint16_type nbPtsPerFace = ( nOrder==0 )?( ( nDim==2 )?1:0 ):mpl::at<points_face_t, mpl::int_<nDim> >::type::value;
    static const uint16_type nbPtsPerVolume = ( nOrder==0 )?( ( nDim==3 )?1:0 ):mpl::at<points_volume_t, mpl::int_<nDim> >::type::value;
    static const uint16_type numPoints = ( numVertices * nbPtsPerVertex +
                                           numEdges * nbPtsPerEdge +
                                           numFaces * nbPtsPerFace +
                                           numVolumes * nbPtsPerVolume );
#else
    static constexpr uint16_type numberOfVertices()  { return numVertices; }
    static constexpr uint16_type numberOfEdges()  { return numEdges; }
    static constexpr uint16_type numberOfFaces()  { return numGeometricFaces; }
    static constexpr uint16_type numberOfTopologicalFaces()   { return numTopologicalFaces; }
    static constexpr uint16_type numberOfVolumes()  { return numVolumes; }
    
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

#endif

    typedef mpl::vector<boost::none_t, details::line, details::quad, details::hexa > map_entity_to_point_t;
    
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

    template<uint16_type shape_dim, uint16_type R=nDim>
    using shape =  Hypercube<shape_dim, R>;

    constexpr Hypercube() : Convex( Dim, 1, RDim ) {}
    constexpr Hypercube( uint16_type O ) : Convex( Dim, O, RDim ) {}
 

#if 0
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
#endif
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
    static constexpr uint32_type polyDims( int n )
    {
        return uint32_type( math::pow( double( n+1 ), double( nDim ) ) );
    }

    /**
     * Given an edge \p e in the element and the local index \p p (0
     * or 1) of a point in the edge \p e , \return the index in the
     * element of the point.
     */
    constexpr int e2p( int e,  int p ) const
    {
        return edge_to_point_t::e2p( this->order(), e, p );
    }

    /**
     * Given a face \p f in the element and the local index \p e of an
     * edge in the face \p f, \return the index in the element of the
     * edge.
     */
    constexpr int f2e( int f,  int e ) const
    {
        return face_to_edge_t::f2e( f, e );
    }

    /**
     * Given a face \p f in the element and the local index \p p of a
     * point in the face \p f , \return the index in the element of
     * the point.
     */
    constexpr int f2p( int f,  int p ) const
    {
        return face_to_point_t::f2p( this->order(), f, p );
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
