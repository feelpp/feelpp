/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-03-06

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
   \file gausslobatto.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-03-06
 */
#ifndef __GaussLobatto_H
#define __GaussLobatto_H 1

#include <feel/feelmesh/refentity.hpp>

#include <feel/feelcore/visitor.hpp>

#include <stdexcept>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


#include <feel/feelcore/traits.hpp>

#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/quadpoint.hpp>
#include <feel/feelpoly/pointsetinterpolation.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelalg/glas.hpp>



namespace Feel
{
namespace ublas = boost::numeric::ublas;


template<class Convex, uint16_type Integration_Degree, typename T> class PointSetQuadrature;
template<int Dim, int Order, int RealDim, template<uint16_type,uint16_type,uint16_type> class Entity, typename T> struct GT_Lagrange;

/*!
 * \class GaussLobatto
 * \brief Gauss quadrature points
 *
 * \code
 * // generate a Gauss point set that would integrate exactly linear
 * // functions using double precision numerical type
 * Gauss<Simplex,1,double> gauss;
 * \endcode
 *
 * \ingroup Polynomial
 * @author Gilles Steiner
 * @author Christophe Prud'homme
 */
template<class Convex, uint16_type Integration_Degree, typename T>
class GaussLobatto : public PointSetQuadrature<Convex, Integration_Degree, T>  {};

/// \cond detail
template< uint16_type Integration_Degree, typename T>
class GaussLobatto<Simplex<1,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<1,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<1,1> , Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const uint16_type Degree = ( Integration_Degree+3 )/2+1;
    static const uint32_type Npoints = Degree;

    GaussLobatto()
        :
        super( Npoints )
    {
        ublas::vector<T> px( Npoints );

        details::gausslobattojacobi<Npoints, T, ublas::vector<T>, ublas::vector<T> >( this->M_w, px );
        ublas::row( this->M_points, 0 ) = px;
    }

    ~GaussLobatto() {}

    FEELPP_DEFINE_VISITABLE();
};


/** Gauss-Lobatto x Left-Radau Quadrature on a triangle **/

template< uint16_type Integration_Degree, typename T>
class GaussLobatto<Simplex<2,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<2,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<2,1> , Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    typedef GaussLobatto<Simplex<1,1>,Integration_Degree, T> face_quad_type;


    static const uint16_type DegreeX = ( Integration_Degree+3 )/2+1;
    static const uint16_type DegreeY = ( Integration_Degree+2 )/2+1;
    static const uint32_type Npoints = DegreeX*DegreeY;

    GaussLobatto()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( DegreeX );
        weights_type px( DegreeX );
        details::gausslobattojacobi<DegreeX,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        weights_type wy( DegreeY );
        weights_type py( DegreeY );
        details::left_gaussradaujacobi<DegreeY,T, ublas::vector<T>, ublas::vector<T> >( wy, py, 1.0, 0.0 );

        // coordinate in cartesian space
        node_type eta( 2 );
        details::xi<TRIANGLE, value_type> to_xi;

        for ( int i = 0,  k = 0; i < DegreeX; ++i )
        {
            for ( int j = 0; j < DegreeY; ++j, ++k )
            {
                // computes the weight of the k-th node
                this->M_w( k ) = 0.5 * wx( i ) * wy( j );
                // use expansion for the collapsed triangle to compute the points
                // coordinates (from cartesian to collapsed coordinates)
                eta( 0 ) = px( i );
                eta( 1 ) = py( j );
                ublas::column( this->M_points, k ) = to_xi( eta );
            }
        }

        boost::shared_ptr<GT_Lagrange<2,1, 2 ,Simplex, T> > gm( new GT_Lagrange<2, 1, 2, Simplex,T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<2, 1, 2>,2,1>(), gm, face_qr );

    }

    ~GaussLobatto() {}

    FEELPP_DEFINE_VISITABLE();
};

/** Gauss-Lobatto x Left-Radau x Left-Radau Quadrature on a tetrahedra **/

template< uint16_type Integration_Degree, typename T>
class GaussLobatto<Simplex<3,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<3,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<3,1> , Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    typedef GaussLobatto<Simplex<2,1>,Integration_Degree, T> face_quad_type;


    static const uint16_type DegreeX = ( Integration_Degree+3 )/2+1;
    static const uint16_type DegreeY = ( Integration_Degree+2 )/2+1;
    static const uint16_type DegreeZ = ( Integration_Degree+2 )/2+1;
    static const uint32_type Npoints = DegreeX*DegreeY*DegreeZ;

    GaussLobatto()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( DegreeX );
        weights_type px( DegreeX );
        details::gausslobattojacobi<DegreeX,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        weights_type wy( DegreeY );
        weights_type py( DegreeY );
        details::left_gaussradaujacobi<DegreeY,T, ublas::vector<T>, ublas::vector<T> >( wy, py, 1.0, 0.0 );

        weights_type wz( DegreeZ );
        weights_type pz( DegreeZ );
        details::left_gaussradaujacobi<DegreeZ,T, ublas::vector<T>, ublas::vector<T> >( wz, pz, 2.0, 0.0 );

        // coordinate in cartesian space
        node_type eta( 3 );
        details::xi<TETRAHEDRON, value_type> to_xi;

        for ( int i = 0,  k = 0; i < DegreeX; ++i )
        {
            for ( int j = 0; j < DegreeY; ++j )
            {
                for ( int l = 0; l < DegreeZ; ++l, ++k )
                {
                    // computes the weight of the k-th node
                    this->M_w( k ) = 0.125 * wx( i ) * wy( j ) * wz( l );
                    // use expansion for the collapsed triangle to compute the points
                    // coordinates (from cartesian to collapsed coordinates)
                    eta( 0 ) = px( i );
                    eta( 1 ) = py( j );
                    eta( 2 ) = pz( l );
                    ublas::column( this->M_points, k ) = to_xi( eta );
                }
            }
        }

        boost::shared_ptr<GT_Lagrange<3, 1, 3, Simplex, T> > gm( new GT_Lagrange<3, 1, 3, Simplex, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<3, 1, 3>,3,1>(), gm, face_qr );
    }

    ~GaussLobatto() {}

    FEELPP_DEFINE_VISITABLE();
};




/** GaussLobatto Quadrature on Simplex Product **/

/** GaussLobatto Quadrature on the quadrangle [-1,1]x[-1,1] **/

template< uint16_type Integration_Degree, typename T>
class GaussLobatto<Hypercube<2,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<2,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<2,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;
    typedef GaussLobatto<Hypercube<1,1>,Integration_Degree, T> face_quad_type;
    static const uint16_type Degree = ( Integration_Degree+3 )/2+1;
    static const uint32_type Npoints = Degree*Degree;

    GaussLobatto()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gausslobattojacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j, ++k )
            {
                // computes the weight of the k-th node
                this->M_w( k ) = wx( i ) * wx( j );
                this->M_points( 0, k ) = px( i );
                this->M_points( 1, k ) = px( j );
            }
        }

        boost::shared_ptr<GT_Lagrange<2, 1, 2, Hypercube, T> > gm( new GT_Lagrange<2, 1, 2, Hypercube, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<2, 1, 2>,2,1>(), gm, face_qr );
    }

    ~GaussLobatto() {}

    FEELPP_DEFINE_VISITABLE();
};

/** GaussLobatto Quadrature on the hexahedra [-1,1]x[-1,1]x[-1,1] **/

template< uint16_type Integration_Degree, typename T>
class GaussLobatto<Hypercube<3,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<3,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<3,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;
    typedef GaussLobatto<Hypercube<2,1>,Integration_Degree, T> face_quad_type;
    static const uint16_type Degree = ( Integration_Degree+3 )/2+1;
    static const uint32_type Npoints = Degree*Degree*Degree;

    GaussLobatto()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gausslobattojacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j )
            {
                for ( int l = 0; l < Degree ; ++l, ++k )
                {
                    // computes the weight of the k-th node
                    this->M_w( k ) = wx( i ) * wx( j ) * wx( l );
                    this->M_points( 0, k ) = px( i );
                    this->M_points( 1, k ) = px( j );
                    this->M_points( 2, k ) = px( l );
                }
            }
        }

        boost::shared_ptr<GT_Lagrange<3, 1, 3, Hypercube, T> > gm( new GT_Lagrange<3, 1, 3, Hypercube, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<3, 1, 3>,3,1>(), gm, face_qr );
    }

    ~GaussLobatto() {}
    FEELPP_DEFINE_VISITABLE();
};
/// \endcond


/**
 * \class GaussLobatto
 * \brief Gauss Lobatto pointset
 *
 * \ingroup Polynomial
 * @author Goncalo Pena
 * @see
 */
template< class Convex,
          uint16_type Order,
          typename T = double >
class PointSetGaussLobatto : public PointSetInterpolation<Convex::nDim, Order, T, Hypercube>
{
public:

    typedef PointSetInterpolation< Convex::nDim, Order, T, Hypercube> super;

    typedef typename super::return_type return_type;

    typedef T value_type;

    static const uint32_type Dim = Convex::nDim;
    static const uint32_type nPoints = Order+1;
    static const uint32_type topological_dimension = Convex::topological_dimension;
    static const uint32_type nRealDim = Convex::nRealDim;

    static const size_type Shape = Convex::Shape;

    static const bool is_simplex = Convex::is_simplex;
    static const bool is_hypercube = Convex::is_hypercube;

    typedef Reference<Convex, Dim, Convex::nOrder, Convex::nDim/*Convex::nRealDim*/, value_type> reference_convex_type;

    typedef ublas::vector<value_type> vector_type;

    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;

    typedef typename reference_convex_type::points_type points_type;

    typedef typename Convex::edge_to_point_t edge_to_point_t;
    typedef typename Convex::face_to_point_t face_to_point_t;
    typedef typename Convex::face_to_edge_t face_to_edge_t;

    typedef Hypercube<Dim, Order, Dim> conv_order_type;
    static const uint32_type numPoints = conv_order_type::numPoints;

    reference_convex_type RefConv;

    typedef typename super::range_type range_type;
    typedef typename super::index_map_type index_map_type;

    PointSetGaussLobatto( int interior = 0 )
    {
        FEELPP_ASSERT( is_hypercube || ( Dim == 1 )  ).error( "gauss lobatto points are just defined in simplex products" );

        nodes_type pts( Dim, numPoints );

        calculate_gl_points( mpl::bool_< ( Order > 1 )>() );

        if ( interior == 0 && Order > 0 )
        {
            for ( uint16_type d = 0, p = 0; d < topological_dimension+1; ++d )
            {
                for ( int e = RefConv.entityRange( d ).begin();
                        e < RefConv.entityRange( d ).end();
                        ++e )
                {
                    nodes_type Gt ( makePoints( d, e ) );

                    if ( Gt.size2() )
                    {
                        ublas::subrange( pts, 0, Dim, p, p+Gt.size2() ) = Gt;

                        for ( size_type j = 0; j < Gt.size2(); ++j )
                        {
                            this->addToEid( d, p+j );
                            this->addToPtE( p+j, std::make_pair( d, e ) );
                        }

                        p+=Gt.size2();
                    }
                }
            }

            this->setPoints( pts );
        }

        else if ( interior == 1 && Order > 0 )
            this->setPoints( makePoints( Dim, 0 ) );

        else if ( Order == 0 )
            this->setPoints( glas::average( RefConv.vertices() ) );

        this->setName( "gausslobatto", Order );

        //std::cout << "end constructor...\n";
    }

    ~PointSetGaussLobatto() {}

private:

    vector_type gl_pts;

    void calculate_gl_points ( mpl::bool_<true> )
    {
        vector_type wx( Order+1 );
        vector_type px( Order+1 );

        details::gausslobattojacobi<Order+1, T, vector_type, vector_type >( wx, px );

        ublas::vector_range<vector_type> inner_pts ( px, ublas::range( 1, px.size()-1 ) );

        gl_pts = inner_pts;
    }

    void calculate_gl_points ( mpl::bool_<false> )
    {}

    points_type makePoints( uint16_type topo_dim, uint16_type __id )
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
                return makeLattice<Shape>();

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
                G = makeLattice<SHAPE_LINE>();
                Gret.resize( nRealDim, G.size2() );

                pt_to_entity_hexahedron<Shape, 1> p_to_e( __id );

                for ( size_type i = 0; i < G.size2(); ++i )
                    ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );

                return Gret;
            }

            else if ( topo_dim == 2 )
            {
                G = makeLattice<SHAPE_QUAD>();
                Gret.resize( nRealDim, G.size2() );
                pt_to_entity_hexahedron<Shape, 2> p_to_e( __id );

                for ( size_type i = 0; i < G.size2(); ++i )
                    ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );

                return Gret;
            }
        }

        return points_type();
    }

    template<size_type shape>
    points_type makeLattice()
    {
        points_type G;

        if ( Order > 0 )
        {
            if ( shape == SHAPE_LINE )
                G = make_line_points();

            else if ( shape == SHAPE_QUAD )
                return make_quad_points();

            else if ( shape == SHAPE_HEXA )
                return make_hexa_points();
        }

        else if ( Order == 0 )
            G = glas::average( RefConv.vertices() );

        return G;
    }

    int n_line_points()
    {
        return std::max( 0, int( Order ) - 1 );
    }

    int n_quad_points() const
    {
        return std::max( 0, ( int( Order ) - 1 )*( int( Order ) - 1 ) );

    }

    int n_hexa_points() const
    {
        return std::max( 0, ( int( Order )-1 )*( int( Order )-1 )*( int( Order )-1 ) );
    }

    points_type
    make_line_points()
    {
        points_type G ( Dim, n_line_points() );

        if ( Order > 1 )
        {
            vector_type ones ( ublas::scalar_vector<value_type>( G.size2(), value_type( 1 ) ) );

            ublas::row( G, 0 ) = gl_pts;

            for ( uint16_type i=1; i<Dim; i++ )
                ublas::row( G, i ) = - ones;
        }

        return G;
    }

    points_type
    make_quad_points()
    {
        points_type G( Dim, n_quad_points() );

        if ( Order > 1 )
        {
            uint16_type numInterior = Order - 1;

            for ( uint16_type i = 0,  k = 0; i < numInterior; ++i )
            {
                for ( uint16_type j = 0; j < numInterior; ++j, ++k )
                {
                    G( 0, k ) = gl_pts( i );
                    G( 1, k ) = gl_pts( j );
                }
            }

            vector_type ones ( ublas::scalar_vector<value_type>( G.size2(), value_type( 1 ) ) );

            if ( Dim == 3 )
                ublas::row( G, 2 ) = - ones;
        }

        return G;
    }

    points_type
    make_hexa_points()
    {
        points_type G( Dim, n_hexa_points() );

        if ( Order > 1 )
        {
            uint16_type numInterior = Order - 1;

            for ( uint16_type i = 0,  k = 0; i < numInterior; ++i )
            {
                for ( uint16_type j = 0; j < numInterior; ++j )
                {
                    for ( uint16_type l = 0; l < numInterior; ++l, ++k )
                    {
                        G( 0, k ) = gl_pts( i );
                        G( 1, k ) = gl_pts( j );
                        G( 2, k ) = gl_pts( l );
                    }
                }
            }
        }

        return G;
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
#endif /* __GaussLobatto_H */
