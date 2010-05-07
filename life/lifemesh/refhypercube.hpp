/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-05-06

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file refhypercube.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-05-06
 */
#ifndef __refhypercube_H
#define __refhypercube_H 1

namespace Life
{
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
class Reference<SimplexProduct<Dim, Order, RDim>, Dim, Order, RDim, T>
    :
    public SimplexProduct<Dim, Order, RDim>
{
public:


    /** @name Typedefs
     */
    //@{

    typedef SimplexProduct<Dim, Order, RDim> super;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    static const uint16_type topological_dimension = super::topological_dimension;
    static const uint16_type real_dimension = super::real_dimension;

    static const size_type Shape = super::Shape;
    static const size_type Geometry = super::Geometry;

    typedef T value_type;
    typedef Reference<SimplexProduct<Dim, Order, RDim>, Dim, Order, RDim, T> self_type;

    typedef typename mpl::if_<boost::is_same<typename super::element_type, boost::none_t>,
                              mpl::identity<boost::none_t>,
                              mpl::identity<Reference<typename super::element_type, nDim+1, nOrder, nRealDim, T> > >::type::type element_type;
    typedef Reference<typename super::topological_face_type, super::topological_face_type::nDim, nOrder, nRealDim, T> topological_face_type;

    typedef typename super::edge_to_point_t edge_to_point_t;
    typedef typename super::face_to_point_t face_to_point_t;
    typedef typename super::face_to_edge_t face_to_edge_t;

    static const uint16_type numVertices = super::numVertices;
    static const uint16_type numFaces = super::numFaces;
    static const uint16_type numGeometricFaces = super::numGeometricFaces;
    static const uint16_type numTopologicalFaces = super::numTopologicalFaces;
    static const uint16_type numEdges = super::numEdges;
    static const uint16_type numNormals = super::numNormals;

    static const uint16_type numPoints = super::numPoints;
    static const uint16_type nbPtsPerVertex = super::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = super::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = super::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = super::nbPtsPerVolume;

    typedef typename node<value_type>::type node_type;
    typedef typename matrix_node<value_type>::type points_type;

    typedef typename node<value_type>::type normal_type;
    typedef ublas::vector<normal_type> normals_type;
    typedef typename normals_type::const_iterator normal_const_iterator;
    typedef typename super::permutation_type permutation_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Reference()
        :
        super(),
        _M_id( 0 ),
        _M_vertices( nDim, numVertices ),
        _M_points( nDim, numPoints ),
        _M_normals( numNormals )
        {
            if ( nDim == 1 )
            {
                _M_vertices( 0, 0 ) = -1.0;
                _M_vertices( 0, 1 ) =  1.0;

                _M_points = make_line_points();
            }
            if ( nDim == 2 )
            {
                _M_vertices( 0, 0 ) = -1.0;
                _M_vertices( 1, 0 ) = -1.0;

                _M_vertices( 0, 1 ) =  1.0;
                _M_vertices( 1, 1 ) = -1.0;

                _M_vertices( 0, 2 ) =  1.0;
                _M_vertices( 1, 2 ) =  1.0;

                _M_vertices( 0, 3 ) = -1.0;
                _M_vertices( 1, 3 ) =  1.0;

                _M_points = make_quad_points();
            }
            if ( nDim == 3 )
            {
                // z = -1
                _M_vertices( 0, 0 ) = -1.0;
                _M_vertices( 1, 0 ) = -1.0;
                _M_vertices( 2, 0 ) = -1.0;

                _M_vertices( 0, 1 ) =  1.0;
                _M_vertices( 1, 1 ) = -1.0;
                _M_vertices( 2, 1 ) = -1.0;

                _M_vertices( 0, 2 ) =  1.0;
                _M_vertices( 1, 2 ) =  1.0;
                _M_vertices( 2, 2 ) = -1.0;

                _M_vertices( 0, 3 ) = -1.0;
                _M_vertices( 1, 3 ) =  1.0;
                _M_vertices( 2, 3 ) = -1.0;

                // z = 1
                _M_vertices( 0, 4 ) = -1.0;
                _M_vertices( 1, 4 ) = -1.0;
                _M_vertices( 2, 4 ) =  1.0;

                _M_vertices( 0, 5 ) =  1.0;
                _M_vertices( 1, 5 ) = -1.0;
                _M_vertices( 2, 5 ) =  1.0;

                _M_vertices( 0, 6 ) =  1.0;
                _M_vertices( 1, 6 ) =  1.0;
                _M_vertices( 2, 6 ) =  1.0;

                _M_vertices( 0, 7 ) = -1.0;
                _M_vertices( 1, 7 ) =  1.0;
                _M_vertices( 2, 7 ) =  1.0;

                _M_points = make_hexa_points();
            }
            //std::cout << "P = " << _M_points << "\n";
            make_normals();
        }

    Reference( element_type const& e, uint16_type __f )
        :
        super(),
        _M_id( __f ),
        _M_vertices( nRealDim, numVertices ),
        _M_points( nRealDim, numPoints ),
        _M_normals( numNormals )
        {
          if ( __f >= element_type::numTopologicalFaces )
            {
                std::ostringstream str;
                str << "invalid face number " << __f << "\n"
                    << "must be 0 <= f < " << element_type::numTopologicalFaces << "\n";
                throw std::invalid_argument( str.str() );
            }

            for ( int i = 0;i < numVertices; ++i )
            {
                if ( real_dimension == 3 )
                    ublas::column( _M_vertices, i ) = e.vertex( element_type::f2p( __f, i ) );
                else
                    ublas::column( _M_vertices, i ) = e.vertex( element_type::e2p( __f, i ) );
            }
            for ( int i = 0;i < numPoints; ++i )
            {
                if ( real_dimension == 3 )
                    ublas::column( _M_points, i ) = e.point( element_type::f2p( __f, i ) );
                else
                    ublas::column( _M_points, i ) = e.point( element_type::e2p( __f, i ) );
            }
            _M_points = _M_vertices;
            make_normals();
        }

    Reference( Reference const & r )
        :
        super( r ),
        _M_id( r._M_id ),
        _M_vertices( r._M_vertices ),
        _M_points( r._M_points ),
        _M_normals( r._M_normals )
        {

        }

    ~Reference() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    Reference& operator=( Reference const& r )
        {
            if ( this != &r )
            {
                _M_id = r._M_id;
                _M_vertices = r._M_vertices;
                _M_points = r._M_points;
                _M_normals = r._M_normals;
            }
            return *this;
        }

    //@}

    /** @name Accessors
     */
    //@{

    uint16_type topologicalDimension() const { return topological_dimension; }
    uint16_type dimension() const { return real_dimension; }

    uint16_type nVertices() const { return numVertices; }
    uint16_type nPoints() const { return numPoints; }
    uint16_type nEdges() const { return numEdges; }
    uint16_type nFaces() const { return numFaces; }

    points_type const& vertices() const { return _M_vertices; }

    ublas::matrix_column<points_type const> vertex( uint16_type __i ) const { return ublas::column( _M_vertices, __i ); }

    ublas::matrix_column<points_type const> edgeVertex( uint16_type __e, uint16_type __p ) const
        {
            return ublas::column( _M_vertices, edge_to_point_t::e2p(__e,__p) );
        }

    points_type const& points() const { return _M_points; }

    ublas::matrix_column<points_type const> point( uint16_type __i ) const { return ublas::column( _M_points, __i ); }

    /**
     * get the normals array
     *
     *
     * @return the normals
     */
    normals_type const& normals() const { return _M_normals; }

    /**
     * get the n-th normal
     *
     * @param __n the index of the normal
     *
     * @return the n-th normal of the quad
     */
    node_type const& normal( uint16_type __n ) const { return _M_normals[__n]; }

    /**
     * the first iterator of the normal vector
     *
     *
     * @return the begin() iterator of the normal vector
     */
    normal_const_iterator beginNormal() const { return _M_normals.begin(); }

    /**
     * the end() iterator
     *
     *
     * @return
     */
    normal_const_iterator endNormal() const { return _M_normals.end(); }

    topological_face_type topologicalFace( uint16_type __f ) const
        {
            topological_face_type ref( *this, __f );
            return ref;
        }

    points_type const& G() const { return _M_points; }

    size_type id() const { return 0; }
    size_type marker() const { return 0; }
    flag_type marker2() const { return 0; }
    flag_type marker3() const { return 0; }
    uint16_type ad_first() const { return -1; }
    uint16_type ad_second() const { return -1; }
    uint16_type pos_first() const { return -1; }
    uint16_type pos_second() const { return -1; }
    permutation_type permutation( uint16_type /*f*/ ) const { return permutation_type(); }
    double h() const
        {
            // FIXME: should be computed once for all in constructor
            double __max = 0.0;
            for ( int __e = 0; __e < numEdges;++ __e )
            {
                double __len = ublas::norm_2( edgeVertex( __e, 1 ) - edgeVertex( __e, 0 ) );
                __max = ( __max < __len )?__len:__max;
            }
            return __max;
        }


    double hFace( int /*__f*/ ) const
        {
#if 0
            // FIXME: should be computed once for all in constructor
            if ( nDim == 1 )
                return 0.0;
            else
                return face( __f ).h();
#else
            return 0.0;
#endif
        }

    EntityRange<self_type> entityRange( uint16_type d ) const { return EntityRange<self_type>( d ); }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    /**
     *
     * \return a positive or null number if pt is in the convex
     */
    boost::tuple<bool, value_type>
    isIn( typename node<value_type>::type const& pt ) const
        {
            return isIn( pt, mpl::int_<nDim>() );
        }

    boost::tuple<bool, value_type>
    isIn( typename node<value_type>::type const& pt, mpl::int_<1> ) const
    {
            typename matrix_node<value_type>::type M( nDim+1,nDim+1);
            typename matrix_node<value_type>::type P( nDim,nDim+1);
            std::fill( M.data().begin(), M.data().end(), value_type( 1.0 ) );

            // compute the measure(times 2) of the segments, triangles or
            // tetras formed by the point pt and the vertices of the
            // simplex. The simplices are formed such that when the
            // point is in the simplex the measure are all
            // positive. In the case the point is outside the simplex
            // the measure is negative
            double meas_times = 1e30;
            for( int n = 0; n < numVertices; ++n )
                {
                    ublas::column( P, 0 ) = pt;
                    ublas::column( P, 1 ) = this->vertex( n );
                    ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P;
                    double res = details::det( M, mpl::int_<nDim+1>() )/2;
                    meas_times = ( meas_times < res )?meas_times:res;
                }
            return boost::make_tuple( meas_times>=-1e10, meas_times);
        }
    boost::tuple<bool, value_type>
    isIn( typename node<value_type>::type const& pt, mpl::int_<2> ) const
        {
            typename matrix_node<value_type>::type M( nDim+1,nDim+1);
            typename matrix_node<value_type>::type P( nDim,nDim+1);
            std::fill( M.data().begin(), M.data().end(), value_type( 1.0 ) );

            // compute the measure(times 2) of the segments, triangles or
            // tetras formed by the point pt and the vertices of the
            // simplex. The simplices are formed such that when the
            // point is in the simplex the measure are all
            // positive. In the case the point is outside the simplex
            // the measure is negative
            double meas_times = 1e30;
            for( int n = 0; n < numVertices; ++n )
                {
                    ublas::column( P, 0 ) = pt;
                    ublas::column( P, 1 ) = this->vertex( n );
                    ublas::column( P, 2 ) = this->vertex( (n+1)%numVertices );
                    ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P;
                    double res = details::det( M, mpl::int_<nDim+1>() )/2;
                    meas_times = ( meas_times < res )?meas_times:res;
                }
            return boost::make_tuple( meas_times>=-1e10, meas_times);
        }
    boost::tuple<bool, value_type>
    isIn( typename node<value_type>::type const& pt, mpl::int_<3> ) const
        {
            typename matrix_node<value_type>::type M( nDim+1,nDim+1);
            typename matrix_node<value_type>::type P( nDim,nDim+1);
            std::fill( M.data().begin(), M.data().end(), value_type( 1.0 ) );

            // compute the measure(times 2) of the segments, triangles or
            // tetras formed by the point pt and the vertices of the
            // simplex. The simplices are formed such that when the
            // point is in the simplex the measure are all
            // positive. In the case the point is outside the simplex
            // the measure is negative
            double meas_times = 1e30;
            for( int n = 0; n < numVertices; ++n )
                {
                    ublas::column( P, 0 ) = pt;
                    ublas::column( P, 1 ) = this->vertex( n );
                    ublas::column( P, 2 ) = this->vertex( (n+1)%numVertices );
                    ublas::column( P, 3 ) = this->vertex( (n+2)%numVertices );
                    ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P;
                    double res = details::det( M, mpl::int_<nDim+1>() )/2;
                    meas_times = ( meas_times < res )?meas_times:res;
                }
            return boost::make_tuple( meas_times>=-1e10, meas_times);
        }

    points_type makePoints( uint16_type topo_dim, uint16_type __id, int interior = 1 ) const
        {
            // vertices
            if ( topo_dim == 0 )
                {

                    //std::cerr << "coucpu 2 " << topo_dim << " " << topological_dimension << " " << __id << "\n";
                    points_type G( _M_vertices.size1(), 1 );
                    ublas::column( G, 0 ) = ublas::column( _M_vertices, __id );
                    return G;
                }
            // interior points of the convex
            else if ( topo_dim == topological_dimension )
                {
                    //std::cerr << "coucpu 2 " << topo_dim << " " << topological_dimension << " " << __id << "\n";
                    if ( __id == 0 )
                        return makeLattice<Shape>( interior );
                    throw std::logic_error( "cannot make those points" );
                    return points_type();
                }
            // all the other points
            else
                {
                    points_type G;
                    if ( topo_dim == 1 )
                        {
                            G = Reference<SimplexProduct<1, Order, 1>, 1, Order, 1, T>().makeLattice<SHAPE_LINE>( interior );
                            pt_to_entity<Shape,1> p_to_e( __id );
                            points_type Gret( nRealDim, G.size2() );
                            for ( size_type i = 0; i < G.size2(); ++i )
                                ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );

                            return Gret;
                        }
                    else if ( topo_dim == 2 )
                        {
                            G = Reference<SimplexProduct<2, Order, 2>, 2, Order, 2, T>().makeLattice<SHAPE_QUAD>( interior );
                            pt_to_entity<Shape,2> p_to_e( __id );
                            points_type Gret( nRealDim, G.size2() );
                            for ( size_type i = 0; i < G.size2(); ++i )
                                ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );

                            return Gret;
                        }
                }
            return points_type();
        }

    template<size_type shape>
    points_type
    makeLattice( uint16_type interior = 0 ) const
        {
            if ( nOrder > 0 )
                {
                    if ( shape == SHAPE_LINE )
                        return make_line_points( interior );
                    else if ( shape == SHAPE_QUAD )
                        return make_quad_points( interior );
                    else if ( shape == SHAPE_HEXA )
                        return make_hexa_points( interior );
                }
            else if ( nOrder == 0 )
                return glas::average( _M_vertices );
        }
    //@}

private:

    int n_line_points( int interior = 0 ) const
        {
            return std::max( 0, int( Order )+1-2*interior );
        }
    int n_quad_points( int interior = 0 ) const
        {
            if ( interior == 1 )
                return std::max( 0, ( int( Order )+1-2*interior )*( int( Order )+1-2*interior ) );
            return ( Order+1 )*( Order+1 );
        }
    int n_hexa_points( int interior = 0 ) const
        {
            if ( interior == 1 )
                return std::max( 0, ( int( Order )+1-2*interior )*( int( Order )+1-2*interior )*( int( Order )+1-2*interior ) );
            return ( Order+1 )*( Order+1 )*( Order+1 );
        }


    points_type
    make_line_points( int interior = 0 ) const
    {
        if ( nOrder > 0 )
            {
                ublas::vector<node_type> h ( 1 );
                h( 0 ) = vertex( 1 ) - vertex( 0 );

                points_type p( 1, n_line_points( interior ) );
                for ( int i = interior, indp = 0; i < int( Order )+1-interior;++i, ++indp )
                    {
                        ublas::column( p, indp ) = vertex( 0 ) + ( h(0) * value_type(i) )/value_type(Order);
                    }
                return p;
            }
        else
            return glas::average( _M_vertices );
    }

    points_type
    make_quad_points( int interior = 0 ) const
        {
            if ( nOrder > 0 )
                {
                    ublas::vector<node_type> h ( 2 );
                    h( 0 ) = vertex( 1 ) - vertex( 0 );
                    h( 1 ) = vertex( 3 ) - vertex( 0 );
                    //Debug( 4005 ) << "h = " << h << "\n";
                    //Debug( 4005 ) << "n quad pts = " << n_quad_points( interior ) << "\n";
                    points_type G( 2, n_quad_points( interior ) );
                    for ( int i = interior, p = 0; i < int( Order )+1-interior;++i )
                        {
                            for ( int j = interior; j < int( Order ) + 1 -interior;++j, ++p )
                                {
                                    ublas::column( G, p ) = vertex( 0 ) + ( value_type(i) * h( 1 )  +
                                                                            value_type(j) * h( 0 ) )/ value_type(Order);
                                }
                        }
                    return G;
                }
            else
                return glas::average( _M_vertices );
        }

    points_type
    make_hexa_points( int interior = 0 ) const
        {
            if ( nOrder > 0 )
                {
                    ublas::vector<node_type> h ( 3 );
                    h( 0 ) = vertex( 1 ) - vertex( 0 );
                    h( 1 ) = vertex( 3 ) - vertex( 0 );
                    h( 2 ) = vertex( 4 ) - vertex( 0 );
                    points_type G( 3, n_hexa_points( interior ) );
                    //Debug( 4005 ) << "n hexa pts = " << n_hexa_points( interior ) << "\n";
                    for ( int i = interior, p = 0; i < int( Order )+1-interior;++i )
                        {
                            for ( int j = interior; j < int( Order ) + 1 - interior;++j )
                                {
                                    for ( int k = interior; k < int( Order ) + 1 - interior;++k, ++p )
                                        {
                                            ublas::column( G, p ) = vertex( 0 ) + ( value_type(i) * h( 2 ) +
                                                                                    value_type(j) * h( 1 ) +
                                                                                    value_type(k) * h( 0 ) ) / value_type(Order);

                                        }
                                }
                        }
                    //std::cout << "G = " << G << "\n";
                    return G;
                }
            else
                return glas::average( _M_vertices );
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
    struct pt_to_face
    {
        pt_to_face( std::vector<uint16_type> vert_ids )
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
                return u + 0.5*(x[ 0 ]+1.0) * diff[ 0 ] + 0.5*(x[ 1 ]+1.0) * diff[ 1 ];
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
    struct pt_to_entity
    {
        typedef typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_LINE> >,
                                  mpl::identity<mpl::vector<boost::none_t,pt_to_edge<shape>,pt_to_edge<shape> > >,
                                  typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_QUAD> >,
                                                    mpl::identity<mpl::vector<boost::none_t, pt_to_edge<shape>, pt_to_element<shape> > >,
                                                    mpl::identity<mpl::vector<boost::none_t, pt_to_edge<shape>, pt_to_face<shape>, pt_to_element<shape> > >
        >::type // 2
        >::type::type _type;
        typedef typename mpl::at<_type, mpl::int_<topo_dim> >::type mapping_type;
        typedef mpl::vector<boost::none_t, edge_to_point_t, face_to_point_t> list_v;

        pt_to_entity( uint16_type entity_id )
            :
            mapping( typename mpl::at<list_v, mpl::int_<topo_dim> >::type().entity( topo_dim, entity_id ) )
            {}

        node_type operator()( node_type const& x ) const
            {
                return mapping( x );
            }
        mapping_type mapping;
    };

    void make_normals()
        {
            _M_normals.resize( numNormals );

            for ( int n = 0; n < numNormals; ++n )
                {
                    _M_normals[n].resize( nDim );
                    _M_normals[n].clear();
                }
            if ( nDim == 1 )
                {
                    _M_normals[0][0] = -1;
                    _M_normals[1][0] =  1;
                }
            if ( nDim == 2 )
                {
                    _M_normals[0][1] = -1;
                    _M_normals[1][0] =  1;
                    _M_normals[2][1] =  1;
                    _M_normals[3][0] = -1;
                }
            if ( nDim == 3 )
                {
                    _M_normals[0][2] = -1;
                    _M_normals[1][1] = -1;
                    _M_normals[2][0] =  1;
                    _M_normals[3][1] =  1;
                    _M_normals[4][0] = -1;
                    _M_normals[5][2] =  1;
                }
        }
private:

    uint16_type _M_id;

    points_type _M_vertices;

    points_type _M_points;

    normals_type _M_normals;


    points_type M_barycenter;

    points_type M_barycenterfaces;

};

template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<SimplexProduct<Dim, Order, RDim>, Dim, Order, RDim, T>::nbPtsPerVertex;
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<SimplexProduct<Dim, Order, RDim>, Dim, Order, RDim, T>::nbPtsPerEdge;
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<SimplexProduct<Dim, Order, RDim>, Dim, Order, RDim, T>::nbPtsPerFace;
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<SimplexProduct<Dim, Order, RDim>, Dim, Order, RDim, T>::numGeometricFaces;



template<typename T> class Entity<SHAPE_QUAD, T>: public Reference<SimplexProduct<2, 1, 2>, 2, 1, 2, T> {};
template<typename T> class Entity<SHAPE_HEXA, T>: public Reference<SimplexProduct<3, 1, 3>, 3, 1, 3, T> {};

}
#endif /* __refhypercube_H */
