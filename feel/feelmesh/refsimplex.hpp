/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \file refsimplex.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-05-06
 */
#ifndef __refsimplex_H
#define __refsimplex_H 1

namespace Feel
{
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
class Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>
    :
public Simplex<Dim, Order, RDim>
{
public:
    typedef Simplex<Dim, Order, RDim> super;

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    static const uint16_type topological_dimension = super::topological_dimension;
    static const uint16_type real_dimension = super::real_dimension;

    typedef super GeoShape;
    static const size_type Shape = super::Shape;
    static const size_type Geometry = super::Geometry;

    typedef T value_type;
    typedef Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T> self_type;

    typedef typename mpl::if_<boost::is_same<typename super::element_type, boost::none_t>,
            mpl::identity<boost::none_t>,
            mpl::identity<Reference<typename super::element_type, nDim+1, nOrder, nRealDim, T> > >::type::type element_type;
    typedef Reference<typename super::topological_face_type,
            super::topological_face_type::nDim, nOrder, nRealDim, T> topological_face_type;

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
    typedef points_type matrix_node_type;

    typedef node_type normal_type;
    typedef ublas::vector<normal_type> normals_type;
    typedef typename normals_type::const_iterator normal_const_iterator;

    typedef node_type edge_tangent_type;
    typedef ublas::vector<edge_tangent_type> edge_tangents_type;
    typedef typename edge_tangents_type::const_iterator edge_tangent_const_iterator;

    typedef typename super::permutation_type permutation_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Reference()
        :
        super(),
        M_id( 0 ),
        M_vertices( nRealDim, numVertices ),
        M_points( nRealDim, numPoints ),
        M_normals( numNormals ),
        M_edge_tangents( numEdges ),
        M_barycenter( nRealDim ),
        M_barycenterfaces( nRealDim, numTopologicalFaces ),
        M_meas( 0 )
    {
        M_vertices *= 0;
        M_points *= 0;

        if ( nDim == 1 )
        {
            M_vertices( 0, 0 ) = -1.0;
            M_vertices( 0, 1 ) =  1.0;

            M_points = make_line_points();
        }

        if ( nDim == 2 )
        {
            M_vertices( 0, 0 ) = -1.0;
            M_vertices( 1, 0 ) = -1.0;
            M_vertices( 0, 1 ) =  1.0;
            M_vertices( 1, 1 ) = -1.0;
            M_vertices( 0, 2 ) = -1.0;
            M_vertices( 1, 2 ) =  1.0;

            M_points = make_triangle_points();
        }

        if ( nDim == 3 )
        {
            M_vertices( 0, 0 ) = -1.0;
            M_vertices( 1, 0 ) = -1.0;
            M_vertices( 2, 0 ) = -1.0;
            M_vertices( 0, 1 ) =  1.0;
            M_vertices( 1, 1 ) = -1.0;
            M_vertices( 2, 1 ) = -1.0;
            M_vertices( 0, 2 ) = -1.0;
            M_vertices( 1, 2 ) =  1.0;
            M_vertices( 2, 2 ) = -1.0;
            M_vertices( 0, 3 ) = -1.0;
            M_vertices( 1, 3 ) = -1.0;
            M_vertices( 2, 3 ) =  1.0;

            M_points = make_tetrahedron_points();
        }

        //std::cout << "P = " << M_points << "\n";
        make_normals();

        make_edge_tangents();

        computeBarycenters();
        computeMeasure();
    }

    Reference( element_type const& e, uint16_type __f, uint16_type __p = permutation_type::IDENTITY )
        :
        super(),
        M_id( __f ),
        M_vertices( nRealDim, numVertices ),
        M_points( nRealDim, numPoints ),
        M_normals( numNormals ),
        M_edge_tangents( numNormals ),
        M_barycenter( nRealDim ),
        M_barycenterfaces( nRealDim, numTopologicalFaces ),
        M_meas( 0 )
    {
        if ( __f >= element_type::numTopologicalFaces )
        {
            std::ostringstream str;
            str << "invalid face number " << __f << "\n"
                << "must be 0 <= f < " << element_type::numTopologicalFaces << "\n";
            throw std::invalid_argument( str.str() );
        }



        CHECK( nDim <3 ) << "nDim must be less than 3 here\n";
        if ( nDim == 2 )
        {

            std::map<uint16_type, std::vector<uint16_type> > permTriangles = {
                {triangular_faces_type::IDENTITY, {0,1,2}},
                {triangular_faces_type::ROTATION_ANTICLOCK, {1,2,0}},
                {triangular_faces_type::ROTATION_CLOCKWISE, {2,0,1}},
                {triangular_faces_type::REVERSE_BASE,{0,2,1}},
                {triangular_faces_type::REVERSE_HYPOTENUSE,{1,0,2}},
                {triangular_faces_type::REVERSE_HEIGHT,{2,1,0}}
            };

            DCHECK( permTriangles.find( __p )!=permTriangles.end() ) << "invalid permutation :" << __p << "\n";

            for ( int i = 0; i < numVertices; ++i )
            {
                const int iperm = permTriangles.find( __p )->second[ i ];
                ublas::column( M_vertices, iperm ) = e.vertex( element_type::f2p( __f, i ) );
            }

            M_points = make_triangle_points();
        }
        else if ( nDim == 1 )
        {
            std::map<uint16_type, std::vector<uint16_type> > permLines{
                {line_permutations::IDENTITY, {0,1} },
                {line_permutations::REVERSE_PERMUTATION, {1,0} } };

            DCHECK( permLines.find( __p )!=permLines.end() ) << "invalid permutation :" << __p << "\n";

            for ( int i = 0; i < numVertices; ++i )
            {
                const int iperm = permLines.find( __p )->second[ i ];
                ublas::column( M_vertices, iperm ) = e.vertex( element_type::e2p( __f, i ) );
            }

            M_points = make_line_points();
        }

        make_normals();
        computeBarycenters();
        computeMeasure();
    }

    Reference( Reference const & r )
        :
        super( r ),
        M_id( r.M_id ),
        M_vertices( r.M_vertices ),
        M_points( r.M_points ),
        M_normals( r.M_normals ),
        M_edge_tangents( r.M_edge_tangents ),
        M_barycenter( r.M_barycenter ),
        M_barycenterfaces( r.M_barycenterfaces ),
        M_meas( r.M_meas )
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
            M_id = r.M_id;
            M_vertices = r.M_vertices;
            M_points = r.M_points;
            M_normals = r.M_normals;
            M_edge_tangents = r.M_edge_tangents;
            M_barycenter = r.M_barycenter;
            M_barycenterfaces = r.M_barycenterfaces;
            M_meas = r.M_meas;
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    uint16_type topologicalDimension() const
    {
        return topological_dimension;
    }
    uint16_type dimension() const
    {
        return real_dimension;
    }

    uint16_type nVertices() const
    {
        return numVertices;
    }
    uint16_type nPoints() const
    {
        return numPoints;
    }
    uint16_type nEdges() const
    {
        return numEdges;
    }
    uint16_type nFaces() const
    {
        return numFaces;
    }

    points_type const& vertices() const
    {
        return M_vertices;
    }

    ublas::matrix_column<points_type const> vertex( uint16_type __i ) const
    {
        return ublas::column( M_vertices, __i );
    }

    ublas::matrix_column<points_type const> edgeVertex( uint16_type __e, uint16_type __p ) const
    {
        return ublas::column( M_vertices, edge_to_point_t::e2p( __e,__p ) );
    }

    /**
     * \return the vertex \p p of the face \f
     */
    ublas::matrix_column<points_type const> faceVertex( uint16_type f, uint16_type p ) const
    {
        return ublas::column( M_vertices, face_to_point_t::f2p( f,p ) );
    }

    /**
     * \return the vertices of the face \p f
     */
    matrix_node_type faceVertices( uint16_type f ) const
    {
        matrix_node_type v( nRealDim, nDim  );

        // there is exactely nDim vertices on each face on a d-simplex
        for ( int p = 0; p < nDim; ++p )
        {
            switch ( nDim )
            {
            case 1:
            case 3:
                ublas::column( v, p ) = ublas::column( M_vertices, face_to_point_t::f2p( f,p ) );
                break;

            case 2:
                ublas::column( v, p ) = ublas::column( M_vertices, face_to_point_t::e2p( f,p ) );
                break;
            }
        }

        return v;
    }

    points_type const& points() const
    {
        return M_points;
    }

    ublas::matrix_column<points_type const> point( uint16_type __i ) const
    {
        return ublas::column( M_points, __i );
    }

    /**
     * \return the barycenter of the reference simplex
     */
    node_type barycenter() const
    {
        return M_barycenter;
    }

    /**
     * \return the barycenter of the faces of the reference simplex
     */
    points_type barycenterFaces() const
    {
        return M_barycenterfaces;
    }

    /**
     * \return the barycenter of the face \p f of the reference simplex
     */
    ublas::matrix_column<matrix_node_type const> faceBarycenter( uint16_type f ) const
    {
        return ublas::column( M_barycenterfaces, f );
    }

    /**
     * get the normals array
     *
     *
     * @return the normals
     */
    normals_type const& normals() const
    {
        return M_normals;
    }

    /**
     * get the n-th normal
     *
     * @param __n the index of the normal
     *
     * @return the n-th normal of the triangle
     */
    node_type const& normal( uint16_type __n ) const
    {
        return M_normals[__n];
    }

    /**
     * get the n-th unit edge tangent
     *
     * @param __n the index of the normal
     *
     * @return the n-th normal of the triangle
     */
    edge_tangent_type const& tangent( uint16_type __n ) const
    {
        return M_edge_tangents[__n];
    }

    /**
     * get the n-th unit edge tangent
     *
     * @param __n the index of the normal
     *
     * @return the n-th normal of the triangle
     */
    edge_tangent_type const& edgeTangent( uint16_type __n ) const
    {
        return M_edge_tangents[__n];
    }


    /**
     * the first iterator of the normal vector
     *
     *
     * @return the begin() iterator of the normal vector
     */
    normal_const_iterator beginNormal() const
    {
        return M_normals.begin();
    }

    /**
     * the end() iterator
     *
     *
     * @return
     */
    normal_const_iterator endNormal() const
    {
        return M_normals.end();
    }

    topological_face_type topologicalFace( uint16_type __f, uint16_type __p = permutation_type::IDENTITY ) const
    {
        topological_face_type ref( *this, __f, __p );
        return ref;
    }

    points_type const& G() const
    {
        return M_points;
    }

    /**
     * \return the measure of the reference element
     */
    double measure() const
    {
        return M_meas;
    }

    size_type id() const
    {
        return 0;
    }
    flag_type marker() const
    {
        return 0;
    }
    flag_type marker2() const
    {
        return 0;
    }
    flag_type marker3() const
    {
        return 0;
    }
    uint16_type ad_first() const
    {
        return -1;
    }
    uint16_type ad_second() const
    {
        return -1;
    }
    uint16_type pos_first() const
    {
        return -1;
    }
    uint16_type pos_second() const
    {
        return -1;
    }
    permutation_type permutation( uint16_type /*f*/ ) const
    {
        return permutation_type();
    }
    double h() const
    {
        // FIXME: should be computed once for all in constructor
        double __max = 0.0;

        for ( int __e = 0; __e < numEdges; ++ __e )
        {
            double __len = ublas::norm_2( edgeVertex( __e, 1 ) - edgeVertex( __e, 0 ) );
            __max = ( __max < __len )?__len:__max;
        }

        return __max;
    }
    double hMin() const
    {
        // FIXME: should be computed once for all in constructor
        double __min = 0.0;

        for ( int __e = 0; __e < numEdges; ++ __e )
        {
            double __len = ublas::norm_2( edgeVertex( __e, 1 ) - edgeVertex( __e, 0 ) );
            __min = ( __min > __len )?__len:__min;
        }

        return __min;
    }
    double h( int e ) const
    {
        return ublas::norm_2( edgeVertex( e, 1 ) - edgeVertex( e, 0 ) );
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

    EntityRange<self_type> entityRange( uint16_type d ) const
    {
        return EntityRange<self_type>( d );
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    points_type makePoints( uint16_type topo_dim, uint16_type __id, int interior = 1 ) const
    {
        // vertices
        if ( topo_dim == 0 )
        {

            //std::cerr << "coucpu 2 " << topo_dim << " " << topological_dimension << " " << __id << "\n";
            points_type G( M_vertices.size1(), 1 );
            ublas::column( G, 0 ) = ublas::column( M_vertices, __id );
            return G;
        }

        // interior points of the convex
        else if ( topo_dim == topological_dimension )
        {
            //std::cerr << "coucpu 2 " << topo_dim << " " << topological_dimension << " " << __id << "\n";
            if ( __id == 0 )
                return this->template makeLattice<Shape>( interior );

            throw std::logic_error( "cannot make those points" );
            return points_type();
        }

        // all the other points
        else
        {
            points_type G;

            if ( topo_dim == 1 )
            {
                Reference<Simplex<1, Order, 1>, 1, Order, 1, T> refline;
                G = refline.template makeLattice<SHAPE_LINE>( interior );
                pt_to_entity<Shape,1> p_to_e( __id );
                points_type Gret( nRealDim, G.size2() );

                for ( size_type i = 0; i < G.size2(); ++i )
                    ublas::column( Gret, i ) = p_to_e( ublas::column( G, i ) );

                return Gret;
            }

            else if ( topo_dim == 2 )
            {
                Reference<Simplex<2, Order, 2>, 2, Order, 2, T> refface;
                G = refface.template makeLattice<SHAPE_TRIANGLE>( interior );
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

            else if ( shape == SHAPE_TRIANGLE )
                return make_triangle_points( interior );

            else if ( shape == SHAPE_TETRA )
                return make_tetrahedron_points( interior );
        }

        DCHECK ( nOrder == 0 ) << "Invalid polynomial order";
        return glas::average( M_vertices );
    }

    /**
     *Let the tetrahedron have vertices
     V1 = (x1, y1, z1)
     V2 = (x2, y2, z2)
     V3 = (x3, y3, z3)
     V4 = (x4, y4, z4)

     and your test point be
     P = (x, y, z).

     Then the point P is in the tetrahedron if following five determinants all have the same sign.
             |x1 y1 z1 1|
        D0 = |x2 y2 z2 1|
             |x3 y3 z3 1|
             |x4 y4 z4 1|

             |x  y  z  1|
        D1 = |x2 y2 z2 1|
             |x3 y3 z3 1|
             |x4 y4 z4 1|

             |x1 y1 z1 1|
        D2 = |x  y  z  1|
             |x3 y3 z3 1|
             |x4 y4 z4 1|

             |x1 y1 z1 1|
        D3 = |x2 y2 z2 1|
             |x  y  z  1|
             |x4 y4 z4 1|

             |x1 y1 z1 1|
        D4 = |x2 y2 z2 1|
             |x3 y3 z3 1|
             |x  y  z  1|

     Some additional notes:
      -# If by chance the D0=0, then your tetrahedron is degenerate (the points are coplanar).
      -# If any other Di=0, then P lies on boundary i (boundary i being that boundary formed by the three points other than Vi).
      -# If the sign of any Di differs from that of D0 then P is outside boundary i.
      -# If the sign of any Di equals that of D0 then P is inside boundary i.
      -# If P is inside all 4 boundaries, then it is inside the tetrahedron.
      -# As a check, it must be that D0 = D1+D2+D3+D4.
      -# The pattern here should be clear; the computations can be extended
      to simplicies of any dimension. (The 2D and 3D case are the
      triangle and the tetrahedron).
      -# the quantities bi = Di/D0 are the usual barycentric
      coordinates.  Comparing signs of Di and D0 is only a check that
      P and Vi are on the same side of boundary i.
      *
     * \return a positive or null number if pt is in the convex
     */
    boost::tuple<bool, value_type>
    isIn( typename node<value_type>::type const& pt ) const
    {
        //std::cout << "[isIn] pt =" << pt << "\n";
        ublas::vector<double> D( nDim + 1, 0 );
        typename matrix_node<value_type>::type M( nDim+1,nDim+1 );
        typename matrix_node<value_type>::type P( nDim,nDim+1 );
        typename matrix_node<value_type>::type P0( nDim,nDim+1 );
        std::fill( M.data().begin(), M.data().end(), value_type( 1.0 ) );

        for ( int p = 0; p < numVertices; ++p )
        {
            ublas::column( P0, p ) = this->vertex( p );
        }

        ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P0;
        //VLOG(1) << "M=" << M << "\n";



        // compute the measure(times 2) of the segments, triangles or
        // tetras formed by the point pt and the vertices of the
        // simplex. The simplices are formed such that when the
        // point is in the simplex the measure are all
        // positive. In the case the point is outside the simplex
        // the measure is negative
        double meas_times = details::det( M, mpl::int_<nDim+1>() );

        for ( int n = 0; n < numVertices; ++n )
        {
            ublas::noalias( P ) = P0;
            ublas::column( P, n ) = pt;
            ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P;
            //VLOG(1) << "M=" << M << "\n";

            // multiply by -1 such that the volume of the
            // reference is > 0
            D( n )= meas_times*details::det( M, mpl::int_<nDim+1>() );

            //std::cout.setf( std::ios::scientific );
            //std::cout.precision( 10 );
            //std::cout   << "D " << n << "=" << D(n) << "\n";
            //VLOG(1) << "sign " << n << "=" << sign << "\n";
            //meas_times = ( meas_times < res )?meas_times:res;
        }

        double dmin = *std::min_element( D.begin(), D.end() );
        ublas::vector<double>::const_iterator Dit = std::find_if( D.begin(), D.end(), lambda::_1 < -5e-7 );

        //VLOG(1) << "meas=" << meas_times << "\n";
        //return meas_times;
        // sign must be > 0 if the point is inside the simplex
        return boost::make_tuple( Dit == D.end(), dmin );
    }
    //@}

private:

    int n_line_points( int interior = 0 ) const
    {
        return std::max( 0, int( Order )+1-2*interior );
    }
    int n_triangle_points( int interior = 0 ) const
    {
        if ( interior == 1 )
            return std::max( 0, ( int( Order )+1-2*interior )*( int( Order )-2*interior )/2 );

        return ( Order+1 )*( Order+2 )/2;
    }
    int n_tetrahedron_points( int interior = 0 ) const
    {
        if ( interior == 1 )
            return std::max( 0, ( int( Order )+1-2*interior )*( int( Order )-2*interior )*( int( Order )-1-2*interior )/6 );

        return ( Order+1 )*( Order+2 )*( Order+3 )/6;
    }

    points_type
    make_line_points( int interior = 0 ) const
    {
        if ( nOrder > 0 )
        {
            ublas::vector<node_type> h ( 1 );
            h( 0 ) = vertex( 1 ) - vertex( 0 );

            points_type p( nRealDim, n_line_points( interior ) );

            for ( int i = interior, indp = 0; i < int( Order )+1-interior; ++i, ++indp )
            {
                ublas::column( p, indp ) = vertex( 0 ) + ( h( 0 ) * value_type( i ) )/value_type( Order );
            }

            return p;
        }

        else
            return glas::average( M_vertices );

    }

    points_type
    make_triangle_points( int interior = 0 ) const
    {
        if ( nOrder > 0 )
        {
            ublas::vector<node_type> h ( 2 );
            h *= 0;
            h( 0 ) = vertex( 1 ) - vertex( 0 );
            h( 1 ) = vertex( 2 ) - vertex( 0 );
            //std::cout << "h = " << h << "\n";
            //DVLOG(2) << "n triangle pts = " << n_triangle_points( interior ) << "\n";
            points_type G( nRealDim, n_triangle_points( interior ) );

            for ( int i = interior, p = 0; i < int( Order )+1-interior; ++i )
            {
                for ( int j = interior; j < int( Order ) + 1 - i-interior; ++j, ++p )
                {
                    ublas::column( G, p ) = vertex( 0 ) + ( value_type( j ) * h( 0 )+
                                                            value_type( i ) * h( 1 ) )/ value_type( Order );
                }
            }

            return G;
        }

        else
            return glas::average( M_vertices );
    }

    points_type
    make_tetrahedron_points( int interior = 0 ) const
    {
        if ( nOrder > 0 )
        {
            ublas::vector<node_type> h ( 3 );
            h( 0 ) = vertex( 1 ) - vertex( 0 );
            h( 1 ) = vertex( 2 ) - vertex( 0 );
            h( 2 ) = vertex( 3 ) - vertex( 0 );
            points_type G( 3, n_tetrahedron_points( interior ) );

            //DVLOG(2) << "n tetra pts = " << n_tetrahedron_points( interior ) << "\n";
            for ( int i = interior, p = 0; i < int( Order )+1-interior; ++i )
            {
                for ( int j = interior; j < int( Order ) + 1 - i - interior; ++j )
                {
                    for ( int k = interior; k < int( Order ) + 1 - i - j - interior; ++k, ++p )
                    {
                        ublas::column( G, p ) = vertex( 0 ) + ( value_type( i ) * h( 2 ) +
                                                                value_type( j ) * h( 1 ) +
                                                                value_type( k ) * h( 0 ) ) / value_type( Order );

                    }
                }
            }

            //std::cout << "G = " << G << "\n";
            return G;
        }

        else
            return glas::average( M_vertices );
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
    struct pt_to_face
    {
        pt_to_face( std::vector<uint16_type> vert_ids )
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
                typename mpl::if_<mpl::equal_to<mpl::size_t<shape>, mpl::size_t<SHAPE_TRIANGLE> >,
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

    void make_edge_tangents()
    {
        uint16_type reindex1[][4] = { { 1, 0, 0, 0},
            { 0, 1, 2, 0},
            { 0, 1, 2, 3}
        };

        for ( uint16_type __n = 0; __n < numEdges; ++__n )
        {
            //const int ind_normal = reindex1[nDim-1][__n];
            M_edge_tangents[__n].resize( nDim );
        }
        for(int edge = 0; edge < numEdges; ++edge )
        {
            M_edge_tangents[edge] = edgeVertex(edge,1)-edgeVertex(edge,0);
            //M_edge_tangents[edge] /= ublas::norm_2( M_edge_tangents[edge] );
        }
#if 0
[0]=-1/math::sqrt( 2. );
            M_edge_tangents[0][1]= 1/math::sqrt( 2. );
        M_edge_tangents[1][0]=  0;
        M_edge_tangents[1][1]= -1;
        M_edge_tangents[2][0]= 1;
        M_edge_tangents[2][1]= 0;
#endif
    }
    void make_normals()
    {
        /*uint16_type reindex[][4] = { { 1, 0, 0, 0},
                                     { 1, 2, 0, 0},
                                     { 2, 3, 1, 0} };
        */
        uint16_type reindex1[][4] = { { 1, 0, 0, 0},
            { 0, 1, 2, 0},
            { 0, 1, 2, 3}
        };


#if 0
        node_type zero = ublas::zero_vector<value_type>( nDim );
        uint16_type d = nDim;
        std::for_each( M_normals.begin(), M_normals.end(),
                       ( lambda::bind( &node_type::resize, lambda::_1,
                                       lambda::constant( d ), false ),
                         lambda::_1 = lambda::constant( zero ) ) );
#endif

        for ( uint16_type __n = 0; __n < numNormals; ++__n )
        {
            //const uint16_type ind_normal = reindex[nDim-1][__n];
            const int ind_normal = reindex1[nDim-1][__n];
            M_normals[ind_normal].resize( nDim );
            M_normals[ind_normal].clear();// = ublas::zero_vector<value_type>( nDim );

            if ( __n > 0 )
            {
                M_normals[ind_normal][__n-1] = -1;
            }

            else
            {
                typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<0> >,
                        mpl::int_<1>,
                        mpl::int_<nDim> >::type denom;
                M_normals[ind_normal] = ublas::scalar_vector<value_type>( nDim, math::sqrt( value_type( 1.0 )/value_type( denom::value ) ) );
            }

            //DVLOG(2) << "normal[" << ind_normal << "]=" << M_normals[ind_normal] << "\n";
        }
    }

    void computeBarycenters();
    void computeMeasure();
private:

    uint16_type M_id;

    points_type M_vertices;

    points_type M_points;

    normals_type M_normals;
    edge_tangents_type M_edge_tangents;

    node_type M_barycenter;

    points_type M_barycenterfaces;

    value_type M_meas;
};

template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>::nbPtsPerVertex;
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>::nbPtsPerEdge;
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>::nbPtsPerFace;
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
const uint16_type Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>::numGeometricFaces;



template<typename T> class Entity<SHAPE_LINE, T>: public Reference<Simplex<1, 1, 1>,1,1, 1, T> {};
template<typename T> class Entity<SHAPE_TRIANGLE, T>: public Reference<Simplex<2, 1, 2>,2,1, 2, T> {};
template<typename T> class Entity<SHAPE_TETRA, T>: public Reference<Simplex<3, 1, 3>,3,1, 3, T> {};


template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
void
Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>::computeBarycenters()
{
    M_barycenter = ublas::column( glas::average( M_vertices ), 0 );

    for ( int f = 0; f < numTopologicalFaces; ++f )
    {
        //std::cout << "face " << f << " vertices " << faceVertices( f ) << "\n";
        ublas::column( M_barycenterfaces, f ) = ublas::column( glas::average( faceVertices( f ) ), 0 );
    }
}
template<uint16_type Dim, uint16_type Order, uint16_type RDim,  typename T>
void
Reference<Simplex<Dim, Order, RDim>, Dim, Order, RDim, T>::computeMeasure()
{
    if ( nDim == nRealDim )
    {

        typename matrix_node<value_type>::type M( nDim,nDim );
        value_type factor( 1 );

        switch ( nDim )
        {
        case 1:
            ublas::column( M, 0 ) = this->vertex( 1 )-this->vertex( 0 );
            factor = 1;
            break;

        case 2:
            ublas::column( M, 0 ) = this->vertex( 0 )-this->vertex( 1 );
            ublas::column( M, 1 ) = this->vertex( 1 )-this->vertex( 2 );
            factor = 2;
            break;

        case 3:
            /**
             * tetrahedron with vertices
             * a = (a1, a2, a3), b = (b1, b2, b3), c = (c1, c2, c3), and d = (d1, d2, d3),
             * the volume is (1/6)·|det(a−b, b−c, c−d)|
             */
            ublas::column( M, 0 ) = this->vertex( 0 )-this->vertex( 1 );
            ublas::column( M, 1 ) = this->vertex( 1 )-this->vertex( 2 );
            ublas::column( M, 2 ) = this->vertex( 2 )-this->vertex( 3 );
            factor = 6;
            break;
        }

        M_meas = math::abs( details::det( M, mpl::int_<nDim>() ) )/factor;
    }

    else
    {
#if 0
        /**
           In three dimensions, the area of a general triangle {A = (xA, yA,
           zA), B = (xB, yB, zB) and C = (xC, yC, zC)} is the Pythagorean sum of
           the areas of the respective projections on the three principal planes
           (i.e. x = 0, y = 0 and z = 0):
        */
        typename matrix_node<value_type>::type M( nRealDim,nRealDim );
        value_type factor( 1 );

        switch ( nDim )
        {
        case 1:
            M_meas = ublas::norm2( this->vertex( 1 )-this->vertex( 0 ) );
            break;

        case 2:
            ublas::column( M, 0 ) = this->vertex( 0 )-this->vertex( 1 );
            ublas::column( M, 1 ) = this->vertex( 1 )-this->vertex( 2 );
            factor = 2;
            break;
        }

#endif
    }
}

}
#endif /* __refsimplex_H */
