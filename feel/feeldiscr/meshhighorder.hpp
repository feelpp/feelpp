/*
  This file is part of the Feel library

  Copyright (C) 2007,2008 EPFL

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

#ifndef __MeshHighOrder
#define __MeshHighOrder 1


#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feelvf/vf.hpp>


namespace Feel
{
/**
 * \class MeshHighOrder
 * Class to handle high order meshes (only works for 2D meshes)
 */
template< class Convex >
class MeshHighOrder
{
    static const bool is_simplex = Convex::is_simplex;
    static const uint16_type Dim = Convex::nDim;
    static const uint16_type Order = Convex::nOrder;


    typedef typename mpl::if_< mpl::bool_< is_simplex >, Simplex<Dim, 1>, Hypercube<Dim, 1> >::type convex_type;
    typedef Mesh< convex_type > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    typedef typename mpl::if_< mpl::bool_< is_simplex >, Simplex<2, Order>, Hypercube<2, Order> >::type new_convex_type;
    typedef Mesh< new_convex_type > new_mesh_type;

    typedef boost::shared_ptr<new_mesh_type> new_mesh_ptrtype;

    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::element_const_iterator element_const_iterator;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::point_type point_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::node_type node_type;

    typedef typename new_mesh_type::edge_type new_edge_type;
    typedef typename new_mesh_type::element_type new_element_type;
    typedef typename new_mesh_type::element_const_iterator new_element_const_iterator;

    typedef PointSetEquiSpaced<Hypercube<1,1>, Order, double> oned_pointset_type;
    typedef PointSetEquiSpaced<convex_type, Order, double> interior_pointset_type;

    typedef typename oned_pointset_type::points_type node_points_type;
    typedef typename interior_pointset_type::points_type points_type;

public:

    typedef new_mesh_type ho_mesh_type;

    MeshHighOrder( mesh_ptrtype& mesh,
                   std::string projectionStrategy = "axis" );

    MeshHighOrder( MeshHighOrder const& tc );

    ~MeshHighOrder();

    void clearMesh();

    void updateP1Mesh ( mesh_ptrtype const& mesh );

    new_mesh_ptrtype getMesh() const;

    template< typename elem_type >
    void generateMesh( std::vector<flag_type> const& flags, std::vector<elem_type> const& polyBoundary )
    {
        new_mesh->clear();

        //Create pointset
        oned_pointset_type pts_space;
        node_points_type span = pts_space.pointsBySubEntity( 1,0 );

        // stores a boolean indicating if the vertex is in the new mesh
        std::vector<bool> vertexAdd( old_mesh->numVertices() );

        // stores a boolean indicating if the edge is in the new mesh
        std::vector<bool> edgeAdd( old_mesh->numEdges() );

        // connectivity map for edges
        std::vector<uint16_type> edgesMap ( old_mesh->numEdges() );

        std::fill( vertexAdd.begin(), vertexAdd.end(), 0 );
        std::fill( edgeAdd.begin(), edgeAdd.end(), 0 );
        std::fill( edgesMap.begin(), edgesMap.end(), 0 );

        uint16_type nodesCount = old_mesh->numVertices();

        // for in the elements
        for ( element_const_iterator elt = old_mesh->beginElement();
                elt != old_mesh->endElement(); ++elt )
        {
            new_element_type new_element;
            new_element.setId( elt->id() );
            new_element.setOnBoundary( elt->isOnBoundary() );

            // add vertices to mesh
#if !defined ( NDEBUG )
            DVLOG(2) << "Add vertices...\n";
#endif
            this->addVertices( *elt, new_element, new_mesh, vertexAdd );


            // add edges and nodes in edges to mesh
#if !defined ( NDEBUG )
            DVLOG(2) << "Add edges...\n";
#endif

            for ( uint16_type i = 0; i< element_type::numEdges; ++i )
            {
                edge_type old_edge = elt->edge( i );

                bool hasEdge = edgeAdd[ old_edge.id() ];

                new_edge_type new_edge;
                new_edge.setId( old_edge.id() );
                new_edge.marker().assign( old_edge.marker().value() );
                new_edge.setOnBoundary( old_edge.isOnBoundary() );

                // if the edge is not found in the list, create it and add the points in the boundary
                if ( hasEdge == 0 )
                {
                    // add points to edge:
                    for ( uint16_type j = 0; j<2; ++j )
                    {
                        uint16_type localId = j;

                        if ( !swapEdge( elt->id(), i ) )
                            new_edge.setPoint( j, new_mesh->point( old_edge.point( localId ).id() ) );

                        else
                        {
                            localId = 1-j;
                            new_edge.setPoint( j, new_mesh->point( old_edge.point( localId ).id() ) );
                        }

#if !defined ( NDEBUG )
                        DVLOG(2) << "[AddPointToEdge] Point localId=" << j
                                      << ": globalToMeshId=" << new_mesh->point( old_edge.point( localId ).id() ).id()
                                      << "; " << new_mesh->point( old_edge.point( localId ).id() ).node();
#endif
                    }

                    edgesMap[ new_edge.id() ] = nodesCount;
                }


                node_type node( Dim );

                for ( uint16_type j=2; j<=Order; ++j )
                {
                    double t = span( 0,j-2 );

                    if ( swapEdge( elt->id(), i ) )
                        node = this->affineEdge( old_edge.point( 1 ).node(), old_edge.point( 0 ).node(), t );

                    else
                        node = this->affineEdge( old_edge.point( 0 ).node(), old_edge.point( 1 ).node(), t );

                    //strategies for curved edges:
                    if ( old_edge.isOnBoundary() )
                    {
                        this->shiftCoordinates( node, old_edge.point( 0 ).node(), old_edge.point( 1 ).node(),
                                                old_edge.marker(), flags, polyBoundary );
                    }

                    point_type new_point( nodesCount, node, old_edge.isOnBoundary() );
                    new_point.marker().assign( old_edge.marker().value() );
                    //uint16_type localId = j;

                    if ( hasEdge == 0 )
                    {
                        new_mesh->addPoint( new_point );

                        new_edge.setPoint( j, new_mesh->point( nodesCount ) );

#if !defined ( NDEBUG )
                        DVLOG(2) << "[AddPointToMesh] Point " << j
                                      << ": id=" << new_mesh->point( nodesCount ).id()
                                      << "; " << new_mesh->point( nodesCount ).node();
#endif
                        nodesCount++;
                    }

                    else
                    {
                        if ( !swapEdge( elt->id(), i ) )
                            new_point.setId( edgesMap[ new_edge.id() ] + j - 2 );

                        else
                            new_point.setId( edgesMap[ new_edge.id() ] + Order - j );
                    }

                    new_element.setPoint( element_type::numVertices +  ( Order-1 )*i +j-2, new_mesh->point( new_point.id() ) );

#if !defined ( NDEBUG )
                    DVLOG(2) << "[AddPointToElement] Add point with local id "
                                  << element_type::numVertices + ( Order-1 )*i +j-2
                                  << " with global id " << new_mesh->point( new_point.id() ).id()
                                  << " to element " << new_element.id() << "\n";
#endif
                }

                if ( hasEdge == 0 )
                {
                    new_mesh->addFace( new_edge );

                    edgeAdd[ new_edge.id() ] = 1;
                }

#if !defined ( NDEBUG )
                DVLOG(2) << "[AddToMesh] Edge of id " << new_edge.id() << " has been added\n";
                DVLOG(2) << "\n";
#endif
            } // end for in edges

            // for in the points in the interior of the element
            if ( Order > ( 1 + is_simplex ) )
            {
                interior_pointset_type pts_interior( 1 );

                points_type pts = pts_interior.points();

                this->GordonHall( *elt, pts, flags, polyBoundary );

                // add interior points given by the GordonHall transformation to element and mesh
                for ( uint16_type i=0; i < pts.size2(); ++i )
                {
                    // add point to mesh
                    node_type node ( Dim );
                    node[0] = pts( 0,i );
                    node[1] = pts( 1,i );

                    point_type new_point( nodesCount, node, false );
                    new_point.marker().assign( 0 );

                    new_mesh->addPoint( new_point );

                    new_element.setPoint( element_type::numVertices + ( Order-1 )*element_type::numEdges +i,
                                          new_mesh->point( new_point.id() ) );

#if !defined ( NDEBUG )
                    DVLOG(2) << "Added Point "
                                  << element_type::numVertices + ( Order-1 )*element_type::numEdges +i
                                  << ": id=" << new_mesh->point( new_point.id() ).id() << "; "
                                  << new_mesh->point( new_point.id() ).node() << "\n";
#endif

                    nodesCount++;
                }

            }

            new_mesh->addElement( new_element );

#if !defined ( NDEBUG )
            DVLOG(2) << "[AddToMesh] Element of id " << new_element.id() << " has been added\n";
            DVLOG(2) << "-------------------------------------------------------\n\n";
#endif

        } // end for in elements

        new_mesh->setNumVertices( old_mesh->numVertices() );

#if !defined ( NDEBUG )
        DVLOG(2) << "Number of elements in the new mesh: " << new_mesh->numElements() << "\n";
#endif

        new_mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        new_mesh->updateForUse();
    }


    template< typename elem_type >
    void GordonHall( element_type const& elt, points_type& pts,
                     std::vector<flag_type> const& flags, std::vector<elem_type> const& polyBoundary )
    {
        points_type final_pts = 0*pts;

        node_type node1 ( Dim );
        node_type node2 ( Dim );
        ublas::vector<double> ones( ublas::scalar_vector<double>( pts.size2(), 1.0 ) );

        ublas::vector<double> xi( pts.size2(), 0.0 );
        ublas::vector<double> eta( pts.size2(), 0.0 );
        xi = ublas::row( pts, 0 );
        eta = ublas::row( pts, 1 );

        for ( uint16_type i = 0; i< element_type::numEdges; ++i )
        {
            std::vector<flag_type>::const_iterator result;
            result = find( flags.begin(), flags.end(), elt.edge( i ).marker().value() );

            uint16_type pos = distance( flags.begin(), result );

            node1 = elt.edge( i ).point( 0 ).node();
            node2 = elt.edge( i ).point( 1 ).node();

            if ( this->swapEdge( elt.id(), i ) )
                std::swap( node1,node2 );

            bool isAffine = ( result == flags.end() );

            this->applyDeformation( i, node1, node2,
                                    xi, eta, ones, polyBoundary[pos], isAffine,
                                    final_pts, mpl::bool_<is_simplex>() );
        }

        pts = final_pts;
    }


private:

    mesh_ptrtype old_mesh;

    new_mesh_ptrtype new_mesh;

    std::string strategy;

    std::map< size_type, std::map< size_type, bool > > M_swap_edges;

    // affine transformation from [-1,1] to the edge defined by node1 and node2
    node_type affineEdge( node_type const& node1, node_type const& node2, double const& x ) const;

    void createSwapEdgesMap ( mpl::bool_<true>  );

    void createSwapEdgesMap ( mpl::bool_<false>  );

    void updatePts( ublas::vector<double> const& x, points_type const& pts, points_type& final_pts ) const;

    void addVertices( element_type const& elt, new_element_type& new_element,
                      new_mesh_ptrtype& new_mesh, std::vector<bool>& vertexAdd ) const;

    bool swapEdge ( size_type eltId, size_type edgeId );

#if 0
    template< typename elem_type >
    class ProjectPoints
    {
    public:

        typedef typename node<double>::type node_t_type;
        typedef typename matrix_node<double>::type matrix_node_t_type;
        typedef typename elem_type::functionspace_type::node_type oned_node_type;

        ProjectPoints( elem_type& _p, node_type& _node, node_type const& node0, node_type const& node1 )
            :
            p( _p ),
            m( ( node1[1] - node0[1] )/( node1[0] - node0[0] ) ),
            node( _node )

        {
            using namespace Feel::vf;
            //interpolate(p.functionSpace(), gradv(p), der );
            //der = project( p.functionSpace(), p.functionSpace()->mesh(), gradv(p) );
            //std::cout << der << "\n";

        }

        double operator()( const node_t_type& x ) const
        {
            //oned_node_type pt(1);
            //pt[0] = x[0];

            //double feval = p(pt)(0,0,0) - node[1] + (1/m)*(pt[0] - node[0]);
            std::cout << "Eval f( " << x[0] << " ) = " << p( x )( 0,0,0 ) - node[1] + ( 1/m )*( x[0] - node[0] ) << "\n";
            return p( x )( 0,0,0 ) - node[1] + ( 1/m )*( x[0] - node[0] );
        }

        void operator()( const node_t_type& x, node_t_type& gr ) const
        {
            typename elem_type::grad_type gradient( p.grad( x ) );
            double g_v1_x = gradient( 0,0,0 );
            std::cout << "Gradient " << g_v1_x << "\n";

            /*oned_node_type pt(1);
              pt[0] = x[0];*/
            std::cout << "Eval der( " << x[0] << " ) = " << -x[0]/math::sqrt( 1-x[0]*x[0] ) + 1/m << "\n";
            gr[0] = -x[0]/math::sqrt( 1-x[0]*x[0] ) + 1/m;
        }

    private:

        elem_type p, der;
        double m;
        node_type node;

    };

#endif







    template< typename elem_type >
    void shiftCoordinates( node_type& node,
                           node_type const& node0,
                           node_type const& node1,
                           Marker1& marker,
                           std::vector<flag_type> const& flags,
                           std::vector<elem_type> const& polys ) const
    {
        std::vector<flag_type>::const_iterator result;
        result = find( flags.begin(), flags.end(), marker.value() );

        if ( result != flags.end() )
        {
            // determine the position of the flag found in vector flags
            uint16_type pos = distance( flags.begin(), result );

            typedef typename elem_type::functionspace_type::node_type oned_node_type;

            oned_node_type pt( 1 );

            pt[0] = node[0];

            elem_type p = polys[pos];

            if ( strategy == "axis" )
            {

#if !defined ( NDEBUG )
                DVLOG(2) << "Point (" << node[0] << "," << node[1]
                              << ") moves to (" << node[0] << "," << p( pt )( 0,0,0 )
                              << ")" << "\n";
#endif

                node[1] = p( pt )( 0,0,0 );
            }

            else
            {
                if ( strategy == "normal" )
                {
#if 0
                    ProjectPoints<elem_type> Projector( p, node, node0, node1 );

                    iteration_ptrtype iter( Iteration<double>::New() );
                    iter->setMaximumNumberOfIterations( 50 );
                    iter->setRelativePrecision( 1e-8 );

                    bfgs( Projector, Projector, pt, 10, *iter );

                    std::cout << "Result: " << pt << "\n";
#endif

                    double n0 = node0[0];
                    double n1 = node1[0];

                    if ( n0 > n1 )
                        std::swap( n0,n1 );

                    DVLOG(2) << "Initial interval: [" << n0 << "," << n1 << "]\n";
                    std::pair<double, double> interval = std::make_pair( n0,n1 );

                    double m = ( node1[1] - node0[1] )/( node1[0] - node0[0] );
                    DVLOG(2) << "Slope: " << m << "\n";

                    double error = 1;
                    uint16_type iter = 0;

                    while ( error > 1e-14 && iter<100  )
                    {
                        double midpoint = ( interval.first+interval.second )/2;

                        pt[0] = midpoint;
                        double feval = p( pt )( 0,0,0 ) - node[1] + ( 1/m )*( pt[0] - node[0] );

                        pt[0] = interval.first;
                        double feval0 = p( pt )( 0,0,0 ) - node[1] + ( 1/m )*( pt[0] - node[0] );

                        if ( feval*feval0 > 0 )
                            interval.first = midpoint;

                        else
                            interval.second = midpoint;

                        error = math::abs( feval );

                        DVLOG(2) << "Abs(residual) = " << error << "\n";

                        ++iter;
                    }

                    DVLOG(2) << "Finished in " << iter << " iterations...\n";

                    DVLOG(2) << "Point (" << node[0] << "," << node[1]
                                  << ") moves to (" << pt[0] << "," << p( pt )( 0,0,0 )
                                  << ")" << "\n";

                    node[1] = p( pt )( 0,0,0 );
                    node[0] = pt[0];

                }
            }
        }
    }


    template< typename elem_type >
    void deformEdge( node_type const& node1, node_type const& node2,
                     elem_type const& p, ublas::vector<double> const& x,
                     points_type& pts, bool const& isAffine = 1 ) const
    {
        if ( isAffine )
        {
            ublas::vector<double> one( ublas::scalar_vector<double>( x.size(), 1.0 ) );

            ublas::row( pts, 0 ) = 0.5*( node1[0]*( one - x ) + node2[0]*( one+x ) );
            ublas::row( pts, 1 ) = 0.5*( node1[1]*( one - x ) + node2[1]*( one+x ) );
        }

        else
        {
            double a = node1[0];
            double b = node2[0];

            typedef typename elem_type::functionspace_type::node_type oned_node_type;
            oned_node_type pt( 1 );

            for ( uint16_type i = 0; i < x.size(); i++ )
            {
                pt[0] = 0.5*( ( b-a )*x( i ) + b+a );
                pts( 0,i ) = pt[0];
                pts( 1,i ) = p( pt )( 0,0,0 );
            }
        }
    }

    template< typename elem_type >
    void applyDeformation( uint16_type i, node_type const& node1, node_type const& node2,
                           ublas::vector<double> const& xi, ublas::vector<double> const& eta,
                           ublas::vector<double> const& ones,
                           elem_type const& p, bool const& isAffine,
                           points_type& final_pts, mpl::bool_<true> ) const
    {
        points_type nodes3 = points_type( Dim, final_pts.size2() );

        switch ( i )
        {
        case 0:
        {
            // calculates nodes3 = \pi_0 ( -\xi )
            deformEdge( node1, node2, p, -xi, nodes3, isAffine );

            // calculates ( 1 + 0.5*(\xi+\eta) )*nodes3
            updatePts( ones + 0.5*( xi+eta ), nodes3, final_pts );

            deformEdge( node1, node2, p, -( ones+xi+eta ), nodes3, isAffine );
            updatePts( - 0.5*( ones + xi ), nodes3, final_pts );
            break;
        }

        case 1:
        {
            deformEdge( node1, node2, p, -eta, nodes3, isAffine );
            updatePts( 0.5*( ones-xi ), nodes3, final_pts );

            deformEdge( node1, node2, p, xi, nodes3, isAffine );
            updatePts( - 0.5*( ones + eta ), nodes3, final_pts );

            deformEdge( node1, node2, p, ones, nodes3, isAffine );
            updatePts( 0.5*( eta+xi ), nodes3, final_pts );
            break;
        }

        default:
        {
            deformEdge( node1, node2, p, xi, nodes3, isAffine );
            updatePts( 0.5*( ones-eta ), nodes3, final_pts );

            deformEdge( node1, node2, p, -eta, nodes3, isAffine );
            updatePts( - 0.5*( ones + xi ), nodes3, final_pts );

            deformEdge( node1, node2, p, ones, nodes3, isAffine );
            updatePts( 0.5*( ones+xi ), nodes3, final_pts );
        }
        }
    }



    template< typename elem_type >
    void applyDeformation( uint16_type i, node_type const& node1, node_type const& node2,
                           ublas::vector<double> const& xi, ublas::vector<double> const& eta,
                           ublas::vector<double> const& ones,
                           elem_type const& p, bool const& isAffine,
                           points_type& final_pts, mpl::bool_<false> ) const
    {
        points_type nodes3 = points_type( Dim, final_pts.size2() );

        switch ( i )
        {
        case 0:
        {
            // calculates nodes3 = \pi_0 ( xi )
            deformEdge( node1, node2, p, xi, nodes3, isAffine );

            // calculates ( 0.5*(1-eta) )*nodes3
            updatePts( 0.5*( ones-eta ), nodes3, final_pts );

            // calculates nodes3 = \pi_0 ( 1 )
            deformEdge( node1, node2, p, ones, nodes3, isAffine );

            // calculates ( -0.25*( (1+xi)*(1-eta)) )*nodes3
            updatePts( -0.25*ublas::element_prod( ones+xi, ones - eta ), nodes3, final_pts );
        }

        case 1:
        {
            // calculates nodes3 = \pi_1 ( eta )
            deformEdge( node1, node2, p, eta, nodes3, isAffine );

            // calculates ( 0.5*(1+xi) )*nodes3
            updatePts( 0.5*( ones+xi ), nodes3, final_pts );

            // calculates nodes3 = \pi_1 ( 1 )
            deformEdge( node1, node2, p, ones, nodes3, isAffine );

            // calculates ( -0.25*( (1+xi)*(1+eta)) )*nodes3
            updatePts( -0.25*ublas::element_prod( ones+xi, ones + eta ), nodes3, final_pts );
        }

        case 2:
        {
            // calculates nodes3 = \pi_2 ( -xi )
            deformEdge( node1, node2, p, -xi, nodes3, isAffine );

            // calculates ( 0.5*(1+eta) )*nodes3
            updatePts( 0.5*( ones+eta ), nodes3, final_pts );

            // calculates nodes3 = \pi_2 ( 1 )
            deformEdge( node1, node2, p, ones, nodes3, isAffine );

            // calculates ( -0.25*( (1-xi)*(1+eta)) )*nodes3
            updatePts( -0.25*ublas::element_prod( ones-xi, ones + eta ), nodes3, final_pts );
        }

        default:
        {
            // calculates nodes3 = \pi_3 ( -eta )
            deformEdge( node1, node2, p, -eta, nodes3, isAffine );

            // calculates ( 0.5*(1-xi) )*nodes3
            updatePts( 0.5*( ones-xi ), nodes3, final_pts );

            // calculates nodes3 = \pi_3 ( 1 )
            deformEdge( node1, node2, p, ones, nodes3, isAffine );

            // calculates ( -0.25*( (1-xi)*(1-eta)) )*nodes3
            updatePts( -0.25*ublas::element_prod( ones-xi, ones - eta ), nodes3, final_pts );
        }
        }
    }

};
}

#if !defined(FEELPP_INSTANTIATION_MODE)
# include <feel/feeldiscr/meshhighorderimpl.hpp>
#endif //

#endif // __MeshHighOrder
