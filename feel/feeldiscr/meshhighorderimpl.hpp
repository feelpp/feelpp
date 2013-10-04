/*
  This file is part of the Feel library

  Copyright (C) 2007,2008 University of Coimbra
  Copyright (C) 2010 Université de Grenoble 1

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
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <feel/feeldiscr/meshhighorder.hpp>

namespace Feel
{
template < class Convex >
MeshHighOrder<Convex>::MeshHighOrder( mesh_ptrtype& mesh,
                                      std::string projectionStrategy )
    :
    old_mesh( mesh ),
    new_mesh( new new_mesh_type ),
    strategy ( projectionStrategy )
{
    this->createSwapEdgesMap( mpl::bool_< is_simplex >() );
}

template < class Convex >
MeshHighOrder<Convex>::MeshHighOrder( MeshHighOrder const& tc )
    :
    old_mesh( tc.old_mesh ),
    new_mesh( tc.new_mesh ),
    strategy( tc.strategy ),
    M_swap_edges( tc.M_swap_edges )
{
}

template < class Convex >
MeshHighOrder<Convex>::~MeshHighOrder()
{
}


template < class Convex >
void
MeshHighOrder<Convex>::clearMesh()
{
    new_mesh->clear();
}

template < class Convex >
void
MeshHighOrder<Convex>::updateP1Mesh ( mesh_ptrtype const& mesh )
{
    old_mesh = mesh;
}

template < class Convex >
typename MeshHighOrder<Convex>::new_mesh_ptrtype
MeshHighOrder<Convex>::getMesh() const
{
    return new_mesh;
}


// affine transformation from [-1,1] to the edge defined by node1 and node2
template < class Convex >
typename MeshHighOrder<Convex>::node_type
MeshHighOrder<Convex>::affineEdge( node_type const& node1, node_type const& node2, double const& x ) const
{
    return ( 0.5*( ( 1-x )*node1 + ( x+1 )*node2 ) );
}


template < class Convex >
void
MeshHighOrder<Convex>::createSwapEdgesMap ( mpl::bool_<true>  )
{
    // Construct a map that returns if for a given element and an edge,
    // the edge points should be inverted (for triangles only)
    for ( element_const_iterator elt = old_mesh->beginElement();
            elt != old_mesh->endElement(); ++elt )
    {
        //check if point(0) is elt->point(1)
        // if they are not the same, points need to be swapped
        if ( elt->edge( 0 ).point( 0 ).id() != elt->point( 1 ).id() )
            M_swap_edges[elt->id()][0] = 1;

        //check if point(0) is elt->point(2)
        // if they are not the same, points need to be swapped
        if ( elt->edge( 1 ).point( 0 ).id() != elt->point( 2 ).id() )
            M_swap_edges[elt->id()][1] = 1;

        //check if point(0) is elt->point(0)
        // if they are not the same, points need to be swapped
        if ( elt->edge( 2 ).point( 0 ).id() != elt->point( 0 ).id() )
            M_swap_edges[elt->id()][2] = 1;
    }
}

template < class Convex >
void
MeshHighOrder<Convex>::createSwapEdgesMap ( mpl::bool_<false>  )
{
    // Construct a map that returns if for a given element and an edge,
    // the edge points should be inverted (for triangles only)
    for ( element_const_iterator elt = old_mesh->beginElement();
            elt != old_mesh->endElement(); ++elt )
    {
#if 0

        //check if point(0) is elt->point(1)
        // if they are not the same, points need to be swapped
        if ( elt->edge( 0 ).point( 0 ).id() != elt->point( 1 ).id() )
            M_swap_edges[elt->id()][0] = 1;

        //check if point(0) is elt->point(2)
        // if they are not the same, points need to be swapped
        if ( elt->edge( 1 ).point( 0 ).id() != elt->point( 2 ).id() )
            M_swap_edges[elt->id()][1] = 1;

        //check if point(0) is elt->point(0)
        // if they are not the same, points need to be swapped
        if ( elt->edge( 2 ).point( 0 ).id() != elt->point( 0 ).id() )
            M_swap_edges[elt->id()][2] = 1;

#endif
    }
}



template < class Convex >
bool
MeshHighOrder<Convex>::swapEdge ( size_type eltId, size_type edgeId )
{
    return M_swap_edges[eltId][edgeId];
}

template < class Convex >
void
MeshHighOrder<Convex>::updatePts( ublas::vector<double> const& x, points_type const& pts,
                                  points_type& final_pts ) const
{
    ublas::row( final_pts, 0 ) += element_prod( ublas::row( pts, 0 ), x );
    ublas::row( final_pts, 1 ) += element_prod( ublas::row( pts, 1 ), x );
}


template < class Convex >
void
MeshHighOrder<Convex>::addVertices( element_type const& elt, new_element_type& new_element,
                                    new_mesh_ptrtype& new_mesh, std::vector<bool>& vertexAdd ) const
{
    for ( uint16_type i = 0; i< element_type::numVertices; ++i )
    {
        point_type old_point = elt.point( i );

        point_type new_point( old_point.id(), old_point.node(), old_point.isOnBoundary() );
        new_point.marker().assign( old_point.marker().value() );

        // if the created point has not been added to the mesh then add it
        if ( vertexAdd[ old_point.id() ] == 0 )
        {
            new_mesh->addPoint( new_point );

#if !defined ( NDEBUG )
            DVLOG(2) << "[AddPointToMesh] Vertex of id: " << new_point.id() << " has coordinates "
                          << new_point.node();
#endif

            vertexAdd[ new_point.id() ] = 1;
        }

        // add point to the element it belongs
        new_element.setPoint( i, new_mesh->point( new_point.id() ) );

#if !defined ( NDEBUG )
        DVLOG(2) << "[AddVertexElement] Add point: localId=" << i
                      << " globalToMeshId=" << new_mesh->point( new_point.id() ).id()
                      << " to element " << new_element.id() << "\n";
#endif

    }
}

}



