/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-03-03

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
   \file pointsetinterpolation.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-03-03
 */
#ifndef __PointSetInterpolation_H
#define __PointSetInterpolation_H 1

#include <feel/feelmesh/pointset.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

/**
 * \class PointSetInterpolation
 * \brief Good interpolation point sets on a convex
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type Dim,
         uint16_type Order,
         typename T,
         template<uint16_type,uint16_type,uint16_type> class Convex = Simplex>
class PointSetInterpolation : public PointSet<Convex<Dim,Order,Dim>,T>
{
    typedef PointSet<Convex<Dim,Order,Dim>,T> super;

public:

    /** @name Typedefs
     */
    //@{

    typedef typename super::return_type return_type;
    typedef T value_type;
    typedef Convex<Dim,Order,Dim> convex_type;

    typedef typename super::nodes_type nodes_type;
    typedef typename matrix_node<value_type>::type points_type;

    static const bool is_simplex = convex_type::is_simplex;
    static const bool is_hypercube = convex_type::is_hypercube;

    static const uint32_type convexOrder = convex_type::nOrder;
    static const uint32_type topological_dimension = convex_type::topological_dimension;

    typedef mpl::if_< mpl::bool_< is_simplex >,
            Simplex<Dim, Order, Dim> ,
            Hypercube<Dim, Order, Dim> > conv_order_type;

    typedef Reference<convex_type, Dim, convexOrder, Dim, value_type> RefElem;

    static const uint32_type numPoints = conv_order_type::type::numPoints;

    static const uint32_type nbPtsPerVertex = conv_order_type::type::nbPtsPerVertex;
    static const uint32_type nbPtsPerEdge = conv_order_type::type::nbPtsPerEdge;
    static const uint32_type nbPtsPerFace = conv_order_type::type::nbPtsPerFace;
    static const uint32_type nbPtsPerVolume = conv_order_type::type::nbPtsPerVolume;

    RefElem RefConv;

    typedef std::pair<uint16_type, uint16_type> range_type;
    typedef std::vector< std::vector<size_type> > index_map_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    PointSetInterpolation()
        :
        super()
    {
        M_eid.resize( topological_dimension + 1 );
        M_pt_to_entity.resize( numPoints );
    }
    PointSetInterpolation( size_type np )
        :
        super( np )
    {
        M_eid.resize( topological_dimension + 1 );
        M_pt_to_entity.resize( numPoints );
    }
    PointSetInterpolation( PointSetInterpolation & psi )
        :
        super( psi ),
        M_eid( psi.getEid() ),
        M_pt_to_entity( psi.getPtE() )

    {}

    ~PointSetInterpolation() {}

    ublas::matrix_range<nodes_type const> pointsByEntity( uint16_type e ) const
    {
        FEELPP_ASSERT( M_eid[e].size() )( e ).error( "no points defined on this entity" );

        return ublas::project( this->points(),
                               ublas::range( 0,Dim ),
                               ublas::range( *M_eid[e].begin(), *M_eid[e].rbegin() + 1 ) );
    }

    range_type interiorRangeById( uint16_type e, uint16_type id ) const
    {
        uint16_type numEntities = 1;

        if ( e == 0 )
            numEntities = convex_type::numVertices;

        else if ( e == 1 )
            numEntities = convex_type::numEdges;

        else if ( e == 2 )
            numEntities = convex_type::numFaces;

        uint16_type N = M_eid[e].size()/numEntities;

        return std::make_pair( *M_eid[e].begin() + id*N, *M_eid[e].begin() + ( id+1 )*N );
    }

    ublas::matrix_range<nodes_type const> interiorPointsById( uint16_type e, uint16_type id ) const
    {
        range_type position = interiorRangeById( e, id );

        ublas::matrix_range<nodes_type const> G = ublas::project( this->points(),
                ublas::range( 0, Dim ),
                ublas::range( position.first, position.second ) );

        return G;
    }

    uint32_type entityIds( int i, int j ) const
    {
        return M_eid[i][j];
    }

    uint32_type numEntities( int i ) const
    {
        return M_eid[i].size();
    }

    range_type const& pointToEntity( int p ) const
    {
        return M_pt_to_entity[p];
    }

    //Returns the local indices of all the subentities that compose the entity
    index_map_type entityToLocal ( uint16_type top_dim, uint16_type local_id, bool boundary = 0 )
    {
        index_map_type indices( top_dim+1 );

        if ( top_dim == 0 && boundary )
        {
            range_type pair = interiorRangeById( top_dim, local_id );

            indices[0].push_back( pair.first );
        }

        else
        {
            if ( !boundary && top_dim > 0 )
            {
                if ( numEntities( top_dim ) != 0 )
                    indices[top_dim].push_back( local_id );
            }

            else
            {
                if ( top_dim > 0 )
                {
                    //For the time being, this definition works, but in a more general framework,
                    //a entity should be created here with the respective top_dim
                    //in order to retrieve the number of vertices, edges, faces.

                    //number of vertices, number of edges and number of faces in the volume
                    std::map<uint16_type, uint16_type> numPointsInEntity;

                    numPointsInEntity[0] = ( is_simplex )?( top_dim+1 ):( 2 + ( top_dim-1 )*( 2+3*( top_dim-2 ) ) );
                    numPointsInEntity[1] = ( 2 + ( top_dim-1 )*( 1+2*is_simplex + ( 5-4*is_simplex )*( top_dim-1 ) ) )/2;
                    numPointsInEntity[2] = 6-2*is_simplex;

                    for ( uint16_type i=0; i < numPointsInEntity[0] ; i++ )
                    {
                        if ( top_dim == 1 )
                            indices[0].push_back( RefConv.e2p( local_id, i ) );

                        else if ( top_dim == 2 )
                            indices[0].push_back( RefConv.f2p( local_id, i ) );
                    }

                    if ( ( top_dim == 2 ) && ( M_eid[1].size() != 0 ) )
                    {
                        for ( uint16_type i=0; i < numPointsInEntity[1] ; i++ )
                            indices[1].push_back( RefConv.f2e( local_id, i ) );
                    }

                    if ( top_dim == 3 )
                    {
                        for ( uint16_type k=0; k<numPointsInEntity.size(); k++ )
                        {
                            for ( uint16_type i=0; i < numPointsInEntity[k]; i++ )
                            {
                                indices[k].push_back( i );
                            }
                        }
                    }

                    if ( M_eid[top_dim].size() )
                    {
                        indices[top_dim].push_back( local_id );
                    }
                }
            }
        }

        return indices;
    }


    points_type pointsBySubEntity( uint16_type top_dim, uint16_type local_id, bool boundary = 0 )
    {
        index_map_type index_list = entityToLocal( top_dim, local_id, boundary );

        uint16_type matrix_size = 0;

        if ( index_list[0].size() != 0 )
            matrix_size +=index_list[0].size()*nbPtsPerVertex;

        if ( ( top_dim>=1 ) && ( index_list[1].size() != 0 ) )
            matrix_size +=index_list[1].size()*nbPtsPerEdge;

        if ( ( top_dim >=2 ) && ( index_list[2].size() != 0 ) )
            matrix_size +=index_list[2].size()*nbPtsPerFace;

        if ( ( top_dim ==3 ) && ( index_list[3].size() != 0 ) )
            matrix_size +=nbPtsPerVolume;

        points_type G ( Dim, matrix_size );

        for ( uint16_type i=0, p=0; i < top_dim+1; i++ )
        {
            if ( index_list[i].size() )
            {
                for ( uint16_type j=0; j < index_list[i].size(); j++ )
                {
                    points_type aux = interiorPointsById( i, index_list[i][j] );

                    ublas::subrange( G, 0, Dim, p, p+aux.size2() ) = aux;

                    p+=aux.size2();
                }
            }
        }

        return G;
    }

    index_map_type getEid ()
    {
        return M_eid;
    }

    std::vector<range_type> getPtE()
    {
        return M_pt_to_entity;
    }

    void setEid ( index_map_type eid )
    {
        M_eid = eid;
    }

    void setPtE ( std::vector<range_type> pt_ent )
    {
        M_pt_to_entity = pt_ent;
    }

    void addToEid ( uint16_type p, uint16_type q )
    {
        M_eid[p].push_back( q );
    }

    void addToPtE ( uint16_type p, range_type q )
    {
        M_pt_to_entity[p] = q;
    }

    FEELPP_DEFINE_VISITABLE();

private:

    index_map_type M_eid;
    std::vector<range_type> M_pt_to_entity;

};
}
#endif /* __PointSetInterpolation_H */
