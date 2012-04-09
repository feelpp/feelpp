/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-06-30

  Copyright (C) 2005,2006 EPFL

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
   \file importergambit.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-06-30
 */
#ifndef __GAMBIT_HPP
#define __GAMBIT_HPP 1

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/importer.hpp>

namespace Feel
{
/// \cond detail
namespace gambit
{
namespace Element
{
enum
{
    EDGE = 1,
    QUADRILATERAL = 2,
    TRIANGLE = 3,
    BRICK = 4,
    WEDGE = 5,
    TETRAHEDRON = 6,
    PYRAMID = 7
};
}
namespace Material
{
enum
{
    UNDEFINED = 0,
    CONJUGATE = 1,
    FLUID = 2,
    POROUS = 3,
    SOLID = 4,
    DEFORMABLE = 5
};
}
namespace IBCode
{
enum
{
    UNSPECIFIED = 0,
    AXIS = 1,
    CONJUGATE = 2,
    CONVECTION = 3,
    CYCLIC = 4,
    DEAD = 5,
    ELEMENT_SIDE = 6,
    ESPECIES = 7,
    EXHAUST_FAN = 8,
    FAN = 9,
    FREE_SURFACE = 10,
    GAP = 11,
    INFLOW = 12,
    INLET = 13,
    INLET_VENT = 14,
    INTAKE_FAN = 15,
    INTERFACE = 16,
    INTERIOR = 17,
    INTERNAL = 18,
    LIVE = 19,
    MASS_FLOW_INLET = 20,
    MELT = 21,
    MELT_INTERFACE = 22,
    MOVING_BOUNDARY = 23,
    NODE = 24,
    OUTFLOW = 25,
    OUTLET = 26,
    OUTLET_VENT = 27,
    PERIODIC = 28,
    PLOT = 29,
    POROUS = 30,
    POROUS_JUMP = 31,
    PRESSURE = 32,
    PRESSURE_FAR_FIELD = 33,
    PRESSURE_INFLOW = 34,
    PRESSURE_INLET = 35,
    PRESSURE_OUTFLOW = 36,
    PRESSURE_OUTLET = 37,
    RADIATION = 38,
    RADIATOR = 39,
    RECIRCULATION_INLET = 40,
    RECIRCULATION_OUTLET = 41,
    SLIP = 42,
    SREACTION = 43,
    SURFACE = 44,
    SYMMETRY = 45,
    TRACTION = 46,
    TRAJECTORY = 47,
    VELOCITY = 48,
    VELOCITY_INLET = 49,
    VENT = 50,
    WALL = 51,
    SPRING = 52
};
}

struct quad4
{
    static const int face[][2];
};


struct quad8
{
    static const int face[][3];
};

struct tria3
{
    static const int face[][2];
};

struct tria6
{
    static const int face[][3];
};


struct tetra4
{
    static const int face[][3];
};

struct tetra10
{
    static const int face[][6];
};

typedef std::vector<double> nodes_type;
typedef std::vector<boost::tuple<bool, int> >  nodes_boundary_type;

typedef std::vector<boost::tuple<int, int, std::vector<int>, boost::tuple<int, int> > > elements_type;

bool read( std::string const& filename,
           nodes_type& nodes,
           nodes_boundary_type& boundary,
           elements_type& __elements );


} // gambit namespace
/// \endcond detail
/**
 * \class ImporterGambit
 * \brief gambit(fluent mesh generator) importer class
 *
 * the importer concept follows the visitor pattern
 *
 * \code
 * ImporterGambit<mesh_type> import( "mesh.msh");
 * mesh.accept( import );
 * \endcode
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 * @see Importer, gmsh
 */
template<typename MeshType>
class ImporterGambit
    :
public Importer<MeshType>

{
    typedef Importer<MeshType> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef typename super::mesh_type mesh_type;
    typedef typename super::point_type point_type;
    typedef typename super::node_type node_type;
    typedef typename super::edge_type edge_type;
    typedef typename super::face_type face_type;
    typedef typename super::element_type element_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    ImporterGambit()
        :
        super( GAMBIT )
    {}

    ImporterGambit( std::string const& fname )
        :
        super( fname, GAMBIT )
    {}
    ImporterGambit( ImporterGambit const & i )
        :
        super( i )
    {}
    ~ImporterGambit()
    {}


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh );

    //@}
};

template <typename MeshType>
void
ImporterGambit<MeshType>::visit( mesh_type* mesh )
{
    gambit::nodes_type nodes;
    gambit::nodes_boundary_type boundary;
    gambit::elements_type __elements;

    // read the file and store the info
    gambit::read( this->filename(), nodes, boundary, __elements );


    //mesh->setMaxNumPoints( nodes.size()/2, false );
    //mesh->numPoints() = nodes.size()/2;
    //mesh->setNumBPoints( std::count( __isonboundary.begin(), __isonboundary.end(), true ) );

    uint16_type ncoord = mesh_type::nDim;

    // add the points to the mesh
    for ( uint __i = 0; __i < nodes.size()/ncoord; ++__i )
    {
        node_type __n( ncoord );

        for ( uint16_type j = 0; j < ncoord; ++j )
            __n[j] = nodes[ncoord*__i+j];

        point_type __pt( __i, __n, boundary[ __i ].get<0>() );
        __pt.setMarker( boundary[__i].get<1>() );
        mesh->addPoint( __pt );
    }

    std::vector<int> __n_vertices( nodes.size()/ncoord );
    __n_vertices.assign( nodes.size()/ncoord, 0 );

    std::vector<int> __n_b_vertices(  nodes.size()/ncoord  );
    __n_b_vertices.assign(  nodes.size()/ncoord , 0 );

    // add the element to the mesh
    for ( size_type __i = 0; __i < __elements.size(); ++__i )
    {
        //pe = &( mesh->addElement() );
        element_type pe;

        pe.marker().assign( ( __elements[__i].get<1>() ) );

        // this handle linear and quadratic elements
        for ( uint n = 0; n < __elements[__i].get<2>().size(); ++n )
        {
            pe.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[n] ) );
            __n_vertices[ boost::get<2>( __elements[__i] )[n] ] = 1;
        }

        mesh->addElement( pe );

        // add face of the element if it is on a boundary of
        // the domain
        if ( __elements[__i].get<3>().get<0>() != -1 )
        {
            face_type pf;
            pf.setOnBoundary( true );

            // local id of the face in the element
            uint16_type numface = __elements[__i].get<3>().get<0>();
            Debug( 8012 ) << "adding face " << numface
                          << " on boundary : " << __elements[__i].get<3>().get<1>() << "\n";

            pf.marker().assign( ( __elements[__i].get<3>().get<1>() ) );

            if ( mesh_type::nDim == 2 && __elements[__i].get<2>().size() == 3 ) // linear triangle 3-node
            {
                Debug( 8012 ) << "setting for edge from triangle  : " << __i << " face : " << numface << "\n";

                for ( int n = 0; n < 3; ++n )
                {
                    pf.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[gambit::tria3::face[numface][n]] ) );
                    __n_vertices[boost::get<2>( __elements[__i] )[gambit::tria3::face[numface][n]] ] = 1;
                    __n_b_vertices[boost::get<2>( __elements[__i] )[gambit::tria3::face[numface][n]] ] = 1;
                }

                Debug( 8012 ) << "added face on boundary ("
                              << pf.isOnBoundary() << ") with id :" << pf.id()
                              << " n1: " << mesh->point( boost::get<2>( __elements[__i] )[gambit::tria3::face[numface][0]]   ).node()
                              << " n2: " << mesh->point( boost::get<2>( __elements[__i] )[gambit::tria3::face[numface][1]]  ).node() << "\n";
            }

            else if ( mesh_type::nDim == 2 && __elements[__i].get<2>().size() == 6 ) // quadratic triangle 6-node
                for ( int n = 0; n < 3; ++n )
                    pf.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[gambit::tria6::face[numface][n]] ) );

            else if ( mesh_type::nDim == 2 && __elements[__i].get<2>().size() == 4 ) // linear quadrangle 4-node
                for ( int n = 0; n < 2; ++n )
                    pf.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[gambit::quad4::face[numface][n]] ) );

            else if ( mesh_type::nDim == 2 && __elements[__i].get<2>().size() == 8 ) // quadratic triangle 8-node
                for ( int n = 0; n < 3; ++n )
                    pf.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[gambit::quad8::face[numface][n]] ) );

            else if ( mesh_type::nDim == 3 && __elements[__i].get<2>().size() == 4 ) // linear tetra 4-node
            {
                for ( int n = 0; n < 3; ++n )
                {
                    int g2l[]= { 3, 2, 1, 0 };
                    int local_index = Feel::details::tetra<1/*mesh_type::nOrder*/>::f2p( g2l[numface], n );//gambit::tetra4::face[numface][n];
                    int index = boost::get<2>( __elements[__i] )[local_index];
                    pf.setPoint( n, mesh->point( index ) );
                    __n_vertices[index] = 1;
                    __n_b_vertices[index] = 1;
                }

                Debug( 8012 ) << "added face on boundary ("
                              << pf.isOnBoundary() << ") with id :" << pf.id()
                              << " n1: " << mesh->point( boost::get<2>( __elements[__i] )[gambit::tetra4::face[numface][0]]   ).node()
                              << " n2: " << mesh->point( boost::get<2>( __elements[__i] )[gambit::tetra4::face[numface][1]]  ).node() << "\n";
            }

            else if ( mesh_type::nDim == 3 && __elements[__i].get<2>().size() == 10 ) // quadratic tetra 10-node
                for ( int n = 0; n < 6; ++n )
                    pf.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[gambit::tetra10::face[numface][n]] ) );


            mesh->addFace( pf );


        }

    }
}

#if 0
template <typename MeshType>
void
ImporterGambit<MeshType>::visit( mesh_type* mesh,
                                 mpl::int_<3> )
{
    gambit::nodes_type nodes;
    gambit::nodes_boundary_type boundary;
    gambit::elements_type __elements;

    // read the file and store the info
    gambit::read( this->filename(), nodes, boundary, __elements );

    Debug( 8012 ) << " o- inserting " << nodes.size()/3 << " points to mesh" << "\n";

    // add the points to the mesh
    for ( uint __i = 0; __i < nodes.size()/3; ++__i )
    {
        node_type __n( 3 );
        __n[0] = nodes[3*__i];
        __n[1] = nodes[3*__i+1];
        __n[2] = nodes[3*__i+2];
        point_type __pt( __i, __n, boundary[ __i ].get<0>() );
        __pt.marker().assign( boundary[__i].get<1>() );
        mesh->addPoint( __pt );
    }

    std::vector<int> __n_vertices( nodes.size()/3 );
    __n_vertices.assign( nodes.size()/3, 0 );

    std::vector<int> __n_b_vertices(  nodes.size()/3  );
    __n_b_vertices.assign(  nodes.size()/3 , 0 );

    Debug( 8012 ) << " o- inserting " << __elements.size() << " elements to mesh" << "\n";

    // add the element to the mesh
    for ( uint __i = 0; __i < __elements.size(); ++__i )
    {
        element_type pv;
        pv.setId( __i );

        pv.marker().assign( ( __elements[__i].get<1>() ) );

        // this handle linear and quadratic elements
        for ( uint n = 0; n < __elements[__i].get<2>().size(); ++n )
        {
            pv.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[n] ) );
            __n_vertices[ boost::get<2>( __elements[__i] )[n] ] = 1;
        }

        Debug( 8013 ) << " o- inserted element " << pv.id() << " in mesh ("  << pv.marker() << ")\n";

        mesh->addElement( pv );


        if ( __elements[__i].get<3>().get<0>() != -1 )
        {
            Debug( 8013 ) << "adding face on boundary : " << __elements[__i].get<3>().get<1>() << "\n";
            face_type pf;

            // local id of the face in the element
            int numface = __elements[__i].get<3>().get<0>();

            pf.marker().assign( ( __elements[__i].get<3>().get<1>() ) );


            if ( __elements[__i].get<2>().size() == 4 ) // linear tetra 4-node
                for ( int n = 0; n < 3; ++n )
                {
                    int g2l[]= { 3, 2, 1, 0 };
                    int local_index = Feel::details::tetra::f2p( g2l[numface], n );//gambit::tetra4::face[numface][n];
                    int index = boost::get<2>( __elements[__i] )[local_index];
                    pf.setPoint( n, mesh->point( index ) );
                    __n_vertices[index] = 1;
                    __n_b_vertices[index] = 1;
                }

            else if ( __elements[__i].get<2>().size() == 10 ) // quadratic tetra 10-node
                for ( int n = 0; n < 6; ++n )
                    pf.setPoint( n, mesh->point( boost::get<2>( __elements[__i] )[gambit::tetra10::face[numface][n]] ) );

            pf.setOnBoundary( true );
            mesh->addFace( pf );
        }
    }

    mesh->setNumVertices( std::accumulate( __n_vertices.begin(), __n_vertices.end(), 0 ) );
    Debug( 8012 ) << " o- number of vertices : " << mesh->numVertices() << "\n";
}
#endif

} // Feel namespace
#endif /* __GAMBIT_HPP */
