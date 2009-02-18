/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-12

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file importergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-12
 */
#if 0
#include <life/lifefilters/importergmsh.hpp>
#endif

#include <boost/algorithm/string/trim.hpp>
#include <life/lifediscr/mesh.hpp>


namespace Life
{
template<typename MeshType>
void
ImporterGmsh<MeshType>::visit( mesh_type* mesh )
{
    if ( this->version() != "1.0" && this->version() != "2.0" )
        throw std::logic_error( "invalid gmsh file format version" );

    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit()] starts\n";

    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit( "  << mesh_type::nDim << "D )] starts\n";
    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit( "  << mesh_type::nDim << "D )] filename = " << this->filename() << "\n";

    std::ifstream __is ( this->filename().c_str() );

    char __buf[256];
    __is >> __buf;

    if ( this->version() == "2.0" && std::string( __buf ) == "$MeshFormat" )
        {
            std::string theversion;
            // version file-type(0=ASCII,1=BINARY) data-size(sizeof(double))
            __is >> theversion >> __buf >> __buf;
            LIFE_ASSERT( theversion == "2" )( theversion )( this->version() ).error( "invalid gmsh file format version ");
            // should be $EndMeshFormat
            __is >> __buf;
            LIFE_ASSERT( std::string( __buf ) == "$EndMeshFormat" )
                ( __buf )
                ( "$EndMeshFormat").error ( "invalid file format entry" );
            __is >> __buf;
            if ( std::string( __buf ) == "$PhysicalNames" )
                {
                    int nnames;
                    __is >> nnames;
                    for( int n = 0; n < nnames; ++n )
                        {
                            int id;
                            std::string name;
                            __is >> id >> name;
                            boost::trim( name );
                            boost::trim_if(name,boost::is_any_of("\""));

                            mesh->addMarkerName( std::make_pair( name, id ) );
                        }
                    __is >> __buf;
                    LIFE_ASSERT( std::string( __buf ) == "$EndPhysicalNames" )
                        ( __buf )
                        ( "$EndPhysicalNames").error ( "invalid file format entry" );
                    __is >> __buf;
                }
        }

    //
    // Read NODES
    //
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    LIFE_ASSERT( std::string( __buf ) == "$NOD" ||
                 std::string( __buf ) == "$Nodes" )
        ( __buf )
        ( "$NOD" )( "$Nodes" )
        .error("invalid nodes string in gmsh importer");

    uint __n;
    __is >> __n;
    Debug( 8011 ) << "number of nodes: " << __n;

    std::vector<double> __x( 3*__n );
    std::vector<bool> __isonboundary(__n);
    std::vector<uint> __whichboundary(__n);
    Debug( 8011 ) << "reading "<< __n << " nodes\n";
    std::map<int,int> itoii;
    for( uint __i = 0; __i < __n;++__i )
        {
            uint __ni;
            __is >> __ni
                 >> __x[3*__i]
                 >> __x[3*__i+1]
                 >> __x[3*__i+2];

            // stores mapping to be able to reorder the indices
            // so that they are contiguous
            itoii[__ni] = __i;
        }

    __is >> __buf;
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    // make sure that we have read all the points
    LIFE_ASSERT( std::string( __buf ) == "$ENDNOD" ||
                 std::string( __buf ) == "$EndNodes" )
        ( __buf )
        ( "$ENDNOD" )( "$EndNodes" ).error("invalid end nodes string in gmsh importer");


    //
    // Read ELEMENTS
    //
    __is >> __buf;
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    LIFE_ASSERT( std::string( __buf ) == "$ELM" ||
                 std::string( __buf ) == "$Elements" )
        ( __buf )
        ( "$ELM" )( "$Elements" )
        .error("invalid elements string in gmsh importer");

    uint __nele;
    __is >> __nele;

    Debug( 8011 ) << "number of elements: " << __nele << "\n";
    std::vector<std::vector<int> > __e(__nele);
    std::vector<int> __et(__nele);
    std::vector<GMSH_ENTITY> __etype( __nele );
    std::vector<int> __gt(16);
    __gt.assign( 16, 0 );

    int nptable[20];
    nptable[1]=2; // node line.
    nptable[2]=3; // node triangle.
    nptable[3]=4; // node quadrangle.
    nptable[4]=4; // node tetrahedron.
    nptable[5]=8; // node hexahedron.
    nptable[6]=6; // node prism.
    nptable[7]=5; // node pyramid.
    nptable[8]=3; // node second order line (2 nodes associated with the vertices and 1 with the edge).
    nptable[9]=6; // node second order triangle (3 nodes associated with the vertices and 3 with the edges).
    nptable[10]=9; // node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
    nptable[11]=10; // node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
    nptable[12]=27; // node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
    nptable[13]=18; //node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
    nptable[14]=14; // node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
    nptable[15]=1; // node point.
    nptable[16]=8; // node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
    nptable[17]=20; // node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
    nptable[18]=15; // node second order prism (6 nodes associated with the vertices and 9 with the edges).
    nptable[19]=13; // node second order pyramid (5 nodes associated with the vertices and 8 with the edges).

    for( uint __i = 0; __i < __nele;++__i )
        {
            int __ne, __t, __tag, __np, __dummy;

            if ( this->version() == "1.0" )
                {
                    __is >> __ne  // elm-number
                         >> __t // elm-type
                         >> __tag // reg-phys
                         >> __dummy // reg-elem
                         >> __np; // number-of-nodes
                    LIFE_ASSERT( __np == nptable[__t] )( __np )( __t )( nptable[__t] ).error( "invalid number of nodes" );
                }
            else if ( this->version() == "2.0" )
                {
                    int __ntag;
                    __is >> __ne  // elm-number
                         >> __t // elm-type
                         >> __ntag // number-of-tags
                         >> __tag; // reg-phys
                    LIFE_ASSERT( __ntag >= 2 )( __ntag )( __tag )( __dummy ).warn( "invalid number of tags" );
                    for( int nt = 1; nt < __ntag; ++nt )
                        __is >> __dummy;

                    __np = nptable[__t];
                }

            ++__gt[ __t ];
            __etype[__i] = GMSH_ENTITY(__t);
            __et[__i] = __tag;
            __e[__i].resize( __np );
            int __p = 0;
            while ( __p != __np )
                {
                    __is >> __e[__i][__p];

                    // reorder the nodes since they may not have had a contiguous ordering
                    __e[__i][__p] = itoii[ __e[__i][__p]];

                    ++__p;
                }
        }

    // make sure that we have read everything
    __is >> __buf;
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    LIFE_ASSERT( std::string( __buf ) == "$ENDELM" ||
                 std::string( __buf ) == "$EndElements" )
        ( __buf )
        ( "$ENDELM" )( "$EndElements" ).error("invalid end elements string in gmsh importer");

    // read physical names
    if ( this->version() == "2.0" )
        {

        }

    //
    // FILL Mesh Data Structure
    //
    Debug( 8011 ) << "number of edges= " << __gt[1] << "\n";

    __isonboundary.assign( __n, false );
    __whichboundary.assign( __n, 0 );
    for( uint __i = 0; __i < __nele;++__i )
        {
            switch( __etype[__i] )
                {
                case GMSH_POINT:
                    if ( mesh_type::nDim == 1 )
                        {
                            __isonboundary[ __e[__i][0] ] = true;

                            __whichboundary[__e[__i][0] ] = __et[__i];
                        }
                    break;
                case GMSH_LINE:
                case GMSH_LINE_2:
                case GMSH_LINE_3:
                case GMSH_LINE_4:
                case GMSH_LINE_5:
                    if ( mesh_type::nDim == 2 )
                        {
                            for( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
                                {
                                    __isonboundary[ __e[__i][jj] ] = true;
                                    __whichboundary[__e[__i][jj] ] = __et[__i];
                                }
                        }
                    break;
                case GMSH_TRIANGLE:
                case GMSH_TRIANGLE_2:
                case GMSH_TRIANGLE_3:
                case GMSH_TRIANGLE_4:
                case GMSH_TRIANGLE_5:
                    if ( mesh_type::nDim == 3 )
                        {
                            for( uint16_type jj = 0; jj < npoints_per_face; ++jj )
                                {
                                    __isonboundary[ __e[__i][jj] ] = true;
                                    __whichboundary[__e[__i][jj] ] = __et[__i];
                                }
                        }
                    break;
                default:
                    break;
                }
        }

    // add the points to the mesh
    for( uint __i = 0; __i < __n;++__i )
        {
            node_type __n( mesh_type::nRealDim );
            for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
                __n[j] = __x[3*__i+j];
            point_type __pt( __i,__n, __isonboundary[ __i ] );
            __pt.setOnBoundary( __isonboundary[ __i ] );
            __pt.marker().value( __whichboundary[__i] );
            mesh->addPoint( __pt );
        }

    _M_n_vertices.resize( __n );
    _M_n_vertices.assign( __n, 0 );
    _M_n_b_vertices.resize( __n );
    _M_n_b_vertices.assign( __n, 0 );

    // add the element to the mesh
    for( uint __i = 0; __i < __nele;++__i )
        {
            switch( __etype[__i] )
                {
                    // Points
                case GMSH_POINT:
                    addPoint( mesh, __e[__i], __et[__i], __etype[__i] );
                    break;

                    // Edges
                case GMSH_LINE:
                case GMSH_LINE_2:
                case GMSH_LINE_3:
                case GMSH_LINE_4:
                case GMSH_LINE_5:
                    addEdge( mesh, __e[__i], __et[__i], __etype[__i] );
                    break;

                    // Faces
                case GMSH_TRIANGLE:
                case GMSH_TRIANGLE_2:
                case GMSH_TRIANGLE_3:
                case GMSH_TRIANGLE_4:
                case GMSH_TRIANGLE_5:
                case GMSH_QUADRANGLE:
                case GMSH_QUADRANGLE_2:
                    addFace( mesh, __e[__i], __et[__i], __etype[__i] );
                    break;

                    // Volumes
                case GMSH_TETRAHEDRON:
                case GMSH_TETRAHEDRON_2:
                case GMSH_TETRAHEDRON_3:
                case GMSH_TETRAHEDRON_4:
                case GMSH_TETRAHEDRON_5:
                case GMSH_HEXAHEDRON:
                case GMSH_HEXAHEDRON_2:

                    addVolume( mesh, __e[__i], __et[__i], __etype[__i] );
                    break;
                default:
                    break;
                }
        } // loop over geometric entities in gmsh file (can be elements or faces)

    mesh->setNumVertices( std::accumulate( _M_n_vertices.begin(), _M_n_vertices.end(), 0 ) );

    Debug( 8011 ) << "done with reading and creating mesh from gmsh file\n";
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type )
{
    addPoint( mesh, __e, tag, type, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY /*type*/, mpl::int_<1> )
{
    face_type pf;
    pf.setId( mesh->numFaces() );
    pf.marker().value(  tag  );
    pf.setPoint( 0, mesh->point( __e[0] ) );

    _M_n_vertices[ __e[0] ] = 1;

    _M_n_b_vertices[ __e[0] ] = 1;

    pf.setOnBoundary( true );
    face_iterator fit;
    bool inserted;
    boost::tie( fit, inserted ) = mesh->addFace( pf );
    Debug( 8011 ) << "added point on boundary ("
                  << fit->isOnBoundary() << ") with id :" << fit->id() << " and marker " << pf.marker()
                  << " n1: " << mesh->point( __e[0] ).node() << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*/*mesh*/, std::vector<int> const& /*__e*/, int /*tag*/, GMSH_ENTITY /*type*/, mpl::int_<2> )
{
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*/*mesh*/, std::vector<int> const& /*__e*/, int /*tag*/, GMSH_ENTITY /*type*/, mpl::int_<3> )
{
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type )
{
    addEdge( mesh, __e, tag, type, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type*mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<1> )
{
    element_type e;

    e.marker().value(  tag  );
    if ( type == GMSH_LINE ||
         type == GMSH_LINE_2 ||
         type == GMSH_LINE_3 ||
         type == GMSH_LINE_4 ||
         type == GMSH_LINE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_edge; ++j )
                e.setPoint( jj, mesh->point( __e[jj] ) );
        }

    mesh->addElement( e );

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    Debug( 8011 ) << "added edge with id :" << e.id()
                  << " n1: " << mesh->point( __e[0] ).node()
                  << " n2: " << mesh->point( __e[1] ).node() << "\n";
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<2> )
{
    face_type pf;
    pf.setId( mesh->numFaces() );
    pf.marker().value(  tag  );

    if ( type == GMSH_LINE ||
         type == GMSH_LINE_2 ||
         type == GMSH_LINE_3 ||
         type == GMSH_LINE_4 ||
         type == GMSH_LINE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_edge; ++j )
                pf.setPoint( jj, mesh->point( __e[jj] ) );
        }

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;

    _M_n_b_vertices[ __e[0] ] = 1;
    _M_n_b_vertices[ __e[1] ] = 1;

    pf.setOnBoundary( true );
    bool inserted;
    face_iterator fit;
    boost::tie( fit, inserted ) = mesh->addFace( pf );
    Debug( 8011 ) << "added edge on boundary ("
                  << fit->isOnBoundary() << ") with id :" << fit->id()
                  << " n1: " << mesh->point( __e[0] ).node()
                  << " n2: " << mesh->point( __e[1] ).node() << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<3> )
{}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type )
{
    addFace( mesh, __e, tag, type, mpl::int_<mesh_type::nDim>() );
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<1> )
{}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<2> )
{
    element_type pf;
    pf.marker().value(  tag  );

    if ( type == GMSH_QUADRANGLE ||
         type == GMSH_TRIANGLE ||
         type == GMSH_TRIANGLE_2 ||
         type == GMSH_TRIANGLE_3 ||
         type == GMSH_TRIANGLE_4 ||
         type == GMSH_TRIANGLE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_element; ++j )
                pf.setPoint( jj, mesh->point( __e[jj] ) );
        }

    mesh->addElement( pf );

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    _M_n_vertices[ __e[2] ] = 1;
    if ( type == GMSH_QUADRANGLE )
        _M_n_vertices[ __e[3] ] = 1;
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<3> )
{
    face_type pf;
    pf.setId( mesh->numFaces() );
    pf.marker().value(  tag  );

    if ( type == GMSH_QUADRANGLE ||
         type == GMSH_TRIANGLE ||
         type == GMSH_TRIANGLE_2 ||
         type == GMSH_TRIANGLE_3 ||
         type == GMSH_TRIANGLE_4 ||
         type == GMSH_TRIANGLE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_face; ++j )
                pf.setPoint( jj, mesh->point( __e[jj] ) );
        }

    pf.setOnBoundary( true );
    bool inserted;
    face_iterator fit;
    boost::tie( fit, inserted ) = mesh->addFace( pf );
    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    _M_n_vertices[ __e[2] ] = 1;

    if ( type == GMSH_QUADRANGLE )
        _M_n_vertices[ __e[3] ] = 1;
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type )
{
    addVolume( mesh, __e, tag, type, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<1> )
{}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<2> )
{}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<3> )
{
    element_type pv;

    pv.marker().value(  tag  );
    pv.setPoint( 0, mesh->point( __e[0] ) );
    pv.setPoint( 1, mesh->point( __e[1] ) );
    pv.setPoint( 2, mesh->point( __e[2] ) );
    pv.setPoint( 3, mesh->point( __e[3] ) );
    if ( type == GMSH_HEXAHEDRON )
        {
            pv.setPoint( 4, mesh->point( __e[4] ) );
            pv.setPoint( 5, mesh->point( __e[5] ) );
            pv.setPoint( 6, mesh->point( __e[6] ) );
            pv.setPoint( 7, mesh->point( __e[7] ) );
        }
    if ( type == GMSH_TETRAHEDRON_2 )
        {
            // warning the numerotation differs between life and gmsh
            // for the next three points
            pv.setPoint( 4, mesh->point( __e[5] ) );
            pv.setPoint( 5, mesh->point( __e[6] ) );
            pv.setPoint( 6, mesh->point( __e[4] ) );

            pv.setPoint( 7, mesh->point( __e[7] ) );
            pv.setPoint( 8, mesh->point( __e[8] ) );
            pv.setPoint( 9, mesh->point( __e[9] ) );
        }
    mesh->addElement( pv );

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    _M_n_vertices[ __e[2] ] = 1;
    _M_n_vertices[ __e[3] ] = 1;
    if ( type == GMSH_HEXAHEDRON )
        {
            _M_n_vertices[ __e[4] ] = 1;
            _M_n_vertices[ __e[5] ] = 1;
            _M_n_vertices[ __e[6] ] = 1;
            _M_n_vertices[ __e[7] ] = 1;
        }

}

#if 0
//
// Explicit instantiations
//
template class ImporterGmsh<Mesh<Simplex<1,1> > >;
template class ImporterGmsh<Mesh<Simplex<1,1,2> > >;
template class ImporterGmsh<Mesh<Simplex<2,1> > >;
template class ImporterGmsh<Mesh<Simplex<2,1,3> > >;
template class ImporterGmsh<Mesh<Simplex<3,1> > >;
template class ImporterGmsh<Mesh<Simplex<2,2> > >;
template class ImporterGmsh<Mesh<Simplex<3,2> > >;
template class ImporterGmsh<Mesh<SimplexProduct<1,1> > >;
template class ImporterGmsh<Mesh<SimplexProduct<2,1> > >;
template class ImporterGmsh<Mesh<SimplexProduct<3,1> > >;
}
#endif // 0
