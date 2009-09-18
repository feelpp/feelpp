/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-16

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file importergmsh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-16
 */

#ifndef __ImporterGmsh_H
#define __ImporterGmsh_H 1

#include <iostream>
#include <fstream>
#include <iomanip>

#include <life/lifefilters/importer.hpp>
#include <boost/algorithm/string/trim.hpp>


namespace Life
{
/**
 * \enum GMSH_ENTITY
 *
 * enumerate the various elements available in gmsh
 */
enum GMSH_ENTITY
    {
        GMSH_LINE = 1, //!< Line (2 nodes).
        GMSH_TRIANGLE = 2, //!< Triangle (3 nodes).
        GMSH_QUADRANGLE = 3,  //!< Quadrangle (4 nodes).
        GMSH_TETRAHEDRON = 4, //!< Tetrahedron (4 nodes).
        GMSH_HEXAHEDRON = 5, //!< Hexahedron (8 nodes).
        GMSH_PRISM = 6,  //!< Prism (6 nodes).
        GMSH_PYRAMID = 7, //!< Pyramid (5 nodes).
        GMSH_LINE_2 = 8,  //!< Second order line (3 nodes: 2 associated
        //with the vertices and 1 with the edge).
        GMSH_TRIANGLE_2 = 9, //!< Second order triangle (6 nodes: 3
        //associated with the vertices and 3 with the
        //edges).
        GMSH_QUADRANGLE_2 = 10, //!<Second order quadrangle (9 nodes: 4
        //associated with the vertices, 4 with the
        //edges and 1 with the face).
        GMSH_TETRAHEDRON_2 = 11, //!< Second order tetrahedron (10 nodes:
        //4 associated with the vertices and 6
        //with the edges).
        GMSH_HEXAHEDRON_2 = 12,  //!< Second order hexahedron (27 nodes: 8
        //associated with the vertices, 12 with
        //the edges, 6 with the faces and 1 with
        //the volume).
        GMSH_PRISM_2 = 13, //!<Second order prism (18 nodes: 6 associated
        //with the vertices, 9 with the edges and 3 with
        //the quadrangular faces).
        GMSH_PYRAMID_2 = 14, //!<Second order pyramid (14 nodes: 5
        //associated with the vertices, 8 with the
        //edges and 1 with the quadrangular face).
        GMSH_POINT = 15, //!< Point (1 node).

        GMSH_TRIANGLE_INCOMPLETE_3=20, //!< triangle of order 3
        GMSH_TRIANGLE_3=21, //!< triangle of order 3
        GMSH_TRIANGLE_INCOMPLETE_4=22, //!< triangle of order 4
        GMSH_TRIANGLE_4=23, //!< triangle of order 4
        GMSH_TRIANGLE_INCOMPLETE_5=24, //!< triangle of order 5
        GMSH_TRIANGLE_5=25, //!< triangle of order 5

        GMSH_LINE_3=26, //!< line of order 3
        GMSH_LINE_4=27, //!< line of order 4
        GMSH_LINE_5=28, //!< line of order 5

        GMSH_TETRAHEDRON_3=29, //!< tetra of order 3
        GMSH_TETRAHEDRON_4=30, //!< tetra of order 4
        GMSH_TETRAHEDRON_5=31, //!< tetra of order 5
    };

/**
 * \class ImporterGmsh
 * \brief gmsh importer class
 *
 * the importer concept follows the visitor pattern
 *
 * \code
 * typename Mesh2D<LinearTetra> mesh_type;
 * mesh_type mesh;
 *
 * ImporterGmsh<mesh_type> import( "mesh.msh");
 * mesh.accept( import );
 * \endcode
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 */
template<typename MeshType>
class ImporterGmsh
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

    typedef typename mesh_type::face_iterator face_iterator;
    //@}

    //@{
    /** @name Constants
     */
    static const uint16_type npoints_per_edge = ( edge_type::numVertices*edge_type::nbPtsPerVertex+
                                                  edge_type::numEdges*edge_type::nbPtsPerEdge+
                                                  edge_type::numFaces*edge_type::nbPtsPerFace );

    static const uint16_type npoints_per_face = ( face_type::numVertices*face_type::nbPtsPerVertex+
                                                  face_type::numEdges*face_type::nbPtsPerEdge+
                                                  face_type::numFaces*face_type::nbPtsPerFace );

    static const uint16_type npoints_per_element = element_type::numPoints;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    ImporterGmsh()
        :
        super( GMSH ),
        _M_version( "2.1" )
    {
        showMe();
    }

    explicit ImporterGmsh( std::string const& fname, std::string version = "2.1" )
        :
        super( fname, GMSH ),
        _M_version( version )
    {
        showMe();
    }
    ImporterGmsh( ImporterGmsh const & i )
        :
        super( i ),
        _M_version( i._M_version )
    {
        showMe();
    }
    ~ImporterGmsh()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @return the file format version
     */
    std::string version() const { return _M_version; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setVersion( std::string const& version ) { _M_version = version; }

    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh );

    void showMe() const;

    //@}



protected:

private:

    void addPoint( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type );
    void addPoint( mesh_type*mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY /*type*/, mpl::int_<1> );
    void addPoint( mesh_type*/*mesh*/, std::vector<int> const& /*__e*/, int /*tag*/, GMSH_ENTITY /*type*/, mpl::int_<2> );
    void addPoint( mesh_type*/*mesh*/, std::vector<int> const& /*__e*/, int /*tag*/, GMSH_ENTITY /*type*/, mpl::int_<3> );

    void addEdge( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type );
    void addEdge( mesh_type*mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<1> );
    void addEdge( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<2> );

    void addEdge( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<3> );

    void addFace( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type );
    void addFace( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<1> );
    void addFace( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<2> );
    void addFace( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<3> );

    void addVolume( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type );
    void addVolume( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<1> );
    void addVolume( mesh_type*, std::vector<int> const&, int, GMSH_ENTITY, mpl::int_<2> );
    void addVolume( mesh_type* mesh, std::vector<int> const& __e, int tag, GMSH_ENTITY type, mpl::int_<3> );

private:

    std::string _M_version;
    std::vector<int> _M_n_vertices;
    std::vector<int> _M_n_b_vertices;

};



template<typename MeshType>
void
ImporterGmsh<MeshType>::showMe() const
{
    Log() << "[ImporterGmsh::showMe] npoints_per_element = " << npoints_per_element << "\n";
    Log() << "[ImporterGmsh::showMe]    npoints_per_face = " << npoints_per_face << "\n";
    Log() << "[ImporterGmsh::showMe]    npoints_per_edge = " << npoints_per_edge << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::visit( mesh_type* mesh )
{
    if ( this->version() != "1.0" &&
         this->version() != "2.0" &&
         this->version() != "2.1" )
        throw std::logic_error( "invalid gmsh file format version" );

    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit()] starts\n";

    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit( "  << mesh_type::nDim << "D )] starts\n";
    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit( "  << mesh_type::nDim << "D )] filename = " << this->filename() << "\n";

    std::ifstream __is ( this->filename().c_str() );

    char __buf[256];
    __is >> __buf;

    if ( ( (this->version() == "2.0") ||
           (this->version() == "2.1") ) &&
         std::string( __buf ) == "$MeshFormat" )
        {
            std::string theversion;
            // version file-type(0=ASCII,1=BINARY) data-size(sizeof(double))
            __is >> theversion >> __buf >> __buf;
            LIFE_ASSERT( boost::lexical_cast<double>( theversion ) >= 2 )( theversion )( this->version() ).warn( "invalid gmsh file format version ");
            // should be $EndMeshFormat
            __is >> __buf;
            LIFE_ASSERT( std::string( __buf ) == "$EndMeshFormat" )
                ( __buf )
                ( "$EndMeshFormat").error ( "invalid file format entry" );
            __is >> __buf;
            Debug() << "[importergmsh] " << __buf << " (expect $PhysicalNames)\n";
            if ( std::string( __buf ) == "$PhysicalNames" )
                {
                    int nnames;
                    __is >> nnames;
                    for( int n = 0; n < nnames; ++n )
                        {
                            int id, topodim;
                            std::string name;
                            if ( this->version() == "2.1" )
                                __is >> topodim >> id >> name;
                            else if ( this->version() == "2.0" )
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
    std::vector<int> __gt(32);
    __gt.assign( 32, 0 );

    int nptable[32];
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

    // ho triangle
    nptable[20]=10;
    nptable[21]=10;
    nptable[22]=15;
    nptable[23]=15;
    nptable[24]=21;
    nptable[25]=21;

    // ho line
    nptable[26]=4;
    nptable[27]=5;
    nptable[28]=6;

    // ho tetra
    nptable[29]=20;
    nptable[30]=35;
    nptable[31]=56;

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
            else if ( boost::lexical_cast<double>( this->version()) >= 2  )
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

                    Debug( 8011 ) << "element type: " << __t << " nb pts: " << __np << "\n";
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
    if ( boost::lexical_cast<double>( this->version()) >= 2  )
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
            __pt.marker().assign( __whichboundary[__i] );
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
    pf.marker().assign(  tag  );
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

    e.marker().assign(  tag  );
    if ( type == GMSH_LINE ||
         type == GMSH_LINE_2 ||
         type == GMSH_LINE_3 ||
         type == GMSH_LINE_4 ||
         type == GMSH_LINE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_element; ++jj )
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
    pf.marker().assign(  tag  );

    if ( type == GMSH_LINE ||
         type == GMSH_LINE_2 ||
         type == GMSH_LINE_3 ||
         type == GMSH_LINE_4 ||
         type == GMSH_LINE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
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
    pf.marker().assign(  tag  );

    if ( type == GMSH_QUADRANGLE ||
         type == GMSH_TRIANGLE ||
         type == GMSH_TRIANGLE_2 ||
         type == GMSH_TRIANGLE_3 ||
         type == GMSH_TRIANGLE_4 ||
         type == GMSH_TRIANGLE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_element; ++jj )
                {
                    //std::cout << "jj = " << jj << "\n";
                    if ( jj < element_type::numVertices*element_type::nbPtsPerVertex )
                        pf.setPoint( jj, mesh->point( __e[jj] ) );
                    else if ( (jj >= element_type::numVertices*element_type::nbPtsPerVertex) &&
                              (jj < (element_type::numVertices*element_type::nbPtsPerVertex + element_type::numEdges*element_type::nbPtsPerEdge ) ))
                        {
                            const uint16_type nbPtsPerEdge = (element_type::nbPtsPerEdge==0)?1:element_type::nbPtsPerEdge;
                            uint16_type edge_id = ( jj - element_type::numVertices*element_type::nbPtsPerVertex ) / nbPtsPerEdge;
                            //std::cout << "edge_id = " << edge_id << "\n";
                            if ( edge_id == 0 )
                                pf.setPoint( jj+(element_type::numEdges-1)*element_type::nbPtsPerEdge, mesh->point( __e[jj] ) );
                            else if ( edge_id == 1 )
                                pf.setPoint( jj-element_type::nbPtsPerEdge, mesh->point( __e[jj] ) );
                            else if ( edge_id == 2 )
                                pf.setPoint( jj-element_type::nbPtsPerEdge, mesh->point( __e[jj] ) );
                        }
                    // face interior pts when order <= 4
                    else if ( mesh_type::nOrder < 5 )
                        pf.setPoint( jj, mesh->point( __e[jj] ) );
                    // face interior pts when order == 5 (needs re-ordering)
                    else
                        {
                            uint16_type pt_id_start = (element_type::numVertices*element_type::nbPtsPerVertex +
                                                       element_type::numEdges*element_type::nbPtsPerEdge );
                            uint16_type pt_id = ( jj - pt_id_start );
                            //uint16_type pt_id_mod_3 = pt_id % 3;
                            if ( pt_id == 0 )
                                pf.setPoint( jj, mesh->point( __e[jj] ) );
                            if ( pt_id == 1 )
                                pf.setPoint( pt_id_start+2, mesh->point( __e[jj] ) );
                            if ( pt_id == 2 )
                                pf.setPoint( pt_id_start+5, mesh->point( __e[jj] ) );
                            if ( pt_id == 3 )
                                pf.setPoint( pt_id_start+1, mesh->point( __e[jj] ) );
                            if ( pt_id == 4 )
                                pf.setPoint( pt_id_start+4, mesh->point( __e[jj] ) );
                            if ( pt_id == 5 )
                                pf.setPoint( pt_id_start+3, mesh->point( __e[jj] ) );
                        }
                }
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
    pf.marker().assign(  tag  );

    if ( type == GMSH_QUADRANGLE ||
         type == GMSH_TRIANGLE ||
         type == GMSH_TRIANGLE_2 ||
         type == GMSH_TRIANGLE_3 ||
         type == GMSH_TRIANGLE_4 ||
         type == GMSH_TRIANGLE_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_face; ++jj )
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

    pv.marker().assign(  tag  );

    if ( type == GMSH_HEXAHEDRON ||
         type == GMSH_HEXAHEDRON_2 ||
         type == GMSH_TETRAHEDRON ||
         type == GMSH_TETRAHEDRON_2 ||
         type == GMSH_TETRAHEDRON_3 ||
         type == GMSH_TETRAHEDRON_4 ||
         type == GMSH_TETRAHEDRON_5 )
        {
            for( uint16_type jj = 0; jj < npoints_per_element; ++jj )
                pv.setPoint( jj, mesh->point( __e[jj] ) );
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


} // Life


#endif /* __ImporterGmsh_H */
