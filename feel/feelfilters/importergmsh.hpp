/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

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
#include <algorithm>

#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelfilters/importer.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <boost/algorithm/string/trim.hpp>


namespace Feel
{
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

    /**
     * setElementRegionAsPhysicalRegion(bool parameter)
     *
     * change reading for importing specific meshes for which the gmsh reader
     * consider the Physical flag as null
     */
    void setElementRegionAsPhysicalRegion( const bool param )
    {
        M_use_elementary_region_as_physical_region=param;
    }

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

    ImporterGmsh( WorldComm const& _worldcomm = Environment::worldComm() )
        :
        super( GMSH, _worldcomm ),
        _M_version( FEELPP_GMSH_FORMAT_VERSION ),
        M_use_elementary_region_as_physical_region( false )
    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
    }

    explicit ImporterGmsh( std::string const& _fname, std::string _version = FEELPP_GMSH_FORMAT_VERSION,
                           WorldComm const& _worldcomm = Environment::worldComm() )
        :
        super( _fname, GMSH, _worldcomm ),
        _M_version( _version ),
        M_use_elementary_region_as_physical_region( false )
    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
    }
    ImporterGmsh( ImporterGmsh const & i )
        :
        super( i ),
        _M_version( i._M_version ),
        M_use_elementary_region_as_physical_region( false ),
        _M_ignorePhysicalGroup( i._M_ignorePhysicalGroup ),
        _M_ignorePhysicalName( i._M_ignorePhysicalName )
    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
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
    std::string version() const
    {
        return _M_version;
    }

    /**
     * \return true if the element is on processor or is a ghost cell
     *         true if a ghost cell
     */
    boost::tuple<bool,boost::tuple<bool,int> > isElementOnProcessor( std::vector<int> const& tags ) const;

    //@}

    /** @name  Mutators
     */
    //@{

    void setVersion( std::string const& version )
    {
        _M_version = version;
    }

    void setIgnorePhysicalGroup( int i )
    {
        _M_ignorePhysicalGroup.insert( i );
    }
    void setIgnorePhysicalName( std::string s )
    {
        _M_ignorePhysicalName.insert( s );
    }

    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh );

    void showMe() const;

    //@}



protected:

private:

    void addPoint( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel );
    void addPoint( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<1> );
    void addPoint( mesh_type* /*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> /*tag*/, GMSH_ENTITY /*type*/, int & /*__idGmshToFeel*/, mpl::int_<2> );
    void addPoint( mesh_type* /*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> /*tag*/, GMSH_ENTITY /*type*/, int & /*__idGmshToFeel*/, mpl::int_<3> );

    void addEdge( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel );
    void addEdge( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<1> );
    void addEdge( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<2> );
    void addEdge( mesh_type* /*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> const& /*tag*/, GMSH_ENTITY /*type*/, int & /*__idGmshToFeel*/, mpl::int_<3> );

    void addFace( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel );
    void addFace( mesh_type* /*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> const& /*tag*/, GMSH_ENTITY /*type*/, int & /*__idGmshToFeel*/, mpl::int_<1> );
    void addFace( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<2> );
    void addFace( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<3> );

    void addVolume( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel );
    void addVolume( mesh_type* /*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> const& /*tag*/, GMSH_ENTITY,int & /*__idGmshToFeel*/ , mpl::int_<1> );
    void addVolume( mesh_type* /*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> const& /*tag*/, GMSH_ENTITY,int & /*__idGmshToFeel*/ , mpl::int_<2> );
    void addVolume( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & /*__idGmshToFeel*/, mpl::int_<3> );

    void updateGhostCellInfo( mesh_type* mesh, std::vector<int> const& __idGmshToFeel, std::map<int,boost::tuple<int,int> > const& __mapGhostElt );


private:

    std::string _M_version;
    std::vector<int> _M_n_vertices;
    std::vector<int> _M_n_b_vertices;

    std::set<int> _M_ignorePhysicalGroup;
    std::set<std::string> _M_ignorePhysicalName;
    bool M_use_elementary_region_as_physical_region;


};



template<typename MeshType>
void
ImporterGmsh<MeshType>::showMe() const
{
    Debug( 8011 ) << "[ImporterGmsh::showMe] npoints_per_element = " << npoints_per_element << "\n";
    Debug( 8011 ) << "[ImporterGmsh::showMe]    npoints_per_face = " << npoints_per_face << "\n";
    Debug( 8011 ) << "[ImporterGmsh::showMe]    npoints_per_edge = " << npoints_per_edge << "\n";
}
template<typename MeshType>
boost::tuple<bool,boost::tuple<bool,int> >
ImporterGmsh<MeshType>::isElementOnProcessor( std::vector<int> const& tag ) const
{
    bool is_element_on_processor = false;

    // if there is no partition tag (only 2 tags) then consider all elements on
    // current processor
    if ( tag.size() == 2 )
        return boost::make_tuple( true,false );

    // tag[2] is the number of partition tags (1 in general 2 if cell share an
    // face with 2 processors)
    int rank = this->worldComm().localRank();
    auto it = std::find_if( tag.begin()+3, tag.end(), [&rank] ( int i )
    {
        return i == rank;
    } );

    if ( it != tag.end() )
        is_element_on_processor = true;

    bool isGhostCell=false;
    int idOfGhostCell=0;

    if ( tag[3]!=rank )
    {
        isGhostCell=true;
        idOfGhostCell=tag[3];
    }

    return boost::make_tuple( is_element_on_processor, boost::make_tuple( isGhostCell,idOfGhostCell ) );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::visit( mesh_type* mesh )
{

    if ( this->version() != "1.0" &&
            this->version() != "2.0" &&
            this->version() != "2.1" &&
            this->version() != "2.2" &&
            this->version() != FEELPP_GMSH_FORMAT_VERSION )
        throw std::logic_error( "invalid gmsh file format version" );

    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit()] starts\n";

    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit( "  << mesh_type::nDim << "D )] starts\n";
    Debug( 8011 ) << "[ImporterGmsh<" << typeid( *mesh ).name() << ">::visit( "  << mesh_type::nDim << "D )] filename = " << this->filename() << "\n";

    std::ifstream __is ( this->filename().c_str() );

    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        ostr << "Invalid file name " << this->filename() << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    char __buf[256];
    __is >> __buf;

    std::string theversion;

    if ( ( ( this->version() == "2.0" ) ||
            ( this->version() == "2.1" ) ||
            ( this->version() == "2.2" ) ||
            ( this->version() == FEELPP_GMSH_FORMAT_VERSION ) )  &&
            std::string( __buf ) == "$MeshFormat" )
    {

        // version file-type(0=ASCII,1=BINARY) data-size(sizeof(double))
        __is >> theversion >> __buf >> __buf;
        FEELPP_ASSERT( boost::lexical_cast<double>( theversion ) >= 2 )( theversion )( this->version() ).warn( "invalid gmsh file format version " );
        // should be $EndMeshFormat
        __is >> __buf;
        FEELPP_ASSERT( std::string( __buf ) == "$EndMeshFormat" )
        ( __buf )
        ( "$EndMeshFormat" ).error ( "invalid file format entry" );
        __is >> __buf;
        Debug( 8011 ) << "[importergmsh] " << __buf << " (expect $PhysicalNames)\n";

        if ( std::string( __buf ) == "$PhysicalNames" )
        {
            int nnames;
            __is >> nnames;

            for ( int n = 0; n < nnames; ++n )
            {
                int id, topodim;
                std::string name;

                if ( boost::lexical_cast<double>( this->version() ) >= 2.1  )
                {
                    __is >> topodim >> id >> name;
                    Debug( 8011 ) << "[importergmsh] reading topodim: "  << topodim << " id: " << id << " name: " << name << "\n";
                }

                else if ( this->version() == "2.0" )
                    __is >> id >> name;

                boost::trim( name );
                boost::trim_if( name,boost::is_any_of( "\"" ) );

                std::vector<int> data = {id, topodim};
                mesh->addMarkerName( name, id, topodim );
                if ( _M_ignorePhysicalName.find( name )!=_M_ignorePhysicalName.end() ) this->setIgnorePhysicalGroup( id );
            }

            FEELPP_ASSERT( mesh->markerNames().size() == ( size_type )nnames )( mesh->markerNames().size() )( nnames ).error( "invalid number of physical names" );
            __is >> __buf;
            FEELPP_ASSERT( std::string( __buf ) == "$EndPhysicalNames" )
            ( __buf )
            ( "$EndPhysicalNames" ).error ( "invalid file format entry" );
            __is >> __buf;
        }
    }

    //
    // Read NODES
    //
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    FEELPP_ASSERT( std::string( __buf ) == "$NOD" ||
                   std::string( __buf ) == "$Nodes" ||
                   std::string( __buf ) == "$ParametricNodes" )
    ( __buf )
    ( "$NOD" )( "$Nodes" )( "$ParametricNodes" )
    .error( "invalid nodes string in gmsh importer" );
    bool has_parametric_nodes = ( std::string( __buf ) == "$ParametricNodes" );
    uint __n;
    __is >> __n;
    Debug( 8011 ) << "number of nodes: " << __n;

    std::vector<double> __x( 3*__n );
    std::vector<int> __gdim( __n );
    std::vector<int> __gtag( __n );
    std::vector<double> __uv( 2*__n );
    std::fill( __gdim.begin(), __gdim.end(), 0 );
    std::fill( __gtag.begin(), __gtag.end(), 0 );
    std::fill( __uv.begin(), __uv.end(), 0. );

    std::vector<bool> __isonboundary( __n );
    std::vector<std::vector<int> > __whichboundary( __n );
    Debug( 8011 ) << "reading "<< __n << " nodes\n";
    std::map<int,int> itoii;

    for ( uint __i = 0; __i < __n; ++__i )
    {
        uint __ni;
        __is >> __ni
             >> __x[3*__i]
             >> __x[3*__i+1]
             >> __x[3*__i+2];

        if ( has_parametric_nodes )
        {
            __is >> __gdim[__i] >> __gtag[__i];

            // if gdim == 0 then u = v = 0
            // if gdim == 3 then no need for a parametric point
            // this logic is done later when filling the mesh data structure
            if ( __gdim[__i] == 1 )
                __is >> __uv[2*__i];

            else if ( __gdim[__i] == 2 )
                __is >> __uv[2*__i] >> __uv[2*__i+1];
        }

        // stores mapping to be able to reorder the indices
        // so that they are contiguous
        itoii[__ni] = __i;
    }

    __is >> __buf;
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    // make sure that we have read all the points
    FEELPP_ASSERT( std::string( __buf ) == "$ENDNOD" ||
                   std::string( __buf ) == "$EndNodes" ||
                   std::string( __buf ) == "$EndParametricNodes"
                 )
    ( __buf )
    ( "$ENDNOD" )( "$EndNodes" )( "$EndParametricNodes" ).error( "invalid end nodes string in gmsh importer" );


    //
    // Read ELEMENTS
    //
    __is >> __buf;
    Debug( 8011 ) << "buf: "<< __buf << "\n";
    FEELPP_ASSERT( std::string( __buf ) == "$ELM" ||
                   std::string( __buf ) == "$Elements" )
    ( __buf )
    ( "$ELM" )( "$Elements" )
    .error( "invalid elements string in gmsh importer" );

    uint __nele;
    __is >> __nele;

    Debug( 8011 ) << "number of elements: " << __nele << "\n";
    std::vector<std::vector<int> > __e( __nele ); // nodes in each element
    std::vector<std::vector<int> > __et( __nele ); // tags in each element
    std::vector<int> __idGmshToFeel( __nele ); // id Gmsh to id Feel
    std::vector<GMSH_ENTITY> __etype( __nele );
    std::vector<int> __gt( 32 );
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
    nptable[20]=9;// node order 3 triangle (9 nodes)
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

    for ( uint __i = 0; __i < __nele; ++__i )
    {
        int __ne, __t, __physical_region, __np, __elementary_region = 1;
        std::vector<int> partitions;

        if ( this->version() == "1.0" )
        {
            __is >> __ne  // elm-number
                 >> __t // elm-type
                 >> __physical_region // reg-phys
                 >> __elementary_region // reg-elem
                 >> __np; // number-of-nodes
            FEELPP_ASSERT( __np == nptable[__t] )( __np )( __t )( nptable[__t] ).error( "invalid number of nodes" );
        }

        else if ( boost::lexical_cast<double>( this->version() ) >= 2  )
        {
            int __ntag;
            __is >> __ne  // elm-number
                 >> __t // elm-type
                 >> __ntag; // number-of-tags

            for ( int t = 0; t < __ntag; ++t )
            {
                int tag;
                // tag=1 physical region
                // tag=2 elementary region
                // tag=3 n partition tags
                // tag=4.. partition ids
                __is >> tag;

                if ( tag < 0 )
                    __et[__i].push_back( -tag );

                else
                    __et[__i].push_back( tag );
            }

            if ( ( theversion=="2.0" || theversion=="2.1" ) &&  __ntag==3 )
            {
                __et[__i].push_back( __et[__i][2] );
                __et[__i][2]=1;
                ++__ntag;
            }

            // shift partition id according to processor ids
            for ( int ii = 3; ii < __ntag; ++ ii )
                __et[__i][ii] -= 1;

            for ( int jj=0; jj < ( int )__et[__i].size(); ++jj )
                Debug( 8011 ) << __et[__i][jj] << " ";

            Debug( 8011 ) << " is on proc:" << isElementOnProcessor( __et[__i] ).template get<0>() << "\n";
            __np = nptable[__t];

            Debug( 8011 ) << "element type: " << __t << " nb pts: " << __np << "\n";
        }

        ++__gt[ __t ];
        __etype[__i] = GMSH_ENTITY( __t );

        if ( M_use_elementary_region_as_physical_region )
        {
            __et[__i][0]=__et[__i][1];
        }

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
    FEELPP_ASSERT( std::string( __buf ) == "$ENDELM" ||
                   std::string( __buf ) == "$EndElements" )
    ( __buf )
    ( "$ENDELM" )( "$EndElements" ).error( "invalid end elements string in gmsh importer" );

    // read physical names
    if ( boost::lexical_cast<double>( this->version() ) >= 2  )
    {

    }

    //
    // FILL Mesh Data Structure
    //
    Debug( 8011 ) << "number of edges= " << __gt[1] << "\n";

    __isonboundary.assign( __n, false );
    std::vector<int> vv( 2, 0 );
    __whichboundary.assign( __n, vv );

    for ( uint __i = 0; __i < __nele; ++__i )
    {
        if ( isElementOnProcessor( __et[__i] ).template get<0>() == false ||
                _M_ignorePhysicalGroup.find( __et[__i][0] ) != _M_ignorePhysicalGroup.end() )
            continue;

        switch ( __etype[__i] )
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
                for ( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
                {
                    __isonboundary[ __e[__i][jj] ] = true;
                    __whichboundary[__e[__i][jj] ] = __et[__i];
                }
            }

            break;

        case GMSH_QUADRANGLE:
        case GMSH_QUADRANGLE_2:
        case GMSH_TRIANGLE:
        case GMSH_TRIANGLE_2:
        case GMSH_TRIANGLE_3:
        case GMSH_TRIANGLE_4:
        case GMSH_TRIANGLE_5:
            if ( mesh_type::nDim == 3 )
            {
                for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
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

    std::map<int,boost::tuple<int,int> > mapGhostElt;

    // add the points to the mesh
    for ( uint __i = 0; __i < __n; ++__i )
    {
        node_type __n( mesh_type::nRealDim );

        for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
            __n[j] = __x[3*__i+j];

        point_type __pt( __i,__n, __isonboundary[ __i ] );
        __pt.setOnBoundary( __isonboundary[ __i ] );
        __pt.setTags( __whichboundary[__i] );

        if ( has_parametric_nodes )
        {
            __pt.setGDim( __gdim[__i] );
            __pt.setGTag( __gtag[__i] );

            if ( __gdim[__i] < 3 )
            {
                __pt.setParametricCoordinates( __uv[2*__i], __uv[2*__i+1] );
                mesh->setParametric( true );
            }
        }

        mesh->addPoint( __pt );
    }

    _M_n_vertices.resize( __n );
    _M_n_vertices.assign( __n, 0 );
    _M_n_b_vertices.resize( __n );
    _M_n_b_vertices.assign( __n, 0 );

    // add the element to the mesh
    for ( uint __i = 0; __i < __nele; ++__i )
    {
        if ( isElementOnProcessor( __et[__i] ).template get<0>() == false ||
                _M_ignorePhysicalGroup.find( __et[__i][0] ) != _M_ignorePhysicalGroup.end() )
            continue;

        switch ( __etype[__i] )
        {
            // Points
        case GMSH_POINT:
        {
            addPoint( mesh, __e[__i], __et[__i], __etype[__i], __idGmshToFeel[__i] );

            if ( isElementOnProcessor( __et[__i] ).template get<1>().template get<0>() == true ) // a ghost cell
            {
                auto idProc=isElementOnProcessor( __et[__i] ).template get<1>().template get<1>();
                mapGhostElt.insert( std::make_pair( __i,boost::make_tuple( __idGmshToFeel[__i], idProc ) ) );
            }

            break;
        }

        // Edges
        case GMSH_LINE:
        case GMSH_LINE_2:
        case GMSH_LINE_3:
        case GMSH_LINE_4:
        case GMSH_LINE_5:
        {
            addEdge( mesh, __e[__i], __et[__i], __etype[__i], __idGmshToFeel[__i] );

            if ( isElementOnProcessor( __et[__i] ).template get<1>().template get<0>() == true ) // a ghost cell
            {
                auto idProc=isElementOnProcessor( __et[__i] ).template get<1>().template get<1>();
                mapGhostElt.insert( std::make_pair( __i,boost::make_tuple( __idGmshToFeel[__i], idProc ) ) );
            }

            break;
        }

        // Faces
        case GMSH_TRIANGLE:
        case GMSH_TRIANGLE_2:
        case GMSH_TRIANGLE_3:
        case GMSH_TRIANGLE_4:
        case GMSH_TRIANGLE_5:
        case GMSH_QUADRANGLE:
        case GMSH_QUADRANGLE_2:
        {
            addFace( mesh, __e[__i], __et[__i], __etype[__i], __idGmshToFeel[__i] );

            if ( isElementOnProcessor( __et[__i] ).template get<1>().template get<0>() == true ) // a ghost cell
            {
                auto idProc=isElementOnProcessor( __et[__i] ).template get<1>().template get<1>();
                mapGhostElt.insert( std::make_pair( __i,boost::make_tuple( __idGmshToFeel[__i], idProc ) ) );
            }

            break;
        }

        // Volumes
        case GMSH_TETRAHEDRON:
        case GMSH_TETRAHEDRON_2:
        case GMSH_TETRAHEDRON_3:
        case GMSH_TETRAHEDRON_4:
        case GMSH_TETRAHEDRON_5:
        case GMSH_HEXAHEDRON:
        case GMSH_HEXAHEDRON_2:
        {
            addVolume( mesh, __e[__i], __et[__i], __etype[__i], __idGmshToFeel[__i] );

            if ( isElementOnProcessor( __et[__i] ).template get<1>().template get<0>() == true ) // a ghost cell
            {
                auto idProc=isElementOnProcessor( __et[__i] ).template get<1>().template get<1>();
                mapGhostElt.insert( std::make_pair( __i,boost::make_tuple( __idGmshToFeel[__i], idProc ) ) );
            }

            break;
        }

        default:
            break;
        }

    } // loop over geometric entities in gmsh file (can be elements or faces)

    if ( this->worldComm().localSize()>1 )
        updateGhostCellInfo( mesh, __idGmshToFeel,  mapGhostElt );

    mesh->setNumVertices( std::accumulate( _M_n_vertices.begin(), _M_n_vertices.end(), 0 ) );

    Debug( 8011 ) << "done with reading and creating mesh from gmsh file\n";
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel )
{
    addPoint( mesh, __e, tag, type, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<1> )
{
    face_type pf;
    //pf.setWorldComm(this->worldComm());
    pf.setProcessIdInPartition( this->worldComm().localRank() );
    pf.setId( mesh->numFaces() );
    pf.setTags(  tag  );
    pf.setPoint( 0, mesh->point( __e[0] ) );

    _M_n_vertices[ __e[0] ] = 1;

    _M_n_b_vertices[ __e[0] ] = 1;

    pf.setOnBoundary( true );
    face_iterator fit;
    bool inserted;
    boost::tie( fit, inserted ) = mesh->addFace( pf );
    __idGmshToFeel=pf.id();

    auto theface = mesh->faceIterator( pf.id() );
    mesh->faces().modify( theface, detail::update_id_in_partition_type( this->worldComm().localRank(), pf.id() ) );

    Debug( 8011 ) << "added point on boundary ("
                  << fit->isOnBoundary() << ") with id :" << fit->id() << " and marker " << pf.marker()
                  << " n1: " << mesh->point( __e[0] ).node() << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*/*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> /*tag*/, GMSH_ENTITY /*type*/, int & /*__idGmshToFeel*/, mpl::int_<2> )
{
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*/*mesh*/, std::vector<int> const& /*__e*/, std::vector<int> /*tag*/, GMSH_ENTITY /*type*/, int & /*__idGmshToFeel*/, mpl::int_<3> )
{
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel )
{
    addEdge( mesh, __e, tag, type, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type*mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<1> )
{
    element_type e;
    //e.setWorldComm(this->worldComm());
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setTags(  tag  );

    if ( type == GMSH_LINE ||
            type == GMSH_LINE_2 ||
            type == GMSH_LINE_3 ||
            type == GMSH_LINE_4 ||
            type == GMSH_LINE_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
            e.setPoint( jj, mesh->point( __e[jj] ) );
    }

    mesh->addElement( e );
    __idGmshToFeel=e.id();

    auto theelt = mesh->elementIterator( e.id(), e.partitionId() );
    mesh->elements().modify( theelt, detail::update_id_in_partition_type( this->worldComm().localRank(), e.id() ) );

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    Debug( 8011 ) << "added edge with id :" << e.id()
                  << " n1: " << mesh->point( __e[0] ).node()
                  << " n2: " << mesh->point( __e[1] ).node() << "\n";
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<2> )
{
    face_type pf;
    pf.setProcessIdInPartition( this->worldComm().localRank() );
    pf.setId( mesh->numFaces() );
    pf.setTags(  tag  );

    if ( type == GMSH_LINE ||
            type == GMSH_LINE_2 ||
            type == GMSH_LINE_3 ||
            type == GMSH_LINE_4 ||
            type == GMSH_LINE_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
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
    __idGmshToFeel=pf.id();

    auto theface = mesh->faceIterator( pf.id() );
    mesh->faces().modify( theface, detail::update_id_in_partition_type( this->worldComm().localRank(), pf.id() ) );

    Debug( 8011 ) << "added edge on boundary ("
                  << fit->isOnBoundary() << ") with id :" << fit->id()
                  << " n1: " << mesh->point( __e[0] ).node()
                  << " n2: " << mesh->point( __e[1] ).node() << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type*, std::vector<int> const&, std::vector<int> const&, GMSH_ENTITY, int & /*__idGmshToFeel*/, mpl::int_<3> )
{}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel )
{
    addFace( mesh, __e, tag, type, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type*, std::vector<int> const&, std::vector<int> const&, GMSH_ENTITY, int & /*__idGmshToFeel*/, mpl::int_<1> )
{}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<2> )
{
    GmshOrdering<element_type> ordering;

    element_type pf;
    pf.setProcessIdInPartition( this->worldComm().localRank() );
    pf.setTags(  tag  );

    if ( type == GMSH_QUADRANGLE ||
            type == GMSH_QUADRANGLE_2 ||
            type == GMSH_TRIANGLE ||
            type == GMSH_TRIANGLE_2 ||
            type == GMSH_TRIANGLE_3 ||
            type == GMSH_TRIANGLE_4 ||
            type == GMSH_TRIANGLE_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            //std::cout << "gmsh index " << jj << " -> " << ordering.fromGmshId(jj) << " -> " << mesh->point( __e[jj] ).id()+1 << " : " << mesh->point( __e[jj] ).node() << "\n";
            pf.setPoint( ordering.fromGmshId( jj ), mesh->point( __e[jj] ) );
        }
    }

    mesh->addElement( pf );
    __idGmshToFeel=pf.id();

    auto theelt = mesh->elementIterator( pf.id(), pf.partitionId() );
    mesh->elements().modify( theelt, detail::update_id_in_partition_type( this->worldComm().localRank(), pf.id() ) );

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    _M_n_vertices[ __e[2] ] = 1;

    if ( type == GMSH_QUADRANGLE ||
            type == GMSH_QUADRANGLE_2 )
        _M_n_vertices[ __e[3] ] = 1;
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<3> )
{
    GmshOrdering<face_type> ordering;

    face_type pf;
    pf.setProcessIdInPartition( this->worldComm().localRank() );
    pf.setId( mesh->numFaces() );
    pf.setTags(  tag  );

    if ( type == GMSH_QUADRANGLE ||
            type == GMSH_QUADRANGLE_2 ||
            type == GMSH_TRIANGLE ||
            type == GMSH_TRIANGLE_2 ||
            type == GMSH_TRIANGLE_3 ||
            type == GMSH_TRIANGLE_4 ||
            type == GMSH_TRIANGLE_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
            pf.setPoint( ordering.fromGmshId( jj ), mesh->point( __e[jj] ) );

        //pf.setPoint( jj, mesh->point( __e[jj] ) );
    }

    pf.setOnBoundary( true );
    bool inserted;
    face_iterator fit;
    boost::tie( fit, inserted ) = mesh->addFace( pf );

    __idGmshToFeel=pf.id();

    auto theface = mesh->faceIterator( pf.id() );
    mesh->faces().modify( theface, detail::update_id_in_partition_type( this->worldComm().localRank(), pf.id() ) );

    _M_n_vertices[ __e[0] ] = 1;
    _M_n_vertices[ __e[1] ] = 1;
    _M_n_vertices[ __e[2] ] = 1;

    if ( type == GMSH_QUADRANGLE ||
            type == GMSH_QUADRANGLE_2 )
        _M_n_vertices[ __e[3] ] = 1;
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel )
{
    addVolume( mesh, __e, tag, type, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type*, std::vector<int> const&, std::vector<int> const&, GMSH_ENTITY, int & /*__idGmshToFeel*/, mpl::int_<1> )
{}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type*, std::vector<int> const&, std::vector<int> const&, GMSH_ENTITY, int & /*__idGmshToFeel*/, mpl::int_<2> )
{}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, std::vector<int> const& __e, std::vector<int> const& tag, GMSH_ENTITY type, int & __idGmshToFeel, mpl::int_<3> )
{
    element_type pv;
    pv.setProcessIdInPartition( this->worldComm().localRank() );
    GmshOrdering<element_type> ordering;
    pv.setTags(  tag  );

    //
    // WARNING: not yet done for high order in 3D !!!
    //
    if ( type == GMSH_HEXAHEDRON ||
            type == GMSH_HEXAHEDRON_2 ||
            type == GMSH_TETRAHEDRON ||
            type == GMSH_TETRAHEDRON_2 ||
            type == GMSH_TETRAHEDRON_3 ||
            type == GMSH_TETRAHEDRON_4 ||
            type == GMSH_TETRAHEDRON_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            //std::cout << "gmsh index " << jj << " -> " << ordering.fromGmshId(jj) << " -> " << mesh->point( __e[jj] ).id()+1 << " : " << mesh->point( __e[jj] ).node() << "\n";
            pv.setPoint( ordering.fromGmshId( jj ), mesh->point( __e[jj] ) );
        }
    }

    mesh->addElement( pv );
    __idGmshToFeel=pv.id();

    auto theelt = mesh->elementIterator( pv.id(), pv.partitionId() );
    mesh->elements().modify( theelt, detail::update_id_in_partition_type( this->worldComm().localRank(), pv.id() ) );

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

template<typename MeshType>
void
ImporterGmsh<MeshType>::updateGhostCellInfo( mesh_type* mesh, std::vector<int> const& __idGmshToFeel, std::map<int,boost::tuple<int,int> > const& __mapGhostElt )
{
    // counter of msg sent for each process
    std::vector<int> nbMsgToSend( this->worldComm().localSize() );
    std::fill( nbMsgToSend.begin(),nbMsgToSend.end(),0 );

    // map usefull to get final result
    std::vector< std::map<int,int> > mapMsg( this->worldComm().localSize() );

    // iterate over ghost elt
    auto it_map = __mapGhostElt.begin();
    auto en_map = __mapGhostElt.end();

    for ( int cpt=0; it_map!=en_map; ++it_map,++cpt )
    {
        auto idGmsh = it_map->first;
        auto idProc = it_map->second.template get<1>();
#if 0
        std::cout << "[updateGhostCellInfo]----1---\n"
                  << "I am the proc " << this->worldComm().globalRank()
                  << " local proc " << this->worldComm().localRank()
                  << " I send to the proc " << idProc << " for idGmsh " << idGmsh+1
                  << " with tag "<< nbMsgToSend[idProc]
                  << " the G " << mesh->element( it_map->second.template get<0>(),idProc ).G()
                  << std::endl;
#endif
        // send
        this->worldComm().localComm().send( idProc , nbMsgToSend[idProc], idGmsh );
        // save tag of request
        mapMsg[idProc].insert( std::make_pair( nbMsgToSend[idProc],it_map->second.template get<0>() ) );
        // update nb send
        nbMsgToSend[idProc]++;
    }

    // counter of msg received for each process
    std::vector<int> nbMsgToRecv;
    mpi::all_to_all( this->worldComm().localComm(),
                     nbMsgToSend,
                     nbMsgToRecv );

    // get gmsh id asked and re-send the correspond id Feel
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0 ; cpt<nbMsgToRecv[proc] ; ++cpt )
        {
            int idGmsh;
            //reception idGmsh
            this->worldComm().localComm().recv( proc, cpt, idGmsh );
#if 0
            std::cout << "[updateGhostCellInfo]----2---\n"
                      << "I am the proc" << this->worldComm().localRank()
                      << " I receive of the proc " << proc << " for idGmsh " << idGmsh+1
                      << " with tag "<< cpt
                      << " the G " << mesh->element( __idGmshToFeel[idGmsh] ).G()
                      << " with idFeel Classic " << __idGmshToFeel[idGmsh]
                      << std::endl;
#endif
            //re-send idFeel
            this->worldComm().localComm().send( proc, cpt, __idGmshToFeel[ idGmsh ] );
        }
    }

    // get response to initial request and update Feel::Mesh data
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        {
            int idFeel;
            // receive idFeel
            this->worldComm().localComm().recv( proc, cpt, idFeel );
            // update data
            auto elttt = mesh->elementIterator( mapMsg[proc][cpt],proc );
            mesh->elements().modify( elttt, detail::update_id_in_partition_type( proc, idFeel ) );
#if 0
            std::cout << "[updateGhostCellInfo]----3---\n"
                      << "END! I am the proc" << this->worldComm().localRank()<<" I receive of the proc " << proc
                      <<" with tag "<< cpt
                      << " for the G " << mesh->element( mapMsg[proc][cpt], proc ).G()
                      << " with idFeel Classic " << mapMsg[proc][cpt]
                      << " with idFeel " << idFeel
                      << " and modif " << mesh->element( mapMsg[proc][cpt] , proc ).idInPartition( proc )
                      << std::endl;
#endif
        }
    }


    //check
#if 0
    std::vector<int> nbMsgToSendCheck( this->worldComm().localSize() );
    std::fill( nbMsgToSendCheck.begin(),nbMsgToSendCheck.end(),0 );

    std::vector< std::map<int,int> > mapMsgCheck( this->worldComm().localSize() );
    auto it_ghost=mesh->beginGhostElement();
    auto en_ghost=mesh->endGhostElement();

    for ( int cpt=0 ; it_ghost != en_ghost; ++it_ghost,++cpt )
    {

        //ATTENTION WARNING int OK, but auto no!!!!!!!!!!!!!!!!!!!!!!!
        int/*auto*/ IdProcessOfGhost = it_ghost->processId();
        int/*auto*/ RealIdOfGhost = it_ghost->idInPartition( IdProcessOfGhost );
        this->worldComm().localComm().send( IdProcessOfGhost, nbMsgToSendCheck[IdProcessOfGhost], RealIdOfGhost );
#if 0
        std::cout<< "I am the proc" << this->worldComm().localRank() << " I send to the proc " << IdProcessOfGhost
                 <<" with tag "<< nbMsgToSendCheck[IdProcessOfGhost]
                 << " idSend " << RealIdOfGhost
                 << " it_ghost->G() " << it_ghost->G()
                 << std::endl;
#endif
        mapMsgCheck[IdProcessOfGhost].insert( std::make_pair( cpt,nbMsgToSendCheck[IdProcessOfGhost] ) );

        ++nbMsgToSendCheck[IdProcessOfGhost];
    }

    // counter of msg received for each process
    std::vector<int> nbMsgToRecvCheck;
    mpi::all_to_all( this->worldComm().localComm(),
                     nbMsgToSendCheck,
                     nbMsgToRecvCheck );

    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToRecvCheck[proc]; ++cpt )
        {
            int idRecv;
            this->worldComm().localComm().recv( proc, cpt,idRecv );
#if 0
            std::cout<< "I am the proc" << this->worldComm().localRank()<<" I receive to proc " << proc
                     <<" with tag "<< cpt
                     << " idRecv " << idRecv
                     << " it_ghost->G() " << mesh->element( idRecv ).G()
                     << std::endl;
#endif
            this->worldComm().localComm().send( proc, cpt, mesh->element( idRecv ).G() );
        }
    }

    typedef typename matrix_node<double>::type matrix_node_type;
    std::vector< std::vector<matrix_node_type> > getFinalInfoCheck( this->worldComm().localSize() );

    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        getFinalInfoCheck[proc].resize( nbMsgToSendCheck[proc] );

        for ( int cpt=0; cpt<nbMsgToSendCheck[proc]; ++cpt )
        {
            this->worldComm().localComm().recv( proc, cpt, getFinalInfoCheck[proc][cpt] );
#if 0
            std::cout<< "I am the proc " << this->worldComm().localRank()<<" I receive to proc " << proc
                     <<" with tag "<< cpt
                     << " points G " << getFinalInfoCheck[proc][cpt]
                     << std::endl;
#endif
        }
    }

    // check final data
    it_ghost=mesh->beginGhostElement();
    en_ghost=mesh->endGhostElement();

    for ( int cpt=0; it_ghost!=en_ghost; ++it_ghost,++cpt )
    {
        int/*auto*/ IdProcessOfGhost = it_ghost->processId();
        auto indic = mapMsgCheck[IdProcessOfGhost][cpt];
#if 0
        std::cout<< "FINAL CHECK! I am the proc" << this->worldComm().localRank()
                 << " G() initial " << it_ghost->G()
                 << " G() returned" << getFinalInfoCheck[IdProcessOfGhost/*idProc*/][indic]
                 << std::endl;
#endif
        // check test
        auto eltG=it_ghost->G();
        auto elt2G = getFinalInfoCheck[IdProcessOfGhost/*idProc*/][indic];

        for ( int i =0 ; i<eltG.size2(); i++ )
        {
            for ( int j =0 ; j<eltG.size1(); j++ )
                if ( eltG( i,j )!=elt2G( i,j ) ) std::cout << "\n BAD!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
        }
    }

#endif

}

} // Feel



#endif /* __ImporterGmsh_H */
