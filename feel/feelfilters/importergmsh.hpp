/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-16

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2009-2014 Feel++ Consortium

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-16
 */

#ifndef __ImporterGmsh_H
#define __ImporterGmsh_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelcore/feelgmsh.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/importer.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feeltiming/tic.hpp>
#include <boost/algorithm/string/trim.hpp>

#if defined( FEELPP_HAS_GMSH_H )
#include <GModel.h>
#include <MElement.h>
#include <MVertex.h>
#include <MEdge.h>
#include <MPoint.h>
#include <MLine.h>
#include <MTriangle.h>
#include <MQuadrangle.h>
#include <MTetrahedron.h>
#include <MHexahedron.h>
#endif
// there is a macro called sign in Gmsh that conflicts with
// at least one member function sign() from DofTable.
// hence we undefine the macro sign after including Gmsh headers
#undef sign



// from Gmsh
void SwapBytes(char *array, int size, int n);

namespace Feel
{
namespace detail
{
class GMSHPoint
{
public:
    Eigen::Vector3d x;
    Eigen::Vector2d uv;
    int id;
    bool onbdy;
    bool parametric;
    int gdim,gtag;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    GMSHPoint()
        :
        x( Eigen::Vector3d::Zero() ),
        uv( Eigen::Vector2d::Zero() ),
        id( 0 ),
        onbdy( false ),
        parametric( false ),
        gdim( 0 ),
        gtag( 0 )
        {}
    GMSHPoint(GMSHPoint const& p )
        :
        x( p.x ),
        uv( p.uv ),
        id ( p.id ),
        onbdy( p.onbdy ),
        parametric( p.parametric ),
        gdim( p.gdim ),
        gtag( p.gtag )
        {}
};
struct GMSHElement
{
    GMSHElement()
        :
        num( 0 ),
        type( MSH_PNT ),
        physical( 0 ),
        elementary( 0 ),
        numPartitions( 1 ),
        partition( 0 ),
        ghosts(),
        is_on_processor( false ),
        is_ghost( false ),
        ghost_partition_id( -1 ),
        parent( 0 ),
        dom1(0), dom2(0),
        numVertices(0),
        indices()
        {}

    GMSHElement( int n,
                 int t,
                 int p,
                 int e,
                 rank_type _numPartitions,
                 rank_type _partition,
                 std::vector<rank_type> const& _ghosts,
                 int _parent,
                 int _dom1, int _dom2,
                 int _numVertices,
                 std::vector<int> const& _indices,
                 rank_type worldcommrank,
                 rank_type worldcommsize,
                 bool use_partitioning = false )
    :
        num( n ),
        type( t ),
        physical( p ),
        elementary( e ),
        numPartitions( _numPartitions ),
        partition( use_partitioning?_partition:(_partition % worldcommsize) ),
        ghosts( _ghosts ),
        is_on_processor( false ),
        is_ghost( false ),
        ghost_partition_id( partition ),
        parent( _parent ),
        dom1(_dom1), dom2(_dom2),
        numVertices( _numVertices ),
        indices( _indices )
        {
            setPartition(worldcommrank,worldcommsize);
        }
    template<typename T>
    GMSHElement( T* ele, int n, int elementary, int physical )
        :
        GMSHElement(n,
                    ele->getTypeForMSH(),
                    physical,
                    elementary,
                    1, // numpartition
                    ele->getPartition(),
                    std::vector<rank_type>(),
                    0, 0, 0, // parent, dom1, dom2
                    ele->getNumVerticesForMSH(),
                    std::vector<int>( ele->getNumVerticesForMSH(), 0 ),
                    0, 1
                    )
        {
#if GMSH_VERSION_GREATER_OR_EQUAL_THAN(2,6,1)
          std::vector<int> verts;
          ele->getVerticesIdForMSH(verts);
          indices = verts;
#else
          int* verts = ele->getVerticesIdForMSH();
          std::copy( verts, verts+indices.size(), indices.begin() );
#endif
#if 0
          std::cout << "Adding element index " << n << " type " << type << " phys " << physical
                    << " elementary " << elementary << " np " << numPartitions << " part " << partition << " nverts " << verts.size() << std::endl;
          std::cout << "   - ";
          std::for_each( verts.begin(), verts.end(), []( int v ) { std::cout << v << " "; } );
          std::cout << "\n";
#endif
        }
    void setPartition(rank_type worldcommrank, rank_type worldcommsize)
        {
            // maybe proc id not start to 0
            for ( auto _itghost=ghosts.begin(),_enghost=ghosts.end() ; _itghost!=_enghost ; ++_itghost )
                *_itghost = ( (*_itghost) % worldcommsize);


            if ( worldcommsize == 1 )
            {
                is_on_processor = true;
                is_ghost = false;
            }
            else if ( worldcommrank == partition )
            {
                is_on_processor = true;
                is_ghost = false;
            }
            else
            {
                // is the element a ghost cell
                // look into ghosts if 'partition' is present
                auto it = std::find( ghosts.begin(), ghosts.end(), worldcommrank );
                if ( it != ghosts.end() )
                {
                    is_on_processor = true;
                    is_ghost = true;
                    ghost_partition_id = partition;
                }
            }
        }

   GMSHElement( GMSHElement const& g )
    :
        num( g.num ),
        type( g.type ),
        physical( g.physical ),
        elementary( g.elementary ),
        numPartitions( g.numPartitions ),
        partition( g.partition ),
        ghosts( g.ghosts ),
        is_on_processor( g.is_on_processor ),
        is_ghost( g.is_ghost ),
        ghost_partition_id( g.ghost_partition_id ),
        parent( g.parent ),
        dom1( g.dom1 ), dom2( g.dom2 ),
        numVertices( g.numVertices ),
        indices( g.indices )
        {}

    bool isOnProcessor() const { return is_on_processor; }
    bool isGhost() const { return is_ghost; }
    rank_type ghostPartitionId() const { return ghost_partition_id; }

    template<typename IteratorBegin,typename IteratorEnd>
    bool isIgnored(IteratorBegin it, IteratorEnd en ) const { return std::find( it, en, physical ) != en; }

    void updatePartition( std::map<rank_type,rank_type> const& p2e, rank_type worldcommrank, rank_type worldcommsize )
        {
            partition = num;
            for( auto& g : ghosts )
                g = p2e.at(g)-1;
            setPartition( worldcommrank, worldcommsize );
        }
    int num;
    int type;
    int physical;
    int elementary;

    //! partitioning info
    rank_type numPartitions;
    rank_type partition;
    std::vector<rank_type> ghosts;
    bool is_on_processor;
    bool is_ghost;
    rank_type ghost_partition_id;

    int parent;
    int dom1, dom2;

    // vertices
    int numVertices;
    std::vector<int> indices;


};
}
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
        M_version( FEELPP_GMSH_FORMAT_VERSION ),
        M_in_memory( false ),
        M_use_elementary_region_as_physical_region( false ),
        M_respect_partition( false )
    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
    }

    explicit ImporterGmsh( std::string const& _fname, std::string _version = FEELPP_GMSH_FORMAT_VERSION,
                           WorldComm const& _worldcomm = Environment::worldComm() )
        :
        super( _fname, GMSH, _worldcomm ),
        M_version( _version ),
        M_in_memory( false ),
        M_use_elementary_region_as_physical_region( false ),
        M_respect_partition( false )
    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
    }
    ImporterGmsh( ImporterGmsh const & i )
        :
        super( i ),
        M_version( i.M_version ),
        M_in_memory( i.M_in_memory ),
        M_use_elementary_region_as_physical_region( false ),
        M_ignorePhysicalGroup( i.M_ignorePhysicalGroup ),
        M_ignorePhysicalName( i.M_ignorePhysicalName ),
        M_respect_partition( i.M_respect_partition )
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
        return M_version;
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
        M_version = version;
    }
    void setInMemory( bool in )
    {
        M_in_memory = in;
    }
    void setIgnorePhysicalGroup( int i )
    {
        M_ignorePhysicalGroup.insert( i );
    }
    void setIgnorePhysicalName( std::string s )
    {
        M_ignorePhysicalName.insert( s );
    }

    void setRespectPartition( bool r )
        {
            M_respect_partition = r;
        }

#if defined( FEELPP_HAS_GMSH_H)
    /**
     * @brief set the GMsh GModel object
     * @param gmodel GMsh GModel object
     */
    void setGModel( GModel* gmodel ) { M_gmodel = gmodel; }
#endif

    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh );

    void showMe() const;

    //@}



protected:

private:
    void readFromMemory( mesh_type* mesh );
    void readFromFile( mesh_type* mesh );

    void addVertices( mesh_type* mesh, Feel::detail::GMSHElement const& elt,
                      std::map<int, Feel::detail::GMSHPoint > const& gmshpts );

    void addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel );
    void addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<1> );
    void addPoint( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, int & /*__idGmshToFeel*/, mpl::int_<2> );
    void addPoint( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, int & /*__idGmshToFeel*/, mpl::int_<3> );

    void addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel );
    void addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<1> );
    void addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<2> );
    void addEdge( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, int & /*__idGmshToFeel*/, mpl::int_<3> );

    void addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel );
    void addFace( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, int & /*__idGmshToFeel*/, mpl::int_<1> );
    void addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<2> );
    void addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<3> );

    void addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel );
    void addVolume( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, int & /*__idGmshToFeel*/ , mpl::int_<1> );
    void addVolume( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, int & /*__idGmshToFeel*/ , mpl::int_<2> );
    void addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & /*__idGmshToFeel*/, mpl::int_<3> );

    void updateGhostCellInfoByUsingBlockingComm( mesh_type* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
                                                 std::vector<int> const& nbMsgToRecv );

    void updateGhostCellInfoByUsingNonBlockingComm( mesh_type* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
                                                    std::vector<int> const& nbMsgToRecv );


private:

    std::string M_version;
    bool M_in_memory;
    std::map<int,int> M_n_vertices;
    //std::vector<int> M_n_b_vertices;

    std::set<int> M_ignorePhysicalGroup;
    std::set<std::string> M_ignorePhysicalName;
    bool M_use_elementary_region_as_physical_region;
    bool M_respect_partition;
    //std::map<int,int> itoii;
    //std::vector<int> ptseen;
    GModel* M_gmodel;

};



template<typename MeshType>
void
ImporterGmsh<MeshType>::showMe() const
{
    DVLOG(2) << "[ImporterGmsh::showMe] npoints_per_element = " << npoints_per_element << "\n";
    DVLOG(2) << "[ImporterGmsh::showMe]    npoints_per_face = " << npoints_per_face << "\n";
    DVLOG(2) << "[ImporterGmsh::showMe]    npoints_per_edge = " << npoints_per_edge << "\n";
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
ImporterGmsh<MeshType>::addVertices( mesh_type* mesh, Feel::detail::GMSHElement const& elt,
                                     std::map<int, Feel::detail::GMSHPoint > const& gmshpts )
{
    node_type coords( mesh_type::nRealDim );
    // add the points associates to the element on the processor
    for ( uint16_type p = 0; p < elt.numVertices; ++p )
    {
        int ptid = elt.indices[p];
        // don't do anything if the point is already registered
        if ( mesh->hasPoint( ptid ) )
            continue;

        auto const& gmshpt = gmshpts.find(ptid)->second;
        for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
            coords[j] = gmshpt.x[j];

        point_type pt( ptid, coords, gmshpt.onbdy );
        pt.setProcessIdInPartition( this->worldComm().localRank() );
        if ( gmshpt.parametric )
        {
            pt.setGDim( gmshpt.gdim );
            pt.setGTag( gmshpt.gtag );

            if ( gmshpt.gdim < 3 )
            {
                pt.setParametricCoordinates( gmshpt.uv[0], gmshpt.uv[1] );
                mesh->setParametric( true );
            }
        }
        mesh->addPoint( pt );
    } // loop over local points
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::visit( mesh_type* mesh )
{
    DVLOG(2) << "visit("  << mesh_type::nDim << "D ) starts\n";
    if ( ( M_gmodel != 0 ) && ( M_in_memory == true ) )
        readFromMemory( mesh );
    else
        readFromFile( mesh );

    mesh->setNumVertices( std::accumulate( M_n_vertices.begin(), M_n_vertices.end(), 0,
                                           []( int lhs, std::pair<int,int> const& rhs )
                                           {
                                               return lhs+rhs.second;
                                           } ) );
    M_n_vertices.clear();
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::readFromMemory( mesh_type* mesh )
{
    tic();
    LOG(INFO) << "Reading Msh from memory ";
    // get the number of vertices and index the vertices in a continuous
    // sequence
    bool saveAll=true;
    bool saveSinglePartition=false;
    int numVertices = M_gmodel->indexMeshVertices(saveAll, saveSinglePartition);

    if( M_gmodel->numPhysicalNames() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "  +- number of physicals: " << M_gmodel->numPhysicalNames() << "\n";
        for(GModel::piter it = M_gmodel->firstPhysicalName(); it != M_gmodel->lastPhysicalName(); it++)
        {
            int id = it->first.second;
            int topodim = it->first.first;
            std::vector<int> data = {id, topodim};
            mesh->addMarkerName( it->second.c_str(), id, topodim );
        }
    }

    // get the number of elements we need to save
    int numElements = M_gmodel->getNumMeshElements();
    if ( Environment::isMasterRank() )
        std::cout << "  +- number of vertices: " << numVertices << "\n"
                  << "  +- number of elements: " << numElements << "\n";

    std::map<int, Feel::detail::GMSHPoint > gmshpts;
    std::vector<GEntity*> entities;
    M_gmodel->getEntities(entities);
    for(unsigned int i = 0; i < entities.size(); i++)
        for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
        {
            auto const& pt = entities[i]->mesh_vertices[j];
            int id = pt->getIndex();
            Eigen::Vector3d x(pt->x(), pt->y(), pt->z());
            gmshpts[id].id = id;
            gmshpts[id].x = x;
            gmshpts[id].parametric = false;
        }
    std::map<int,int> __idGmshToFeel; // id Gmsh to id Feel
    int num = 0;
    // read vertices
    for( GModel::viter it = M_gmodel->firstVertex(); it != M_gmodel->lastVertex(); ++it )
    {
        int elementary = (*it)->tag();
        auto const& physicals = (*it)->physicals;
        for( auto const& p : (*it)->points )
        {
            if ( physicals.size() )
            {
                Feel::detail::GMSHElement ele(p, ++num, elementary, physicals.size()?physicals[0]:0);
                this->addVertices( mesh, ele, gmshpts );
                addPoint( mesh, ele, __idGmshToFeel[p->getNum()] );
            }
        }
    }
    // read edges
    for( GModel::eiter it = M_gmodel->firstEdge(); it != M_gmodel->lastEdge(); ++it )
    {
        int elementary = (*it)->tag();
        auto const& physicals = (*it)->physicals;
        for( auto const& edge : (*it)->lines )
        {
            Feel::detail::GMSHElement ele(edge, ++num, elementary, physicals.size()?physicals[0]:0);
            this->addVertices( mesh, ele, gmshpts );
            addEdge( mesh, ele, __idGmshToFeel[edge->getNum()] );
        }
    }
    // read triangles
    for( GModel::fiter it = M_gmodel->firstFace(); it != M_gmodel->lastFace(); ++it)
    {
        int elementary = (*it)->tag();
        auto const& physicals = (*it)->physicals;
        for( auto const& triangle : (*it)->triangles )
        {
            Feel::detail::GMSHElement ele(triangle, ++num, elementary, physicals.size()?physicals[0]:0);
            this->addVertices( mesh, ele, gmshpts );
            addFace( mesh, ele, __idGmshToFeel[triangle->getNum()] );
        }
    }
    // read quadrangles
    for( GModel::fiter it = M_gmodel->firstFace(); it != M_gmodel->lastFace(); ++it)
    {
        int elementary = (*it)->tag();
        auto const& physicals = (*it)->physicals;
        for( auto const& quad : (*it)->quadrangles )
        {
            Feel::detail::GMSHElement ele(quad, ++ num, elementary, physicals.size()?physicals[0]:0);
            this->addVertices( mesh, ele, gmshpts );
            addFace( mesh, ele, __idGmshToFeel[quad->getNum()] );
        }
    }
    // read tetras
    for( GModel::riter it = M_gmodel->firstRegion(); it != M_gmodel->lastRegion(); ++it)
    {
        int elementary = (*it)->tag();
        auto const& physicals = (*it)->physicals;
        for( auto const& tetra : (*it)->tetrahedra )
        {
            Feel::detail::GMSHElement ele(tetra, ++num, elementary, physicals.size()?physicals[0]:0);
            this->addVertices( mesh, ele, gmshpts );
            addVolume( mesh, ele, __idGmshToFeel[tetra->getNum()] );   }
    }
    // read hexa
    for( GModel::riter it = M_gmodel->firstRegion(); it != M_gmodel->lastRegion(); ++it)
    {
        int elementary = (*it)->tag();
        auto const& physicals = (*it)->physicals;
        for( auto const& hexa : (*it)->hexahedra )
        {
            Feel::detail::GMSHElement ele(hexa, ++num, elementary, physicals.size()?physicals[0]:0);
            this->addVertices( mesh, ele, gmshpts );
            addVolume( mesh, ele, __idGmshToFeel[hexa->getNum()] );
        }
    }
    LOG(INFO) << "Reading Msh from memory done";
    toc("read msh from memory");
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::readFromFile( mesh_type* mesh )
{
    tic();
    LOG(INFO) << "Reading Msh file " << this->filename();
    if ( this->version() != "1.0" &&
         this->version() != "2.0" &&
         this->version() != "2.1" &&
         this->version() != "2.2" &&
         this->version() != FEELPP_GMSH_FORMAT_VERSION )
        throw std::logic_error( "invalid gmsh file format version" );


    std::ifstream __is ( this->filename().c_str() );

    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        LOG(ERROR) << "Invalid file name " << this->filename() << " (file not found)";
        ostr << "Invalid file name " << this->filename() << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    char __buf[256];
    __is >> __buf;

    std::string theversion;
    double version = 2.2;

    bool binary = false, swap = false;

    if ( ( ( this->version() == "2.0" ) ||
           ( this->version() == "2.1" ) ||
           ( this->version() == "2.2" ) ||
           ( this->version() == FEELPP_GMSH_FORMAT_VERSION ) )  &&
         std::string( __buf ) == "$MeshFormat" )
    {

        // version file-type(0=ASCII,1=BINARY) data-size(sizeof(double))
        int format, size;
        __is >> theversion >> format >> size;
        LOG(INFO) << "GMSH mesh file version : " << theversion << " format: " << (format?"binary":"ascii") << " size of double: " << size << "\n";
        CHECK( boost::lexical_cast<double>( theversion ) >= 2  )
            << "Feel++ supports only Gmsh version >= 2.\n";
        version =  boost::lexical_cast<double>( theversion );
        if(format)
        {
            char c;
            CHECK( (c=__is.get()) == '\n' ) << "Invalid character " << c << " should be newline\n";
            binary = true;
            LOG(INFO) << "GMSH mesh is in binary format\n";
            int one;
            __is.read( (char*)&one, sizeof(int) );
            if(one != 1)
            {
                swap = true;
                LOG(INFO) << "one before swap : " << one << "\n";
                if(swap) SwapBytes((char*)&one, sizeof(int), 1);
                LOG(INFO) << "one after swap : " << one << "\n";
                LOG(INFO) <<"Swapping bytes from binary file (to be done)\n";
            }
            //__is >> c;
            //LOG(INFO) << "character= " << c << "\n";
        }
        // should be $EndMeshFormat
        __is >> __buf;
        CHECK( std::string( __buf ) == "$EndMeshFormat" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndMeshFormat\n";
        __is >> __buf;
        DVLOG(2) << "[importergmsh] " << __buf << " (expect $PhysicalNames)\n";

        std::vector<MeshMarkerName> meshMarkerNameMap = markerMap(MeshType::nDim);
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
                    DVLOG(2) << "[importergmsh] reading topodim: "  << topodim << " id: " << id << " name: " << name << "\n";
                }

                else if ( this->version() == "2.0" )
                    __is >> id >> name;

                boost::trim( name );
                boost::trim_if( name,boost::is_any_of( "\"" ) );

                if ( meshMarkerNameMap.empty() )
                {
                    std::vector<int> data = {id, topodim};
                    mesh->addMarkerName( name, id, topodim );
                }
                if ( M_ignorePhysicalName.find( name )!=M_ignorePhysicalName.end() ) this->setIgnorePhysicalGroup( id );
            }
            if ( meshMarkerNameMap.empty() )
            {
                FEELPP_ASSERT( mesh->markerNames().size() == ( size_type )nnames )( mesh->markerNames().size() )( nnames ).error( "invalid number of physical names" );
            }
            __is >> __buf;
            FEELPP_ASSERT( std::string( __buf ) == "$EndPhysicalNames" )
                ( __buf )
                ( "$EndPhysicalNames" ).error ( "invalid file format entry" );
            __is >> __buf;
        }

        for( auto const& it : meshMarkerNameMap )
        {
            mesh->addMarkerName( it.name, it.ids[0], it.ids[1] );
        }
    }

    //
    // Read NODES
    //
    DVLOG(2) << "buf: "<< __buf << "\n";
    CHECK( std::string( __buf ) == "$NOD" ||
           std::string( __buf ) == "$Nodes" ||
           std::string( __buf ) == "$ParametricNodes" )
        << "invalid nodes string '" << __buf << "' in gmsh importer. It should be either $NOD or $Nodes or $ParametricNodes.\n";
    bool has_parametric_nodes = ( std::string( __buf ) == "$ParametricNodes" );
    uint __n;
    __is >> __n;

    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    if ( binary )
        __is.get();


    std::map<int, Feel::detail::GMSHPoint > gmshpts;
    LOG(INFO) << "Reading "<< __n << " nodes\n";

    Eigen::Vector3d x;
    Eigen::Vector2d uv;
    for ( uint __i = 0; __i < __n; ++__i )
    {
        int id = 0;
        if ( !binary )
        {
            __is >> id
                 >> x[0]
                 >> x[1]
                 >> x[2];
        }
        else
        {
            __is.read( (char*)&id, sizeof(int) );
            if(swap) SwapBytes((char*)&id, sizeof(int), 1);
            __is.read( (char*)&x[0], 3*sizeof(double) );
            if(swap) SwapBytes((char*)&x[0], sizeof(double), 3);

        }
        gmshpts[id].id = id;
        gmshpts[id].x = x;

        if ( has_parametric_nodes )
        {

            gmshpts[id].parametric = true;
            CHECK( !binary ) << "GMSH Binary format not yet supported for parametric nodes\n";
            int gdim = 0, gtag = 0;
            if ( !binary )
            {
                __is >> gdim >> gtag;

                // if gdim == 0 then u = v = 0
                // if gdim == 3 then no need for a parametric point
                // this logic is done later when filling the mesh data structure
                if ( gdim == 1 )
                    __is >> uv[0];

                else if ( gdim == 2 )
                    __is >> uv[0] >> uv[1];
            }
            gmshpts[id].gdim = gdim;
            gmshpts[id].gtag = gtag;
            gmshpts[id].uv = uv;
        }

        // stores mapping to be able to reorder the indices
        // so that they are contiguous
        //itoii[idpts[__i]] = __i;
    }
    //ptseen.resize( __n );
    //std::fill( ptseen.begin(), ptseen.end(), -1 );
    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    if ( binary )
        __is.get();

    __is >> __buf;
    DVLOG(2) << "buf: "<< __buf << "\n";
    // make sure that we have read all the points
    CHECK( std::string( __buf ) == "$ENDNOD" ||
           std::string( __buf ) == "$EndNodes" ||
           std::string( __buf ) == "$EndParametricNodes" )
        << "invalid end nodes string '" << __buf
        << "' in gmsh importer. It should be either $ENDNOD or $EndNodes or $EndParametricNodes\n";

    //
    // Read ELEMENTS
    //
    __is >> __buf;
    CHECK( std::string( __buf ) == "$ELM" ||
           std::string( __buf ) == "$Elements" )
        << "invalid elements string " << __buf << " in gmsh importer, it should be either $ELM or $Elements\n";

    int numElements;
    __is >> numElements;

    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    if ( binary )
        __is.get();

    LOG(INFO) << "Reading " << numElements << " elements...\n";
    std::list<Feel::detail::GMSHElement> __et; // tags in each element
    std::map<int,int> __idGmshToFeel; // id Gmsh to id Feel
    std::map<int,int> __gt;

    if ( !binary )
    {
        for(int i = 0; i < numElements; i++)
        {
          int num, type, physical = 0, elementary = 0, parent = 0;
          int dom1 = 0, dom2 = 0, numVertices;
          std::vector<rank_type> ghosts;
          int numTags;
          // some faces may not be associated to a partition in the mesh file,
          // hence will be read given the partition id 0 and will be discarded
          rank_type partition = (this->worldComm().globalSize()>1)?this->worldComm().localRank():0;
          __is >> num  // elm-number
               >> type // elm-type
               >> numTags; // number-of-tags

          rank_type numPartitions = 1;

          for(int j = 0; j < numTags; j++)
          {
              int tag;
              __is >> tag;
              if(j == 0) physical = tag;
              else if(j == 1) elementary = tag;
              else if(version < 2.2 && j == 2) partition = tag;
              else if(version >= 2.2 && j == 2 && numTags > 3) numPartitions = tag;
              else if(version >= 2.2 && j == 3) partition = tag-1;
              else if(j >= 4 && j < 4 + numPartitions - 1) ghosts.push_back((-tag)-1);
              else if(j == 3 + numPartitions && (numTags == 4 + numPartitions))
                  parent = tag;
              else if(j == 3 + numPartitions && (numTags == 5 + numPartitions)) {
                  dom1 = tag; j++;
                  __is >> dom2;
              }
          }

          CHECK(type != MSH_POLYG_ && type != MSH_POLYH_ && type != MSH_POLYG_B)
              << "GMSH Element type " << type << " not supported by Feel++\n";
          numVertices = MElement::getInfoMSH(type);
          CHECK(numVertices!=0) << "Unknown number of vertices for element type " << type << "\n";

          std::vector<int> indices(numVertices);
          for(int j = 0; j < numVertices; j++)
          {
              __is >> indices[j];
              // we do not renumber anymore
#if 0
              indices[j] = itoii[ indices[j] ];
#endif
          }
          if ( M_use_elementary_region_as_physical_region )
          {
              physical = elementary;
          }

          Feel::detail::GMSHElement gmshElt( num, type, physical, elementary,
                                             numPartitions, partition, ghosts,
                                             parent, dom1, dom2,
                                             numVertices, indices,
                                             this->worldComm().localRank(),
                                             this->worldComm().localSize(),
                                             M_respect_partition );

          // WARNING: we had another condition if the number of processors and
          // elements is the same, in that case we store everything for now
          if ( (( gmshElt.isOnProcessor() == false) ||
                gmshElt.isIgnored(M_ignorePhysicalGroup.begin(), M_ignorePhysicalGroup.end()) ) )
              continue;
          __et.push_back( gmshElt );

          if ( __gt.find( type ) != __gt.end() )
              ++__gt[ type ];
          else
              __gt[type]=1;


        }  // element description loop
    } // !binary
    else // binary case
    {
        int numElementsPartial = 0;
        while(numElementsPartial < numElements)
        {
            int header[3];

            __is.read( (char*)&header, 3*sizeof(int) );
            if(swap) SwapBytes((char*)header, sizeof(int), 3);

            int type = header[0];
            int numElems = header[1];
            int numTags = header[2];
            char const* name;
            CHECK( type < MSH_NUM_TYPE ) << "Invalid GMSH element type " << type << "\n";
            int numVertices = MElement::getInfoMSH(type,&name);
            CHECK( numVertices > 0 ) << "Unsupported element type " << name << "\n";

            unsigned int n = 1 + numTags + numVertices;
            std::vector<int> data(n);
            std::vector<int> indices( numVertices );
            std::vector<rank_type> ghosts;

            for(int i = 0; i < numElems; i++)
            {
                ghosts.clear();
                __is.read( (char*)data.data(), sizeof(int)*n );
                if(swap) SwapBytes((char*)data.data(), sizeof(int), n);

                //if(swap) SwapBytes((char*)data, sizeof(int), n);
                int num = data[0];
                int physical = (numTags > 0) ? data[1] : 0;
                int elementary = (numTags > 1) ? data[2] : 0;
                rank_type numPartitions = (version >= 2.2 && numTags > 3) ? data[3] : 1;
                rank_type partition = (version < 2.2 && numTags > 2) ? data[3]-1 :
                    (version >= 2.2 && numTags > 3) ? data[4]-1 : 0;
                if(numPartitions > 1)
                    for(int j = 0; j < numPartitions - 1; j++)
                        ghosts.push_back( (-data[5 + j]) -1 );
                int parent = (version < 2.2 && numTags > 3) ||
                    (version >= 2.2 && numPartitions && numTags > 3 + numPartitions) ||
                    (version >= 2.2 && !numPartitions && numTags > 2) ?
                    data[numTags] : 0;
                int dom1 = 0, dom2 = 0;

                std::copy( &data[numTags + 1], &data[numTags + 1]+numVertices, indices.begin() );

                // we do not renumber anymore
#if 0
                for(int j = 0; j < numVertices; j++)
                {
                    indices[j] = itoii[ indices[j] ];
                }
#endif
                if ( M_use_elementary_region_as_physical_region )
                {
                    physical = elementary;
                }

                Feel::detail::GMSHElement gmshElt( num, type, physical, elementary,
                                                   numPartitions, partition, ghosts,
                                                   parent, dom1, dom2,
                                                   numVertices, indices,
                                                   this->worldComm().localRank(),
                                                   this->worldComm().localSize(),
                                                   M_respect_partition );

                if ( ( gmshElt.isOnProcessor() == false ||
                       gmshElt.isIgnored(M_ignorePhysicalGroup.begin(), M_ignorePhysicalGroup.end()) ) )
                    continue;
                __et.push_back( gmshElt );

                if ( __gt.find( type ) != __gt.end() )
                    ++__gt[ type ];
                else
                    __gt[type]=1;

            }
            numElementsPartial += numElems;

        } // while
        CHECK( numElementsPartial == numElements ) << "Invalid number of elements read from GMSH, read " << numElementsPartial << " element but expected " << numElements << "\n";
    } // binary

    for ( auto const& it : __gt )
    {
        const char* name;
        MElement::getInfoMSH( it.first, &name );
        LOG(INFO) << "Read " << it.second << " " << name << " elements\n";
    }

    if ( binary )
        __is >> __buf;

    // make sure that we have read everything
    __is >> __buf;
    CHECK( std::string( __buf ) == "$ENDELM" ||
           std::string( __buf ) == "$EndElements" )
        << "invalid end elements string " << __buf
        << " in gmsh importer. It should be either $ENDELM or $EndElements\n";

    // read periodic data if present
    __is >> __buf;
    std::vector<PeriodicEntity> periodic_entities;
    if ( std::string( __buf ) == "$Periodic" )
    {
        int count;
        __is >> count;
        LOG(INFO) << "Reading " << count << " periodic entities\n";
        for(int i = 0; i < count; i++)
        {
            int dim,slave,master;
            __is >> dim >> slave >> master;
            PeriodicEntity e( dim, slave, master );
            int numv;
            __is >> numv;
            for(int j = 0; j < numv; j++)
            {
                int v1,v2;
                __is >> v1 >> v2;
                e.correspondingVertices[v1] = v2;
            }
            CHECK( e.correspondingVertices.size() == numv ) << "Invalid number of vertices in periodic entity"
                                                            << " dim: " << e.dim
                                                            << " slave: " << e.slave
                                                            << " master: " << e.master
                                                            << " got: " << e.correspondingVertices.size()
                                                            << " expected : " << numv << "\n";
            periodic_entities.push_back( e );
        }
        __is >> __buf;
        CHECK( std::string( __buf ) == "$EndPeriodic" )
            << "invalid end $Periodic string " << __buf
            << " in gmsh importer. It should be either $EndPeriodic\n";
    }
    // we are done reading the MSH file


    std::map<int,boost::tuple<int,rank_type> > mapGhostElt;
    std::vector<int> nbMsgToRecv( this->worldComm().localSize(),0 );

    node_type coords( mesh_type::nRealDim );
    //M_n_b_vertices.resize( __n );
    //M_n_b_vertices.assign( __n, 0 );
    for ( auto const& it_gmshElt : __et )
    // add the element to the mesh
    {

        //std::cout<<"isOnProcessor= "<< __et[__i].isOnProcessor()  <<"\n";
        // if the element is not associated to the processor (in partition or ghost) or
        // if the physical entity is ignored
        if ( it_gmshElt.isOnProcessor() == false ||
             it_gmshElt.isIgnored(M_ignorePhysicalGroup.begin(), M_ignorePhysicalGroup.end()) )
            continue;
        // add the points associates to the element on the processor
        for ( uint16_type p = 0; p < it_gmshElt.numVertices; ++p )
        {
            int ptid = it_gmshElt.indices[p];
            // don't do anything if the point is already registered
            if ( mesh->hasPoint( ptid ) )
                continue;

            auto const& gmshpt = gmshpts.find(ptid)->second;
            for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
                coords[j] = gmshpt.x[j];

            point_type pt( ptid, coords, gmshpt.onbdy );
            pt.setProcessIdInPartition( this->worldComm().localRank() );
            if ( gmshpt.parametric )
            {
                pt.setGDim( gmshpt.gdim );
                pt.setGTag( gmshpt.gtag );

                if ( gmshpt.gdim < 3 )
                {
                    pt.setParametricCoordinates( gmshpt.uv[0], gmshpt.uv[1] );
                    mesh->setParametric( true );
                }
            }
            mesh->addPoint( pt );
        } // loop over local points


        switch ( it_gmshElt.type )
        {
            // Points
        case GMSH_POINT:
        {
            addPoint( mesh, it_gmshElt, __idGmshToFeel[it_gmshElt.num] );


            break;
        }

        // Edges
        case GMSH_LINE:
        case GMSH_LINE_2:
        case GMSH_LINE_3:
        case GMSH_LINE_4:
        case GMSH_LINE_5:
        {
            addEdge( mesh, it_gmshElt, __idGmshToFeel[it_gmshElt.num] );

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
            addFace( mesh, it_gmshElt, __idGmshToFeel[it_gmshElt.num] );

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
            addVolume( mesh, it_gmshElt, __idGmshToFeel[it_gmshElt.num] );

            break;
        }

        default:
            break;
        }

        if ( it_gmshElt.isGhost() )
        {
            mapGhostElt.insert( std::make_pair( it_gmshElt.num,boost::make_tuple( __idGmshToFeel[it_gmshElt.num], it_gmshElt.ghostPartitionId() ) ) );
        }
        else if( it_gmshElt.ghosts.size() > 0 )
        {
            for ( auto const& itg : it_gmshElt.ghosts )
                nbMsgToRecv[itg]++;
        }

    } // loop over geometric entities in gmsh file (can be elements or faces)

    // treat periodic entities if any
    for ( auto const& eit : periodic_entities )
    {
        for ( auto const& vit : eit.correspondingVertices )
        {
            auto pit1 = mesh->pointIterator( vit.first );
            auto pit2 = mesh->pointIterator( vit.second );
            CHECK( pit1 != mesh->endPoint() &&
                   pit2 != mesh->endPoint() )
                << "Periodic points data is screwd in periodic entity slave " << eit.slave
                << " master : " << eit.master << " dimension: " << eit.dim;
            auto p1 = *pit1;
            auto p2 = *pit2;
            p1.setMasterId( p2.id() );
            p1.setMasterVertex( boost::addressof( *pit2 ) );

        }
    }
    mesh->setPeriodicEntities( periodic_entities );

    if (VLOG_IS_ON(4))
    {
        //for(int i = 0; i < ptseen.size();  ++i )
        //if ( ptseen[i] == -1 )
        //LOG(WARNING) << "Point with id " << i << " not in element connectivity";
    }

#if !defined( NDEBUG )
    int ne = mesh->numElements();
    int gne = 0;
    mpi::all_reduce( this->worldComm(), ne, gne, [] ( int x, int y )
                     {
                         return x + y;
                     } );

    CHECK( gne > 0 ) << "The mesh does not have any elements.\n"
                     << "something was not right with GMSH mesh importation.\n"
                     << "please check that there are elements of topological dimension "
                     << mesh_type::nDim << "  in the mesh\n";
#endif

    if ( this->worldComm().localSize()>1 )
    {
        if ( false )
            updateGhostCellInfoByUsingBlockingComm( mesh, __idGmshToFeel,  mapGhostElt, nbMsgToRecv );
        else
            updateGhostCellInfoByUsingNonBlockingComm( mesh, __idGmshToFeel,  mapGhostElt, nbMsgToRecv );
    }


    if ( !mesh->markerNames().empty() &&
         ( mesh->markerNames().find("CrossPoints") != mesh->markerNames().end() ) &&
         ( mesh->markerNames().find("WireBasket") != mesh->markerNames().end() ) )
    {
        //LOG(INFO) << "[substructuring] marker cp" << mesh->markerName("CrossPoints")  << "\n";
        //LOG(INFO) << "[substructuring] marker wb" << mesh->markerName("WireBasket")  << "\n";
        //LOG(INFO) << "[substructuring] n cp: " << std::distance( mesh->beginPointWithMarker( mesh->markerName("CrossPoints") ), mesh->endPointWithMarker( mesh->markerName("CrossPoints") ) ) << "\n";
    }
    DVLOG(2) << "done with reading and creating mesh from gmsh file\n";

    toc("read msh from file", FLAGS_v > 0);
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel )
{
    addPoint( mesh, __e, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type*mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<1> )
{
    face_type pf;
    pf.setProcessIdInPartition( this->worldComm().localRank() );
    pf.setId( mesh->numFaces() );
    pf.setMarker( __e.physical );
    pf.setMarker2( __e.elementary );
    pf.setNumberOfPartitions( __e.numPartitions );
    pf.setProcessId( __e.partition );
    pf.setNeighborPartitionIds( __e.ghosts );

    pf.setPoint( 0, mesh->point( __e.indices[0] ) );
//    ptseen[mesh->point( __e.indices[0] ).id()]=1;
    if ( mesh->point( __e.indices[0] ).isOnBoundary() )
    {
        pf.setOnBoundary( true );
    }
    M_n_vertices[ __e.indices[0] ] = 1;

    //M_n_b_vertices[ __e.indices[0] ] = 1;

    face_iterator fit;
    bool inserted;
    boost::tie( fit, inserted ) = mesh->addFace( pf );
    __idGmshToFeel=pf.id();

    DVLOG(2) << "added point on boundary ("
                  << fit->isOnBoundary() << ") with id :" << fit->id() << " and marker " << pf.marker()
                  << " n1: " << mesh->point( __e.indices[0] ).node() << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<2> )
{
    auto pit = mesh->points().modify( mesh->pointIterator(__e.indices[0]),
                                      [&__e]( point_type& e )
                                      {
                                          e.setMarker( __e.physical );
                                          e.setMarker2( __e.elementary );
                                          e.setNumberOfPartitions( __e.numPartitions );
                                          e.setProcessId( __e.partition );
                                          e.setNeighborPartitionIds( __e.ghosts );
                                      } );
    DVLOG(2) << "added point with id :" << mesh->pointIterator(__e.indices[0])->id() << " and marker " << mesh->pointIterator(__e.indices[0])->marker()
                  << " n1: " << mesh->point( __e.indices[0] ).node() << "\n";

}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<3> )
{
    auto pit = mesh->pointIterator(__e.indices[0]);
    bool mod = mesh->points().modify( pit,
                                      [&__e]( point_type& e )
                                      {
                                          e.setMarker( __e.physical );
                                          e.setMarker2( __e.elementary );
                                          e.setNumberOfPartitions( __e.numPartitions );
                                          e.setProcessId( __e.partition );
                                          e.setNeighborPartitionIds( __e.ghosts );
                                      } );
    DVLOG(2) << "added point (modified: " << mod << ")with id :" << pit->id() << " and marker " << pit->marker()
                  << " n1: " << mesh->point( __e.indices[0] ).node() << "\n";
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel )
{
    addEdge( mesh, __e, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type*mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<1> )
{
    element_type e;
    //e.setWorldComm(this->worldComm());
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setMarker( __e.physical );
    e.setMarker2( __e.elementary );
    e.setNumberOfPartitions( __e.numPartitions );
    e.setProcessId( __e.partition );
    e.setNeighborPartitionIds( __e.ghosts );

    if ( __e.type == GMSH_LINE ||
         __e.type == GMSH_LINE_2 ||
         __e.type == GMSH_LINE_3 ||
         __e.type == GMSH_LINE_4 ||
         __e.type == GMSH_LINE_5 )
    {
        int count_pt_on_boundary = 0;
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            if (!e.isGhostCell()) mesh->points().modify( mesh->pointIterator( __e.indices[jj] ), Feel::detail::UpdateProcessId(e.processId()) );
            e.setPoint( jj, mesh->point( __e.indices[jj] ) );
            //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
            if ( mesh->point( __e.indices[jj] ).isOnBoundary() )
                ++count_pt_on_boundary;
        }
        if ( count_pt_on_boundary >= 1 )
        {
            e.setOnBoundary( true );
        }

    }

    mesh->addElement( e );
    __idGmshToFeel=e.id();

    M_n_vertices[ __e.indices[0] ] = 1;
    M_n_vertices[ __e.indices[1] ] = 1;
    DVLOG(2) << "added edge with id :" << e.id()
                  << " n1: " << mesh->point( __e.indices[0] ).node()
                  << " n2: " << mesh->point( __e.indices[1] ).node() << "\n";
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<2> )
{
    face_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setId( mesh->numFaces() );
    e.setMarker( __e.physical );
    e.setMarker2( __e.elementary );
    e.setNumberOfPartitions( __e.numPartitions );
    e.setProcessId( __e.partition );
    e.setNeighborPartitionIds( __e.ghosts );

    if ( __e.type == GMSH_LINE ||
         __e.type == GMSH_LINE_2 ||
         __e.type == GMSH_LINE_3 ||
         __e.type == GMSH_LINE_4 ||
         __e.type == GMSH_LINE_5 )
    {
        DCHECK( __e.indices.size() == npoints_per_edge )
            << "Invalid element indices, "
            << "indices size : " << __e.indices.size() << " points per edge : " << npoints_per_edge;
        int count_pt_on_boundary = 0;
        for ( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
        {
            e.setPoint( jj, mesh->point( __e.indices[jj] ) );
            if ( mesh->point( __e.indices[jj] ).isOnBoundary() )
                ++count_pt_on_boundary;
        }
        if ( count_pt_on_boundary >= 2 )
            e.setOnBoundary( true );

    }

    M_n_vertices[ __e.indices[0] ] = 1;
    M_n_vertices[ __e.indices[1] ] = 1;

    //M_n_b_vertices[ __e.indices[0] ] = 1;
    //M_n_b_vertices[ __e.indices[1] ] = 1;

    bool inserted;
    face_iterator fit;
    boost::tie( fit, inserted ) = mesh->addFace( e );
    __idGmshToFeel=e.id();

    DVLOG(2) << "added edge on boundary ("
                  << fit->isOnBoundary() << ") with id :" << fit->id()
                  << " n1: " << mesh->point( __e.indices[0] ).node()
                  << " n2: " << mesh->point( __e.indices[1] ).node() << "\n";
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addEdge( mesh_type*mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<3> )
{
    edge_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setId( mesh->numEdges() );
    e.setMarker( __e.physical );
    e.setMarker2( __e.elementary );
    e.setNumberOfPartitions( __e.numPartitions );
    // warning : process id is define after (when call mesh->updateForUse()
    // and only for edges which belong to an active 3d element )
    e.setProcessId( invalid_rank_type_value/*__e.partition*/ );
    e.setNeighborPartitionIds( __e.ghosts );

    if ( __e.type == GMSH_LINE ||
         __e.type == GMSH_LINE_2 ||
         __e.type == GMSH_LINE_3 ||
         __e.type == GMSH_LINE_4 ||
         __e.type == GMSH_LINE_5 )
    {
        int count_pt_on_boundary = 0;
        for ( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
        {
            e.setPoint( jj, mesh->point( __e.indices[jj] ) );
            if ( mesh->point( __e.indices[jj] ).isOnBoundary() )
                ++count_pt_on_boundary;
        }
        if ( count_pt_on_boundary >= 2 )
        {
            e.setOnBoundary( true );
        }
    }

    M_n_vertices[ __e.indices[0] ] = 1;
    M_n_vertices[ __e.indices[1] ] = 1;

    //M_n_b_vertices[ __e.indices[0] ] = 1;
    //M_n_b_vertices[ __e.indices[1] ] = 1;

    auto eit = mesh->addEdge( e );
    __idGmshToFeel=eit.id();

    if ( npoints_per_edge == 2 )
        DVLOG(2) << "added edge on boundary ("
                      << eit.isOnBoundary() << ") with id :" << eit.id()
                      << " n1: " << eit.point( 0 ).node()
                      << " n2: " << eit.point( 1 ).node() << "\n";

}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel )
{
    addFace( mesh, __e, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type*, Feel::detail::GMSHElement const&, int & /*__idGmshToFeel*/, mpl::int_<1> )
{}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<2> )
{
    GmshOrdering<element_type> ordering;

    element_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setMarker( __e.physical );
    e.setMarker2( __e.elementary );
    e.setNumberOfPartitions( __e.numPartitions );
    e.setProcessId( __e.partition );
    e.setNeighborPartitionIds( __e.ghosts );

    if ( __e.type == GMSH_QUADRANGLE ||
         __e.type == GMSH_QUADRANGLE_2 ||
         __e.type == GMSH_TRIANGLE ||
         __e.type == GMSH_TRIANGLE_2 ||
         __e.type == GMSH_TRIANGLE_3 ||
         __e.type == GMSH_TRIANGLE_4 ||
         __e.type == GMSH_TRIANGLE_5 )
    {
        int count_pt_on_boundary = 0;
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
            if (!e.isGhostCell()) mesh->points().modify( mesh->pointIterator( __e.indices[jj] ), Feel::detail::UpdateProcessId(e.processId()) );
            e.setPoint( ordering.fromGmshId( jj ), mesh->point( __e.indices[jj] ) );
            if ( mesh->point( __e.indices[jj] ).isOnBoundary() )
                ++count_pt_on_boundary;
        }
        if ( ( __e.type == GMSH_TRIANGLE ||
               __e.type == GMSH_TRIANGLE_2 ||
               __e.type == GMSH_TRIANGLE_3 ||
               __e.type == GMSH_TRIANGLE_4 ||
               __e.type == GMSH_TRIANGLE_5 ) &&
             count_pt_on_boundary >= 2 )
        {
            e.setOnBoundary( true );
        }
        if ( (__e.type == GMSH_QUADRANGLE ||
              __e.type == GMSH_QUADRANGLE_2) &&
             count_pt_on_boundary >= 2 )
            e.setOnBoundary( true );

    }

    mesh->addElement( e );
    __idGmshToFeel=e.id();

    M_n_vertices[ __e.indices[0] ] = 1;
    M_n_vertices[ __e.indices[1] ] = 1;
    M_n_vertices[ __e.indices[2] ] = 1;

    if ( __e.type == GMSH_QUADRANGLE ||
         __e.type == GMSH_QUADRANGLE_2 )
        M_n_vertices[ __e.indices[3] ] = 1;
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<3> )
{
    GmshOrdering<face_type> ordering;

    face_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setId( mesh->numFaces() );
    e.setMarker( __e.physical );
    e.setMarker2( __e.elementary );
    e.setNumberOfPartitions( __e.numPartitions );
    e.setProcessId( __e.partition );
    e.setNeighborPartitionIds( __e.ghosts );

    // we consider in 3D that all faces provided by Gmsh are on the boundary
    // which might not be true
    // \warning there might be a bug here
    e.setOnBoundary( true );

    if ( __e.type == GMSH_QUADRANGLE ||
         __e.type == GMSH_QUADRANGLE_2 ||
         __e.type == GMSH_TRIANGLE ||
         __e.type == GMSH_TRIANGLE_2 ||
         __e.type == GMSH_TRIANGLE_3 ||
         __e.type == GMSH_TRIANGLE_4 ||
         __e.type == GMSH_TRIANGLE_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
        {
            //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
            e.setPoint( ordering.fromGmshId( jj ), mesh->point( __e.indices[jj] ) );
        }

        //e.setPoint( jj, mesh->point( __e[jj] ) );
    }

    bool inserted;
    face_iterator fit;
    boost::tie( fit, inserted ) = mesh->addFace( e );

    __idGmshToFeel=e.id();

    M_n_vertices[ __e.indices[0] ] = 1;
    M_n_vertices[ __e.indices[1] ] = 1;
    M_n_vertices[ __e.indices[2] ] = 1;

    if ( __e.type == GMSH_QUADRANGLE ||
         __e.type == GMSH_QUADRANGLE_2 )
        M_n_vertices[ __e.indices[3] ] = 1;

}

template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel )
{
    addVolume( mesh, __e, __idGmshToFeel, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type*, Feel::detail::GMSHElement const&, int & /*__idGmshToFeel*/, mpl::int_<1> )
{}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type*, Feel::detail::GMSHElement const&, int & /*__idGmshToFeel*/, mpl::int_<2> )
{}
template<typename MeshType>
void
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e, int & __idGmshToFeel, mpl::int_<3> )
{
    element_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    GmshOrdering<element_type> ordering;
    e.setMarker( __e.physical );
    e.setMarker2( __e.elementary );
    e.setNumberOfPartitions( __e.numPartitions );
    e.setProcessId( __e.partition );
    e.setNeighborPartitionIds( __e.ghosts );


    //
    // WARNING: not yet done for high order in 3D !!!
    //
    if ( __e.type == GMSH_HEXAHEDRON ||
         __e.type == GMSH_HEXAHEDRON_2 ||
         __e.type == GMSH_TETRAHEDRON ||
         __e.type == GMSH_TETRAHEDRON_2 ||
         __e.type == GMSH_TETRAHEDRON_3 ||
         __e.type == GMSH_TETRAHEDRON_4 ||
         __e.type == GMSH_TETRAHEDRON_5 )
    {
        int count_pt_on_boundary = 0;
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
            if (!e.isGhostCell()) mesh->points().modify( mesh->pointIterator( __e.indices[jj] ), Feel::detail::UpdateProcessId(e.processId()) );
            //std::cout << "gmsh index " << jj << " -> " << ordering.fromGmshId(jj) << " -> " << mesh->point( __e[jj] ).id()+1 << " : " << mesh->point( __e[jj] ).node() << "\n";
            e.setPoint( ordering.fromGmshId( jj ), mesh->point( __e.indices[jj] ) );

            if ( mesh->point( __e.indices[jj] ).isOnBoundary() )
                ++count_pt_on_boundary;
        }
        // the tet share a face with the boundary
        if ( ( __e.type == GMSH_TETRAHEDRON ||
               __e.type == GMSH_TETRAHEDRON_2 ||
               __e.type == GMSH_TETRAHEDRON_3 ||
               __e.type == GMSH_TETRAHEDRON_4 ||
               __e.type == GMSH_TETRAHEDRON_5 )
             && count_pt_on_boundary >= 3)
            e.setOnBoundary( true );
        // the hex share a face with the boundary
        if ( ( __e.type == GMSH_HEXAHEDRON ||
               __e.type == GMSH_HEXAHEDRON_2 )
             && count_pt_on_boundary >= 4)
            e.setOnBoundary( true );
    }

    mesh->addElement( e );
    __idGmshToFeel=e.id();

    M_n_vertices[ __e.indices[0] ] = 1;
    M_n_vertices[ __e.indices[1] ] = 1;
    M_n_vertices[ __e.indices[2] ] = 1;
    M_n_vertices[ __e.indices[3] ] = 1;

    if ( __e.type == GMSH_HEXAHEDRON )
    {
        M_n_vertices[ __e.indices[4] ] = 1;
        M_n_vertices[ __e.indices[5] ] = 1;
        M_n_vertices[ __e.indices[6] ] = 1;
        M_n_vertices[ __e.indices[7] ] = 1;
    }
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::updateGhostCellInfoByUsingBlockingComm( mesh_type* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
                                                                std::vector<int> const& nbMsgToRecv )
{
    // counter of msg sent for each process
    std::vector<int> nbMsgToSend( this->worldComm().localSize(),0 );
    // map usefull to get final result
    std::vector< std::map<int,int> > mapMsg( this->worldComm().localSize() );

    // iterate over ghost elt
    auto it_map = __mapGhostElt.begin();
    auto const en_map = __mapGhostElt.end();
    for ( ; it_map!=en_map ; ++it_map )
    {
        auto const idGmsh = it_map->first;
        auto const idProc = it_map->second.template get<1>();
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

#if !defined( NDEBUG )
    // check nbMsgToRecv computation
    std::vector<int> nbMsgToRecv2;
    mpi::all_to_all( this->worldComm().localComm(),
                     nbMsgToSend,
                     nbMsgToRecv2 );
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        CHECK( nbMsgToRecv[proc]==nbMsgToRecv2[proc] ) << "paritioning data incorect "
                                                       << "myrank " << this->worldComm().localRank() << " proc " << proc
                                                       << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
                                                       << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc] << "\n";
    }
#endif

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
                      << " the G " << mesh->element( __idGmshToFeel.find(idGmsh)->second ).G()
                      << " with idFeel Classic " << __idGmshToFeel.find(idGmsh)->second
                      << std::endl;
#endif
            //re-send idFeel
            this->worldComm().localComm().send( proc, cpt, __idGmshToFeel.find( idGmsh )->second );
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
            mesh->elements().modify( elttt, Feel::detail::updateIdInOthersPartitions( proc, idFeel ) );
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

}

template<typename MeshType>
void
ImporterGmsh<MeshType>::updateGhostCellInfoByUsingNonBlockingComm( mesh_type* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
                                                                   std::vector<int> const& nbMsgToRecv )
{

    DVLOG(1) << "updateGhostCellInfoNonBlockingComm : start on rank "<< this->worldComm().localRank() << "\n";

    const int nProc = this->worldComm().localSize();

    //-----------------------------------------------------------//
    // compute size of container to send
    std::map< rank_type, int > nDataInVecToSend;
    auto it_map = __mapGhostElt.begin();
    auto const en_map = __mapGhostElt.end();
    for ( ; it_map!=en_map ; ++it_map )
    {
        const rank_type idProc = it_map->second.template get<1>();
        if ( nDataInVecToSend.find(idProc) == nDataInVecToSend.end() )
            nDataInVecToSend[idProc]=0;
        nDataInVecToSend[idProc]++;
    }
    //-----------------------------------------------------------//
    // init and resize the container to send
    std::map< rank_type, std::vector<int> > dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //-----------------------------------------------------------//
    // prepare container to send
    std::map< rank_type, std::map<int,int> > memoryMsgToSend;
    std::map< rank_type, int > nDataInVecToSendBis;
    it_map = __mapGhostElt.begin();
    for ( ; it_map!=en_map ; ++it_map )
    {
        const int idGmsh = it_map->first;
        const int idFeel = it_map->second.template get<0>();
        const rank_type idProc = it_map->second.template get<1>();

        if ( nDataInVecToSendBis.find(idProc) == nDataInVecToSendBis.end() )
            nDataInVecToSendBis[idProc]=0;
        // save request
        memoryMsgToSend[idProc][nDataInVecToSendBis[idProc]] = idFeel;
        // update container
        dataToSend[idProc][nDataInVecToSendBis[idProc]] = idGmsh;
        // update counter
        nDataInVecToSendBis[idProc]++;
    }
    //-----------------------------------------------------------//
    // counter of request
    int nbRequest=0;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( dataToSend.find(proc) != dataToSend.end() )
            ++nbRequest;
        if ( nbMsgToRecv[proc] > 0 )
            ++nbRequest;
    }
    if ( nbRequest ==0 ) return;

    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    //-----------------------------------------------------------//
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToSend->first , 0, itDataToSend->second );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // first recv
    std::map<rank_type,std::vector<int> > dataToRecv;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( nbMsgToRecv[proc] > 0 )
        {
            reqs[cptRequest] = this->worldComm().localComm().irecv( proc , 0, dataToRecv[proc] );
            ++cptRequest;
        }
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //-----------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector<int> > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        dataToReSend[idProc].resize( nDataRecv );
        //store the idFeel corresponding
        for ( int k=0; k<nDataRecv; ++k )
            dataToReSend[idProc][k] = __idGmshToFeel.find( itDataRecv->second[k] )->second;
    }
    //-----------------------------------------------------------//
    // send respond to the request
    cptRequest=0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend!=enDataToReSend ; ++itDataToReSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToReSend->first , 0, itDataToReSend->second );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // recv the initial request
    std::map<rank_type, std::vector<int> > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc, 0, finalDataToRecv[idProc] );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //-----------------------------------------------------------//
    // update mesh : id in other partitions for the ghost cells
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const rank_type idProc = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        for ( int k=0; k<nDataRecv; ++k )
        {
            auto eltToUpdate = mesh->elementIterator( memoryMsgToSend[idProc][k],idProc );
            mesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( idProc, itFinalDataToRecv->second[k] ) );
        }
    }
    //-----------------------------------------------------------//
    DVLOG(1) << "updateGhostCellInfoNonBlockingComm : finish on rank "<< this->worldComm().localRank() << "\n";
}


} // Feel



#endif /* __ImporterGmsh_H */
