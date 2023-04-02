/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-16

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2009-2017 Feel++ Consortium

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
#ifndef FEELPP_IMPORTERGMSH_HPP
#define FEELPP_IMPORTERGMSH_HPP 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <unordered_map>

#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelcore/feelgmsh.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/importer.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feeltiming/tic.hpp>
#include <boost/algorithm/string/trim.hpp>

#if !defined( FEELPP_HAS_GMSH_API )
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

#endif

// from Gmsh
#if defined( FEELPP_HAS_GMSH_API )
int getInfoMSH(const int typeMSH, std::string & elementName);
#else
int getInfoMSH(const int typeMSH, const char **const name);
#endif


namespace Feel
{
void SwapBytes( char* array, int size, int n );
namespace detail
{
#pragma GCC visibility push(hidden)
class FEELPP_NO_EXPORT  GMSHPoint
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
    GMSHPoint( int _id, Eigen::Vector3d const& _x, bool _p, int _gdim, int _gtag, Eigen::Vector2d _uv )
        :
        x(_x),
        uv( _uv ),
        id( _id ),
        onbdy( false ),
        parametric( _p ),
        gdim( _gdim ),
        gtag( _gtag )
        {}
    GMSHPoint( int _id, Eigen::Vector3d const& _x )
        :
        x(_x),
        uv(),
        id( _id ),
        onbdy( false ),
        parametric( false ),
        gdim( 0 ),
        gtag( 0 )
        {}
    GMSHPoint(GMSHPoint const& p ) = default;
    GMSHPoint(GMSHPoint && p ) = default;
    GMSHPoint& operator=( GMSHPoint const& ) = default;
    GMSHPoint& operator=( GMSHPoint && ) = default;
};
struct FEELPP_NO_EXPORT GMSHElement
{
    GMSHElement()
        :
        num( 0 ),
#if defined( FEELPP_HAS_GMSH_API )
        type( 1 ),
#else
        type( MSH_PNT ),
#endif
        physical(),
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
                 std::vector<int> const& p,
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
#if defined( FEELPP_HAS_GMSH_H )
    template<typename T>
    GMSHElement( T* ele, int n, int elementary, std::vector<int> const& _physicals )
        :
        GMSHElement(n,
                    ele->getTypeForMSH(),
                    _physicals,
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
#endif // FEELPP_HAS_GMSH_H
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

    GMSHElement( GMSHElement const& g ) = default;
    GMSHElement( GMSHElement && g ) = default;
    GMSHElement& operator=( GMSHElement const& ) = default;
    GMSHElement& operator=( GMSHElement && ) = default;

    bool isOnProcessor() const { return is_on_processor; }
    bool isGhost() const { return is_ghost; }
    rank_type ghostPartitionId() const { return ghost_partition_id; }

    template<typename IteratorBegin,typename IteratorEnd>
    bool isIgnored(IteratorBegin it, IteratorEnd en ) const
    {
        if ( physical.empty() )
            return false;
        bool res = true;
        for ( int p : physical )
            res = res && ( std::find( it, en, p ) != en );
        return res;
    }

    void updatePartition( std::map<rank_type,rank_type> const& p2e, rank_type worldcommrank, rank_type worldcommsize )
        {
            partition = num;
            for( auto& g : ghosts )
                g = p2e.at(g)-1;
            setPartition( worldcommrank, worldcommsize );
        }
    int num;
    int type;
    std::vector<int> physical;
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
#pragma GCC visibility pop
//!
//! @return true if element should be store on process, false otherwise
//!
FEELPP_EXPORT bool isOnProcessor( std::vector<rank_type> ghosts, rank_type partition, rank_type worldcommrank, rank_type worldcommsize);

//!
//! @return true if \p physical is found in container, false otherwise
//!
template<typename Iterator>
FEELPP_EXPORT bool isFound( Iterator beg, Iterator end, int physical )
{
    return std::find( beg, end, physical ) != end;
}


template<typename Iterator>
FEELPP_EXPORT std::vector<int> gmshPhysicalTags( Iterator const& it, int elementary, bool use_elementary_region_as_physical_region )
{
    if ( use_elementary_region_as_physical_region )
        return std::vector<int>( { elementary } );
    else
        return (*it)->physicals;
}

FEELPP_STRONG_INLINE rank_type gmshPartitionTagToProcessId( int partitionTag, rank_type worldcommsize )
{
    return partitionTag % worldcommsize;
}

} // detail

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
class FEELPP_EXPORT ImporterGmsh
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

    explicit ImporterGmsh( worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() )
        :
        super( GMSH, UNDEFINED,  _worldcomm ),
        M_version( FEELPP_GMSH_FORMAT_VERSION ),
        M_in_memory( false ),
        M_use_elementary_region_as_physical_region( false ),
        M_respect_partition( false ),
        M_scale( 1 ),
        M_deleteGModelAfterUse( false )
    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
    }

    explicit ImporterGmsh( std::string const& _fname, std::string _version = FEELPP_GMSH_FORMAT_VERSION,
                           worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() )
        :
        super( _fname, GMSH, UNDEFINED, _worldcomm ),
        M_version( _version ),
        M_in_memory( false ),
        M_use_elementary_region_as_physical_region( false ),
        M_respect_partition( false ),
        M_scale( 1 ),
        M_deleteGModelAfterUse( false )
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
        M_respect_partition( i.M_respect_partition ),
        M_scale( i.M_scale ),
        M_gmodelName( i.M_gmodelName ),
        M_deleteGModelAfterUse( i.M_deleteGModelAfterUse )

    {
        this->setIgnorePhysicalName( "FEELPP_GMSH_PHYSICALNAME_IGNORED" );
        //showMe();
    }
    ~ImporterGmsh() override
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
    void setScaling( double scale )
    {
        M_scale = scale;
    }

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

    void setGModelName( std::string const& name ) { M_gmodelName = name; }

    void setDeleteGModelAfterUse( bool b ) { M_deleteGModelAfterUse = b; }


    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh ) override;

    void showMe() const;

    //@}



protected:

private:
    FEELPP_NO_EXPORT void readFromMemory( mesh_type* mesh );
    FEELPP_NO_EXPORT void readFromFile( mesh_type* mesh );
    FEELPP_NO_EXPORT void readFromFileVersion2( mesh_type* mesh, std::ifstream & __is, char __buf[],
                                                double version, bool binary, bool swap );
    template <typename gmsh_size_type,typename gmsh_size_partition_type,typename gmsh_size_periodiclink_type,typename gmsh_elttag_type>
    FEELPP_NO_EXPORT void readFromFileVersion4( mesh_type* mesh, std::ifstream & __is, char __buf[],
                                                double version, bool binary, bool swap );


    FEELPP_NO_EXPORT void addVertices( mesh_type* mesh, Feel::detail::GMSHElement const& elt,
                                       std::map<int, Feel::detail::GMSHPoint > const& gmshpts );

    FEELPP_NO_EXPORT int addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e );

    FEELPP_NO_EXPORT int addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e );
    FEELPP_NO_EXPORT int addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<1> );
    FEELPP_NO_EXPORT int addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<2> );
    FEELPP_NO_EXPORT int addEdge( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, mpl::int_<3> );

    FEELPP_NO_EXPORT int addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e );
    FEELPP_NO_EXPORT int addFace( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, mpl::int_<1> );
    FEELPP_NO_EXPORT int addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<2> );
    FEELPP_NO_EXPORT int addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<3> );

    FEELPP_NO_EXPORT int addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e );
    FEELPP_NO_EXPORT int addVolume( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, mpl::int_<1> );
    FEELPP_NO_EXPORT int addVolume( mesh_type* /*mesh*/, Feel::detail::GMSHElement const& /*__e*/, mpl::int_<2> );
    FEELPP_NO_EXPORT int addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<3> );

    FEELPP_NO_EXPORT void updateGhostCellInfo( mesh_type* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
                                                    std::vector<int> const& nbMsgToRecv );


private:

    std::string M_version;

    bool M_in_memory;

    std::set<int> M_ignorePhysicalGroup;
    std::set<std::string> M_ignorePhysicalName;
    bool M_use_elementary_region_as_physical_region;
    bool M_respect_partition;
    double M_scale;
    //std::map<int,int> itoii;
    //std::vector<int> ptseen;

    std::string M_gmodelName;
    bool M_deleteGModelAfterUse;
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
            coords[j] = gmshpt.x[j]*M_scale;

        point_type pt( ptid, coords );
        pt.setProcessIdInPartition( this->worldComm().localRank() );
        if ( gmshpt.parametric )
        {
            if ( gmshpt.gdim < 3 )
            {
                pt.setParametricCoordinates( gmshpt.gdim, gmshpt.gtag, gmshpt.uv[0], gmshpt.uv[1] );
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

    if ( !M_gmodelName.empty() && ( M_in_memory == true ) )
    {
#if defined( FEELPP_HAS_GMSH_H )
        readFromMemory( mesh );
#else
        LOG(WARNING) << "Gmsh is not available. Cannot read from memory Gmsh directly. Use file instead.";
        readFromFile( mesh );
#endif
    }
    else
        readFromFile( mesh );
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::readFromMemory( mesh_type* mesh )
{
#if defined( FEELPP_HAS_GMSH_H )
    tic();
    LOG(INFO) << "Reading Msh from memory ";

#if defined( FEELPP_HAS_GMSH_API )
    rank_type procId = this->worldComm().localRank();
    rank_type worldSize = this->worldComm().localSize();

    std::map<int,int> __idGmshToFeel; // id Gmsh to id Feel

    gmsh::vectorpair dimTagsEntities;
    gmsh::model::getEntities( dimTagsEntities );

    if ( M_use_elementary_region_as_physical_region )
    {
        for ( auto const& [ entityDim, entityTag ] : dimTagsEntities )
        {
            std::string entityName;
#if GMSH_VERSION_GREATER_OR_EQUAL_THAN(4,2,0)
            gmsh::model::getEntityName( entityDim,entityTag,entityName );
#endif
            boost::trim(entityName); // get rid of eventual trailing spaces
            if ( !entityName.empty() )
                mesh->addMarkerName( entityName, entityTag, entityDim );
        }
    }
    else
    {
        gmsh::vectorpair dimTagsPhysicalGroups;
        gmsh::model::getPhysicalGroups( dimTagsPhysicalGroups );
        for ( auto const& [ dim, tag ] : dimTagsPhysicalGroups )
        {
            std::string physicalName;
            gmsh::model::getPhysicalName( dim,tag,physicalName );
            boost::trim(physicalName); // get rid of eventual trailing spaces
            if ( !physicalName.empty() )
                mesh->addMarkerName( physicalName, tag, dim );
        }
    }

    for ( auto const& [ entityDim, entityTag ] : dimTagsEntities )
    {

#if 0
        std::vector<int> partitions;
        gmsh::model::getPartitions(entityDim,entityTag,partitions);
        std::cout << "partitions="<<partitions<<"\n";
#endif
#if GMSH_VERSION_GREATER_OR_EQUAL_THAN(4,2,0)
        std::vector<std::size_t> nodeTags;
#else
        std::vector<int> nodeTags;
#endif
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes( nodeTags,coord,parametricCoord,
                                     entityDim, entityTag
                                     /*const bool includeBoundary = false*/);
        std::size_t nNodesInEntity = nodeTags.size();
        for ( std::size_t k=0;k<nNodesInEntity;++k )
        {
            std::size_t ptId = nodeTags[k];
            node_type coordFeel( mesh_type::nRealDim );
            for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
                coordFeel[j] = coord[k*3+j]*M_scale;

            point_type pt( ptId, coordFeel );
            pt.setProcessIdInPartition( procId/*this->worldComm().localRank()*/ );
            pt.setProcessId( procId );
            mesh->addPoint( pt );
        }

        // physical tags
        std::vector<int> physicalTags;
        if ( M_use_elementary_region_as_physical_region )
            physicalTags.push_back( entityTag );
        else
            gmsh::model::getPhysicalGroupsForEntity( entityDim, entityTag, physicalTags );

        std::vector<int> elementTypes;
#if GMSH_VERSION_GREATER_OR_EQUAL_THAN(4,2,0)
        std::vector<std::vector<std::size_t> > elementTags;
        std::vector<std::vector<std::size_t> > nodeTagsInElements;
#else
        std::vector<std::vector<int> > elementTags;
        std::vector<std::vector<int> > nodeTagsInElements;
#endif
        gmsh::model::mesh::getElements( elementTypes, elementTags, nodeTagsInElements, entityDim, entityTag );
        for ( int et=0;et<elementTypes.size();++et )
        {
            size_type eltId = 0;
            int elementType = elementTypes[et];
            // partitioning data
            int numPartitions = 1, partition = procId;
            int parent =0 , dom1 = 0, dom2 = 0;
            std::vector<rank_type> ghosts;

            std::string ename;
            int numVertices = getInfoMSH( elementType,ename );

            std::vector<int> indices( numVertices );
            Feel::detail::GMSHElement it_gmshElt( eltId, elementType, physicalTags, entityTag,    //   num, type, physical, elementary,
                                                  numPartitions, partition, ghosts,
                                                  parent, dom1, dom2,
                                                  numVertices, indices,
                                                  this->worldComm().localRank(),
                                                  this->worldComm().localSize(),
                                                  M_respect_partition );

            std::size_t nElementsOfThisType = elementTags[et].size();

            for ( std::size_t k=0;k<nElementsOfThisType;++k )
            {
                it_gmshElt.num = elementTags[et][k];
                for ( int i=0;i<numVertices;++i )
                    indices[i] = nodeTagsInElements[et][k*numVertices+i];
                it_gmshElt.indices = indices;


                switch ( elementType )
                {
                // Points
                case GMSH_POINT:
                {
                    __idGmshToFeel[it_gmshElt.num] = addPoint( mesh, it_gmshElt );
                    break;
                }
                // Edges
                case GMSH_LINE:
                case GMSH_LINE_2:
                case GMSH_LINE_3:
                case GMSH_LINE_4:
                case GMSH_LINE_5:
                {
                    __idGmshToFeel[it_gmshElt.num] = addEdge( mesh, it_gmshElt );
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
                    __idGmshToFeel[it_gmshElt.num] = addFace( mesh, it_gmshElt );
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
                    __idGmshToFeel[it_gmshElt.num] = addVolume( mesh, it_gmshElt );
                    break;
                }
                default:
                    break;
                }

            }
        }
    }
    if ( M_deleteGModelAfterUse )
        gmsh::model::remove();

    // update ordered points in mesh data structure
    mesh->updateOrderedPoints();
#else

    GModel* gmodel = GModel::findByName( M_gmodelName );
    CHECK( gmodel != nullptr ) << "gmodel not found with name " << M_gmodelName;


    // get the number of vertices and index the vertices in a continuous
    // sequence
    bool saveAll=true;
    bool saveSinglePartition=false;
    int numVertices = gmodel->indexMeshVertices(saveAll, saveSinglePartition);


    if( gmodel->numPhysicalNames() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "  +- number of physicals: " << gmodel->numPhysicalNames() << "\n";
        for(GModel::piter it = gmodel->firstPhysicalName(); it != gmodel->lastPhysicalName(); it++)
        {
            int id = it->first.second;
            int topodim = it->first.first;
            std::string marker = it->second.c_str();
            boost::trim(marker); // get rid of eventual trailing spaces
            mesh->addMarkerName( marker, id, topodim );
        }
    }

    // get the number of elements we need to save
    int numElements = gmodel->getNumMeshElements();
    if ( Environment::isMasterRank() )
        std::cout << "  +- number of vertices: " << numVertices << "\n"
                  << "  +- number of elements: " << numElements << "\n";

    std::map<int, Feel::detail::GMSHPoint > gmshpts;
    std::vector<GEntity*> entities;
    gmodel->getEntities(entities);
    for(unsigned int i = 0; i < entities.size(); i++)
        for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
        {
            auto const& pt = entities[i]->mesh_vertices[j];
            int id = pt->getIndex();
            Eigen::Vector3d x(pt->x()*M_scale, pt->y()*M_scale, pt->z()*M_scale);
            gmshpts[id].id = id;
            gmshpts[id].x = x;
            gmshpts[id].parametric = false;
        }
    std::map<int,int> __idGmshToFeel; // id Gmsh to id Feel
    int num = 0;
    // read vertices
    for( GModel::viter it = gmodel->firstVertex(); it != gmodel->lastVertex(); ++it )
    {
        int elementary = (*it)->tag();
        std::vector<int> physicals = Feel::detail::gmshPhysicalTags( it, elementary, M_use_elementary_region_as_physical_region );
        for( auto const& p : (*it)->points )
        {
            if ( physicals.size() )
            {
                Feel::detail::GMSHElement ele(p, ++num, elementary, physicals);
                this->addVertices( mesh, ele, gmshpts );
                __idGmshToFeel[p->getNum()] = addPoint( mesh, ele );
            }
        }
    }
    // read edges
    for( GModel::eiter it = gmodel->firstEdge(); it != gmodel->lastEdge(); ++it )
    {
        int elementary = (*it)->tag();
        std::vector<int> physicals = Feel::detail::gmshPhysicalTags( it, elementary, M_use_elementary_region_as_physical_region );
        for( auto const& edge : (*it)->lines )
        {
            Feel::detail::GMSHElement ele(edge, ++num, elementary, physicals);
            this->addVertices( mesh, ele, gmshpts );
            __idGmshToFeel[edge->getNum()] = addEdge( mesh, ele );
        }
    }
    // read triangles
    for( GModel::fiter it = gmodel->firstFace(); it != gmodel->lastFace(); ++it)
    {
        int elementary = (*it)->tag();
        std::vector<int> physicals = Feel::detail::gmshPhysicalTags( it, elementary, M_use_elementary_region_as_physical_region );
        for( auto const& triangle : (*it)->triangles )
        {
            Feel::detail::GMSHElement ele(triangle, ++num, elementary,physicals);
            this->addVertices( mesh, ele, gmshpts );
            __idGmshToFeel[triangle->getNum()] = addFace( mesh, ele );
        }
    }
    // read quadrangles
    for( GModel::fiter it = gmodel->firstFace(); it != gmodel->lastFace(); ++it)
    {
        int elementary = (*it)->tag();
        std::vector<int> physicals = Feel::detail::gmshPhysicalTags( it, elementary, M_use_elementary_region_as_physical_region );
        for( auto const& quad : (*it)->quadrangles )
        {
            Feel::detail::GMSHElement ele(quad, ++ num, elementary, physicals);
            this->addVertices( mesh, ele, gmshpts );
            __idGmshToFeel[quad->getNum()] = addFace( mesh, ele );
        }
    }
    // read tetras
    for( GModel::riter it = gmodel->firstRegion(); it != gmodel->lastRegion(); ++it)
    {
        int elementary = (*it)->tag();
        std::vector<int> physicals = Feel::detail::gmshPhysicalTags( it, elementary, M_use_elementary_region_as_physical_region );
        for( auto const& tetra : (*it)->tetrahedra )
        {
            Feel::detail::GMSHElement ele(tetra, ++num, elementary, physicals);
            this->addVertices( mesh, ele, gmshpts );
            __idGmshToFeel[tetra->getNum()] = addVolume( mesh, ele );   }
    }
    // read hexa
    for( GModel::riter it = gmodel->firstRegion(); it != gmodel->lastRegion(); ++it)
    {
        int elementary = (*it)->tag();
        std::vector<int> physicals = Feel::detail::gmshPhysicalTags( it, elementary, M_use_elementary_region_as_physical_region );
        for( auto const& hexa : (*it)->hexahedra )
        {
            Feel::detail::GMSHElement ele(hexa, ++num, elementary, physicals);
            this->addVertices( mesh, ele, gmshpts );
            __idGmshToFeel[hexa->getNum()] = addVolume( mesh, ele );
        }
    }

    if ( M_deleteGModelAfterUse )
        delete gmodel;
#endif
    LOG(INFO) << "Reading Msh from memory done";
    toc("read msh from memory");
#else
    LOG(WARNING) << "Gmsh library not available, cannot read mesh in memory ";
#endif
}
template<typename MeshType>
void
ImporterGmsh<MeshType>::readFromFile( mesh_type* mesh )
{
    tic();
    LOG(INFO) << "Reading Msh file " << this->filename();

    std::ifstream __is ( this->filename().c_str() );
    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        LOG(ERROR) << "Invalid file name " << this->filename() << " (file not found)";
        ostr << "Invalid file name " << this->filename() << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    char __buf[256];
    __is.getline(__buf, 256);

    std::string theversion;
    double version = boost::lexical_cast<double>( FEELPP_GMSH_FORMAT_VERSION );
    bool binary = false, swap = false;

    if ( std::string( __buf ) == "$MeshFormat" )
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
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndMeshFormat" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndMeshFormat\n";
        __is.getline(__buf, 256);
    }


    if ( std::string( __buf ) == "$PhysicalNames" )
    {
        std::vector<MeshMarkerName> meshMarkerNameMap = markerMap(MeshType::nDim);

        int nnames;
        __is >> nnames;

        for ( int n = 0; n < nnames; ++n )
        {
            int id, topodim;
            std::string name;

            if ( version >= 2.1 )
            {
                __is >> topodim >> id >> name;
                DVLOG(2) << "[importergmsh] reading topodim: "  << topodim << " id: " << id << " name: " << name << "\n";
            }

            else if ( version == 2.0 )
                __is >> id >> name;

            boost::trim( name );
            boost::trim_if( name,boost::is_any_of( "\"" ) );

            if ( meshMarkerNameMap.empty() )
            {
                mesh->addMarkerName( name, id, topodim );
            }
            if ( M_ignorePhysicalName.find( name )!=M_ignorePhysicalName.end() ) this->setIgnorePhysicalGroup( id );
        }
        if ( meshMarkerNameMap.empty() )
        {
            CHECK( mesh->markerNames().size() == ( size_type )nnames ) << "invalid number of physical names" << mesh->markerNames().size() << " vs " << nnames;
        }
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndPhysicalNames" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndPhysicalNames\n";

        for( auto const& it : meshMarkerNameMap )
        {
            mesh->addMarkerName( it.name, it.ids[0], it.ids[1] );
        }

        __is.getline(__buf, 256);
    }

    if ( version >= 2 && version < 4 )
        this->readFromFileVersion2( mesh, __is, __buf, version, binary, swap );
    else if ( version == 4 )
        this->readFromFileVersion4<unsigned long,int,int,int>( mesh, __is, __buf, version, binary, swap );
    else
        this->readFromFileVersion4<size_t,size_t,size_t,size_t>( mesh, __is, __buf, version, binary, swap );

    toc("read msh from file", FLAGS_v > 0);
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::readFromFileVersion2( mesh_type* mesh, std::ifstream & __is, char __buf[],
                                              double version, bool binary, bool swap )
{
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


    std::unordered_map<int, Feel::detail::GMSHPoint > gmshpts;
    LOG(INFO) << "Reading "<< __n << " nodes\n";

    Eigen::Vector3d x;
    Eigen::Vector2d uv;
    tic();
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
        if ( has_parametric_nodes )
        {
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
            gmshpts.emplace_hint( gmshpts.end(),std::make_pair(id,Feel::detail::GMSHPoint( id, x, true, gdim, gtag, uv ) ) );
        }
        else
            gmshpts.emplace_hint( gmshpts.end(),std::make_pair(id,Feel::detail::GMSHPoint( id, x ) ) );
        
        // stores mapping to be able to reorder the indices
        // so that they are contiguous
        //itoii[idpts[__i]] = __i;
    }
    toc("ImporterGmsh::readFromFile read points", FLAGS_v > 0 );
    //ptseen.resize( __n );
    //std::fill( ptseen.begin(), ptseen.end(), -1 );
    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    if ( binary )
        __is.get();

    __is.getline( __buf, 256 );
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
    __is.getline(__buf, 256);
    CHECK( std::string( __buf ) == "$ELM" ||
           std::string( __buf ) == "$Elements" )
        << "invalid elements string " << __buf << " in gmsh importer, it should be either $ELM or $Elements\n";

    int numElements;
    __is >> numElements;

    // eat  '\n' in binary mode otherwise the next binary read will get screwd
    if ( binary )
        __is.get();

    LOG(INFO) << "Reading " << numElements << " elements...\n";
    std::vector<Feel::detail::GMSHElement> __et; // tags in each element
    __et.reserve( numElements );
    std::map<int,int> __idGmshToFeel; // id Gmsh to id Feel
    std::map<int,int> __gt;

    rank_type theCurrentPartitionId = (this->worldComm().globalSize()>1)?this->worldComm().localRank():0;

    if ( !binary )
    {
        tic();
        for(int i = 0; i < numElements; i++)
        {
          int num, type, physical = 0, elementary = 0, parent = 0;
          int dom1 = 0, dom2 = 0, numVertices;
          std::vector<rank_type> ghosts;
          int numTags;
          // some faces may not be associated to a partition in the mesh file,
          // hence will be read given the partition id 0 and will be discarded
          rank_type partition = theCurrentPartitionId;
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
              else if(version >= 2.2 && j == 2 && numTags > 3) { numPartitions = tag; ghosts.reserve( numPartitions-1 ); }
              else if(version >= 2.2 && j == 3) partition = tag-1;
              else if(j >= 4 && j < 4 + numPartitions - 1) ghosts.push_back((-tag)-1);
              else if(j == 3 + numPartitions && (numTags == 4 + numPartitions))
                  parent = tag;
              else if(j == 3 + numPartitions && (numTags == 5 + numPartitions)) {
                  dom1 = tag; j++;
                  __is >> dom2;
              }
          }
#if defined( FEELPP_HAS_GMSH_API )
          std::string ename;
          numVertices = getInfoMSH(type,ename);
#else
          const char* ename;
          numVertices = getInfoMSH(type,&ename);
#endif

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
          // discard elements not on processor or with an ignored physical
          // WARNING: we had another condition if the number of processors and
          // elements is the same, in that case we store everything for now
          if ( Feel::detail::isOnProcessor( ghosts,
                                            partition,
                                            this->worldComm().localRank(),
                                            this->worldComm().localSize() ) == false ||
               Feel::detail::isFound( M_ignorePhysicalGroup.begin(), M_ignorePhysicalGroup.end(), physical ) )
              continue;
          
          __et.emplace_back( num, type, std::vector<int>({physical}), elementary,
                             numPartitions, partition, ghosts,
                             parent, dom1, dom2,
                             numVertices, indices,
                             this->worldComm().localRank(),
                             this->worldComm().localSize(),
                             M_respect_partition );
          
          if ( __gt.find( type ) != __gt.end() )
              ++__gt[ type ];
          else
              __gt[type]=1;


        }  // element description loop
        toc("ImporterGmsh::readFromFile read and store GMSHElement", FLAGS_v > 0 );
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
#if defined( FEELPP_HAS_GMSH_API )
            std::string name;
            int numVertices = getInfoMSH(type,name);
#else
            char const* name;
            int numVertices = getInfoMSH(type,&name);
#endif
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
                ghosts.reserve( numPartitions-1 );
                rank_type partition = (version < 2.2 && numTags > 2) ? data[3]-1 :
                    (version >= 2.2 && numTags > 3) ? data[4]-1 : theCurrentPartitionId;
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

                Feel::detail::GMSHElement gmshElt( num, type, std::vector<int>({physical}), elementary,
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

#if defined( FEELPP_HAS_GMSH_H )
    for ( auto const& it : __gt )
    {
#if defined( FEELPP_HAS_GMSH_API )
        std::string name;
        getInfoMSH( it.first,name );
#else
        const char* name;
        getInfoMSH( it.first, &name );
#endif
        LOG(INFO) << "Read " << it.second << " " << name << " elements\n";
    }
#endif
    if ( binary )
        __is.getline(__buf, 256);

    // make sure that we have read everything
    __is.getline(__buf, 256);
    CHECK( std::string( __buf ) == "$ENDELM" ||
           std::string( __buf ) == "$EndElements" )
        << "invalid end elements string " << __buf
        << " in gmsh importer. It should be either $ENDELM or $EndElements\n";

    // read periodic data if present
    __is.getline(__buf, 256);
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
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndPeriodic" )
            << "invalid end $Periodic string " << __buf
            << " in gmsh importer. It should be either $EndPeriodic\n";
    }
    // we are done reading the MSH file


    std::map<int,boost::tuple<int,rank_type> > mapGhostElt;
    std::vector<int> nbMsgToRecv( this->worldComm().localSize(),0 );

    node_type coords( mesh_type::nRealDim );
    tic();
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
                coords[j] = gmshpt.x[j]*M_scale;

            point_type pt( ptid, coords, gmshpt.onbdy );
            pt.setProcessIdInPartition( this->worldComm().localRank() );
            if ( gmshpt.parametric )
            {
                if ( gmshpt.gdim < 3 )
                {
                    pt.setParametricCoordinates( gmshpt.gdim, gmshpt.gtag, gmshpt.uv[0], gmshpt.uv[1] );
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
            __idGmshToFeel[it_gmshElt.num] = addPoint( mesh, it_gmshElt );
            break;
        }

        // Edges
        case GMSH_LINE:
        case GMSH_LINE_2:
        case GMSH_LINE_3:
        case GMSH_LINE_4:
        case GMSH_LINE_5:
        {
            __idGmshToFeel[it_gmshElt.num] = addEdge( mesh, it_gmshElt );
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
            __idGmshToFeel[it_gmshElt.num] = addFace( mesh, it_gmshElt );
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
            __idGmshToFeel[it_gmshElt.num] = addVolume( mesh, it_gmshElt );
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
    mesh->updateOrderedPoints();
    toc( "ImporterGmsh::readFromFile store elements in Mesh", FLAGS_v > 0 );
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
            auto & p1 = pit1->second;
            auto & p2 = pit2->second;
            p1.setMasterId( p2.id() );
            p1.setMasterVertex( boost::addressof( p2 ) );

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
    mpi::all_reduce( this->worldComm(), ne, gne,  std::plus<int>() );

    CHECK( gne > 0 ) << "The mesh does not have any elements.\n"
                     << "something was not right with GMSH mesh importation.\n"
                     << "please check that there are elements of topological dimension "
                     << mesh_type::nDim << "  in the mesh\n";
#endif

    if ( this->worldComm().localSize()>1 )
    {
        this->updateGhostCellInfo( mesh, __idGmshToFeel,  mapGhostElt, nbMsgToRecv );
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

}

template<typename MeshType>
template <typename gmsh_size_type,typename gmsh_size_partition_type,typename gmsh_size_periodiclink_type,typename gmsh_elttag_type>
void
ImporterGmsh<MeshType>::readFromFileVersion4( mesh_type* mesh, std::ifstream & __is, char __buf[],
                                              double version, bool binary, bool swap )
{
    rank_type procId = this->worldComm().localRank();
    rank_type worldSize = this->worldComm().localSize();

    //useful for binary read
    std::vector<int> _vectmpint;
    std::vector<gmsh_size_type> _vectmpsizet;
    std::vector<gmsh_elttag_type> _vectmpelttag;

    std::vector<std::map<int,std::set<int>>> entityTagToPhysicalMarkers(4);
    if ( std::string( __buf ) == "$Entities" )
    {
        VLOG(2) << "Reading $Entities ...";
        // eat  '\n' in binary mode otherwise the next binary read will get screwd
        if ( binary )
            __is.get();

        std::vector<gmsh_size_type> nGeoEntities(4); // points, curves, surfaces, volumes
        if ( !binary )
        {
            __is >> nGeoEntities[0] >> nGeoEntities[1] >> nGeoEntities[2] >> nGeoEntities[3];
        }
        else
        {
            __is.read( (char*)&nGeoEntities[0], 4*sizeof(gmsh_size_type) );
            CHECK( !swap ) << "TODO SWAP BINARY";
            //if(swap) SwapBytes((char*)&nGeoEntities[0], sizeof(size_type), 4);
            //__is.read( (char*)&x[0], 3*sizeof(double) );
            //if(swap) SwapBytes((char*)&x[0], sizeof(double), 3);
        }

        int entityTag, physicalTag;
        std::vector<double> coordGeoPoints(3);
        if ( version == 4 )
            coordGeoPoints.resize( 6 );
        gmsh_size_type numPhysicalTags = 0;
        for ( size_type k=0; k<nGeoEntities[0]; ++k )
        {
            if ( !binary )
            {
                 if ( version == 4 )
                     __is >> entityTag >> coordGeoPoints[0] >> coordGeoPoints[1] >> coordGeoPoints[2] >> coordGeoPoints[3] >> coordGeoPoints[4] >> coordGeoPoints[5] >> numPhysicalTags;
                 else
                     __is >> entityTag >> coordGeoPoints[0] >> coordGeoPoints[1] >> coordGeoPoints[2] >> numPhysicalTags;

                for ( int l=0;l<numPhysicalTags;++l )
                {
                    __is  >> physicalTag;
                    entityTagToPhysicalMarkers[0][entityTag].insert( physicalTag );
                }
            }
            else
            {
                __is.read( (char*)&entityTag, sizeof(int) );
                __is.read( (char*)&coordGeoPoints[0], coordGeoPoints.size()*sizeof(double) );
                __is.read( (char*)&numPhysicalTags, sizeof(gmsh_size_type) );
                _vectmpint.resize( numPhysicalTags );
                __is.read( (char*)&_vectmpint[0], numPhysicalTags*sizeof(int) );
                entityTagToPhysicalMarkers[0][entityTag].insert( _vectmpint.begin(), _vectmpint.end() );
            }
        }

        std::vector<double> coordMinMax(6);
        gmsh_size_type nBoundingEntity;
        for (int e = 1;e<4;++e)
        {
            for ( size_type k=0; k<nGeoEntities[e]; ++k )
            {
                if ( !binary )
                {
                    __is >> entityTag
                         >> coordMinMax[0] >> coordMinMax[1] >> coordMinMax[2] // min
                         >> coordMinMax[3] >> coordMinMax[4] >> coordMinMax[5] //max
                         >> numPhysicalTags;
                    for ( int l=0;l<numPhysicalTags;++l )
                    {
                        __is  >> physicalTag;
                        entityTagToPhysicalMarkers[e][entityTag].insert( physicalTag );
                    }
                    __is >> nBoundingEntity;
                    for ( int l=0;l<nBoundingEntity;++l )
                        __is  >> physicalTag;
                }
                else
                {
                    __is.read( (char*)&entityTag, sizeof(int) );
                    __is.read( (char*)&coordMinMax[0], 6*sizeof(double) );
                    __is.read( (char*)&numPhysicalTags, sizeof(gmsh_size_type) );
                    _vectmpint.resize( numPhysicalTags );
                    __is.read( (char*)&_vectmpint[0], numPhysicalTags*sizeof(int) );
                    entityTagToPhysicalMarkers[e][entityTag].insert( _vectmpint.begin(), _vectmpint.end() );
                    __is.read( (char*)&nBoundingEntity, sizeof(gmsh_size_type) );
                    _vectmpint.resize( nBoundingEntity );
                    __is.read( (char*)&_vectmpint[0], nBoundingEntity*sizeof(int) );
                }
            }
        }

        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndEntities" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndEntities\n";
        VLOG(2) << "Reading $Entities done";

        __is.getline(__buf, 256);
    }

    gmsh_size_partition_type numPartitions = 1;
    std::vector< std::map<int,std::set<int>>> entityTagInCurrentPartitionToPartitions;

    if ( std::string( __buf ) == "$PartitionedEntities" )
    {
        VLOG(2) << "Reading $PartitionedEntities ...";
        // eat  '\n' in binary mode otherwise the next binary read will get screwd
        if ( binary )
            __is.get();

        gmsh_size_partition_type numGhostEntities;
        if ( !binary )
            __is >> numPartitions >> numGhostEntities;
        else
        {
            __is.read( (char*)&numPartitions, sizeof(gmsh_size_partition_type) );
            __is.read( (char*)&numGhostEntities, sizeof(gmsh_size_partition_type) );
        }

        // ghost partitions tag
        std::set<int> allPartitionTags;
        if ( !binary )
        {
            int ghostEntityTag, partition;
            for ( size_type k=0; k<numGhostEntities; ++k )
            {
                __is >> ghostEntityTag >> partition;
                allPartitionTags.insert( partition );
            }
        }
        else
        {
            _vectmpint.resize( 2*numGhostEntities );
            __is.read( (char*)&_vectmpint[0], 2*numGhostEntities*sizeof(int) );
            for ( int k=0;k<numGhostEntities;++k )
                allPartitionTags.insert( _vectmpint[2*k+1] );
        }
        if ( worldSize > 1 )
            CHECK( allPartitionTags.size() == numPartitions ) << "the partitioning in Gmsh should be run with Create Ghost Cells options";

        std::vector<gmsh_size_type> nGeoEntities(4); // points, curves, surfaces, volumes
        if ( !binary )
            __is >> nGeoEntities[0] >> nGeoEntities[1] >> nGeoEntities[2] >> nGeoEntities[3];
        else
        {
            __is.read( (char*)&nGeoEntities[0], 4*sizeof(gmsh_size_type) );
        }

        if ( nGeoEntities[0] > 0 || nGeoEntities[1] > 0 || nGeoEntities[2] > 0 || nGeoEntities[3] > 0 )
            entityTagInCurrentPartitionToPartitions.resize( 4 );

        std::vector<double> coordGeoPoints(3);
        int entityTag, parentDim, parentTag, partitionTag, physicalTag, boundingEntityTag;
        gmsh_size_partition_type numPartitionsForThisEntity;
        gmsh_size_type numPhysicalTags;
        std::set<int> partitionTags;
        for ( size_type k=0; k<nGeoEntities[0]; ++k )
        {
            if ( !binary )
            {
                __is >> entityTag >> parentDim >> parentTag >> numPartitionsForThisEntity;
            }
            else
            {
                __is.read( (char*)&entityTag, sizeof(int) );
                __is.read( (char*)&parentDim, sizeof(int) );
                __is.read( (char*)&parentTag, sizeof(int) );
                __is.read( (char*)&numPartitionsForThisEntity, sizeof(gmsh_size_partition_type) );
            }
            partitionTags.clear();
            bool isOnCurrentPartition = false;
            for ( size_type p=0;p<numPartitionsForThisEntity;++p )
            {
                if ( !binary )
                    __is >> partitionTag;
                else
                {
                    __is.read( (char*)&partitionTag, sizeof(int) );
                }

                partitionTags.insert( partitionTag );
                if ( procId == Feel::detail::gmshPartitionTagToProcessId( partitionTag, worldSize ) )
                    isOnCurrentPartition = true;
            }
            if ( isOnCurrentPartition )
                entityTagInCurrentPartitionToPartitions[0][entityTag] = partitionTags;
            else if ( parentDim == 0 )
            {
                // WARNING : partitioning information with points entities is quite strange for some node
                // Consequence some extra node are inserted in the mesh, TODO remove maybe this node
                entityTagInCurrentPartitionToPartitions[0][entityTag] = allPartitionTags;
            }

            // here we dont add this physical markers (generated by Gmsh with the partitioning)
            if ( !binary )
            {
                __is >> coordGeoPoints[0] >> coordGeoPoints[1] >> coordGeoPoints[2] >> numPhysicalTags;
                for ( size_type t=0;t<numPhysicalTags;++t )
                    __is >> physicalTag; // not use this physicalTag (generated by gmsh)
            }
            else
            {
                __is.read( (char*)&coordGeoPoints[0], 3*sizeof(double) );
                __is.read( (char*)&numPhysicalTags, sizeof(gmsh_size_type) );
                _vectmpint.resize( numPhysicalTags );
                __is.read( (char*)&_vectmpint[0], numPhysicalTags*sizeof(int) );
            }

            // add physical markers from parent entity (with same dim)
            if ( parentDim == 0 )
            {
                auto itFindPhysicalMarkersFromParent = entityTagToPhysicalMarkers[parentDim].find( parentTag );
                if ( itFindPhysicalMarkersFromParent != entityTagToPhysicalMarkers[parentDim].end() )
                    entityTagToPhysicalMarkers[0][entityTag] = itFindPhysicalMarkersFromParent->second;
            }
        }

        std::vector<double> coordMinMax(6);
        gmsh_size_type numBoundingSubentity;
        for ( int e=1;e<4;++e )
        {
            for ( size_type k=0; k<nGeoEntities[e]; ++k )
            {
                if ( !binary )
                {
                    __is >> entityTag >> parentDim >> parentTag >> numPartitionsForThisEntity;
                }
                else
                {
                    __is.read( (char*)&entityTag, sizeof(int) );
                    __is.read( (char*)&parentDim, sizeof(int) );
                    __is.read( (char*)&parentTag, sizeof(int) );
                    __is.read( (char*)&numPartitionsForThisEntity, sizeof(gmsh_size_partition_type) );
                }
                partitionTags.clear();
                bool isOnCurrentPartition = false;
                for ( size_type p=0;p<numPartitionsForThisEntity;++p )
                {
                    if ( !binary )
                        __is >> partitionTag;
                    else
                    {
                        __is.read( (char*)&partitionTag, sizeof(int) );
                    }

                    partitionTags.insert( partitionTag );
                    if ( procId == Feel::detail::gmshPartitionTagToProcessId( partitionTag, worldSize ) )
                        isOnCurrentPartition = true;
                }
                if ( isOnCurrentPartition )
                    entityTagInCurrentPartitionToPartitions[e][entityTag] = partitionTags;

                if ( !binary )
                {
                    __is >> coordMinMax[0] >> coordMinMax[1] >> coordMinMax[2] // min
                         >> coordMinMax[3] >> coordMinMax[4] >> coordMinMax[5] //max
                         >> numPhysicalTags;
                    for ( size_type t=0;t<numPhysicalTags;++t )
                        __is >> physicalTag; // not use this physicalTag (generated by gmsh)

                    __is >> numBoundingSubentity;
                    for ( size_type t=0;t<numBoundingSubentity;++t )
                        __is >> boundingEntityTag;
                }
                else
                {
                    __is.read( (char*)&coordMinMax[0], 6*sizeof(double) );
                    __is.read( (char*)&numPhysicalTags, sizeof(gmsh_size_type) );
                    _vectmpint.resize( numPhysicalTags );
                    __is.read( (char*)&_vectmpint[0], numPhysicalTags*sizeof(int) );
                    __is.read( (char*)&numBoundingSubentity, sizeof(gmsh_size_type) );
                    _vectmpint.resize( numBoundingSubentity );
                    __is.read( (char*)&_vectmpint[0], numBoundingSubentity*sizeof(int) );
                }

                // add physical markers from parent entity (with same dim)
                if ( parentDim == e )
                {
                    auto itFindPhysicalMarkersFromParent = entityTagToPhysicalMarkers[parentDim].find( parentTag );
                    if ( itFindPhysicalMarkersFromParent != entityTagToPhysicalMarkers[parentDim].end() )
                        entityTagToPhysicalMarkers[e][entityTag] = itFindPhysicalMarkersFromParent->second;
                }

            }
        }

        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndPartitionedEntities" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndPartitionedEntities\n";
        VLOG(2) << "Reading $PartitionedEntities done";

        __is.getline(__buf, 256);
    }

    //std::map<size_type,Feel::detail::GMSHPoint> additionalNodes;
    std::map<rank_type,std::set<size_type>> interprocessNodes; // neighbour partId -> ( points id )

    if ( std::string( __buf ) == "$Nodes" )
    {
        VLOG(2) << "Reading $Nodes ...";
        // eat  '\n' in binary mode otherwise the next binary read will get screwd
        if ( binary )
            __is.get();

        gmsh_size_type nEntityBlocks, nNodes, minNodeTag, maxNodeTag;
        if ( !binary )
        {
            __is >> nEntityBlocks >> nNodes;
            if ( version >= 4.1 )
                __is >> minNodeTag >> maxNodeTag;
        }
        else
        {
            __is.read( (char*)&nEntityBlocks, sizeof(gmsh_size_type) );
            __is.read( (char*)&nNodes, sizeof(gmsh_size_type) );
            if ( version >= 4.1 )
            {
                __is.read( (char*)&minNodeTag, sizeof(gmsh_size_type) );
                __is.read( (char*)&maxNodeTag, sizeof(gmsh_size_type) );
            }
        }

        int entityDim, entityTag, parametric;
        gmsh_size_type numNodesInBlock;
        Eigen::Vector3d x;
        for ( size_type k=0;k<nEntityBlocks;++k )
        {
            if ( version >= 4.1 )
            {
                if ( !binary )
                    __is >> entityDim >> entityTag >> parametric >> numNodesInBlock;
                else
                {
                    __is.read( (char*)&entityDim, sizeof(int) );
                    __is.read( (char*)&entityTag, sizeof(int) );
                    __is.read( (char*)&parametric, sizeof(int) );
                    __is.read( (char*)&numNodesInBlock, sizeof(gmsh_size_type) );
                }
            }
            else // version == 4.0
            {
                if ( !binary )
                    __is >> entityTag >> entityDim  >> parametric >> numNodesInBlock;
                else
                {
                    __is.read( (char*)&entityTag, sizeof(int) );
                    __is.read( (char*)&entityDim, sizeof(int) );
                    __is.read( (char*)&parametric, sizeof(int) );
                    __is.read( (char*)&numNodesInBlock, sizeof(gmsh_size_type) );
                }

            }

            if ( parametric == 1 )
                CHECK( false ) << "TODO parametric nodes";

            bool useThisNode = true;
            std::set<rank_type> connectedProcessIds;

            if ( !entityTagInCurrentPartitionToPartitions.empty() )
            {
                auto itFindEntityTag = entityTagInCurrentPartitionToPartitions[entityDim].find( entityTag );
                useThisNode = itFindEntityTag != entityTagInCurrentPartitionToPartitions[entityDim].end();
                if ( useThisNode )
                {
                    for ( int onPartitionTag : itFindEntityTag->second )
                    {
                        rank_type otherProcId = Feel::detail::gmshPartitionTagToProcessId( onPartitionTag, worldSize );
                        if ( procId != otherProcId )
                            connectedProcessIds.insert( otherProcId );
                    }
                }
            }

            std::vector<gmsh_elttag_type> nodeIds;
            if ( version >= 4.1 )
            {
                nodeIds.resize( numNodesInBlock );
                if ( !binary )
                {
                    for ( size_type p=0;p<numNodesInBlock;++p )
                        __is >> nodeIds[p];
                }
                else
                {
                    __is.read( (char*)&nodeIds[0], numNodesInBlock*sizeof(gmsh_elttag_type) );
                }
            }
            for ( size_type p=0;p<numNodesInBlock;++p )
            {
                size_type id = invalid_size_type_value;
                if ( version >= 4.1 )
                    id = nodeIds[p];
                else
                {
                    if ( !binary )
                        __is >> id;
                    else
                    {
                        int idInt;
                         __is.read( (char*)&idInt, sizeof(int) );
                         id = idInt;
                    }
                }
                if ( !binary )
                {
                    __is >> x[0] >> x[1] >> x[2];
                    // TODO add parametric node if require
                }
                else
                {
                    __is.read( (char*)&x[0], 3*sizeof(double) );
                    // TODO add parametric node if require
                }

                if ( !useThisNode )
                    continue;

                for ( rank_type cpid : connectedProcessIds )
                    interprocessNodes[cpid].insert( id );


                node_type coords( mesh_type::nRealDim );
                for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
                    coords[j] = x[j]*M_scale;

                point_type pt( id, coords /*, gmshpt.onbdy*/ );
                pt.setProcessIdInPartition( procId/*this->worldComm().localRank()*/ );
                pt.setProcessId( procId );
                if ( parametric == 1 )
                {
#if 0
                    pt.setGDim( gmshpt.gdim );
                    pt.setGTag( gmshpt.gtag );
                    if ( gmshpt.gdim < 3 )
                    {
                        pt.setParametricCoordinates( gmshpt.uv[0], gmshpt.uv[1] );
                        mesh->setParametric( true );
                    }
#endif
                }
                mesh->addPoint( pt );
            }
        }
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndNodes" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndNodes\n";
        VLOG(2) << "Reading $Nodes done";

        __is.getline(__buf, 256);
    }


    std::map<int,int> __idGmshToFeel; // id Gmsh to id Feel

    std::map<std::vector<int>,std::tuple<Feel::detail::GMSHElement,bool>> elementsInWaiting;

    if ( std::string( __buf ) == "$Elements" )
    {
        VLOG(2) << "Reading $Elements ...";
        // eat  '\n' in binary mode otherwise the next binary read will get screwd
        if ( binary )
            __is.get();

        gmsh_size_type numEntityBlocks, numElements, minElementTag, maxElementTag;
        if ( !binary )
        {
            __is >> numEntityBlocks >> numElements;
             if ( version >= 4.1 )
                 __is >> minElementTag >> maxElementTag;
        }
        else
        {
            __is.read( (char*)&numEntityBlocks, sizeof(gmsh_size_type) );
            __is.read( (char*)&numElements, sizeof(gmsh_size_type) );
            if ( version >= 4.1 )
            {
                __is.read( (char*)&minElementTag, sizeof(gmsh_size_type) );
                __is.read( (char*)&maxElementTag, sizeof(gmsh_size_type) );
            }
        }


        int entityDim = 0, entityTag = 0, elementType = 0;//, physicalTag = 0;
        std::vector<int> physicalTag;
        gmsh_size_type numElementsInBlock = 0;
        gmsh_elttag_type elementTag = 0;

        // partitioning data
        int numPartitionsElt = 1, partition = procId;
        int parent =0 , dom1 = 0, dom2 = 0;
        std::vector<rank_type> ghosts;

        int numVertices = 0;
        std::vector<int> indices;
        Feel::detail::GMSHElement it_gmshElt( elementTag, elementType, physicalTag, entityTag,    //   num, type, physical, elementary,
                                              numPartitionsElt, partition, ghosts,
                                              parent, dom1, dom2,
                                              numVertices, indices,
                                              this->worldComm().localRank(),
                                              this->worldComm().localSize(),
                                              M_respect_partition );

        for ( size_type k=0;k<numEntityBlocks;++k )
        {
            if ( version >= 4.1 )
            {
                if ( !binary )
                    __is >> entityDim >> entityTag >> elementType >> numElementsInBlock;
                else
                {
                    __is.read( (char*)&entityDim, sizeof(int) );
                    __is.read( (char*)&entityTag, sizeof(int) );
                    __is.read( (char*)&elementType, sizeof(int) );
                    __is.read( (char*)&numElementsInBlock, sizeof(gmsh_size_type) );
                }
            }
            else // version == 4.0
            {
                if ( !binary )
                    __is >> entityTag >> entityDim >> elementType >> numElementsInBlock;
                else
                {
                    __is.read( (char*)&entityTag, sizeof(int) );
                    __is.read( (char*)&entityDim, sizeof(int) );
                    __is.read( (char*)&elementType, sizeof(int) );
                    __is.read( (char*)&numElementsInBlock, sizeof(gmsh_size_type) );
                }
            }

            bool useThisEntity = true;
            if ( !entityTagInCurrentPartitionToPartitions.empty() )
            {
                auto itFindEntityTag = entityTagInCurrentPartitionToPartitions[entityDim].find( entityTag );
                useThisEntity = itFindEntityTag != entityTagInCurrentPartitionToPartitions[entityDim].end();
            }

            if ( true/*useThisEntity*/ )
            {
                // get physical tag
                physicalTag.clear();
                auto itFindEntityTagToPhysicalMarkers = entityTagToPhysicalMarkers[entityDim].find( entityTag );
                if ( itFindEntityTagToPhysicalMarkers != entityTagToPhysicalMarkers[entityDim].end() )
                {
                    for ( int thePhysicalTag : itFindEntityTagToPhysicalMarkers->second )
                        physicalTag.push_back( thePhysicalTag );
                }
            }

#if defined( FEELPP_HAS_GMSH_API )
            std::string ename;
            numVertices = getInfoMSH( elementType,ename );
#else
            const char* ename;
            numVertices = getInfoMSH( elementType,&ename );
#endif
            CHECK(numVertices!=0) << "Unknown number of vertices for element type " << elementType << "\n";

            // update current gmsh element with read data
            it_gmshElt.type = elementType;
            it_gmshElt.physical = physicalTag;
            it_gmshElt.elementary = entityTag;
            it_gmshElt.numVertices = numVertices;

            indices.resize( numVertices );
            _vectmpelttag.resize( numVertices );

            for ( size_type p=0;p<numElementsInBlock;++p )
            {
                if ( !binary )
                {
                    __is >> elementTag;
                    for(int j = 0; j < numVertices; j++)
                        __is >> indices[j];
                }
                else
                {
                    __is.read( (char*)&elementTag, sizeof(gmsh_elttag_type) );
                    __is.read( (char*)&_vectmpelttag[0], numVertices*sizeof(gmsh_elttag_type) );
                    indices.assign( _vectmpelttag.begin(),_vectmpelttag.end() );
                }

                // update current gmsh element with read data
                it_gmshElt.num = elementTag;

                // WARNING : need to apply a special treatment on mesh partitioned with elements of dim <  mesh_type::nDim
                // the partitioning information with this entity is not very clear (some duplications with global entities)
                if ( numPartitions > 1 && entityDim < mesh_type::nDim )
                {
                    it_gmshElt.physical = physicalTag; // here because can be modified just after

                    std::vector<int> indicesSorted = indices;
                    std::sort( indicesSorted.begin(), indicesSorted.end() );
                    auto itFindElementsInWaiting = elementsInWaiting.find( indicesSorted );
                    if ( itFindElementsInWaiting != elementsInWaiting.end() )
                    {
                        Feel::detail::GMSHElement & gmshEltToUpdate = std::get<0>( itFindElementsInWaiting->second );
                        // get new physical markers
                        bool hasNewPhysicalMarkers = false;
                        for ( int ptag : it_gmshElt.physical )
                        {
                            if ( std::find(gmshEltToUpdate.physical.begin(), gmshEltToUpdate.physical.end(), ptag) == gmshEltToUpdate.physical.end() )
                            {
                                gmshEltToUpdate.physical.push_back( ptag );
                                hasNewPhysicalMarkers = true;
                            }
                        }

                        if ( std::get<1>( itFindElementsInWaiting->second ) ) // already inserted
                        {
                            // update physical marker in Feel++ mesh
                            if ( hasNewPhysicalMarkers )
                            {
                                auto itFindGmhsId = __idGmshToFeel.find( gmshEltToUpdate.num );
                                CHECK( itFindGmhsId != __idGmshToFeel.end() ) << "Gmsh element id not found";
                                //CHECK( gmshEltToUpdate.physical.size() == 1 ) "support only one marker";

                                if ( entityDim == ( mesh_type::nDim -1 ) )
                                {
                                    auto faceIt = mesh->faceIterator( itFindGmhsId->second );
                                    CHECK( faceIt != mesh->endFace() ) << "face not found";
                                    faceIt->second.addMarker( gmshEltToUpdate.physical );
                                }
                                else if ( mesh_type::nDim == 3 && entityDim == ( mesh_type::nDim - 2 ) )
                                {
                                    if constexpr ( mesh_type::nDim == 3 )
                                    {
                                        auto edgeIt = mesh->edgeIterator( itFindGmhsId->second );
                                        CHECK( edgeIt != mesh->endEdge() ) << "face not found";
                                        edgeIt->second.addMarker( gmshEltToUpdate.physical );
                                    }
                                }
                                else if ( entityDim == 0 )
                                {
                                    auto pointIt = mesh->pointIterator( itFindGmhsId->second );
                                    CHECK( pointIt != mesh->endPoint() ) << "point not found";
                                    pointIt->second.addMarker( gmshEltToUpdate.physical );
                                }
                            }
                            useThisEntity = false;
                        }
                        else
                        {
                            std::get<1>( itFindElementsInWaiting->second ) = useThisEntity;
                            if ( useThisEntity )
                                it_gmshElt.physical = gmshEltToUpdate.physical;
                        }
                    }
                    else
                    {
                        elementsInWaiting[indicesSorted] = std::make_tuple(it_gmshElt, useThisEntity );
                    }
                }

                if ( useThisEntity )
                {
                    // update current gmsh element with read data
                    it_gmshElt.indices = indices;

                    switch ( elementType )
                    {
                        // Points
                    case GMSH_POINT:
                    {
                        __idGmshToFeel[it_gmshElt.num] = addPoint( mesh, it_gmshElt );
                        break;
                    }

                    // Edges
                    case GMSH_LINE:
                    case GMSH_LINE_2:
                    case GMSH_LINE_3:
                    case GMSH_LINE_4:
                    case GMSH_LINE_5:
                    {
                        __idGmshToFeel[it_gmshElt.num] = addEdge( mesh, it_gmshElt );
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
                        __idGmshToFeel[it_gmshElt.num] = addFace( mesh, it_gmshElt );
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
                        __idGmshToFeel[it_gmshElt.num] = addVolume( mesh, it_gmshElt );
                        break;
                    }

                    default:
                        break;
                    }


                }

            }
        }
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndElements" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndElements";
        VLOG(2) << "Reading $Elements done";
        __is.getline(__buf, 256);
    }

    std::vector<PeriodicEntity> periodic_entities;
    if ( std::string( __buf ) == "$Periodic" )
    {
        VLOG(2) << "Reading $Periodic ...";
        // eat  '\n' in binary mode otherwise the next binary read will get screwd
        if ( binary )
            __is.get();

        gmsh_size_periodiclink_type numPeriodicLinks;
        int entityDim, entityTag, entityTagMaster;
        if ( !binary )
            __is >> numPeriodicLinks;
        else
            __is.read( (char*)&numPeriodicLinks, sizeof(gmsh_size_periodiclink_type) );
        for(size_type i = 0; i < numPeriodicLinks; i++)
        {
            if ( !binary )
                __is >> entityDim >> entityTag >> entityTagMaster;
            else
            {
                __is.read( (char*)&entityDim, sizeof(int) );
                __is.read( (char*)&entityTag, sizeof(int) );
                __is.read( (char*)&entityTagMaster, sizeof(int) );
            }
            PeriodicEntity e( entityDim, entityTag, entityTagMaster );
            if ( version >= 4.1 )
            {
                gmsh_size_type numAffine;
                if ( !binary )
                    __is >> numAffine;
                else
                    __is.read( (char*)&numAffine, sizeof(gmsh_size_type) );
                std::vector<double> valuesAffineTransfo( numAffine );
                if ( !binary )
                {
                    for (size_type k=0;k<numAffine;++k )
                        __is >> valuesAffineTransfo[k];
                }
                else
                {
                    if ( numAffine > 0 )
                        __is.read( (char*)&valuesAffineTransfo[0], numAffine*sizeof(double) );
                }
            }
            gmsh_size_type numCorrespondingNodes;
            if ( !binary )
                __is >> numCorrespondingNodes;
            else
                __is.read( (char*)&numCorrespondingNodes, sizeof(gmsh_size_type) );
            if ( !binary )
            {
                for(size_type j = 0; j < numCorrespondingNodes; j++)
                {
                    gmsh_elttag_type v1,v2;
                    __is >> v1 >> v2;
                    e.correspondingVertices[v1] = v2;
                }
            }
            else
            {
                _vectmpelttag.resize( 2*numCorrespondingNodes );
                __is.read( (char*)&_vectmpelttag[0], 2*numCorrespondingNodes*sizeof(gmsh_elttag_type) );
                for (int k=0;k<numCorrespondingNodes;++k)
                    e.correspondingVertices[_vectmpelttag[k]] = _vectmpelttag[k+1];
            }
            CHECK( e.correspondingVertices.size() == numCorrespondingNodes ) << "Invalid number of vertices in periodic entity"
                                                                             << " dim: " << e.dim
                                                                             << " slave: " << e.slave
                                                                             << " master: " << e.master
                                                                             << " got: " << e.correspondingVertices.size()
                                                                             << " expected : " << numCorrespondingNodes << "\n";
            periodic_entities.push_back( e );
        }
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndPeriodic" )
            << "invalid end $Periodic string " << __buf
            << " in gmsh importer. It should be either $EndPeriodic\n";
        VLOG(2) << "Reading $Periodic done";
        __is.getline(__buf, 256);
    }


    std::map<rank_type,std::set<size_type>> ghostElementToSendToProcessId;
    std::set<rank_type> ghostElementToRecvFromProcessId;

    if ( std::string( __buf ) == "$GhostElements" )
    {
        VLOG(2) << "Reading $GhostElements ...";
        // eat  '\n' in binary mode otherwise the next binary read will get screwd
        if ( binary )
            __is.get();

        gmsh_size_partition_type numGhostElements, numGhostPartitions;
        gmsh_elttag_type elementTag;
        int partitionTag;
        std::vector<int> ghostPartitionTags;
        if ( !binary )
            __is >> numGhostElements;
        else
        {
            __is.read( (char*)&numGhostElements, sizeof(gmsh_size_partition_type) );
        }
        for ( size_type g=0;g<numGhostElements;++g )
        {
            if ( !binary )
            {
                __is >> elementTag >> partitionTag >> numGhostPartitions;
                ghostPartitionTags.resize( numGhostPartitions );
                for ( size_type p=0;p<numGhostPartitions;++p )
                    __is >> ghostPartitionTags[p];
            }
            else
            {
                __is.read( (char*)&elementTag, sizeof(gmsh_elttag_type) );
                __is.read( (char*)&partitionTag, sizeof(int) );
                __is.read( (char*)&numGhostPartitions, sizeof(gmsh_size_partition_type) );
                ghostPartitionTags.resize( numGhostPartitions );
                __is.read( (char*)&ghostPartitionTags[0], numGhostPartitions*sizeof(int) );
            }

            if ( procId == Feel::detail::gmshPartitionTagToProcessId( partitionTag, worldSize ) )
            {
                for ( int k=0;k<numGhostPartitions;++k )
                {
                    rank_type gproc = Feel::detail::gmshPartitionTagToProcessId( ghostPartitionTags[k], worldSize );
                    ghostElementToSendToProcessId[gproc].insert( elementTag );
                }
            }
            else
            {
                auto itHasGhostElement = std::find_if( ghostPartitionTags.begin(), ghostPartitionTags.end(),
                                                       [&procId,&worldSize](int const& ghostPartitionTag )
                                                           {
                                                               return procId == Feel::detail::gmshPartitionTagToProcessId( ghostPartitionTag, worldSize );
                                                           });

                if ( itHasGhostElement != ghostPartitionTags.end() )
                    ghostElementToRecvFromProcessId.insert( Feel::detail::gmshPartitionTagToProcessId( partitionTag, worldSize ) );
            }
        }
        __is.getline(__buf, 256);
        CHECK( std::string( __buf ) == "$EndGhostElements" )
            << "invalid file format entry "
            << __buf
            << " instead of $EndGhostElements";
        VLOG(2) << "Reading $GhostElements done";
        __is.getline(__buf, 256);
    }

    if ( numPartitions > 1 )
    {
        int nbRequest = ghostElementToSendToProcessId.size() + ghostElementToRecvFromProcessId.size();
        if ( nbRequest > 0 )
        {
            mpi::request * reqs = new mpi::request[nbRequest];
            int cptRequest=0;
            // prepare and send ghost information
            for ( auto const& [ gproc, ghostElts ] : ghostElementToSendToProcessId )
            {
                std::vector<boost::tuple<size_type, std::vector<double> > > dataPointsToSend;
                std::vector<std::vector<size_type>> dataGhostElementToSend;// eltId, ptId1, ptId2,...
                std::set<size_type> pointsRegistered;
                for ( size_type ghostEltGmshId : ghostElts )
                {
                    auto itFindGmshId = __idGmshToFeel.find( ghostEltGmshId );
                    CHECK( itFindGmshId != __idGmshToFeel.end() ) << "no element gmsh with id  " << ghostEltGmshId;
                    size_type ghostEltFeelId = itFindGmshId->second;
                    auto const& elt = mesh->element( ghostEltFeelId );
                    std::vector<size_type> theindices( npoints_per_element+1 );
                    theindices[0] = ghostEltFeelId;
                    for ( int p=0;p<npoints_per_element;++ p )
                    {
                        auto const& pt = elt.point( p );
                        size_type ptId = pt.id();
                        theindices[ p+1 ] = ptId;
                        if ( pointsRegistered.find( ptId ) != pointsRegistered.end() )
                            continue;
                        auto itFindIpN = interprocessNodes.find( gproc );
                        if ( itFindIpN == interprocessNodes.end() ||
                             itFindIpN->second.find ( ptId ) == itFindIpN->second.end() )
                        {
                            std::vector<double> thecoord( mesh_type::nRealDim );
                            for ( uint8_type d=0;d<mesh_type::nRealDim;++d )
                                thecoord[d] = pt.node()[d];
                            dataPointsToSend.push_back( boost::make_tuple( ptId,thecoord ) );
                            pointsRegistered.insert( ptId );
                        }
                    }
                    dataGhostElementToSend.push_back( theindices );
                }
                auto fullDataToSend = boost::make_tuple( dataPointsToSend,dataGhostElementToSend );
                reqs[cptRequest++] = this->worldComm().localComm().isend( gproc , 0, fullDataToSend );
            }

            // recv ghost information
            std::map<rank_type, boost::tuple< std::vector<boost::tuple<size_type, std::vector<double> > >,
                                              std::vector<std::vector<size_type>> > > dataToRecv;
            for ( rank_type aproc : ghostElementToRecvFromProcessId )
            {
                reqs[cptRequest++] = this->worldComm().localComm().irecv( aproc , 0, dataToRecv[aproc] );
            }
            // wait all requests
            mpi::wait_all(reqs, reqs + nbRequest);
            // delete reqs because finish comm
            delete [] reqs;

            // create ghost elements
            for ( auto const& [ aproc, dataRecvByProc ] : dataToRecv )
            {
                auto const& dataRecvPoints = boost::get<0>( dataRecvByProc );
                auto const& dataRecvGhostElements = boost::get<1>( dataRecvByProc );
                for ( auto const& dataRecvPoint : dataRecvPoints )
                {
                    size_type ptId = boost::get<0>( dataRecvPoint );
                    auto const& coordsRecv = boost::get<1> (dataRecvPoint );
                    node_type coords( mesh_type::nRealDim );
                    for ( uint16_type j = 0; j < mesh_type::nRealDim; ++j )
                        coords[j] = coordsRecv[j];
                    point_type pt( ptId, coords );
                    pt.setProcessIdInPartition( procId/*this->worldComm().localRank()*/ );
                    mesh->addPoint( pt );
                }
                for ( auto const& dataRecvGhostElement : dataRecvGhostElements )
                {
                    element_type e;
                    e.setId( mesh->elements().size() );
                    e.setProcessIdInPartition( this->worldComm().localRank() );
                    //e.setMarker( __e.physical );
                    //e.setMarker2( __e.elementary );
                    e.setProcessId( aproc );
                    //e.setNeighborPartitionIds( __e.ghosts );
                    size_type idOtherPartition = dataRecvGhostElement[0];
                    e.setIdInOtherPartitions( aproc, idOtherPartition );

                    for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
                    {
                        //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
                        CHECK( mesh->hasPoint( dataRecvGhostElement[jj+1] ) ) << "point not found with id " << dataRecvGhostElement[jj+1];
                        point_type & pt = mesh->pointIterator( dataRecvGhostElement[jj+1] )->second;
                        e.setPoint(  jj, pt );
                    }
                    mesh->addElement( /*mesh->endElement(), */std::move(e) );
                }
            }
        }
    } // if ( numPartitions > 1 )

    // update ordered points in mesh data structure
    mesh->updateOrderedPoints();
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addPoint( mesh_type* mesh, Feel::detail::GMSHElement const& __e )
{
    auto pit = mesh->pointIterator(__e.indices[0]);
    CHECK( pit != mesh->endPoint() ) <<  "point not register in mesh";
    auto & pt = pit->second;

    if ( !__e.physical.empty() )
        pt.addMarker( __e.physical );
    if ( false )
        pt.setMarker2( __e.elementary );
    //pt.setProcessId( __e.partition );
    pt.setNeighborPartitionIds( __e.ghosts );

    DVLOG(2) << "update point with id :" << pt.id() << " and marker " << pt.marker()
             << " n1: " << pt.node() << "\n";
    return pt.id();
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e )
{
    return addEdge( mesh, __e, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
int
ImporterGmsh<MeshType>::addEdge( mesh_type*mesh, Feel::detail::GMSHElement const& __e, mpl::int_<1> )
{
    element_type e;
    e.setId( ( false )? __e.num : mesh->elements().size() );
    //e.setWorldComm(this->worldComm());
    e.setProcessIdInPartition( this->worldComm().localRank() );
    if ( !__e.physical.empty() )
        e.addMarker( __e.physical );
    if ( false )
        e.setMarker2( __e.elementary );
    e.setProcessId( __e.partition );
    e.setNeighborPartitionIds( __e.ghosts );

    if ( __e.type == GMSH_LINE ||
         __e.type == GMSH_LINE_2 ||
         __e.type == GMSH_LINE_3 ||
         __e.type == GMSH_LINE_4 ||
         __e.type == GMSH_LINE_5 )
    {
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            point_type & pt = mesh->pointIterator( __e.indices[jj] )->second;
            if (!e.isGhostCell())
                pt.setProcessId( e.processId() );
            e.setPoint( jj, pt );
        }
    }
    DVLOG(2) << "added edge with id :" << e.id()
             << " n1: " << mesh->point( __e.indices[0] ).node()
             << " n2: " << mesh->point( __e.indices[1] ).node() << "\n";

    auto [eit,inserted] = mesh->addElement( /*mesh->endElement(), */std::move(e) );
    auto const& [eid,eltInserted] = *eit;
    return eid;
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addEdge( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<2> )
{
    face_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setId( mesh->numFaces() );
    if ( !__e.physical.empty() )
        e.addMarker( __e.physical );
    if ( false )
        e.setMarker2( __e.elementary );
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
        for ( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
        {
            DCHECK( mesh->hasPoint( __e.indices[jj] ) ) << "points not registered in mesh";
            point_type & pt = mesh->pointIterator( __e.indices[jj] )->second;
            e.setPoint( jj, pt );
        }
    }

    auto [fit,inserted] = mesh->addFace( std::move(e) );
    auto const& [fid,faceInserted] = *fit;
    DVLOG(2) << "added edge with id :" << faceInserted.id()
             << " n1: " << mesh->point( __e.indices[0] ).node()
             << " n2: " << mesh->point( __e.indices[1] ).node() << "\n";
    return fid;
}
template<typename MeshType>
int
ImporterGmsh<MeshType>::addEdge( mesh_type*mesh, Feel::detail::GMSHElement const& __e, mpl::int_<3> )
{
    edge_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setId( mesh->numEdges() );
    if ( !__e.physical.empty() )
        e.addMarker( __e.physical );
    if ( false )
        e.setMarker2( __e.elementary );
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
        for ( uint16_type jj = 0; jj < npoints_per_edge; ++jj )
        {
            point_type & pt = mesh->pointIterator( __e.indices[jj] )->second;
            e.setPoint( jj, pt );
        }
    }

    auto [edit,inserted] = mesh->addEdge( std::move(e) );
    auto const& [edid,edgeInserted] = *edit;
    if ( npoints_per_edge == 2 )
        DVLOG(2) << "added edge with id :" << edgeInserted.id()
                 << " n1: " << edgeInserted.point( 0 ).node()
                 << " n2: " << edgeInserted.point( 1 ).node() << "\n";
    return edid;
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e )
{
    return addFace( mesh, __e, mpl::int_<mesh_type::nDim>() );
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addFace( mesh_type*, Feel::detail::GMSHElement const&, mpl::int_<1> )
{
    CHECK( false ) << "ImporterGmsh<MeshType>::addFace with dim=1 not valid";
    return 0;
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<2> )
{
    GmshOrdering<element_type> ordering;

    element_type e;
    e.setId( ( false )? __e.num : mesh->elements().size() );
    e.setProcessIdInPartition( this->worldComm().localRank() );
    if ( !__e.physical.empty() )
        e.addMarker( __e.physical );
    if ( false )
        e.setMarker2( __e.elementary );
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
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            point_type & pt = mesh->pointIterator( __e.indices[jj] )->second;
            if ( !e.isGhostCell() )
                pt.setProcessId( e.processId() );
            e.setPoint( ordering.fromGmshId( jj ), pt );
        }
    }
    auto [eit,inserted] = mesh->addElement( /*mesh->endElement(), */std::move(e) );
    auto const& [eid,eltInserted] = *eit;
    return eid;
}
template<typename MeshType>
int
ImporterGmsh<MeshType>::addFace( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<3> )
{
    GmshOrdering<face_type> ordering;

    face_type e;
    e.setProcessIdInPartition( this->worldComm().localRank() );
    e.setId( mesh->numFaces() );
    if ( !__e.physical.empty() )
        e.addMarker( __e.physical );
    if ( false )
        e.setMarker2( __e.elementary );
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
        for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
        {
            //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
            e.setPoint( ordering.fromGmshId( jj ), mesh->point( __e.indices[jj] ) );
        }
    }
    auto [fit,inserted] = mesh->addFace( std::move(e) );
    auto const& [fid,faceInserted] = *fit;
    return fid;
}

template<typename MeshType>
int
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e )
{
    return addVolume( mesh, __e, mpl::int_<mesh_type::nDim>() );
}
template<typename MeshType>
int
ImporterGmsh<MeshType>::addVolume( mesh_type*, Feel::detail::GMSHElement const&, mpl::int_<1> )
{
    CHECK( false ) << "ImporterGmsh<MeshType>::addVolume with dim=1 not valid";
    return 0;
}
template<typename MeshType>
int
ImporterGmsh<MeshType>::addVolume( mesh_type*, Feel::detail::GMSHElement const&, mpl::int_<2> )
{
    CHECK( false ) << "ImporterGmsh<MeshType>::addVolume with dim=2 not valid";
    return 0;
}
template<typename MeshType>
int
ImporterGmsh<MeshType>::addVolume( mesh_type* mesh, Feel::detail::GMSHElement const& __e, mpl::int_<3> )
{
    element_type e;
    e.setId( ( false )? __e.num : mesh->elements().size() );
    e.setProcessIdInPartition( this->worldComm().localRank() );
    GmshOrdering<element_type> ordering;
    if ( !__e.physical.empty() )
        e.addMarker( __e.physical );
    if ( false )
        e.setMarker2( __e.elementary );
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
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            //ptseen[mesh->point( __e.indices[jj] ).id()]=1;
            point_type & pt = mesh->pointIterator( __e.indices[jj] )->second;
            if (!e.isGhostCell()) pt.setProcessId( e.processId() );
            //std::cout << "gmsh index " << jj << " -> " << ordering.fromGmshId(jj) << " -> " << mesh->point( __e[jj] ).id()+1 << " : " << mesh->point( __e[jj] ).node() << "\n";
            e.setPoint( ordering.fromGmshId( jj ), pt );
        }
    }

    auto [eit,inserted] = mesh->addElement( /*mesh->endElement(), */std::move(e) );
    auto const& [eid,eltInserted] = *eit;
    return eid;
}

template<typename MeshType>
void
ImporterGmsh<MeshType>::updateGhostCellInfo( mesh_type* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
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

    for ( auto const& [proc,thedata] : dataToSend )
    {
        int nSendData = thedata.size();
        reqs[cptRequest++] = this->worldComm().localComm().isend( proc , 0, nSendData );
    }

    std::map<rank_type,size_type> sizeRecv;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( nbMsgToRecv[proc] > 0 )
        {
            reqs[cptRequest++] = this->worldComm().localComm().irecv( proc , 0, sizeRecv[proc] );
        }
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    //-----------------------------------------------------------//
    cptRequest=0;
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        int nSendData = itDataToSend->second.size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = this->worldComm().localComm().isend( itDataToSend->first , 0, &(itDataToSend->second[0]), nSendData );
    }
    //-----------------------------------------------------------//
    // first recv
    std::map<rank_type,std::vector<int> > dataToRecv;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( nbMsgToRecv[proc] > 0 )
        {
            int nRecvData = sizeRecv[proc];
            dataToRecv[proc].resize( nRecvData );
            if ( nRecvData > 0 )
                reqs[cptRequest++] = this->worldComm().localComm().irecv( proc , 0, &(dataToRecv[proc][0]), nRecvData );
        }
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);
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
        int nSendData = itDataToReSend->second.size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = this->worldComm().localComm().isend( itDataToReSend->first , 0, &(itDataToReSend->second[0]), nSendData );
    }
    //-----------------------------------------------------------//
    // recv the initial request
    std::map<rank_type, std::vector<int> > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
        int nRecvData = itDataToSend->second.size();
        finalDataToRecv[idProc].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = this->worldComm().localComm().irecv( idProc, 0, &(finalDataToRecv[idProc][0]), nRecvData );
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);
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
            auto & eltModified = mesh->elementIterator( memoryMsgToSend[idProc][k] )->second;
            eltModified.setIdInOtherPartitions( idProc, itFinalDataToRecv->second[k] );
        }
    }
    //-----------------------------------------------------------//
    DVLOG(1) << "updateGhostCellInfoNonBlockingComm : finish on rank "<< this->worldComm().localRank() << "\n";
}


// Gmsh reader factory
#if defined( FEELPP_HAS_GMSH_H )
#if !defined( FEELPP_HAS_GMSH_API )
using GmshReaderFactory = Feel::Singleton< std::map< std::string, std::function<int(std::string const&,std::string const&)> > >;
#endif
#endif


} // Feel



#endif /* __ImporterGmsh_H */
