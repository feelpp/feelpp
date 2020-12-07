/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
             Christophe Prud'homme  <christophe.prudhomme@feelpp.org>
       Date: 2016-01-15

  Copyright (C) 2009-2016 Feel++ Consortium

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
   \file importerAcusimRawMesh.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \author Christophe Prud'homme  <christophe.prudhomme@feelpp.org>
   \date 2016-01-15
 */

#ifndef FEELPP_IMPORTERACUSIMRAWMESH_HPP
#define FEELPP_IMPORTERACUSIMRAWMESH_HPP 1

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/property_tree/ptree.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/importer.hpp>
#include <map>

#if defined( FEELPP_HAS_ACUSIM )
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <acusim.h>
#include <adb.h>
//#include <h3dpublic_export.h>
#include<sstream>
#include<fstream>
#include <math.h>
//#include <nrutil.h>
#endif

namespace pt = boost::property_tree;

namespace Feel
{
namespace pt = boost::property_tree;
namespace detail
{
inline const pt::ptree&
empty_ptree()
{
    static pt::ptree t;
    return t;
}

inline
std::string
markerNameFromAcusimName( std::string const& acusimName, std::string const& elementTypeName, std::string const& sepChar = " " )
{
    boost::char_separator<char> sep( sepChar.c_str() );
    boost::tokenizer< boost::char_separator<char> > kvlist(acusimName, sep);
    std::vector<std::string> acusimNameSplitted;
    for( const auto& ikvl : kvlist )
        acusimNameSplitted.push_back( ikvl );
    //std::cout << "ikvl " <<ikvl << "\n";

    int nNameSplitted = acusimNameSplitted.size();
    for ( int k=0;k<nNameSplitted;k=k+2 )
    {
        if ( (k+1) < nNameSplitted )
        {
            if ( acusimNameSplitted[k+1] == elementTypeName )
                return acusimNameSplitted[k];
        }
    }
    std::string markerName;
    if ( nNameSplitted == 1 )
        markerName = acusimName;

    return markerName;
}

} // namespace detail

/**
 * Reader for Acusim Raw Mesh Format from Acusim (C) Altair
 *
 * the main file is an XML file describing
 *  - the points coordinates
 *  - the connectivity of the elements
 *  - the connectivity of the marked faces
 * the faces are marked using the third string in the name attribute of the surface_set
 */
template <typename MeshType>
class ImporterAcusimRawMesh : public Importer<MeshType>
{
    typedef Importer<MeshType> super;

  public:
    typedef typename super::mesh_type mesh_type;
    typedef typename super::point_type point_type;
    typedef typename super::node_type node_type;
    typedef typename super::edge_type edge_type;
    typedef typename super::face_type face_type;
    typedef typename super::element_type element_type;

    explicit ImporterAcusimRawMesh( worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() )
        :
        super( ACUSIM, ASCII, _worldcomm )
    {
    }

    ImporterAcusimRawMesh( std::string const& filename, worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() );

    ImporterAcusimRawMesh( ImporterAcusimRawMesh const& ) = default;

    void visit( mesh_type* mesh ) override;

    std::string const& filenameNodes() const { return M_filenameNodes; }

    void setFilenameNodes( std::string const& path ) { M_filenameNodes = path; }

  private:
    void readNodes( mesh_type* mesh ) const;
    void readElements( mesh_type* mesh ) const;
    void readFaces( mesh_type* mesh ) const;
	void readMeshAcusolveAPI( mesh_type* mesh ) const;

  private:
    std::string M_filenameNodes;
    std::multimap<std::string, std::string> M_filenamesElements, M_filenamesFaces;
};

template <typename MeshType>
ImporterAcusimRawMesh<MeshType>::ImporterAcusimRawMesh( std::string const& filename, worldcomm_ptr_t const& wc )
    :
    super( ACUSIM, ASCII, wc )
{
    pt::ptree tree;
    if ( fs::exists( filename ) )
        pt::read_xml( filename, tree );
    else
        LOG( ERROR ) << "Invalid AcusimRawMesh filename " << filename;
    fs::path pp = fs::path( filename ).parent_path();
    {
        auto const& attributes = tree.get_child( "mesh.coordinates.<xmlattr>", Feel::detail::empty_ptree() );
        auto c_filename = attributes.get_optional<std::string>( "file" );
        if ( !c_filename )
            LOG( ERROR ) << "AcusimRawMesh coordinate not available";
        this->setFilenameNodes( ( pp / fs::path( *c_filename ) ).string() );
    }
    if ( wc->isMasterRank() )
        std::cout << ". loading AcusimRawMesh, coordinates " << this->filenameNodes() << "\n";
    std::string e_topology;

    for ( auto const& s : tree.get_child( "mesh", Feel::detail::empty_ptree() ) )
    {
        if ( s.first == "element_set" )
        {
            auto const& attributes = s.second.get_child("<xmlattr>", Feel::detail::empty_ptree() );
            auto name = attributes.get_optional<std::string>("name");
            if ( !name )
                LOG(ERROR) << "AcusimRawMesh surface_set has no name";
            std::string marker = *name;
            boost::erase_all(marker, "tet4");
            boost::erase_all(marker, "tria3");
            boost::trim(marker);
            boost::replace_all(marker, " ", "_");

            auto e_filename = attributes.get_optional<std::string>("file");
            e_topology = attributes.get<std::string>("topology", "");
            if ( !e_filename )
                LOG( ERROR ) << "AcusimRawMesh element_set not available";
            auto p = std::make_pair( marker, ( pp / fs::path( *e_filename ) ).string() );
            M_filenamesElements.insert( p );
            if ( wc->isMasterRank() )
                std::cout << "  .. element_set " << p.first << " - " << p.second << " (" << e_topology << ")\n";
        }
        else if ( s.first == "surface_set" )
        {
            auto const& attributes = s.second.get_child("<xmlattr>", Feel::detail::empty_ptree() );
            auto name = attributes.get_optional<std::string>("name");
            if ( !name )
                LOG(ERROR) << "AcusimRawMesh surface_set has no name";
            std::string marker = *name;

            auto parent = attributes.get_optional<std::string>("parent");
            if ( parent )
                boost::erase_all(marker, *parent);
            boost::erase_all(marker, "tet4");
            boost::erase_all(marker, "tria3");
            boost::trim(marker);
            boost::replace_all(marker, " ", "_");

            auto e_filename = attributes.get_optional<std::string>("file");
            e_topology = attributes.get<std::string>("topology", "");
            if ( !e_filename )
                LOG( ERROR ) << "AcusimRawMesh surface_set not available";
            auto p = std::make_pair( marker, ( pp / fs::path( *e_filename ) ).string() );
            M_filenamesFaces.insert( p );
            if ( wc->isMasterRank() )
                std::cout << "  .. surface_set " << p.first << " - " << p.second << " (" << e_topology << ")\n";
        }
    }
}

template <typename MeshType>
void ImporterAcusimRawMesh<MeshType>::visit( mesh_type* mesh )
{
	if ( this->fileType() == ASCII )
	{
        CHECK( fs::exists( this->filenameNodes() ) ) << "filenameNodes not exist : " << this->filenameNodes();

        if ( mesh->worldComm().isMasterRank() )
            std::cout << ".reading nodes file" << std::endl;
        tic();
        this->readNodes( mesh );
        toc( "ImporterAcusimRawMesh read nodes" );
        if ( mesh->worldComm().isMasterRank() )
        {
            std::cout << ".reading nodes file done" << std::endl;
            std::cout << ".reading elements file..." << std::endl;
        }
        tic();
        this->readElements( mesh );
        toc( "ImporterAcusimRawMesh read elements" );
        if ( mesh->worldComm().isMasterRank() )
        {
            std::cout << ".reading elements file done" << std::endl;
            std::cout << ".reading faces files..." << std::endl;
        }
        tic();
        this->readFaces( mesh );
        toc( "ImporterAcusimRawMesh read faces" );
        if ( mesh->worldComm().isMasterRank() )
        {
            std::cout << ".reading faces files done" << std::endl;
        }
    }
	else
	{
		this->readMeshAcusolveAPI( mesh );
	}

}
template <typename MeshType>
void ImporterAcusimRawMesh<MeshType>::readNodes( mesh_type* mesh ) const
{
    std::ifstream __is( this->filenameNodes().c_str() );

    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        LOG( ERROR ) << "Invalid file name " << this->filenameNodes() << " (file not found)";
        ostr << "Invalid file name " << this->filenameNodes() << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    node_type coords( 3 /*mesh_type::nRealDim*/ );
    int ptid = 0;
    rank_type partId = mesh->worldComm().localRank();

    while ( !__is.eof() )
    {
        __is >> ptid;
        if ( __is.eof() )
            break;
        __is >> coords[0] >> coords[1] >> coords[2];

        point_type pt( ptid, coords, false );
        pt.setProcessIdInPartition( partId );
        pt.setProcessId( partId );
        mesh->addPoint( pt );
    }

    __is.close();
}
template <typename MeshType>
void ImporterAcusimRawMesh<MeshType>::readElements( mesh_type* mesh ) const
{

    int eltid = 0;
    rank_type partId = mesh->worldComm().localRank();
    element_type e;
    int ptids = 0;
    int markerId = 1;
    //std::vector<int> ptids( element_type::numPoints );

    for ( auto const& datafileElements : M_filenamesElements )
    {

        std::string const& marker = datafileElements.first;
        std::string const& filenameElements = datafileElements.second;
        if ( mesh->worldComm().isMasterRank() )
        {
            std::cout << "reading element file with marker " << marker << " : " << filenameElements << std::endl;
        }

        std::ifstream __is( filenameElements.c_str() );

        if ( !__is.is_open() )
        {
            std::ostringstream ostr;
            LOG( ERROR ) << "Invalid file name " << filenameElements << " (file not found)";
            ostr << "Invalid file name " << filenameElements << " (file not found)\n";
            throw std::invalid_argument( ostr.str() );
            continue;
        }

        mesh->addMarkerName( marker, markerId, mesh_type::nDim );

        while ( !__is.eof() )
        {
            __is >> eltid;
            if ( __is.eof() )
                break;

            e.setId( eltid );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );
            e.setMarker( markerId );

            for ( size_type k = 0; k < element_type::numPoints; ++k )
            {
                __is >> ptids;
                e.setPoint( k, mesh->point( ptids ) );
            }
            //mesh->addElement( e, false );
            mesh->addElement( e, true );
        }

        __is.close();
        ++markerId;
    }
}

template <typename MeshType>
void ImporterAcusimRawMesh<MeshType>::readFaces( mesh_type* mesh ) const
{
    if ( M_filenamesFaces.empty() )
        return;

    rank_type partId = mesh->worldComm().localRank();
    int eltidConnected = 0, faceid = 0, faceidTemp = 0, ptids = 0;
    face_type f;

    // get next markerId as maxPreviousMarkerId + 1
    int markerId = 0;
    for ( auto const& _marker : mesh->markerNames() )
        markerId = std::max( (int)_marker.second[0], (int)markerId );
    ++markerId;

    for ( auto const& datafileFaces : M_filenamesFaces )
    {

        std::string const& marker = datafileFaces.first;
        std::string const& filenameFaces = datafileFaces.second;
        if ( mesh->worldComm().isMasterRank() )
        {
            std::cout << "reading faces file with marker " << marker << " : " << filenameFaces << std::endl;
        }

        std::ifstream __is( filenameFaces.c_str() );

        if ( !__is.is_open() )
        {
            std::ostringstream ostr;
            LOG( ERROR ) << "Invalid file name " << filenameFaces << " (file not found)";
            ostr << "Invalid file name " << filenameFaces << " (file not found)\n";
            throw std::invalid_argument( ostr.str() );
            continue;
        }

        mesh->addMarkerName( marker, markerId, mesh_type::nDim - 1 );

        while ( !__is.eof() )
        {
            __is >> eltidConnected;
            if ( __is.eof() )
                break;
            __is >> faceidTemp;

            f.setId( faceid++ );
            f.setProcessIdInPartition( partId );
            f.setProcessId( partId );
            f.setMarker( markerId );

            for ( size_type k = 0; k < face_type::numPoints; ++k )
            {
                __is >> ptids;
                f.setPoint( k, mesh->point( ptids ) );
            }
            mesh->addFace( f );
        }

        __is.close();
        ++markerId;
    }
}


template <typename MeshType>
void ImporterAcusimRawMesh<MeshType>::readMeshAcusolveAPI( mesh_type* mesh ) const
{
#if defined( FEELPP_HAS_ACUSIM )
	AdbHd	adbHd;
	//adbHd	= adbNew( "Marche_GDR", "ACUSIM.DIR", "/calcul/BIG_CALCUL/ACUSOLVE/Yoann/feelpp_test/test_install_PO/marche_GDR/", 0 );
	adbHd	= adbNew( (char*)"test", (char*)"ACUSIM.DIR", (char*)"/home/u2/chabannes/po_smallmesh/DATA_PO/", 0 );
	CHECK( adbHd != NULL ) << "Error opening the data base " << adbGetError();

    // init database loading
	bool success = adbOpenRun (adbHd,1);

    //-------------------------------------------------------------//
    // load mesh points
    //-------------------------------------------------------------//
	int nNodes =0;
	success = adbGetInt0  ( adbHd, (char*)"nNodes",&nNodes);

	CHECK( success ) << "fail to get nNodes\n";
	double crd[nNodes*3];
	std::cout << "PO : load nNodes " << nNodes << "\n";
    success = adbGetReals0 ( adbHd,(char*)"crd", crd, 3 ,nNodes);

	int usrIds[nNodes];
	//usrIds = (int*) malloc(nNodes*sizeof(int));
    success = adbGetInts0( adbHd,(char*)"usrIds", usrIds, nNodes , 1);

	// add points in mesh
	node_type coords( 3 );
    int ptid = 0;
    rank_type partId = mesh->worldComm().localRank();
	for( int i=0;i<nNodes;i++)
      {
		  coords[0]=crd[3*i];
		  coords[1]=crd[3*i+1];
		  coords[2]=crd[3*i+2];
		  // feel mesh point
		 point_type pt( /*ptid*/usrIds[i], coords, false );
         pt.setProcessIdInPartition( partId );
         pt.setProcessId( partId );
         mesh->addPoint( pt );
	  }

    //-------------------------------------------------------------//
    // load mesh elements
    //-------------------------------------------------------------//
    int eltid = 0;
    element_type e;
    e.setProcessIdInPartition( partId );
    e.setProcessId( partId );
    int markerId = 1;

    // get number of group of element
    int nElms = 0;
    success = adbGetInt0(adbHd,(char*)"nElms",&nElms);
    std::cout << "number of domain element " << nElms << "\n";
    int numberOfDomain = nElms; // seems necessary, to understand!!!
    int nElmElems = 0 , nElmElemNodes = 0;
    for (int j=0;j<numberOfDomain/*nElms*/;++j)
    {
        char elmName[1024] ;
        success = adbGetStr1(adbHd,(char*)"elmName", elmName,j);
        std::cout << "  Load element group "<< j << " with name : " << elmName << "\n";
        success = adbGetInt1(adbHd,(char*)"nElmElems",&nElmElems,j);
        std::cout << "  number of element : " << nElmElems << "\n";
        success = adbGetInt1(adbHd,(char*)"nElmElemNodes",&nElmElemNodes,j);
        std::cout << "  number of point in element : " << nElmElemNodes << "\n";
        CHECK( nElmElemNodes == element_type::numPoints ) << "number of point in element wrong " << nElmElemNodes <<" vs " << element_type::numPoints;

        // add marker
        std::string marker = Feel::detail::markerNameFromAcusimName( std::string(elmName), "tet4" );
        std::cout << "  marker element : " << marker << "\n";
        mesh->addMarkerName( marker, markerId, mesh_type::nDim );
        e.setMarker( markerId );

        // load element connectivity
        int * elmCnn = new int[nElmElems*nElmElemNodes];
        success = adbGetInts1(adbHd,(char*)"elmCnn",elmCnn,j,nElmElemNodes,nElmElems);
        // add elements in mesh
        for( int i=0;i<nElmElems;++i )
        {
            for ( uint16_type p = 0; p < element_type::numPoints; ++p )
            {
                int ptid = elmCnn[element_type::numPoints*i+p];
                e.setPoint( p, mesh->point( ptid ) );
            }
            mesh->addElement( e, true );
        }
        delete [] elmCnn;
        ++markerId;
    }


    //-------------------------------------------------------------//
    // load mesh faces
    //-------------------------------------------------------------//
    face_type f;
    f.setProcessIdInPartition( partId );
    f.setProcessId( partId );

    //int eltidConnected = 0, faceid = 0, faceidTemp = 0, ptids = 0;
    int faceid = 0;

    // get number of group of element
	int nOsfs = 0;
	success = adbGetInt0(adbHd,(char*)"nOsfs",&nOsfs);
    std::cout << "number of domain face " << nOsfs << "\n";
    int numberOfDomainFace = nOsfs; // seems necessary, to understand!!!
    int nOsfSrfs = 0, nOsfSrfNodes = 0;
	for( int j=0;j<numberOfDomainFace/*nOsfs*/;++j )
	{
        char osfName[1024] ;
        success = adbGetStr1(adbHd,(char*)"osfName",/*&*/osfName,j);
        std::cout << "  Load face group "<< j << " with name : " << osfName << "\n";
        success = adbGetInt1(adbHd,(char*)"nOsfSrfs",&nOsfSrfs,j);
        std::cout << "  number of face : " << nOsfSrfs << "\n";
        success = adbGetInt1(adbHd,(char*)"nOsfSrfNodes",&nOsfSrfNodes,j);
        std::cout << "  number of point in face : " << nOsfSrfNodes << "\n";
        CHECK( nOsfSrfNodes == face_type::numPoints ) << "number of point in element wrong " << nOsfSrfNodes <<" vs " << face_type::numPoints;

        std::string marker = Feel::detail::markerNameFromAcusimName( std::string(osfName), "tria3" );
        std::cout << "  marker face : " << marker << "\n";

        // add marker
        mesh->addMarkerName( marker, markerId, mesh_type::nDim - 1 );
        f.setMarker( markerId );

        // load face connectivity
        int * osfSrfCnn = new int[nOsfSrfNodes*nOsfSrfs];
        success = adbGetInts1(adbHd,(char*)"osfSrfCnn",osfSrfCnn,j,nOsfSrfNodes,nOsfSrfs);

        // add faces in mesh
        for( int i=0;i<nOsfSrfs;++i)
        {
            f.setId( faceid++ );
            for ( uint16_type p = 0; p < face_type::numPoints; ++p )
            {
                int ptid = osfSrfCnn[face_type::numPoints*i+p];
                f.setPoint( p, mesh->point( ptid ) );
            }
            mesh->addFace( f );
        }
        delete [] osfSrfCnn;
        ++markerId;
	}


    success = adbCloseRun( adbHd );

    adbFree( adbHd );
    //std::cout << "FINISH ACUSIMREAD\n";
#endif
}


} // namespace Feel

#endif // FEELPP_IMPORTERACUSIMRAWMESH_HPP
