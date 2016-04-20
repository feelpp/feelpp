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

#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <feel/feelfilters/importer.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace pt =  boost::property_tree;

namespace Feel
{
namespace pt = boost::property_tree;
namespace detail
{
inline 
const pt::ptree& 
empty_ptree()
{
    static pt::ptree t;
    return t;
}
}

/**
 * Reader for Acusim Raw Mesh Format from Acusim (C) Altair
 *
 * the main file is an XML file describing
 *  - the points coordinates
 *  - the connectivity of the elements
 *  - the connectivity of the marked faces
 * the faces are marked using the third string in the name attribute of the surface_set
 */
template<typename MeshType>
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

    ImporterAcusimRawMesh( WorldComm const& _worldcomm = Environment::worldComm() )
        {}

    ImporterAcusimRawMesh( std::string const& filename, WorldComm const& _worldcomm = Environment::worldComm() );

    ImporterAcusimRawMesh( ImporterAcusimRawMesh const& ) = default;

    void visit( mesh_type* mesh );

    std::string const& filenameNodes() const { return M_filenameNodes; }

    void setFilenameNodes( std::string const& path ) { M_filenameNodes = path; }


private :
    void readNodes( mesh_type* mesh ) const;
    void readElements( mesh_type* mesh ) const;
    void readFaces( mesh_type* mesh ) const;

private :
    std::string M_filenameNodes;
    std::multimap<std::string,std::string> M_filenamesElements, M_filenamesFaces;
};

template<typename MeshType>
ImporterAcusimRawMesh<MeshType>::ImporterAcusimRawMesh( std::string const& filename, WorldComm const& wc )
{
    pt::ptree tree;
    if ( fs::exists( filename ) )
        pt::read_xml(filename, tree);
    else
        LOG(ERROR) << "Invalid AcusimRawMesh filename " << filename;
    fs::path pp = fs::path( filename ).parent_path();
    {
        auto const& attributes = tree.get_child("mesh.coordinates.<xmlattr>", Feel::detail::empty_ptree() );
        auto c_filename = attributes.get_optional<std::string>("file");
        if ( !c_filename )
            LOG(ERROR) << "AcusimRawMesh coordinate not available";
        this->setFilenameNodes( (pp/fs::path(*c_filename)).string() );
    }
    if ( wc.isMasterRank() )
        std::cout << ". loading AcusimRawMesh, coordinates " << this->filenameNodes() << "\n";
    std::string e_topology;

    for( auto const& s : tree.get_child("mesh", Feel::detail::empty_ptree()) )
    {
        if ( s.first == "element_set" )
        {
            auto const& attributes = s.second.get_child("<xmlattr>", Feel::detail::empty_ptree() );
            auto name = attributes.get_optional<std::string>("name");
            typedef std::vector< std::string > split_vector_type;
            split_vector_type split_name; // #2: Search for tokens
            boost::split( split_name, *name, boost::is_any_of("\t "), boost::token_compress_on );
            std::string marker = split_name[0];
            auto e_filename = attributes.get_optional<std::string>("file");
            e_topology = attributes.get<std::string>("topology", "");
            if ( !e_filename )
                LOG(ERROR) << "AcusimRawMesh element_set not available";
            auto p = std::make_pair(marker, (pp/fs::path(*e_filename)).string()  );
            M_filenamesElements.insert( p );
            if ( wc.isMasterRank() )
                std::cout << "  .. element_set " << p.first << " - " << p.second << " (" << e_topology << ")\n";
        }
        else if ( s.first == "surface_set" )
        {
            auto const& attributes = s.second.get_child("<xmlattr>", Feel::detail::empty_ptree() );
            auto name = attributes.get_optional<std::string>("name");

            typedef std::vector< std::string > split_vector_type;

            split_vector_type split_name; // #2: Search for tokens
            boost::split( split_name, *name, boost::is_any_of("\t "), boost::token_compress_on );

            std::string marker = split_name[2];
            auto e_filename = attributes.get_optional<std::string>("file");
            e_topology = attributes.get<std::string>("topology", "");
            if ( !e_filename )
                LOG(ERROR) << "AcusimRawMesh surface_set not available";
            auto p = std::make_pair(marker, (pp/fs::path(*e_filename)).string()  );
            M_filenamesFaces.insert( p );
            if ( wc.isMasterRank() )
                std::cout << "  .. surface_set " << p.first << " - " << p.second << " (" << e_topology << ")\n";
        }
    }

}

template<typename MeshType>
void
ImporterAcusimRawMesh<MeshType>::visit( mesh_type* mesh )
{
    CHECK( fs::exists( this->filenameNodes() ) ) << "filenameNodes not exist : " << this->filenameNodes();

    if ( mesh->worldComm().isMasterRank() )
        std::cout << ".reading nodes file" << std::endl;
    tic();
    this->readNodes( mesh );
    toc("ImporterAcusimRawMesh read nodes");
    if ( mesh->worldComm().isMasterRank() )
    {
        std::cout << ".reading nodes file done" << std::endl;
        std::cout << ".reading elements file..." << std::endl;
    }
    tic();
    this->readElements( mesh );
    toc("ImporterAcusimRawMesh read elements");
    if ( mesh->worldComm().isMasterRank() )
    {
        std::cout << ".reading elements file done" << std::endl;
        std::cout << ".reading faces files..." << std::endl;
    }
    tic();
    this->readFaces( mesh );
    toc("ImporterAcusimRawMesh read faces");
    if ( mesh->worldComm().isMasterRank() )
    {
        std::cout << ".reading faces files done" << std::endl;
    }
}
template<typename MeshType>
void
ImporterAcusimRawMesh<MeshType>::readNodes( mesh_type* mesh ) const
{
    std::ifstream __is( this->filenameNodes().c_str() );

    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        LOG(ERROR) << "Invalid file name " << this->filenameNodes() << " (file not found)";
        ostr << "Invalid file name " << this->filenameNodes() << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    node_type coords( 3/*mesh_type::nRealDim*/ );
    int ptid = 0;
    rank_type partId = mesh->worldComm().localRank();

    while ( !__is.eof() )
    {
        __is >> ptid;
        if (__is.eof() )
            break;
        __is >> coords[0] >> coords[1] >> coords[2];

        point_type pt( ptid, coords, false );
        pt.setProcessIdInPartition( partId );
        pt.setProcessId( partId );
        mesh->addPoint( pt );
    }

    __is.close();
}
template<typename MeshType>
void
ImporterAcusimRawMesh<MeshType>::readElements( mesh_type* mesh ) const
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
            std::cout << "reading element file with marker "<< marker  << " : "<< filenameElements << std::endl;
        }

        std::ifstream __is ( filenameElements.c_str() );

        if ( !__is.is_open() )
        {
            std::ostringstream ostr;
            LOG(ERROR) << "Invalid file name " << filenameElements << " (file not found)";
            ostr << "Invalid file name " << filenameElements << " (file not found)\n";
            throw std::invalid_argument( ostr.str() );
            continue;
        }

        mesh->addMarkerName( marker, markerId, mesh_type::nDim );

        while ( !__is.eof() )
        {
            __is >> eltid;
            if (__is.eof() )
                break;

            e.setId( eltid );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );
            e.setMarker( markerId );
            e.setMarker2( markerId );

            for (size_type k = 0; k < element_type::numPoints; ++k)
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

template<typename MeshType>
void
ImporterAcusimRawMesh<MeshType>::readFaces( mesh_type* mesh ) const
{
    if ( M_filenamesFaces.empty() )
        return;

    rank_type partId = mesh->worldComm().localRank();
    int eltidConnected = 0, faceid = 0, faceidTemp = 0, ptids = 0;
    face_type f;

    // get next markerId as maxPreviousMarkerId + 1
    int markerId = 0;
    for( auto const& _marker: mesh->markerNames() )
        markerId = std::max( (int)_marker.second[0], (int)markerId );
    ++markerId;

    for ( auto const& datafileFaces : M_filenamesFaces )
    {

        std::string const& marker = datafileFaces.first;
        std::string const& filenameFaces = datafileFaces.second;
        if ( mesh->worldComm().isMasterRank() )
        {
            std::cout << "reading faces file with marker "<< marker  << " : "<< filenameFaces << std::endl;
        }

        std::ifstream __is ( filenameFaces.c_str() );

        if ( !__is.is_open() )
        {
            std::ostringstream ostr;
            LOG(ERROR) << "Invalid file name " << filenameFaces << " (file not found)";
            ostr << "Invalid file name " << filenameFaces << " (file not found)\n";
            throw std::invalid_argument( ostr.str() );
            continue;
        }

        mesh->addMarkerName( marker, markerId, mesh_type::nDim-1 );

        while ( !__is.eof() )
        {
            __is >> eltidConnected;
            if (__is.eof() )
                break;
            __is >> faceidTemp;

            f.setId( faceid++ );
            f.setProcessIdInPartition( partId );
            f.setProcessId( partId );
            f.setMarker( markerId );
            f.setMarker2( markerId );

            for (size_type k = 0; k < face_type::numPoints; ++k)
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

} // namespace Feel

#endif // FEELPP_IMPORTERACUSIMRAWMESH_HPP
