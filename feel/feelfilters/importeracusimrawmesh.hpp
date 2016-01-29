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
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <feel/feelfilters/importer.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{
namespace detail
{
const pt::ptree& empty_ptree(){
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
    std::string const& filenameElements() const { return M_filenameElements; }

    void setFilenameNodes( std::string const& path ) { M_filenameNodes = path; }
    void setFilenameElements( std::string const& path ) { M_filenameElements = path; }
    

private :
    void readNodes( mesh_type* mesh ) const;
    void readElements( mesh_type* mesh ) const;

private :
    std::string M_filenameNodes, M_filenameElements;
    std::multimap<std::string,std::string> M_filenamesFaces;
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
    std::string e_topology;
    {
        auto const& attributes = tree.get_child("mesh.element_set.<xmlattr>", Feel::detail::empty_ptree() );
        auto e_filename = attributes.get_optional<std::string>("file");
        e_topology = attributes.get<std::string>("topology", "");
        if ( !e_filename )
            LOG(ERROR) << "AcusimRawMesh element_set not available";
        this->setFilenameElements( (pp/fs::path(*e_filename)).string() );
    }
    if ( wc.globalRank() == 0 ) 
        std::cout << ". loading AcusimRawMesh, coordinates " << this->filenameNodes() << " element_set " << this->filenameElements() << " (" << e_topology << ")\n";

    for( auto const& s : tree.get_child("mesh", Feel::detail::empty_ptree()) )
    {
        if ( s.first != "surface_set" )
            continue;
        auto const& attributes = s.second.get_child("<xmlattr>", Feel::detail::empty_ptree() );
        auto name = attributes.get_optional<std::string>("name");
        
        typedef std::vector< std::string > split_vector_type;
    
        split_vector_type split_name; // #2: Search for tokens
        boost::split( split_name, *name, boost::is_any_of("\t "), boost::token_compress_on ); 

        std::string marker = split_name[2];
        auto e_filename = attributes.get_optional<std::string>("file");
        e_topology = attributes.get<std::string>("topology", "");
        if ( !e_filename )
            LOG(ERROR) << "AcusimRawMesh element_set not available";
        auto p = std::make_pair(marker, (pp/fs::path(*e_filename)).string()  );
        M_filenamesFaces.insert( p ); 
        if ( wc.globalRank() == 0 ) 
            std::cout << "  .. surface_set " << p.first << " - " << p.second << " (" << e_topology << ")\n";
    }

}

template<typename MeshType>
void
ImporterAcusimRawMesh<MeshType>::visit( mesh_type* mesh )
{
    CHECK( fs::exists( this->filenameNodes() ) ) << "filenameNodes not exist : " << this->filenameNodes();
    CHECK( fs::exists( this->filenameElements() ) ) << "filenameElements not exist : " << this->filenameElements();
    if ( mesh->worldComm().globalRank() == 0 ) 
        std::cout << ".reading nodes file" << std::endl;
    tic();
    this->readNodes( mesh );
    toc("ImporterAcusimRawMesh read nodes");
    if ( mesh->worldComm().globalRank() == 0 ) 
    {
        std::cout << ".reading nodes file done" << std::endl;
        std::cout << ".reading elements file..." << std::endl;
    }
    tic();
    this->readElements( mesh );
    toc("ImporterAcusimRawMesh read elements");
    if ( mesh->worldComm().globalRank() == 0 ) 
    {
        std::cout << ".reading elements file" << std::endl;
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
        mesh->addPoint( pt );
    }

    __is.close();
}
template<typename MeshType>
void
ImporterAcusimRawMesh<MeshType>::readElements( mesh_type* mesh ) const
{
    std::ifstream __is ( this->filenameElements().c_str() );

    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        LOG(ERROR) << "Invalid file name " << this->filenameElements() << " (file not found)";
        ostr << "Invalid file name " << this->filenameElements() << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    int eltid = 0;
    rank_type partId = mesh->worldComm().localRank();

    element_type e;
    int ptids = 0;
    //std::vector<int> ptids( element_type::numPoints );
    while ( !__is.eof() )
    {
        __is >> eltid;
        if (__is.eof() )
            break;

        e.setId( eltid );
        e.setProcessIdInPartition( partId );
        e.setProcessId( partId );

        for (size_type k = 0; k < element_type::numPoints; ++k)
        {
            __is >> ptids;
            e.setPoint( k, mesh->point( ptids ) );
        }
        //mesh->addElement( e, false );
        mesh->addElement( e, true );
    }

    __is.close();

}

} // namespace Feel

#endif // FEELPP_IMPORTERACUSIMRAWMESH_HPP
