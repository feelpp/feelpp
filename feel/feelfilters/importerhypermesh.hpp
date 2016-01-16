/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
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
   \file importerhypermesh.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-01-15
 */

#ifndef FEELPP_IMPORTERHYPERMESH_HPP
#define FEELPP_IMPORTERHYPERMESH_HPP 1

#include <feel/feelfilters/importer.hpp>
#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{

template<typename MeshType>
class ImporterHyperMesh : public Importer<MeshType>
{
    typedef Importer<MeshType> super;
public:

    typedef typename super::mesh_type mesh_type;
    typedef typename super::point_type point_type;
    typedef typename super::node_type node_type;
    typedef typename super::edge_type edge_type;
    typedef typename super::face_type face_type;
    typedef typename super::element_type element_type;

    ImporterHyperMesh( WorldComm const& _worldcomm = Environment::worldComm() )
        {}

    ImporterHyperMesh( ImporterHyperMesh const& ) = default;

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
};

template<typename MeshType>
void
ImporterHyperMesh<MeshType>::visit( mesh_type* mesh )
{
    CHECK( fs::exists( this->filenameNodes() ) ) << "filenameNodes not exist : " << this->filenameNodes();
    CHECK( fs::exists( this->filenameElements() ) ) << "filenameElements not exist : " << this->filenameElements();
    this->readNodes( mesh );
    this->readElements( mesh );
}

template<typename MeshType>
void
ImporterHyperMesh<MeshType>::readNodes( mesh_type* mesh ) const
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
ImporterHyperMesh<MeshType>::readElements( mesh_type* mesh ) const
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

#endif // FEELPP_IMPORTERHYPERMESH_HPP
