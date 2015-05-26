/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-08-11

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file reloadmesh.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-08-11
 */

#ifndef __RELOADMESH_H
#define __RELOADMESH_H 1

namespace Feel
{
namespace FeelModels
{

    template <typename MeshType>
    boost::shared_ptr<MeshType>
    reloadMesh(std::string nameFile, WorldComm const& worldComm, int straighten=1 )
    {
        typedef MeshType mesh_type;
        //std::string nameFile = "FluidMechanicsMesh.path";
        std::ifstream file( nameFile.c_str() );
        if ( !file )
        {
            CHECK( false ) << "Fail to open the txt file containing path of msh file : " << nameFile << "\n";
        }

        std::string mshfile;
        if ( ! ( file >> mshfile ) )
        {
            CHECK( false ) << "Fail to read the msh path in file : " << nameFile << "\n";
        }
        file.close();

#if 0
        //auto mshfile=this->application()->vm()["fluid.mshfile"].as< std::string >() ;

        ImporterGmsh<mesh_type> import( mshfile );
        import.setVersion( "2.1" );

        boost::shared_ptr<mesh_type> mesh( new mesh_type );
        mesh->accept( import );
        mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
        mesh->updateForUse();


        if ( straighten && mesh_type::nOrder > 1 )
            return straightenMesh( mesh );
        else
            return mesh;
#else

        return loadGMSHMesh(_mesh=new mesh_type,
                            _filename=mshfile,
                            _worldcomm=worldComm,
                            _straighten=straighten,
                            _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );


#endif

    } // reloadFluidMesh()

} // namespace FeelModels
} // namespace Feel

#endif
