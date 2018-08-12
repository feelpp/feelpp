/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
   \file fsimesh.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-08-11
 */

#ifndef FEELPP_MODELS_FSIMESH_H
#define FEELPP_MODELS_FSIMESH_H 1

#include <feel/feeldiscr/mesh.hpp>


namespace Feel
{
namespace FeelModels
{


template< class ConvexType >
class FSIMesh
{
public :

    typedef Mesh<ConvexType> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type mesh_fluid_type;
    typedef mesh_type mesh_solid_type;
    typedef typename mesh_fluid_type::shape_type shape_fluid_type;
    typedef typename mesh_solid_type::shape_type shape_solid_type;

    FSIMesh( std::string prefix, WorldComm & worldcomm = Environment::worldComm() );

    std::string prefix() const { return M_prefix; }
    WorldComm & worldComm() { return *M_worldComm; }
    WorldComm const& worldComm() const { return *M_worldComm; }

    fs::path const& geoPathFSI() const { return M_geoPathFSI; }
    fs::path const& mshPathFSI() const { return M_mshPathFSI; }
    void setGeoPathFSI( fs::path const& path ) { M_geoPathFSI = path; }
    void setMshPathFSI( fs::path const& path ) { M_mshPathFSI = path; }

    fs::path const& mshPathFluidPartN() const { return M_mshfilepathFluidPartN; }
    fs::path const& mshPathSolidPartN() const { return M_mshfilepathSolidPartN; }
    fs::path const& mshPathFluidPart1() const { return M_mshfilepathFluidPart1; }
    fs::path const& mshPathSolidPart1() const { return M_mshfilepathSolidPart1; }
    void setFluidMshPathPart1( fs::path const& path ) { M_mshfilepathFluidPart1 = path; }
    void setSolidMshPathPart1( fs::path const& path ) { M_mshfilepathSolidPart1 = path; }
    void setFluidMshPathPartN( fs::path const& path ) { M_mshfilepathFluidPartN = path; }
    void setSolidMshPathPartN( fs::path const& path ) { M_mshfilepathSolidPartN = path; }

    double meshSize() const { return M_meshSize; }
    void setMeshSize(double d) { M_meshSize=d; }


    bool forceRebuild() const { return M_forceRebuild; }
    void setForceRebuild(bool b) {M_forceRebuild=b; }

    int nPartitions() const { return M_nPartitions; }
    int partitioner() const { return M_partitioner; }
    void setNumberOfPartitions(int p) { M_nPartitions=p; }
    void setPartitioner(int p) { M_partitioner=p; }

    std::set<std::string> const& markersNameFluidVolume() const { return M_markersNameFluidVolume; }
    std::set<std::string> const& markersNameSolidVolume() const { return M_markersNameSolidVolume; }
    void setMarkersNameFluidVolume( std::set<std::string> const& s ) { M_markersNameFluidVolume=s; }
    void setMarkersNameSolidVolume( std::set<std::string> const& s ) { M_markersNameSolidVolume=s; }


    //-------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------------------------------------------//

    void buildFSIMeshFromMsh();
    void buildFSIMeshFromGeo();
    void buildSubMesh( mesh_ptrtype const& fsimesh );
    void buildMeshesPartitioning();

private :

    std::string M_prefix;
    std::shared_ptr<WorldComm> M_worldComm;
    fs::path M_geoPathFSI,M_mshPathFSI;
    fs::path M_mshfilepathFluidPartN,M_mshfilepathFluidPart1;
    fs::path M_mshfilepathSolidPartN,M_mshfilepathSolidPart1;
    std::set<std::string> M_markersNameFluidVolume, M_markersNameSolidVolume;

    double M_meshSize;
    int M_nPartitions;
    int M_partitioner;

    //mesh_ptrtype M_meshFluid;
    //mesh_ptrtype M_meshSolid;

    bool M_forceRebuild;

};

} // namespace FeelModels
} // namespace Feel


#endif // FEELPP_MODELS_FSIMESH_H
