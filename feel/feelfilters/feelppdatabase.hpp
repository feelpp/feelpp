/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2017-02-09

  Copyright (C) 2017 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
 \file feelppdatabase.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2017-02-09
 */
#ifndef FEELPP_FEELPPDATABASE_HPP
#define FEELPP_FEELPPDATABASE_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/removecomments.hpp>
#include <feel/feelcore/utility.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>

namespace Feel
{

template<typename MeshType>
class FeelppDatabase
{
public :
    using mesh_type = decay_type<MeshType>;
    //using space_ptrtype = boost::shared_ptr<space_type>;
    //using mesh_type = typename space_type::mesh_type;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;
    //using element_type = typename space_type::element_type;
    //using element_ptrtype = typename space_type::element_ptrtype;

    FeelppDatabase( WorldComm const& worldComm = Environment::worldComm() )
        :
        M_worldComm( worldComm ),
        M_name( "field" ),
        M_dbRepository( fs::current_path() ),
#ifdef FEELPP_HAS_HDF5
        M_vectorFileFormat( "hdf5" )
#else
        M_vectorFileFormat( "boost" )
#endif
        {}

    WorldComm const& worldComm() const { return M_worldComm; }

    fs::path const& dbRepository() const { return M_dbRepository; }

    std::string const& meshFilename() const { return M_meshFilename; }

    std::vector<double> timeSet() const { return M_timeSet; }

    std::string const& vectorFileFormat() const { return M_vectorFileFormat; }


    void setDbRepository( std::string const& s ) { M_dbRepository = s; }

    void setMesh( mesh_ptrtype const& mesh )
        {
            M_mesh = mesh;
        }

    void loadInfo()
        {
            std::string jsonFilename = "mydb.dbfeelpp.json";
            std::string jsonPath = (M_dbRepository/jsonFilename).string();

            if ( !fs::exists( jsonPath ) )
            {
                LOG(INFO) << "Could not find " << jsonPath << std::endl;
                return;
            }

            auto json_str_wo_comments = removeComments(readFromFile(jsonPath));

            boost::property_tree::ptree ptreeDb;
            std::istringstream istr( json_str_wo_comments );
            boost::property_tree::read_json( istr, ptreeDb );

            try {
                std::string meshFilename  = ptreeDb.template get<std::string>( "mesh-filename" );
                M_meshFilename = meshFilename;
            }
            catch ( boost::property_tree::ptree_bad_path& e )
            {
                std::cout << "no mesh-filename in json\n";
            }

            M_timeSet.clear();
            for (auto const& item : ptreeDb.get_child("time-set"))
            {
                double time = item.second.template get_value<double>();
                std::cout << "time reload " << time << "\n";
                M_timeSet.push_back( time );
            }

        }

    void saveInfo()
        {
            if ( this->worldComm().isMasterRank() )
            {
                boost::property_tree::ptree ptreeDb;
                if ( !M_meshFilename.empty() )
                    ptreeDb.add( "mesh-filename",M_meshFilename );

                if ( !M_timeSet.empty() )
                {
                    boost::property_tree::ptree ptreeDbTimeSet;
                    //ptreeParameterName2.put( "", this->parameterName(d) );
                    for (int k=0;k<M_timeSet.size();++k )
                    {
                        boost::property_tree::ptree ptreeDbTimeValue;
                        ptreeDbTimeValue.put( "", M_timeSet[k] );
                        ptreeDbTimeSet.push_back( std::make_pair("", ptreeDbTimeValue) );
                    }
                    ptreeDb.add_child( "time-set", ptreeDbTimeSet );
                }

                if ( !M_fieldInfoInDb.empty() )
                {
                    boost::property_tree::ptree ptreeDbFieldInfos;
                    for ( auto const& fieldInfoPair : M_fieldInfoInDb )
                    {
                        std::string const& fieldName = fieldInfoPair.first;
                        std::string const& basisName = std::get<0>( fieldInfoPair.second );
                        uint16_type basisOrder = std::get<1>( fieldInfoPair.second );
                        uint16_type nComp = std::get<2>( fieldInfoPair.second );
                        boost::property_tree::ptree ptreeDbFieldInfo;
                        ptreeDbFieldInfo.add("basis-name",basisName);
                        ptreeDbFieldInfo.add("basis-order",basisOrder);
                        ptreeDbFieldInfo.add("nComp",nComp);
                        ptreeDbFieldInfos.add_child( fieldName, ptreeDbFieldInfo );
                    }
                    ptreeDb.add_child( "fields", ptreeDbFieldInfos );
                }

                std::string jsonFilename = "mydb.dbfeelpp.json";
                std::string jsonOutputPath = (M_dbRepository/jsonFilename).string();
                write_json( jsonOutputPath, ptreeDb );
            }
            this->worldComm().barrier();
        }

    mesh_ptrtype loadMesh()
        {
            CHECK( !M_meshFilename.empty() ) << "no mesh filename in database";
            //int nPartition = Environment::worldComm().size();
            //M_meshFilename = (boost::format("mesh_p%1%.json")%nPartition ).str();
            std::string inputMeshPath = (M_dbRepository/ M_meshFilename).string();
            size_type updateCtx = size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES);
            M_mesh = Feel::loadMesh(_mesh=new mesh_type,_filename=inputMeshPath,
                                    _update=updateCtx );
            return M_mesh;
        }

    void saveMeshParititioned( mesh_ptrtype const& meshSeq, int nPartition )
        {
            std::string filename = (boost::format("mesh_p%1%.json")%nPartition ).str();
            M_meshFilename = filename;

            if ( this->worldComm().isMasterRank() )
            {
                if ( !fs::exists( M_dbRepository ) )
                    fs::create_directories( M_dbRepository );

                CHECK( meshSeq ) << "mesh not init";
                std::string outputMeshPath = (M_dbRepository/ filename).string();
                using io_t = PartitionIO<mesh_type>;
                io_t io( outputMeshPath );
                io.write( partitionMesh( meshSeq, nPartition ) );
            }
            this->saveInfo();
            //this->worldComm().barrier();
        }


    template <typename ElementType>
    void save( double time, std::string const& name, ElementType /*const*/& u )
        {
            std::string nameUsed =  ( name.empty() )? "fefunction" : name;
            // TODO check if time exist
            auto space = u.functionSpace();

            if ( M_fieldInfoInDb.find( nameUsed ) == M_fieldInfoInDb.end() )
                M_fieldInfoInDb[nameUsed] = std::make_tuple( space->basisName(),space->basisOrder()[0],ElementType::nComponents  );

            int timeIndexInDatabase = -1;//M_timeSet.size();
            for ( int k=0;k<M_timeSet.size(); ++ k )
            {
                if ( std::abs( M_timeSet[k]-time ) < 1e-9 )
                {
                    timeIndexInDatabase = k;
                    break;
                }
            }
            if ( timeIndexInDatabase == -1 )
            {
                if ( !M_timeSet.empty() )
                    CHECK ( time > (M_timeSet.back()+1e-9) ) << "invalid new time inserted (must be ordered)\n";
                timeIndexInDatabase = M_timeSet.size();
                M_timeSet.push_back( time );
            }

            if ( this->vectorFileFormat() == "hdf5")
            {
#ifdef FEELPP_HAS_HDF5
                std::string filename = (boost::format("%1%-%2%.h5")%nameUsed %timeIndexInDatabase).str();
                u.saveHDF5( (M_dbRepository / filename ).string() );
#else
                CHECK( false ) << "hdf5 not detected";
#endif
            }
            else if ( this->vectorFileFormat() == "boost")
            {
                std::string filename = (boost::format("%1%-%2%_p%3%")%nameUsed %timeIndexInDatabase %u.worldComm().globalRank() ).str();
                fs::ofstream ofs( M_dbRepository / filename );
                boost::archive::binary_oarchive oa( ofs );
                oa << u;
            }
            this->saveInfo();
        }

    template <typename ElementType>
    void load( size_type timeIndex, std::string const& name, ElementType /*const*/& u )
        {
            std::string nameUsed =  ( name.empty() )? "fefunction" : name;
            if ( this->vectorFileFormat() == "hdf5")
            {
                std::string filename = (boost::format("%1%-%2%.h5")%nameUsed %timeIndex).str();
#ifdef FEELPP_HAS_HDF5
                u.loadHDF5( (M_dbRepository / filename ).string() );
#else
                CHECK( false ) << "hdf5 not detected";
#endif
            }
            else if ( this->vectorFileFormat() == "binary")
            {
                std::string filename = (boost::format("%1%-%2%_p%3%")%nameUsed %timeIndex %u.worldComm().globalRank() ).str();
                fs::ifstream ifs( M_dbRepository / filename );
                boost::archive::binary_iarchive ia( ifs );
                ia >> u;
            }
        }


private :
    WorldComm const& M_worldComm;
    std::string M_name;
    fs::path M_dbRepository;
    std::string M_meshFilename;
    std::string M_vectorFileFormat;
    mesh_ptrtype M_mesh;
    std::vector<double> M_timeSet;
    std::map<std::string,std::tuple<std::string,uint16_type,uint16_type> > M_fieldInfoInDb;

};

} // namespace Feel


#endif // FEELPP_FEELPPDATABASE_HPP
