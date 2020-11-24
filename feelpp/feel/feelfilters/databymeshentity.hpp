/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 22 Oct 2020

 Copyright (C) 2020 Feel++ Consortium

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
#pragma once

#include <feel/feelcore/json.hpp>
#include <feel/feelcore/hdf5.hpp>
#include <feel/feelvf/databymeshentity.hpp>


namespace Feel
{

template <typename IndexType>
class DataByMeshEntity
{
public :
    using index_type = IndexType;
    using mesh_base_type = MeshBase<index_type>;
    using value_type = double;
    using mapping_id2value_type = std::map<index_type,value_type>;

    DataByMeshEntity( std::shared_ptr<mesh_base_type> const& mesh, ElementsType entityType, mapping_id2value_type const& entityIdToValue )
        :
        M_mesh( mesh ),
        M_entityType( entityType ),
        M_entityIdToValue( entityIdToValue )
        {}

    DataByMeshEntity( DataByMeshEntity const& ) = default;
    DataByMeshEntity( DataByMeshEntity && ) = default;

    std::shared_ptr<mesh_base_type> mesh() const { return M_mesh; }
    ElementsType entityType() const { return M_entityType; }
    mapping_id2value_type const& entityIdToValue() const { return M_entityIdToValue; }

    std::optional<value_type> valueAtEntityIdIfExists( index_type id ) const
        {
            auto itFindVal = M_entityIdToValue.find( id );
            if ( itFindVal != M_entityIdToValue.end() )
                return std::make_optional( itFindVal->second );
            else
                return std::optional<value_type>{};
        }

    void set( ElementsType entityType, mapping_id2value_type && entityIdToValue )
        {
            M_entityType = entityType;
            M_entityIdToValue = std::forward<mapping_id2value_type>( entityIdToValue );
        }

    void setInCache( ElementsType entityType, int keyInCache, mapping_id2value_type && entityIdToValue )
        {
            if ( M_cache_entityIdToValue.find( keyInCache ) == M_cache_entityIdToValue.end() )
                M_cache_entityIdToValue.emplace( std::make_pair( keyInCache, std::forward<mapping_id2value_type>( entityIdToValue ) ) );
            else
                M_cache_entityIdToValue[keyInCache] = std::forward<mapping_id2value_type>( entityIdToValue );
        }
    bool hasInCache( int keyInCache ) const { return M_cache_entityIdToValue.find( keyInCache ) != M_cache_entityIdToValue.end(); }
    void eraseCacheEntry( int keyInCache )
        {
            auto itFind = M_cache_entityIdToValue.find( keyInCache );
            if ( itFind != M_cache_entityIdToValue.end() )
                M_cache_entityIdToValue.erase( itFind );
        }
    void keepInCache( std::vector<int> const& keyIds )
        {
            std::set<int> keyIdsToErase;
            for ( auto const& [keyId,data] : M_cache_entityIdToValue )
                if ( std::find( keyIds.begin(), keyIds.end(), keyId ) == keyIds.end() )
                    keyIdsToErase.insert( keyId );

            for ( int keyId : keyIdsToErase )
                this->eraseCacheEntry( keyId );
        }
    void setAtThisTimeFromCache( double time, std::map<int,double> const& mapTimeIdToTimeValue )
        {
            if ( mapTimeIdToTimeValue.size() == 1 )
            {
                auto itFindData = M_cache_entityIdToValue.find( mapTimeIdToTimeValue.begin()->first );
                CHECK( itFindData != M_cache_entityIdToValue.end() ) << "not in cache";
                M_entityIdToValue = itFindData->second;
            }
            else if ( mapTimeIdToTimeValue.size() == 2 )
            {
                double time1 = mapTimeIdToTimeValue.begin()->second;
                double time2 = (++mapTimeIdToTimeValue.begin())->second;
                auto itFindData1 = M_cache_entityIdToValue.find( mapTimeIdToTimeValue.begin()->first );
                CHECK( itFindData1 != M_cache_entityIdToValue.end() ) << "not in cache";
                auto itFindData2 = M_cache_entityIdToValue.find( (++mapTimeIdToTimeValue.begin())->first );
                CHECK( itFindData2 != M_cache_entityIdToValue.end() ) << "not in cache";
                auto const& entityIdToValue1 = itFindData1->second;
                auto const& entityIdToValue2 = itFindData2->second;
                CHECK( entityIdToValue1.size() == entityIdToValue2.size() ) << "something wrong";
                for ( auto it1 = entityIdToValue1.begin(), it2 = entityIdToValue2.begin() ; it1 != entityIdToValue1.end() ; ++it1,++it2 )
                {
                    index_type eid1 = it1->first, eid2 = it2->first;
                    CHECK( eid1 == eid2 ) << "something wrong";
                    value_type v1 = it1->second, v2 = it2->second;
                    double a = (v2 - v1)/(time2-time1);
                    double b = v1 - a*time1;
                    M_entityIdToValue[eid1] = a*time+b;
                }
            }
            else
                CHECK( false ) << "something wrong";
        }

private :
    std::shared_ptr<mesh_base_type> M_mesh;
    ElementsType M_entityType;
    mapping_id2value_type M_entityIdToValue;
    std::map<int,mapping_id2value_type> M_cache_entityIdToValue;
};

template <typename IndexType>
class CollectionOfDataByMeshEntity : protected std::map<std::string,DataByMeshEntity<IndexType>>
{
public :
    using index_type = IndexType;
    using data_by_mesh_entity_type = DataByMeshEntity<index_type>;
    using mapping_id2value_type = typename data_by_mesh_entity_type::mapping_id2value_type;
    using mesh_base_type = typename data_by_mesh_entity_type::mesh_base_type;

    CollectionOfDataByMeshEntity() = default;

    template <typename MeshType>
    CollectionOfDataByMeshEntity( std::shared_ptr<MeshType> const& mesh )
        :
        M_mesh( mesh )
        {}

    CollectionOfDataByMeshEntity( CollectionOfDataByMeshEntity const& ) = default;
    CollectionOfDataByMeshEntity( CollectionOfDataByMeshEntity && ) = default;

    //! setup some infos from json
    void setup( pt::ptree const& p );

    template <typename MeshType>
    void setMesh( std::shared_ptr<MeshType> const& mesh ) { M_mesh = mesh; }

    //! update data if mesh has been setup
    template <typename MeshType>
    void updateForUse();

    //! import data
    template <typename MeshType>
        void import( pt::ptree const& p );

    //! get all data names
    std::set<std::string> dataNames() const
    {
        std::set<std::string> res;
        for ( auto const& [name,data] : *this )
            res.insert( name );
        return res;
    }

    //! return true if has data called \dataName
    bool hasData( std::string const& dataName ) const { return this->find( dataName ) != this->end(); }

    //! return data by mesh entity from a name
    data_by_mesh_entity_type const& get( std::string const& dataName ) const
    {
        auto itFindData = this->find( dataName );
        CHECK( itFindData != this->end() ) << "data not found : " << dataName;
        return itFindData->second;
    }

    //! register data by mesh entity with a nam, entity type and values at entity
    void setData( std::string const& dataName, ElementsType entityType, mapping_id2value_type && data )
    {
        if ( this->hasData( dataName ) )
            this->at( dataName ).set( entityType, std::forward<mapping_id2value_type>( data ) );
        else
            this->emplace( std::make_pair( dataName, data_by_mesh_entity_type( M_mesh, entityType, std::forward<mapping_id2value_type>( data ) ) ) );
    }


    //! return true if data and mesh entity are related by a mapping
    bool hasMappingDataToMesh( std::string const& entity ) const { return M_mappingDataToMesh.find( this->elementsTypeMap( entity ) ) != M_mappingDataToMesh.end(); }

    //! return true if data and mesh entity are related by a mapping
    bool hasEmptyMappingDataToMesh( std::string const& entity ) const
    {
        auto itFindMapping = M_mappingDataToMesh.find( this->elementsTypeMap( entity ) );
        if ( itFindMapping == M_mappingDataToMesh.end() )
            return false;
        return itFindMapping->second.empty();
    }

    //! return symbols expression
    auto symbolsExpr( std::string const& prefix_symbol = "data_by_mesh_entity" ) const
    {
        using _expr_type = std::decay_t< decltype( Feel::vf::dataByMeshEntityExpr( this->begin()->second ) ) >;
        symbol_expression_t<_expr_type> seData;
        for ( auto const& [name,data] : *this )
            seData.add( prefixvm( prefix_symbol, name , "_" ), Feel::vf::dataByMeshEntityExpr( data ) );
        return Feel::vf::symbolsExpr( seData );
    }

    //! return the time set
    std::vector<double> const& timeSet() const { return M_timeSet; }

    //! update data from a time
    void updateTime( double time )
    {
        if ( M_timeSet.empty() )
            return;
        if ( M_previousTimeUpdated && (*M_previousTimeUpdated == time) )
            return;

        this->importData( time );

        M_previousTimeUpdated = time;
    }

  private :

    std::vector<index_type> const& dataIdToMeshIdByMapping( std::string const& entity, index_type id ) const
    {
        auto itFindMapping = M_mappingDataToMesh.find( this->elementsTypeMap( entity ) );
        CHECK( itFindMapping != M_mappingDataToMesh.end() ) << "mapping not found";
        auto const& mappingOnEntity = itFindMapping->second;
        auto itFindId = mappingOnEntity.find( id );
        //CHECK( itFindId != mappingOnEntity.end() ) << "id not found in mapping";
        if ( itFindId != mappingOnEntity.end() )
            return itFindId->second;
        else
            return M_dummyEmptyData;
    }

    static ElementsType elementsTypeMap( std::string const& entity )
    {
        if ( entity == "elements" )
            return ElementsType::MESH_ELEMENTS;
        else if ( entity == "faces" )
            return ElementsType::MESH_FACES;
        else if ( entity == "edges" )
            return ElementsType::MESH_EDGES;
        else if ( entity == "points" )
            return ElementsType::MESH_POINTS;
        else
            CHECK( false ) << "invalid entity : " << entity ;
    }

    template <typename MeshType>
    void updateMapping( std::string const& filename );

    void importData( double time = 0 );

    void setDataInCache( std::string const& dataName, ElementsType entityType, int keyInCache, mapping_id2value_type && data )
    {
        if ( !this->hasData( dataName ) )
            this->emplace( std::make_pair( dataName, data_by_mesh_entity_type( M_mesh, entityType, mapping_id2value_type{} ) ) );
        this->at( dataName ).setInCache( entityType, keyInCache, std::forward<mapping_id2value_type>( data ) );
    }

    bool hasDataInCache( std::string const& dataName, int keyInCache ) const
    {
        if ( !this->hasData( dataName ) )
            return false;
        return this->at( dataName ).hasInCache( keyInCache );
    }

    void keepDataInCache( std::string const& dataName, std::vector<int> const& keyIds )
    {
        if ( !this->hasData( dataName ) )
            return;
        this->at( dataName ).keepInCache( keyIds );
    }

    void setDataAtThisTimeFromCache( std::string const& dataName, double time, std::map<int,double> const& mapTimeIdToTimeValue )
    {
        CHECK( this->hasData( dataName ) ) << "data should be init";
        this->at( dataName ).setAtThisTimeFromCache( time, mapTimeIdToTimeValue );
    }

private :
    std::string M_dataFilename, M_mappingFilename;
    std::set<std::string> M_fieldsRequired;

    std::shared_ptr<mesh_base_type> M_mesh;
    std::map<ElementsType,std::map<index_type,std::vector<index_type> > > M_mappingDataToMesh;

    std::string M_entityOfData;
    std::set<std::string> M_fieldsToLoad;
    std::string M_h5file;
    std::vector<double> M_timeSet;
    std::vector<index_type> M_dummyEmptyData;

    std::optional<double> M_previousTimeUpdated;
};


template <typename IndexType>
void
CollectionOfDataByMeshEntity<IndexType>::setup( pt::ptree const& p )
{
    if ( auto pfields = p.get_child_optional("fields") )
    {
        if ( pfields->empty() ) // value case
            M_fieldsRequired.insert( pfields->get_value<std::string>() );
        else // array case
        {
            for ( auto const& item : *pfields )
            {
                CHECK( item.first.empty() ) << "should be an array, not a subtree";
                std::string const& fieldName = item.second.template get_value<std::string>();
                M_fieldsRequired.insert( fieldName );
            }
        }
    }

    if ( auto mappingFilenameOpt = p.template get_optional<std::string>( "mapping_data_to_mesh_nodes" ) )
        M_mappingFilename = Environment::expand( *mappingFilenameOpt );

    if ( auto filenameOpt = p.template get_optional<std::string>( "filename" ) )
       M_dataFilename = Environment::expand( *filenameOpt );
}

template <typename IndexType>
template <typename MeshType>
void
CollectionOfDataByMeshEntity<IndexType>::import( pt::ptree const& p )
{
    this->setup( p );
    this->updateForUse<MeshType>();
}

template <typename IndexType>
template <typename MeshType>
void
CollectionOfDataByMeshEntity<IndexType>::updateForUse()
{
    if ( !M_mesh || !std::dynamic_pointer_cast<MeshType>( M_mesh ) )
        return;

    if ( !M_mappingFilename.empty() )
    {
        this->updateMapping<MeshType>( M_mappingFilename );
    }

    CHECK( !M_dataFilename.empty() ) << "filename not given";
    //M_dataFilename = Environment::expand( *filenameOpt );
    fs::path filedir = fs::path(M_dataFilename).parent_path();
    //std::cout << "filename=" << M_dataFilename << std::endl;
    std::ifstream ifs( M_dataFilename );
    CHECK( ifs ) << "open file fails : " << M_dataFilename;
    nl::json jsonData = nl::json::parse( ifs );

    //std::string entityOfData;
    if ( jsonData.contains("entity") )
        M_entityOfData = jsonData.at( "entity" ).template get<std::string>();
    else
        CHECK( false ) << "entity not given";

    std::string h5file = jsonData.at("filename").template get<std::string>();
    M_h5file = (filedir/h5file).string();
    //std::cout << "M_h5file=" << M_h5file << std::endl;

    //std::set<std::string> fieldsToLoad;
    if ( jsonData.contains("fields") )
        for ( auto const& el : jsonData.at("fields").items() )
        {
            std::string curField = el.value().template get<std::string>();
            if ( M_fieldsRequired.empty() || M_fieldsRequired.find( curField ) != M_fieldsRequired.end() )
                M_fieldsToLoad.insert( curField );
        }

    if ( jsonData.contains("time-set") )
        for ( auto const& el : jsonData.at("time-set").items() )
            M_timeSet.push_back( el.value().template get<double>() );

    // for ( double t :  M_timeSet )
    //     std::cout << "timeSet t " << t << std::endl;
    // for ( std::string const& field : fieldsToLoad )
    //     std::cout << "field="<< field << std::endl;

    if ( M_timeSet.empty() )
        this->importData();
}



template <typename IndexType>
void
CollectionOfDataByMeshEntity<IndexType>::importData( double time )
{
    if ( M_fieldsToLoad.empty() )
        return;

    std::vector<int> timeIndex;
    if ( !M_timeSet.empty() )
    {
        CHECK( time >= M_timeSet.front() && time <= M_timeSet.back() ) << "time " << time << " is out of range : [" << M_timeSet.front() << "," << M_timeSet.back() <<"]";

        for (int k=0;k<M_timeSet.size();++k)
        {
            if ( time == M_timeSet[k] )
            {
                timeIndex = { k };
                break;
            }
            CHECK( k < M_timeSet.size() ) << "something wrong";
            if ( time > M_timeSet[k] && time < M_timeSet[k+1] )
            {
                timeIndex = { k,k+1 };
                break;
            }
        }
    }
    else
        timeIndex = {-1};

    bool hasMappingDataToMesh = this->hasMappingDataToMesh( M_entityOfData );
    if ( hasMappingDataToMesh && this->hasEmptyMappingDataToMesh( M_entityOfData ) ) // needn't to read hdf5 file
    {
        for ( std::string const& field : M_fieldsToLoad )
            this->setData( field, this->elementsTypeMap( M_entityOfData ), std::map<index_type,double>{} );
    }
    else
    {
        std::vector<int> timeIndexToLoad;
        if ( M_timeSet.empty() )
            timeIndexToLoad = timeIndex;
        else
        {
            for ( std::string const& dataName : M_fieldsToLoad )
                this->keepDataInCache( dataName, timeIndex );
            std::string const& firstDataName = *M_fieldsToLoad.begin();
            for ( int tId : timeIndex )
                if ( !this->hasDataInCache( firstDataName, tId ) )
                    timeIndexToLoad.push_back( tId );
        }

        if ( !timeIndexToLoad.empty() )
        {
            HDF5 hdf5;
            hdf5.openFile( M_h5file, Environment::worldCommSeq(), true );

            std::string tableNameIds = M_entityOfData+"_ids";
            hsize_t dimsTable[2];
            hsize_t offsetElt[2] = { 0,0 };
            hdf5.openTable( tableNameIds, dimsTable );
            //std::cout << "dimsTable=" << dimsTable[0] << " and " << dimsTable[1] << std::endl;
            std::vector<unsigned int> entityIds( dimsTable[0]*dimsTable[1] );
            hdf5.read( tableNameIds, H5T_NATIVE_UINT, dimsTable, offsetElt, entityIds.data() );
            hdf5.closeTable( tableNameIds );

            for ( int tId : timeIndexToLoad )
            {
                std::string groupName =  (tId < 0) ? "/data" : (boost::format("/data_%1%")%tId).str();
                for ( std::string const& field : M_fieldsToLoad )
                {
                    std::string tableNameField =  (boost::format("%1%/%2%")%groupName %field).str();
                    hdf5.openTable( tableNameField, dimsTable );
                    std::vector<double> fieldValues( dimsTable[0]*dimsTable[1] );
                    hdf5.read( tableNameField, H5T_NATIVE_DOUBLE, dimsTable, offsetElt, fieldValues.data() );
                    hdf5.closeTable( tableNameField );

                    std::map<index_type,double> entityIdToValue;
                    if ( hasMappingDataToMesh )
                    {
                        for (int k=0;k<entityIds.size();++k )
                        {
                            for ( auto const& meshEntityId : dataIdToMeshIdByMapping( M_entityOfData,entityIds[k]) )
                                entityIdToValue[meshEntityId] = fieldValues[k];
                        }
                    }
                    else
                    {
                        for (int k=0;k<entityIds.size();++k )
                            entityIdToValue[entityIds[k]] = fieldValues[k];
                    }
                    if ( tId < 0 )
                        this->setData( field, this->elementsTypeMap( M_entityOfData ), std::move(entityIdToValue) );
                    else
                        this->setDataInCache( field, this->elementsTypeMap( M_entityOfData ), tId, std::move(entityIdToValue) );
                }
            }

            hdf5.closeFile();
        } //  if ( !timeIndexToLoad.empty() )

        // set data at this time by apply linear interpolation between values at time index
        if ( !M_timeSet.empty() )
        {
            std::map<int,double> mapTimeIdToTimeValue;
            for ( int tId : timeIndex )
                mapTimeIdToTimeValue[tId] = M_timeSet[tId];
            for ( std::string const& dataName : M_fieldsToLoad )
                this->setDataAtThisTimeFromCache( dataName, time, mapTimeIdToTimeValue );
        }
    }

}


template <typename IndexType>
template <typename MeshType>
void
CollectionOfDataByMeshEntity<IndexType>::updateMapping( std::string const& mappingFilename )
{
    LOG(INFO) << "load mapping from : " << mappingFilename;
    std::ifstream ifs( mappingFilename );
    CHECK( ifs ) << "open file fails : " << mappingFilename;
    nl::json mappingDataToMeshNodes = nl::json::parse( ifs );

    std::set<std::string> entities;
    if ( mappingDataToMeshNodes.contains("entities") )
        for ( auto const& el : mappingDataToMeshNodes.at("entities").items() )
            entities.insert( el.value().template get<std::string>() );

    for ( std::string const& entity : entities )
    {
        LOG(INFO) << "get mapping on topological entity : " << entity;
        std::set<std::string> markers;
        if ( mappingDataToMeshNodes.contains( "markers_"+entity ) )
            for ( auto const& el : mappingDataToMeshNodes.at( "markers_"+entity ).items() )
                markers.insert( el.value().template get<std::string>() );

        std::map<std::vector<index_type>,index_type> mapPointIdsToEdgeId;

        auto mesh = std::dynamic_pointer_cast<MeshType>( M_mesh );
        CHECK( mesh ) << "invalid mesh_type";

        if ( entity == "edges" )
        {
            if constexpr ( MeshType::nDim == 3 )
            {
                auto rangeEdge = markers.empty()? edges(mesh) : markededges(mesh,markers);
                std::vector<index_type> pids(2);
                for ( auto const& eWrap : rangeEdge )
                {
                    auto const& e = unwrap_ref( eWrap );
                    pids[0] = e.point(0).id();
                    pids[1] = e.point(1).id();
                    std::sort( pids.begin(), pids.end() );
                    mapPointIdsToEdgeId[ pids ] = e.id();
                }
            }
            else
                CHECK( false ) << "edges entity only in 3d";
        }
        else
            CHECK( false ) << "TODO :" << entity;

        if ( mapPointIdsToEdgeId.empty() )
        {
            M_mappingDataToMesh[this->elementsTypeMap( entity )];
            continue;
        }

        std::vector<index_type> key;
        std::vector<index_type> pids;
        CHECK( mappingDataToMeshNodes.contains( "mapping_"+entity ) )  << "missing mapping";
        for (auto& [dataEntityId, pointIds] : mappingDataToMeshNodes[ "mapping_"+entity].items())
        {
            int nPointOnCurrentEntity = pointIds.size();
            pids.resize( nPointOnCurrentEntity );
            //CHECK( pointIds.size() == 2 ) << "must be a pair : " << pointIds.size();
            for ( int p=0;p<nPointOnCurrentEntity;++p )
                pids[p] = pointIds[p].get<index_type>();
            std::sort( pids.begin(), pids.end() );
            auto itFindEdgeInMesh = mapPointIdsToEdgeId.find( pids );
            //CHECK( itFindEdgeInMesh != mapPointIdsToEdgeId.end() ) <<  "entity not found : " << dataEntityId << " -> " <<  pids;
            if ( itFindEdgeInMesh == mapPointIdsToEdgeId.end() )
                continue;

            M_mappingDataToMesh[this->elementsTypeMap( entity )][ std::stoi( dataEntityId ) ].push_back(  itFindEdgeInMesh->second );
        }
    }
}

} // namespace Feel
