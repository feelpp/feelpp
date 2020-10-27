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

#include <feel/feelvf/databymeshentity.hpp>

namespace Feel
{

template <typename MeshType>
class DataByMeshEntity
{
public :
    using mesh_type = MeshType;
    using index_type = typename mesh_type::index_type;
    using value_type = double;
    using mapping_id2value_type = std::map<index_type,value_type>;

    DataByMeshEntity( std::shared_ptr<mesh_type> const& mesh, ElementsType entityType, mapping_id2value_type const& entityIdToValue )
        :
        M_mesh( mesh ),
        M_entityType( entityType ),
        M_entityIdToValue( entityIdToValue )
        {}

    DataByMeshEntity( DataByMeshEntity const& ) = default;
    DataByMeshEntity( DataByMeshEntity && ) = default;

    std::shared_ptr<mesh_type> const& mesh() const { return M_mesh; }
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
private :
    std::shared_ptr<mesh_type> M_mesh;
    ElementsType M_entityType;
    mapping_id2value_type M_entityIdToValue;
};

template <typename MeshType>
class CollectionOfDataByMeshEntity : protected std::map<std::string,DataByMeshEntity<MeshType>>
{
public :
    using mesh_type = MeshType;
    using index_type = typename mesh_type::index_type;
    using data_by_mesh_entity_type = DataByMeshEntity<MeshType>;

    CollectionOfDataByMeshEntity( std::shared_ptr<mesh_type> const& mesh )
        :
        M_mesh( mesh )
        {}

    //! import data
    void import( pt::ptree const& p );

    //! get all data names
    std::set<std::string> dataNames() const
    {
        std::set<std::string> res;
        for ( auto const& [name,data] : *this )
            res.insert( name );
        return res;
    }

    //! return data by mesh entity from a name
    data_by_mesh_entity_type const& get( std::string const& dataName ) const
    {
        auto itFindData = this->find( dataName );
        CHECK( itFindData != this->end() ) << "data not found : " << dataName;
        return itFindData->second;
    }

    //! register data by mesh entity with a nam, entity type and values at entity
    void add( std::string const& dataName, ElementsType entityType, typename data_by_mesh_entity_type::mapping_id2value_type && data )
    {
        this->emplace( std::make_pair( dataName, data_by_mesh_entity_type( M_mesh, entityType, std::forward<typename data_by_mesh_entity_type::mapping_id2value_type>( data ) ) ) );
    }

    //! return true if data and mesh entity are related by a mapping
    bool hasMappingDataToMesh( std::string const& entity ) const { return M_mappingDataToMesh.find( this->elementsTypeMap( entity ) ) != M_mappingDataToMesh.end(); }

    //! return symbols expression
    auto symbolsExpr( std::string const& prefix_symbol = "data_by_mesh_entity" ) const
    {
        using _expr_type = std::decay_t< decltype( Feel::vf::dataByMeshEntityExpr( this->begin()->second ) ) >;
        symbol_expression_t<_expr_type> seData;
        for ( auto const& [name,data] : *this )
            seData.add( prefixvm( prefix_symbol, name , "_" ), Feel::vf::dataByMeshEntityExpr( data ) );
        return Feel::vf::symbolsExpr( seData );
    }


  private :
    std::vector<index_type> const& dataIdToMeshIdByMapping( std::string const& entity, index_type id ) const
    {
        auto itFindMapping = M_mappingDataToMesh.find( this->elementsTypeMap( entity ) );
        CHECK( itFindMapping != M_mappingDataToMesh.end() ) << "mapping not found";
        auto const& mappingOnEntity = itFindMapping->second;
        auto itFindId = mappingOnEntity.find( id );
        CHECK( itFindId != mappingOnEntity.end() ) << "id not found in mapping";
        return itFindId->second;
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

    void updateMapping( std::string const& filename );
private :
    std::shared_ptr<mesh_type> M_mesh;
    std::map<ElementsType,std::map<index_type,std::vector<index_type> > > M_mappingDataToMesh;
};

template <typename MeshType>
void
CollectionOfDataByMeshEntity<MeshType>::import( pt::ptree const& p )
{

    std::set<std::string> fieldsRequired;
    if ( auto pfields = p.get_child_optional("fields") )
    {
        if ( pfields->empty() ) // value case
            fieldsRequired.insert( pfields->get_value<std::string>() );
        else // array case
        {
            for ( auto const& item : *pfields )
            {
                CHECK( item.first.empty() ) << "should be an array, not a subtree";
                std::string const& fieldName = item.second.template get_value<std::string>();
                fieldsRequired.insert( fieldName );
            }
        }
    }

    for  ( std::string const& f: fieldsRequired )
        std::cout << "FIELDS : "<< f << std::endl;

   if ( auto mappingFilenameOpt = p.template get_optional<std::string>( "mapping_data_to_mesh_nodes" ) )
   {
       this->updateMapping( Environment::expand( *mappingFilenameOpt ) );
   }

   if ( auto filenameOpt = p.template get_optional<std::string>( "filename" ) )
   {
       std::string filename = Environment::expand( *filenameOpt );
       fs::path filedir = fs::path(filename).parent_path();
       //std::cout << "filename=" << filename << std::endl;
       std::ifstream ifs( filename );
       CHECK( ifs ) << "open file fails : " << filename;
       nl::json jsonData = nl::json::parse( ifs );

       std::string entityOfData;
       if ( jsonData.contains("entity") )
           entityOfData = jsonData.at( "entity" ).template get<std::string>();
       else
           CHECK( false ) << "entity not given";

       std::string h5file = jsonData.at("filename").template get<std::string>();
       //std::cout << "h5file=" << h5file << std::endl;

       std::set<std::string> fieldsToLoad;
       if ( jsonData.contains("fields") )
           for ( auto const& el : jsonData.at("fields").items() )
           {
               std::string curField = el.value().template get<std::string>();
               if ( fieldsRequired.empty() || fieldsRequired.find( curField ) != fieldsRequired.end() )
                   fieldsToLoad.insert( curField );
           }

       for ( std::string const& field : fieldsToLoad )
           std::cout << "field="<< field << std::endl;

       HDF5 hdf5;
       hdf5.openFile( (filedir/h5file).string(), Environment::worldCommSeq(), true );

       std::string tableNameIds = entityOfData+"_ids";
       hsize_t dimsTable[2];
       hsize_t offsetElt[2] = { 0,0 };
       hdf5.openTable( tableNameIds, dimsTable );
       //std::cout << "dimsTable=" << dimsTable[0] << " and " << dimsTable[1] << std::endl;
       std::vector<unsigned int> entityIds( dimsTable[0]*dimsTable[1] );
       hdf5.read( tableNameIds, H5T_NATIVE_UINT, dimsTable, offsetElt, entityIds.data() );
       hdf5.closeTable( tableNameIds );

       bool hasMappingDataToMesh = this->hasMappingDataToMesh( entityOfData );
       std::string groupName = "/data";
       for ( std::string const& field : fieldsToLoad )
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
                   for ( auto const& meshEntityId : dataIdToMeshIdByMapping( entityOfData,entityIds[k]) )
                       entityIdToValue[meshEntityId] = fieldValues[k];
               }
           }
           else
           {
               for (int k=0;k<entityIds.size();++k )
                   entityIdToValue[entityIds[k]] = fieldValues[k];
           }
           this->add( field, this->elementsTypeMap( entityOfData ), std::move(entityIdToValue) );
       }

       hdf5.closeFile();

   }
}

template <typename MeshType>
void
CollectionOfDataByMeshEntity<MeshType>::updateMapping( std::string const& mappingFilename )
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

        if ( entity == "edges" )
        {
            auto markersEdge = std::set<std::string>( markers );
            auto rangeEdge = markededges(M_mesh,markersEdge);
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
            CHECK( false ) << "TODO";

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
            CHECK( itFindEdgeInMesh != mapPointIdsToEdgeId.end() ) <<  "entity not found : " << dataEntityId << " -> " <<  pids;

            M_mappingDataToMesh[this->elementsTypeMap( entity )][ std::stoi( dataEntityId ) ].push_back(  itFindEdgeInMesh->second );
        }
    }
}

} // namespace Feel
