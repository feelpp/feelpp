/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2022-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelmodels/modelcore/modelmeshes.hpp>

#if 1
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelmesh/remesher.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#else
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelmesh/remesh.hpp>
#endif

namespace Feel
{
namespace FeelModels
{

template <typename IndexType>
ModelMesh<IndexType>::MeshAdaptation::Setup::Setup( ModelMesh const* parentModelMesh, ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
    :
    M_parentModelMesh( parentModelMesh ),
    M_eventEachTimeStep_frequency( 1 ),
    M_eventEachTimeStep_lastExecutionIndex( invalid_v<size_type> )
{
    if ( jarg.contains( "metric" ) )
    {
        auto j_metric = jarg.at( "metric" );
        if ( j_metric.is_string() || j_metric.is_number() )
            M_metric.setExpr( j_metric, mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
        else
            throw std::runtime_error(  fmt::format("invalid metric type {}", j_metric.dump() ) );
    }

    for ( std::string const& keepMarkerKey : { "required_markers" } )
    {
        if ( jarg.contains( keepMarkerKey ) )
        {
            ModelMarkers _markers;
            _markers.setup( jarg.at(keepMarkerKey)/*, indexes*/ );
            M_requiredMarkers.insert( _markers.begin(),_markers.end() );
        }
    }


    if ( jarg.contains( "setup" ) )
    {
        M_remesherSetup = jarg.at( "setup" );
    }

    nl::json events;
    for ( std::string const& eventKey : { "event","events" } )
    {
        if ( jarg.contains( eventKey ) )
        {
            auto const& j_event = jarg.at( eventKey );
            if ( j_event.is_string() )
                events[j_event.template get<std::string>()] = nl::json{};
             else if ( j_event.is_array() )
                for ( auto const& [j_eventkey,j_eventval] : j_event.items() )
                    events[j_eventval.template get<std::string>()] = nl::json{};
             else if ( j_event.is_object() )
                for ( auto const& [j_eventkey,j_eventval] : j_event.items() )
                    events[j_eventkey] = j_eventval;
         }
    }

    for ( auto const& [event,j_eventval] : events.items() )
    {
        if ( event == "after_import" )
            M_executionEvents.insert( Event::Type::after_import );
        else if ( event == "after_init" )
            M_executionEvents.insert( Event::Type::after_init );
        else if ( event == "each_time_step" )
        {
            M_executionEvents.insert( Event::Type::each_time_step );
            if ( j_eventval.contains( "frequency" ) )
            {
                auto const& j_ets_frequency = j_eventval.at( "frequency" );
                if ( j_ets_frequency.is_number_unsigned() )
                    M_eventEachTimeStep_frequency =  j_ets_frequency.template get<int>();
                else if ( j_ets_frequency.is_string() )
                    M_eventEachTimeStep_frequency = std::stoi( j_ets_frequency.template get<std::string>() );
                else
                    throw std::runtime_error( fmt::format("invalid frequency type {}", j_ets_frequency.dump() ) );
            }
            if ( j_eventval.contains( "condition" ) )
            {
                auto const& j_ets_condition = j_eventval.at( "condition" );
                M_eventEachTimeStep_condition.setExpr( j_ets_condition, mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
                if ( !M_eventEachTimeStep_condition.template hasExpr<1,1>() )
                    throw std::runtime_error( fmt::format("invalid condition {}, should be a boolean scalar expr", j_ets_condition.dump() ) );
            }
        }
        else
            throw std::runtime_error( fmt::format("invalid event {}", event ) );
    }
    M_tmpDir = (fs::path(mMeshes.repository().root())/"meshes")/"tmp";
}

template <typename IndexType>
void
ModelMesh<IndexType>::MeshAdaptation::Setup::updateInformationObject( nl::json & jsonInfo ) const
{
     auto [exprStr,compInfo] = M_metric.exprInformations();
     jsonInfo["metric"] = exprStr;
     jsonInfo["required_markers"] = M_requiredMarkers;
     if ( M_remesherSetup.is_object() )
         jsonInfo["remesher_setup"] = M_remesherSetup;

     std::map<typename Event::Type,std::string> mapEventEnumToString = { { Event::Type::never,"never" },
                                                                         { Event::Type::after_import,"after_import" },
                                                                         { Event::Type::after_init,"after_init" },
                                                                         { Event::Type::each_time_step,"each_time_step" } };
     std::vector<std::string> eventsStr;
     for ( auto event : M_executionEvents )
         eventsStr.push_back( mapEventEnumToString.at(event) );
     if ( !eventsStr.empty() )
         jsonInfo["events"] = eventsStr;
}

template <typename IndexType>
tabulate_informations_ptr_t
ModelMesh<IndexType>::MeshAdaptation::Setup::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    Feel::Table tableInfoBasic;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tableInfoBasic, jsonInfo, tabInfoProp );

    if ( jsonInfo.contains("required_markers") )
    {
        Feel::Table tabInfoMarkers = TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at("required_markers"), true );
        if ( tabInfoMarkers.nRow() > 0 )
            tableInfoBasic.add_row( { "required_markers", tabInfoMarkers } );
    }

    if ( jsonInfo.contains("events") )
    {
        Feel::Table tabInfoEvents = TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at("events"), true );
        if ( tabInfoEvents.nRow() > 0 )
            tableInfoBasic.add_row( { "events", tabInfoEvents } );
    }

    tableInfoBasic.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );

    if ( tableInfoBasic.nRow() > 0 )
        tabInfo->add( "", TabulateInformations::New( tableInfoBasic, tabInfoProp ) );

    if ( jsonInfo.contains("remesher_setup") )
    {
        auto const& j_remesher_setup = jsonInfo.at("remesher_setup");
        Feel::Table tableInfoRemesherSetup;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tableInfoRemesherSetup, j_remesher_setup, tabInfoProp );
        if ( tableInfoRemesherSetup.nRow() > 0 )
            tabInfo->add( "remesher_setup", TabulateInformations::New( tableInfoRemesherSetup, tabInfoProp.newByIncreasingVerboseLevel() ) );
    }

    return tabInfo;
}

template class ModelMesh<uint32_type>;


} // namespace FeelModel
} // namespace Feel
