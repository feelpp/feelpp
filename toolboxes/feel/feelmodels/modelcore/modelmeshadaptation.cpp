/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
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
ModelMesh<IndexType>::MeshAdaptation::Setup::Setup( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
    :
    M_eventEachTimeStep_frequency( 1 ),
    M_eventEachTimeStep_lastExecutionIndex( invalid_v<size_type> )
{
    if ( jarg.contains( "metric" ) )
    {
        auto j_metric = jarg.at( "metric" );
        if ( j_metric.is_string() || j_metric.is_number() )
            M_metric.setExpr( j_metric, mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
        else
            throw std::runtime_error(  fmt::format("invalid metric type {}", j_metric ) );
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
                    throw std::runtime_error( fmt::format("invalid frequency type {}", j_ets_frequency ) );
            }
            if ( j_eventval.contains( "condition" ) )
            {
                auto const& j_ets_condition = j_eventval.at( "condition" );
                M_eventEachTimeStep_condition.setExpr( j_ets_condition, mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
                if ( !M_eventEachTimeStep_condition.template hasExpr<1,1>() )
                    throw std::runtime_error( fmt::format("invalid condition {}, should be a boolean scalar expr", j_ets_condition ) );
            }
        }
        else
            throw std::runtime_error( fmt::format("invalid event {}", event ) );
    }
    M_tmpDir = (fs::path(mMeshes.repository().root())/"meshes")/"tmp";
}


template <typename IndexType>
template <typename MeshType>
std::shared_ptr<MeshType>
ModelMesh<IndexType>::MeshAdaptation::Execute::executeImpl( std::shared_ptr<MeshType> inputMesh )
{
    std::vector<std::string> requiredMarkersElements, requiredMarkersFaces;
    for ( std::string const&  _m : M_mas.requiredMarkers() )
    {
        if ( inputMesh->hasElementMarker( _m ) )
            requiredMarkersElements.push_back( _m );
        else if ( inputMesh->hasFaceMarker( _m ) )
            requiredMarkersFaces.push_back( _m );
    }

    using scalar_metric_field_type = typename std::decay_t<decltype( unwrap_ptr( Pch<1>( inputMesh ) ) )>::element_type;
    auto scalarMetricField = std::dynamic_pointer_cast<scalar_metric_field_type>( M_scalarMetricField );

    if ( inputMesh->worldComm().localSize() == 1 )
    {
        auto r = remesher( inputMesh, requiredMarkersElements, requiredMarkersFaces );

        r.setMetric( *scalarMetricField );

        // auto metExpr = expr( mas.metricExpr().template expr<1,1>(), se );
        // auto Vh = Pch<1>( oldmesh );
        // auto met = Vh->element();
        // met.on(_range=elements(oldmesh),_expr=metExpr );
        // r.setMetric( met );
        return r.execute();
    }
    else
    {
        if ( inputMesh->worldComm().isMasterRank() )
        {
            if ( !fs::exists( M_mas.M_tmpDir ) )
                fs::create_directories( M_mas.M_tmpDir );
        }
        inputMesh->worldComm().barrier();

        std::string imeshParaPath = (M_mas.M_tmpDir/"mesh_i.json").string();
        std::string omeshParaPath = (M_mas.M_tmpDir/"mesh_o.json").string();
        int nPartition = inputMesh->worldComm().localSize();
        inputMesh->saveHDF5( imeshParaPath );

        std::string scalarMetricFieldPath = (M_mas.M_tmpDir/"scalar_metric_field").string();
        scalarMetricField->save(_path=scalarMetricFieldPath);

        std::string spacefileName = (M_mas.M_tmpDir/"scalar_metric_space").string();
        scalarMetricField->functionSpace()->save( spacefileName );

        if ( inputMesh->worldComm().isMasterRank() )
        {
            auto inputMeshSeq = loadMesh(_mesh=new MeshType{Environment::worldCommSeqPtr()}, _savehdf5=0,
                                         _filename=imeshParaPath,
                                         //_update=update_,  TODO test FACE_MINIMAL
                                         _straighten=false );

            auto r = remesher( inputMeshSeq, requiredMarkersElements, requiredMarkersFaces );
#if 1
            //auto VhSeq = Pch<1>( inputMeshSeq );
            auto VhSeq = scalar_metric_field_type::functionspace_type::New(_mesh=inputMeshSeq );
            auto metFieldSeq = VhSeq->elementPtr();
            metFieldSeq->load(_path=scalarMetricFieldPath,_space_path=spacefileName);
            // //auto metExpr = expr( M_mas.metricExpr().template expr<1,1>(), se );
            // auto metExpr = M_mas.metricExpr().template expr<1,1>();
            // metField->on( _range=elements(inputMeshSeq),_expr=metExpr );
            r.setMetric( *metFieldSeq );
#else
            using scalar_metric_field_type = typename std::decay_t<decltype( unwrap_ptr( Pch<1>( inputMesh ) ) )>::element_type;
            auto m = std::dynamic_pointer_cast<scalar_metric_field_type>( M_scalarMetricField );
            auto Vh = Pch<1>( inputMeshSeq );
            auto opI = opInterpolation(_domainSpace=m->functionSpace(),
                                       _imageSpace=Vh );
            auto matrixInterpolation = opI->matPtr();
            auto metField = Vh->elementPtr();
            matrixInterpolation->multVector( *m, *metField );
#endif
            auto outputMeshSeq = r.execute();
            using io_t = PartitionIO<MeshType>;
            io_t io( omeshParaPath );
            io.write( partitionMesh( outputMeshSeq, nPartition/*, partitionByRange, partconfig*/ ) );
        }

        inputMesh->worldComm().barrier();

        auto out = loadMesh(_mesh=new MeshType{inputMesh->worldCommPtr()}, _savehdf5=0,
                            _filename=omeshParaPath/*, _update=update_*/ );
        return out;


    }
    //return std::shared_ptr<MeshType>{};
}



template class ModelMesh<uint32_type>;

template std::shared_ptr<Mesh<Simplex<2,1>>> ModelMesh<uint32_type>::MeshAdaptation::Execute::executeImpl<Mesh<Simplex<2,1>>>( std::shared_ptr<Mesh<Simplex<2,1>>> );
template std::shared_ptr<Mesh<Simplex<3,1>>> ModelMesh<uint32_type>::MeshAdaptation::Execute::executeImpl<Mesh<Simplex<3,1>>>( std::shared_ptr<Mesh<Simplex<3,1>>> );

} // namespace FeelModel
} // namespace Feel


