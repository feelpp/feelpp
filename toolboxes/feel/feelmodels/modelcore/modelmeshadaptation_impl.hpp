/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2022-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#pragma once

#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelmesh/remesher.hpp>
#include <feel/feelfilters/loadmesh.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename IndexType>
template <typename MeshType>
std::shared_ptr<MeshType>
ModelMesh<IndexType>::MeshAdaptation::Execute::executeImpl( std::shared_ptr<MeshType> inputMesh )
{
    auto makeMesh = [&]( worldcomm_ptr_t const& wc ) {
        if constexpr ( MeshType::is_order_dynamic )
            return new MeshType(RuntimeOrder{inputMesh->order()}, "", wc);
        else
            return new MeshType(wc);
    };
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
#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )
        auto r = remesher( inputMesh, requiredMarkersElements, requiredMarkersFaces, {}, {}, { {"remesh", M_mas.M_remesherSetup } } );
        r.setMetric( *scalarMetricField );
        return r.execute();
#else
        CHECK( false ) << "no mmg/parmmg support";
        return {};
#endif
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
            auto inputMeshSeq = loadMesh(_mesh=makeMesh(Environment::worldCommSeqPtr()), _savehdf5=0,
                                         _filename=imeshParaPath,
                                         //_update=update_,  TODO test FACE_MINIMAL
                                         _straighten=false );
#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )
            auto r = remesher( inputMeshSeq, requiredMarkersElements, requiredMarkersFaces, {}, {}, { {"remesh", M_mas.M_remesherSetup } } );

            auto VhSeq = scalar_metric_field_type::functionspace_type::New(_mesh=inputMeshSeq );
            auto metFieldSeq = VhSeq->elementPtr();
            metFieldSeq->load(_path=scalarMetricFieldPath,_space_path=spacefileName);
            r.setMetric( *metFieldSeq );

            auto outputMeshSeq = r.execute();
            using io_t = PartitionIO<MeshType>;
            io_t io( omeshParaPath );
            nl::json partConfig = {
                { "partitioner", M_mas.M_parentModelMesh->partitioningSetup().json() }
            };
            io.write( partitionMesh( outputMeshSeq, nPartition, {}, partConfig ) );
#else
            CHECK( false ) << "no mmg/parmmg support";
#endif
        }

        inputMesh->worldComm().barrier();

        auto out = loadMesh(_mesh=makeMesh(inputMesh->worldCommPtr()), _savehdf5=0,
                            _filename=omeshParaPath/*, _update=update_*/ );
        return out;


    }
    //return std::shared_ptr<MeshType>{};
}



} // namespace FeelModels
} // namespace Feel
