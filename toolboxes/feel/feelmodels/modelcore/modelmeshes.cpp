/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2020-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelmodels/modelcore/modelmeshes.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/analysis.hpp>
#include <feel/feells/distancetorange.hpp>


namespace Feel
{
namespace FeelModels
{


template <typename IndexType>
ModelMeshCommon<IndexType>::ImportConfig::ImportConfig( ModelMeshes<IndexType> const& mMeshes )
    :
    M_generatePartitioning( boption(_prefix=mMeshes.prefix(),_name="gmsh.partition",_vm=mMeshes.clovm()) ),
    M_numberOfPartition( mMeshes.worldComm().localSize() ),
    M_meshSize( doption(_prefix=mMeshes.prefix(),_name="gmsh.hsize",_vm=mMeshes.clovm()) ),
    M_straightenMesh( boption(_prefix=mMeshes.prefix(),_name="gmsh.straighten",_vm=mMeshes.clovm()) ),
    M_meshComponents( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ),
    M_loadByMasterRankOnly( false )
{
    if ( mMeshes.clovm().count( prefixvm(mMeshes.prefix(),"mesh.filename").c_str() ) )
        M_inputFilename = Environment::expand( soption(_prefix=mMeshes.prefix(),_name="mesh.filename",_vm=mMeshes.clovm()) );
}

template <typename IndexType>
void
ModelMeshCommon<IndexType>::ImportConfig::setup( nl::json const& jarg, ModelMeshes<IndexType> const& mMeshes )
{
    if ( jarg.contains("filename") )
    {
        auto const& j_filename = jarg.at("filename");
        if ( j_filename.is_string() )
            M_inputFilename = Environment::expand( j_filename.template get<std::string>() );
    }
    if ( jarg.contains("partition") )
    {
        auto const& j_partition = jarg.at("partition");
        if ( j_partition.is_boolean() )
            M_generatePartitioning = j_partition.template get<bool>();
        else if ( j_partition.is_number_integer() )
            M_generatePartitioning = j_partition.template get<int>() > 0;
        else if ( j_partition.is_string() )
            M_generatePartitioning = boost::lexical_cast<bool>( j_partition.template get<std::string>() );
    }
    if ( jarg.contains("number-of-partition") )
    {
        auto const& j_nparts = jarg.at("number-of-partition");
        if ( j_nparts.is_number_integer() )
            M_numberOfPartition = j_nparts.template get<int>();
        else if ( j_nparts.is_string() )
            M_numberOfPartition = std::stoi( j_nparts.template get<std::string>() );
        CHECK( M_numberOfPartition > 0 && M_numberOfPartition <= mMeshes.worldComm().localSize() ) << "invalid number of partition : " << M_numberOfPartition;
    }
    if ( jarg.contains("hsize") )
    {
        auto const& j_hsize = jarg.at("hsize");
        if ( j_hsize.is_number() )
            M_meshSize = j_hsize.template get<double>();
        else if ( j_hsize.is_string() )
            M_meshSize = std::stod( j_hsize.template get<std::string>() );
    }
}

template <typename IndexType>
void
ModelMeshCommon<IndexType>::ImportConfig::setupInputMeshFilenameWithoutApplyPartitioning( std::string const& filename )
{
    M_inputFilename = filename;
    M_generatePartitioning = false;
}

template <typename IndexType>
void
ModelMeshCommon<IndexType>::ImportConfig::setupSequentialAndLoadByMasterRankOnly()
{
    M_generatePartitioning = false;
    M_numberOfPartition = 1;
    M_loadByMasterRankOnly = true;
}


template <typename IndexType>
void
ModelMeshCommon<IndexType>::ImportConfig::setupFromSubmesh( mesh_base_ptrtype m, std::set<std::string> const& markers )
{
    int codim = 0;
    // M_inputMesh = m;
    // M_submeshCoDimMarkers = { codim, markers };
    // std::make_optional<std::tuple<mesh_base_ptrtype,int,std::set<std::string>>>( )
    M_importFromSubmeshData = std::make_optional<std::tuple<mesh_base_ptrtype,int,std::set<std::string>>>( m, codim, markers );
        //std::make_tuple< m,  codim, markers >
    // TODO VINCENT
}

template <typename IndexType>
void
ModelMeshCommon<IndexType>::ImportConfig::updateForUse( ModelMeshes<IndexType> const& mMeshes )
{
    if ( M_inputFilename.empty() )
        return;
    std::string meshfile = M_inputFilename;
    RemoteData rdTool( M_inputFilename, mMeshes.worldCommPtr() );
    if ( rdTool.canDownload() )
    {
        auto dowloadedData = rdTool.download( (fs::path(Environment::downloadsRepository())/fs::path(mMeshes.prefix())/fs::path("meshes")).string() );
        CHECK( dowloadedData.size() > 0 ) << "no data download";
        meshfile = dowloadedData[0];
        if ( dowloadedData.size() == 2 )
        {
            if ( fs::path( dowloadedData[0] ).extension() == ".h5" && fs::path( dowloadedData[1] ).extension() == ".json" )
                meshfile = dowloadedData[1];
        }
    }
    meshfile = fs::canonical( fs::path( meshfile ) ).string();
    if ( fs::path( meshfile ).extension() == ".geo" )
        M_geoFilename = meshfile;
    else
        M_meshFilename = meshfile;

    if ( this->hasGeoFilename() && M_numberOfPartition > 1 )
        M_generatePartitioning = true;
}

template <typename IndexType>
void
ModelMeshCommon<IndexType>::ImportConfig::updateInformationObject( nl::json & p ) const
{
    if ( this->hasMeshFilename() )
    {
        p.emplace( "mesh-filename", M_meshFilename );
    }
    else
    {
        p.emplace( "geo-filename", M_geoFilename );
        p.emplace( "hsize", M_meshSize );
    }

    p.emplace( "generate-partitioning", M_generatePartitioning );
    if ( M_generatePartitioning )
        p.emplace( "number-of-partition", M_numberOfPartition );
}

template <typename IndexType>
tabulate_informations_ptr_t
ModelMeshCommon<IndexType>::ImportConfig::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    if ( jsonInfo.contains("mesh-filename") )
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "mesh-filename" } );
    else
        TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "geo-filename","hsize" } );

    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { "generate-partitioning", "number-of-partition" } );

    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );

    return TabulateInformations::New( tabInfo );
}


template <typename IndexType>
ModelMesh<IndexType>::FieldsSetup::FieldsSetup( std::string const& name, ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
    :
    M_name( name )
{
    if ( jarg.contains("basis") )
        M_basis = Environment::expand( jarg.at("basis") );
    else
        throw std::runtime_error( "basis required" );

    if ( jarg.contains("parts") )
    {
        auto const& j_parts = jarg.at("parts");
        if ( j_parts.is_array() )
        {
            for ( auto const& [j_partkey,j_partval] : j_parts.items() )
            {
                auto partSetup = PartSetup::create( mMeshes, j_partval );
                if ( partSetup )
                    M_parts.push_back( std::move(*partSetup) );
            }
        }
        else if ( j_parts.is_object() )
        {
            auto partSetup = PartSetup::create( mMeshes, j_parts );
            if ( partSetup )
                M_parts.push_back( std::move(*partSetup) );
        }
    }
    else
    {
        auto partSetup = PartSetup::create( mMeshes, jarg );
        if ( partSetup )
            M_parts.push_back( std::move(*partSetup) );
    }
}

template <typename IndexType>
std::optional<typename ModelMesh<IndexType>::FieldsSetup::PartSetup>
ModelMesh<IndexType>::FieldsSetup::PartSetup::create( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
{
    PartSetup ret;
    if ( jarg.contains("markers") )
        ret.M_markers.setup( jarg.at("markers") );
    if ( jarg.contains("expr") )
    {
        ModelExpression mexpr;
        mexpr.setExpr( jarg.at("expr"), mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
        if ( mexpr.hasAtLeastOneExpr() )
            ret.M_mexpr = mexpr;
    }
    if ( jarg.contains("filename") )
        ret.M_filename = Environment::expand( jarg.at("filename") );

    if ( jarg.contains("space-filename") )
        ret.M_spaceFilename = Environment::expand( jarg.at("space-filename") );

    return ret;
}


template <typename IndexType>
ModelMesh<IndexType>::DistanceToRangeSetup::DistanceToRangeSetup( std::string const& name, ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
    :
    M_name( name )
{
    auto itFind = jarg.find("markers");
    if ( itFind == jarg.end() )
        itFind = jarg.find("marker");
    if ( itFind != jarg.end() )
    {
        ModelMarkers mm;
        mm.setup( *itFind );
        M_markers = mm;
    }

    if ( jarg.contains( "max_distance" ) )
        M_maxDistance.setExpr( jarg.at( "max_distance" ), mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );

    if ( jarg.contains( "normalization" ) )
    {
        auto const& j_normalization = jarg.at( "normalization" );

        auto createNormalization = [this,&mMeshes]( nl::json const& jargLambda )
                                       {
                                           if ( auto optn = Normalization::create( mMeshes,jargLambda ) )
                                           {
                                               M_normalizations.push_back( std::move( optn.value() ) );
                                           }
                                           else
                                               throw std::runtime_error( "wrong json with normalization" );
                                       };

        if ( j_normalization.is_string() || j_normalization.is_object() )
        {
            createNormalization( j_normalization );
        }
        else if ( j_normalization.is_array() )
        {
            for ( auto const& [j_normalizationkey,j_normalizationval] : j_normalization.items() )
                createNormalization( j_normalizationval );
        }
    }
}

template <typename IndexType>
std::optional<typename ModelMesh<IndexType>::DistanceToRangeSetup::Normalization>
ModelMesh<IndexType>::DistanceToRangeSetup::Normalization::create( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
{
    std::string type, name;
    std::optional< std::pair<ModelExpression,ModelExpression> > range;
    if ( jarg.is_string() )
        type = jarg.template get<std::string>();
    else if ( jarg.is_object() )
    {
        if ( jarg.contains("type") )
            type = jarg.at("type").template get<std::string>();
        if ( jarg.contains("name") )
            name = jarg.at("name").template get<std::string>();
        if ( jarg.contains("range") )
        {
            auto const& j_range = jarg.at("range");
            if ( !j_range.is_array() )
                throw std::runtime_error( "range should be an array" );
            if ( j_range.size() != 2 )
                throw std::runtime_error( "range should be an array of size 2" );
            ModelExpression mexpr_a, mexpr_b;
            mexpr_a.setExpr( j_range[0], mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
            mexpr_b.setExpr( j_range[1], mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
            if ( !mexpr_a.hasExpr<1,1>() )
                throw std::runtime_error( "no scalar expr left bound in range" );
            if ( !mexpr_b.hasExpr<1,1>() )
                throw std::runtime_error( "no scalar expr left bound in range" );

            range = std::make_pair( std::move(mexpr_a),std::move(mexpr_b) );
        }
    }
    if ( name.empty() )
        name = type;

    if ( type == "min_max" )
    {
        if ( !range )
            range = std::make_pair( ModelExpression{0.},ModelExpression{1.} );
        return Normalization{ name,"min_max", *range };
    }
    else if ( type == "mean" )
    {
        return Normalization{ name,"mean" };
    }
    return {};
}

template <typename IndexType>
ModelMesh<IndexType>::MeshMotionSetup::MeshMotionSetup( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg )
{
    if ( jarg.contains("ComputationalDomain") )
    {
        auto const& j_cd = jarg.at("ComputationalDomain");
        if ( j_cd.contains( "markers" ) )
        {
            ModelMarkers markers;
            markers.setup( j_cd.at("markers")/*, indexes*/ );
            M_computationalDomainMarkers = markers;
        }
    }

    if ( jarg.contains("Displacement") )
    {
        auto const& j_disp = jarg.at("Displacement");
        if ( j_disp.contains( "Imposed" ) )
        {
            for ( auto const& [j_dispimposedkey,j_dispimposedval] : j_disp.at( "Imposed" ).items() )
            {
                std::set<std::string> markers = { j_dispimposedkey };
                if ( j_dispimposedval.contains( "markers") )
                {
                    ModelMarkers _markers;
                    _markers.setup( j_dispimposedval.at("markers")/*, indexes*/ );
                    markers = _markers;
                }
                ModelExpression mexpr;
                if ( j_dispimposedval.contains( "expr" ) )
                    mexpr.setExpr( j_dispimposedval.at( "expr" ), mMeshes.worldComm(), mMeshes.repository().expr()/*, indexes*/ );
                M_displacementImposed.emplace( j_dispimposedkey, std::make_tuple( std::move(mexpr),std::move( markers ) ) );
            }
        }
        if ( j_disp.contains( "Zero" ) )
        {
             ModelMarkers _markers;
             _markers.setup( j_disp.at("Zero")/*, indexes*/ );
             M_displacementZeroMarkers = _markers;
        }
        if ( j_disp.contains( "Free" ) )
        {
             ModelMarkers _markers;
             _markers.setup( j_disp.at("Free")/*, indexes*/ );
             M_displacementFreeMarkers = _markers;
        }
    }
}

template <typename IndexType>
ModelMesh<IndexType>::ModelMesh( std::string const& name, ModelMeshes<IndexType> const& mMeshes )
    :
    M_name( name ),
    M_mmeshCommon( std::make_shared<ModelMeshCommon<IndexType>>( mMeshes ) )
{}

template <typename IndexType>
void
ModelMesh<IndexType>::setup( nl::json const& jarg, ModelMeshes<IndexType> const& mMeshes )
{
    if ( jarg.contains("Import") )
        M_mmeshCommon->importConfig().setup( jarg.at("Import"), mMeshes );

    if ( jarg.contains("Partitioning") )
        M_mmeshCommon->partitioningSetup().updateForUse( jarg.at("Partitioning") );

    if ( jarg.contains("Data") )
    {
        for ( auto const& el : jarg.at("Data").items() )
         {
             std::string const& dataName = el.key();
             collection_data_by_mesh_entity_type c;
             c.setup( el.value() );
             M_codbme.emplace( std::make_pair( dataName, std::move( c ) ) );
         }
    }

    if ( jarg.contains("Fields") )
    {
        for ( auto const& el : jarg.at("Fields").items() )
        {
            std::string const& dataName = el.key();
            FieldsSetup fs(dataName,mMeshes,el.value() );
            if ( !fs.parts().empty() )
                M_fieldsSetup.push_back( std::move( fs ) );
        }
    }

    if ( jarg.contains("DistanceToRange") )
    {
        for ( auto const& el : jarg.at("DistanceToRange").items() )
        {
            std::string const& dataName = el.key();
            DistanceToRangeSetup dtrs(dataName,mMeshes,el.value() );
            M_distanceToRangeSetup.push_back( std::move( dtrs ) );
        }
    }

    if ( jarg.contains("MeshMotion") )
    {
        MeshMotionSetup mms( mMeshes, jarg.at("MeshMotion") );
        M_meshMotionSetup.emplace( std::move( mms ) );
    }

    if ( jarg.contains("MeshAdaptation") )
    {
        auto const& j_meshadapt = jarg.at("MeshAdaptation");
        if ( j_meshadapt.is_object() )
        {
            typename MeshAdaptation::Setup mas( this, mMeshes, j_meshadapt );
            M_meshAdaptationSetup.push_back( std::move( mas ) );
        }
        else if ( j_meshadapt.is_array() )
        {
            for ( auto const& [j_meshadaptkey,j_meshadaptval] : j_meshadapt.items() )
            {
                typename MeshAdaptation::Setup mas( this, mMeshes, j_meshadaptval );
                M_meshAdaptationSetup.push_back( std::move( mas ) );
            }
        }
        else
            throw std::runtime_error( "meshadaptation JSON value should be a an object or an array" );
    }
}

template <typename IndexType>
void
ModelMesh<IndexType>::setupRestart( ModelMeshes<IndexType> const& mMeshes )
{
    if ( M_mmeshCommon->hasMesh() )
        return;

    std::string meshFilename;
    if ( mMeshes.worldComm().isMasterRank() )
    {
        fs::path thedir = mMeshes.rootRepository();
        std::string fileNameMeshPath = (thedir/prefixvm(M_name,"mesh.path")).string();
        std::ifstream file( fileNameMeshPath.c_str() );
        if ( !file )
            CHECK( false ) << "Fail to open the txt file containing path of mesh file : " << fileNameMeshPath << "\n";

        if ( ! ( file >> meshFilename ) )
            CHECK( false ) << "Fail to read the msh path in file : " << fileNameMeshPath << "\n";
        file.close();

        CHECK( fs::exists( meshFilename ) ) << "restart mesh file doesn't exists : " << meshFilename;
    }
    mpi::broadcast( mMeshes.worldComm().localComm(), meshFilename, mMeshes.worldComm().masterRank() );

    mMeshes.log("ModelMesh","setupRestart", "mesh file : " + meshFilename);
    M_mmeshCommon->importConfig().setupInputMeshFilenameWithoutApplyPartitioning( meshFilename );
}


template <typename IndexType>
void
ModelMesh<IndexType>::updateInformationObject( nl::json & p, std::string const& prefix_symbol ) const
{
    if ( M_mmeshCommon->hasMesh() )
    {
        p["Discretization"] = this->mesh()->journalSection().to_string();
        M_mmeshCommon->importConfig().updateInformationObject( p["Import configuration"] );

        if ( !M_mmeshCommon->meshFilename().empty() )
            p.emplace( "filename", M_mmeshCommon->meshFilename() );
    }

    nl::json::array_t j_meshAdapArray;
    for ( auto const& mas : M_meshAdaptationSetup )
    {
        nl::json j_meshAdap;
        mas.updateInformationObject( j_meshAdap );
        if ( !j_meshAdap.is_null() )
            j_meshAdapArray.push_back( std::move(j_meshAdap) );
    }
    if ( !j_meshAdapArray.empty() )
        p["MeshAdaptation"] = std::move( j_meshAdapArray );

    if ( M_updateFieldInformations )
        M_updateFieldInformations( p["Fields"], prefix_symbol );
}

template <typename IndexType>
tabulate_informations_ptr_t
ModelMesh<IndexType>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    Feel::Table tabInfoOthers;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoOthers, jsonInfo, tabInfoProp );
    tabInfoOthers.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    if ( tabInfoOthers.nRow() > 0 )
        tabInfo->add( "", TabulateInformations::New( tabInfoOthers, tabInfoProp ) );

    if ( jsonInfo.contains("Import configuration") )
    {
        tabInfo->add( "Import configuration", import_config_type::tabulateInformations( jsonInfo.at("Import configuration"), tabInfoProp ) );
    }

    if ( jsonInfo.contains("Discretization") )
    {
        nl::json::json_pointer jsonPointerInfoDiscr( jsonInfo.at( "Discretization" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerInfoDiscr ) )
        {
            auto const& jsonInfoDiscr = JournalManager::journalData().at( jsonPointerInfoDiscr );
            auto tabInfoDiscr = TabulateInformationsSections::New( tabInfoProp );
            Feel::Table tabInfoDiscrEntries;
            TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoDiscrEntries, jsonInfoDiscr, tabInfoProp );
            tabInfoDiscrEntries.format()
                .setShowAllBorders( false )
                .setColumnSeparator(":")
                .setHasRowSeparator( false );
            tabInfoDiscr->add( "", TabulateInformations::New( tabInfoDiscrEntries, tabInfoProp ) );

            if ( jsonInfoDiscr.contains( "partitioning" ) )
            {
                int dim = jsonInfoDiscr.at("dim").template get<int>();
                auto const& jsonInfoDiscrPartitioning = jsonInfoDiscr.at( "partitioning" );
                Feel::Table tabInfoDiscrEntriesDataByPartition;
                tabInfoDiscrEntriesDataByPartition.format().setFirstRowIsHeader( true );
                if ( dim == 1 )
                    tabInfoDiscrEntriesDataByPartition.add_row({"partition id","n_elements","n_elements_with_ghost","n_points"});
                else if ( dim == 2 )
                    tabInfoDiscrEntriesDataByPartition.add_row({"partition id","n_elements","n_elements_with_ghost","n_faces","n_points"});
                else if ( dim == 3 )
                    tabInfoDiscrEntriesDataByPartition.add_row({"partition id","n_elements","n_elements_with_ghost","n_faces","n_edges","n_points"});
                auto jarray_n_elements = jsonInfoDiscrPartitioning.at("n_elements");
                auto jarray_n_elements_with_ghost = jsonInfoDiscrPartitioning.at("n_elements_with_ghost");
                auto itFind_n_faces = jsonInfoDiscrPartitioning.find( "n_faces" );
                auto itFind_n_edges = jsonInfoDiscrPartitioning.find( "n_edges" );
                auto itFind_n_points = jsonInfoDiscrPartitioning.find( "n_points" );
                for (int p=0;p<jarray_n_elements.size();++p)
                {
                    if ( dim == 1 )
                        tabInfoDiscrEntriesDataByPartition.add_row({ p,
                                    jarray_n_elements[p].template get<int>(),
                                    jarray_n_elements_with_ghost[p].template get<int>(),
                                    itFind_n_points.value()[p].template get<int>() });
                    else if ( dim == 2 )
                    {
                        CHECK( itFind_n_faces != jsonInfoDiscrPartitioning.end() ) << "missing face infos";
                        CHECK( itFind_n_points != jsonInfoDiscrPartitioning.end() ) << "missing points infos";
                        tabInfoDiscrEntriesDataByPartition.add_row({ p,
                                    jarray_n_elements[p].template get<int>(),
                                    jarray_n_elements_with_ghost[p].template get<int>(),
                                    itFind_n_faces.value()[p].template get<int>(),
                                    itFind_n_points.value()[p].template get<int>() });
                    }
                    else if ( dim == 3 )
                    {
                        CHECK( itFind_n_faces != jsonInfoDiscrPartitioning.end() ) << "missing face infos";
                        CHECK( itFind_n_edges != jsonInfoDiscrPartitioning.end() ) << "missing edges infos";
                        CHECK( itFind_n_points != jsonInfoDiscrPartitioning.end() ) << "missing points infos";
                        tabInfoDiscrEntriesDataByPartition.add_row({ p,
                                    jarray_n_elements[p].template get<int>(),
                                    jarray_n_elements_with_ghost[p].template get<int>(),
                                    itFind_n_faces.value()[p].template get<int>(),
                                    itFind_n_edges.value()[p].template get<int>(),
                                    itFind_n_points.value()[p].template get<int>() });
                    }
                }
                tabInfoDiscr->add( "", TabulateInformations::New( tabInfoDiscrEntriesDataByPartition, tabInfoProp.newByIncreasingVerboseLevel() ) );
            }

            tabInfo->add("Discretization", tabInfoDiscr );
        }
    }


    if ( jsonInfo.contains( "MeshAdaptation" ) )
    {
        auto tabInfoMeshAdaptation = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_meshadaptkey,j_meshadaptval ] : jsonInfo.at( "MeshAdaptation" ).items() )
            tabInfoMeshAdaptation->add( "", MeshAdaptation::Setup::tabulateInformations( j_meshadaptval, tabInfoProp/*.newByIncreasingVerboseLevel()*/ ) );
        tabInfo->add("Mesh Adaptation", tabInfoMeshAdaptation );
    }

    // fields
    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );

    return tabInfo;
}

template <typename IndexType>
void
ModelMeshes<IndexType>::setup( nl::json const& jarg, std::set<std::string> const& keywordsToSetup )
{
    for ( auto const& el : jarg.items() )
    {
         std::string const& meshName = el.key();
         if ( keywordsToSetup.find( meshName ) == keywordsToSetup.end() )
             continue;

         if ( this->hasModelMesh( meshName ) )
            this->at( meshName )->setup( el.value(), *this );
        else
        {
            auto me = std::make_shared<ModelMesh<IndexType>>( meshName );
            me->setup( el.value(), *this );
            this->emplace( std::make_pair( meshName, std::move( me ) ) );
        }
    }
}

template <typename IndexType>
void
ModelMeshes<IndexType>::updateInformationObject( nl::json & p, std::string const& prefix_symbol ) const
{
    for ( auto & [meshName,mMesh] : *this )
    {
        mMesh->updateInformationObject( p[meshName], prefixvm( prefix_symbol, meshName, "_" ) );
    }
}

template <typename IndexType>
tabulate_informations_ptr_t
ModelMeshes<IndexType>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
#if 0
    tabulate::Table tabInfo;
    for ( auto & [meshName,mMesh] : *this )
    {
        if ( jsonInfo.contains(meshName) )
        {
            tabulate::Table tabInfoMesh;
            tabInfoMesh.add_row( { (boost::format("Mesh : %1%")%meshName ).str() } );
            tabInfoMesh.add_row( { mMesh->tabulateInformation( jsonInfo.at(meshName), tabInfoProp ) } );
            tabInfo.add_row({tabInfoMesh});
        }
    }
    tabInfo.format().hide_border();
    return tabInfo;
#else
    auto tabInfo = TabulateInformationsSections::New();
    for ( auto & [meshName,mMesh] : *this )
    {
        if ( jsonInfo.contains(meshName) )
            tabInfo->add( (boost::format("Mesh : %1%")%meshName ).str(), mMesh->tabulateInformations( jsonInfo.at(meshName), tabInfoProp ) );
    }
    return tabInfo;
#endif
}


template class ModelMeshCommon<uint32_type>;
template class ModelMesh<uint32_type>;
template class ModelMeshes<uint32_type>;


} // namespace FeelModels
} // namespace Feel
