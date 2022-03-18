/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelmeshes.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
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
            FieldsSetup fs(dataName,el.value() );
            M_fieldsSetup.push_back( std::move( fs ) );
        }
    }

    if ( jarg.contains("DistanceToRange") )
    {
        for ( auto const& el : jarg.at("DistanceToRange").items() )
        {
            std::string const& dataName = el.key();
            DistanceToRangeSetup dtrs(dataName,el.value() );
            M_distanceToRangeSetup.push_back( std::move( dtrs ) );
        }
    }

    if ( jarg.contains("MeshMotion") )
    {
        MeshMotionSetup mms( mMeshes, jarg.at("MeshMotion") );
        M_meshMotionSetup.emplace( std::move( mms ) );
    }
}

template <typename IndexType>
void
ModelMesh<IndexType>::setupRestart( ModelMeshes<IndexType> const& mMeshes )
{
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
template <typename MeshType>
void
ModelMesh<IndexType>::updateForUse( ModelMeshes<IndexType> const& mMeshes )
{
    using mesh_type = MeshType;

    if ( !M_mmeshCommon->hasMesh() )
    {
        std::shared_ptr<mesh_type> meshLoaded;
        std::string meshFilename;
        std::string meshFilenameBase = fmt::format("{}.mesh",mMeshes.keyword());
        auto & importConfig = M_mmeshCommon->importConfig();
        importConfig.updateForUse( mMeshes );

        auto wcPtr = ( importConfig.loadByMasterRankOnly() )? mMeshes.worldCommPtr()->subWorldCommSeqPtr() : mMeshes.worldCommPtr();
        if ( importConfig.hasMeshFilename() )
        {
            std::string const& inputMeshFilename = importConfig.meshFilename();
            mMeshes.log("ModelMesh","updateForUse", "load mesh file : " + inputMeshFilename);

            std::string rootpath = mMeshes.rootRepository();
            std::string meshPartitionedFilename = (fs::path( rootpath ) / (meshFilenameBase + ".json")).string();
            std::string meshFileExt = fs::path( inputMeshFilename ).extension().string();
            bool generatePartitioning = importConfig.generatePartitioning();
            if ( generatePartitioning && meshFileExt != ".msh" )
                CHECK( false ) << "Can not rebuild at this time the mesh partitionining with other format than .msh : TODO";

            if ( !importConfig.loadByMasterRankOnly() || mMeshes.worldCommPtr()->isMasterRank() )
            {
                meshLoaded = loadMesh(_mesh=new mesh_type( M_name, wcPtr/*mMeshes.worldCommPtr()*/ ),
                                  _filename=inputMeshFilename,
                                  _prefix=mMeshes.prefix(),
                                  _vm=mMeshes.clovm(),
                                  _worldcomm=wcPtr/*mMeshes.worldCommPtr()*/,
                                  _straighten=importConfig.straightenMesh(),
                                  _rebuild_partitions=generatePartitioning,
                                  _rebuild_partitions_filename=meshPartitionedFilename,
                                  _partitions=importConfig.numberOfPartition(),
                                  _savehdf5=0,
                                  _update= importConfig.meshComponents()/*MESH_UPDATE_EDGES|MESH_UPDATE_FACES*/);
            }

            meshFilename = (generatePartitioning)? meshPartitionedFilename : importConfig.meshFilename();
        }
        else if ( importConfig.hasGeoFilename() )
        {
            std::string const& inputGeoFilename = importConfig.geoFilename();
            mMeshes.log("ModelMesh","updateForUse", "load geo file : " + inputGeoFilename);
            std::string path = mMeshes.rootRepository();

            std::string mshfile = (fs::path( path ) / meshFilenameBase).string();
            if ( importConfig.numberOfPartition() > 1 )
                mshfile += ".json";
            else
                mshfile += ".msh";
            if ( !importConfig.loadByMasterRankOnly() || mMeshes.worldCommPtr()->isMasterRank() )
            {
                gmsh_ptrtype geodesc = geo( _filename=inputGeoFilename,
                                            _prefix=mMeshes.prefix(),
                                            _vm=mMeshes.clovm(),
                                            _worldcomm=wcPtr/*mMeshes.worldCommPtr()*/,
                                            _h=importConfig.meshSize());
                // allow to have a geo and msh file with a filename equal to prefix
                geodesc->setPrefix(meshFilenameBase);
                meshLoaded = createGMSHMesh(_mesh=new mesh_type( M_name, wcPtr/*mMeshes.worldCommPtr()*/ ),
                                        _desc=geodesc,
                                        _prefix=mMeshes.prefix(),
                                        _vm=mMeshes.clovm(),
                                        _worldcomm=wcPtr/*mMeshes.worldCommPtr()*/,
                                        _h=importConfig.meshSize(),
                                        _straighten=importConfig.straightenMesh(),
                                        _partitions=importConfig.numberOfPartition(),
                                        _update=importConfig.meshComponents(),
                                        _directory=mMeshes.rootRepository() );
            }
            meshFilename = mshfile;
        }

        this->setMesh( meshLoaded, meshFilename );

        if ( meshLoaded && !meshFilename.empty() && mMeshes.worldComm().isMasterRank() )
        {
            fs::path thedir = mMeshes.rootRepository();//fs::path( fileSavePath ).parent_path();
            std::string fileNameMeshPath = (thedir/prefixvm(M_name,"mesh.path")).string();
            if ( !fs::exists(thedir))
                fs::create_directories(thedir);
            std::ofstream file( fileNameMeshPath.c_str(), std::ios::out|std::ios::trunc);
            file << meshFilename;
            file.close();
        }
    } // if ( !hasMesh )

    // update data by mesh entity
    for ( auto & [name,data] : M_codbme )
    {
        data.setMesh( this->mesh<mesh_type>() );
        data.template updateForUse<mesh_type>();
    }


    // fields
    static constexpr auto tuple_t_basis = hana::make_tuple( hana::make_tuple( "Pch1", hana::type_c<Lagrange<1,Scalar,Continuous,PointSetFekete>> ),
                                                            hana::make_tuple( "Pch2", hana::type_c<Lagrange<2,Scalar,Continuous,PointSetFekete>> ),
                                                            hana::make_tuple( "Pchv1", hana::type_c<Lagrange<1,Vectorial,Continuous,PointSetFekete>> ),
                                                            hana::make_tuple( "Pchv2", hana::type_c<Lagrange<2,Vectorial,Continuous,PointSetFekete>> )
                                                            );
    for ( auto const& fs : M_fieldsSetup )
    {
        std::string const& basis = fs.basis();
        hana::for_each( tuple_t_basis, [this,&basis,&fs]( auto const& b )
                        {
                            if ( basis == hana::at_c<0>( b ) )
                            {
                                using basis_type = typename std::decay_t<decltype(hana::at_c<1>( b ) )>::type;
                                using space_type = FunctionSpace<MeshType, bases<basis_type> >;

                                auto Vh = M_mmeshCommon->template createFunctionSpace<space_type>( basis );
                                auto u = Vh->elementPtr();
                                u->load(_path=fs.filename(),_type="default");
                                M_fields[fs.name()] = u;
                            }
                        });
    }

    // distance to range
    for ( auto const& dtrs : M_distanceToRangeSetup )
    {
        auto themesh = this->mesh<MeshType>();
        std::string basis = fmt::format( "Pch{}", MeshType::nOrder );
        //std::cout << "distance to range  create basis " << basis << std::endl;
        using distange_to_range_space_type = FunctionSpace<MeshType, bases<Lagrange<MeshType::nOrder,Scalar,Continuous,PointSetFekete>>>;
        auto Vh = M_mmeshCommon->template createFunctionSpace<distange_to_range_space_type>( basis );
        auto u = Vh->elementPtr();
        auto rangeFaces = dtrs.markers().empty()? boundaryfaces(themesh) : markedfaces( themesh, dtrs.markers() );
        *u = distanceToRange( Vh, rangeFaces );
        M_distanceToRanges[dtrs.name()] = u;
    }

    if constexpr ( mesh_type::nDim>1 )
    {
        if ( M_meshMotionSetup )
        {
            auto themesh = this->mesh<MeshType>();
            auto meshALE = meshale( _mesh=themesh,
                                    _prefix=mMeshes.prefix(),
                                    _keyword=fmt::format("meshes_{}_meshmotion",M_name),
                                    _directory=mMeshes.repository() );

            auto const& cdm = M_meshMotionSetup->computationalDomainMarkers();
            if ( cdm.find( "@elements@" ) != cdm.end() )
                meshALE->setWholeMeshAsComputationalDomain( M_name/*mMeshes.keyword()*/ );
            else
                meshALE->setComputationalDomain( M_name/*mMeshes.keyword()*/, markedelements(themesh, cdm) );
            meshALE->init();

            std::set<std::string> markersDispImposedOverFaces;
            for ( auto const& [name,dispData] : M_meshMotionSetup->displacementImposed() )
                markersDispImposedOverFaces.insert(  std::get<1>( dispData ).begin(),  std::get<1>( dispData ).end() );
            if ( !markersDispImposedOverFaces.empty() )
                meshALE->setDisplacementImposedOnInitialDomainOverFaces( M_name/*this->keyword()*/, markersDispImposedOverFaces );

            meshALE->addMarkersInBoundaryCondition( "fixed", M_meshMotionSetup->displacementZeroMarkers() );
            meshALE->addMarkersInBoundaryCondition( "free", M_meshMotionSetup->displacementFreeMarkers() );
            for ( auto const& [name,dispData] : M_meshMotionSetup->displacementImposed() )
                meshALE->addMarkersInBoundaryCondition( "moving", std::get<1>( dispData ) );

            M_mmeshCommon->setMeshMotionTool( meshALE );
        }
    }

}
template <typename IndexType>
void
ModelMesh<IndexType>::updateInformationObject( nl::json & p ) const
{
    if ( M_mmeshCommon->hasMesh() )
    {
        p["Discretization"] = this->mesh()->journalSection().to_string();
        M_mmeshCommon->importConfig().updateInformationObject( p["Import configuration"] );

        if ( !M_mmeshCommon->meshFilename().empty() )
            p.emplace( "filename", M_mmeshCommon->meshFilename() );
    }
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
ModelMeshes<IndexType>::updateInformationObject( nl::json & p ) const
{
    for ( auto & [meshName,mMesh] : *this )
    {
        mMesh->updateInformationObject( p[meshName] );
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

template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,1>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,1>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<1,1,2>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<1,1,3>>>( ModelMeshes<uint32_type> const& );
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,2>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,2>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<1,2,2>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<1,2,3>>>( ModelMeshes<uint32_type> const& );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,3>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,3>>>( ModelMeshes<uint32_type> const& );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,4>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,4>>>( ModelMeshes<uint32_type> const& );
#endif
} // namespace FeelModels
} // namespace Feel
