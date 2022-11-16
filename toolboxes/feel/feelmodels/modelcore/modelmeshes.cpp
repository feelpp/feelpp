/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelmeshes.hpp>

#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <boost/preprocessor/array/to_list.hpp>
#include <boost/preprocessor/list/append.hpp>

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
        else if ( j_partition.is_number_unsigned() )
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
            typename MeshAdaptation::Setup mas( mMeshes, j_meshadapt );
            M_meshAdaptationSetup.push_back( std::move( mas ) );
        }
        else if ( j_meshadapt.is_array() )
        {
            for ( auto const& [j_meshadaptkey,j_meshadaptval] : j_meshadapt.items() )
            {
                typename MeshAdaptation::Setup mas( mMeshes, j_meshadaptval );
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

        if ( meshLoaded )
        {
            if ( !M_metadata.contains("preprocess") )
                M_metadata["preprocess"] = json::array();

            nl::json importMeta = { { "event","import" }, { "filename",meshFilename } };
            M_metadata["preprocess"].push_back( std::move( importMeta ) );
        }

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
    this->updateDistanceToRange<MeshType>();

    // mesh adaptation
    if ( std::find_if( M_meshAdaptationSetup.begin(), M_meshAdaptationSetup.end(),
                       []( auto const& mas ){ return mas.isExecutedWhen( MeshAdaptation::template createEvent<MeshAdaptation::Event::Type::after_import>() ); } ) != M_meshAdaptationSetup.end() )
    {
        const_cast<ModelMeshes<IndexType>&>( mMeshes ).modelProperties().parameters().updateParameterValues();
        auto paramValues = mMeshes.modelProperties().parameters().toParameterValues();
        for ( typename MeshAdaptation::Setup & mas : M_meshAdaptationSetup )
            mas.setParameterValues( paramValues );
        this->updateMeshAdaptation<MeshType>( MeshAdaptation::template createEvent<MeshAdaptation::Event::Type::after_import>(),
                                              Feel::vf::symbolsExpr( mMeshes.symbolsExprParameter(), mMeshes.template symbolsExpr<MeshType>() ) );
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
template <typename MeshType>
void
ModelMesh<IndexType>::updateDistanceToRange()
{
    for ( auto const& dtrs : M_distanceToRangeSetup )
    {
        auto themesh = this->mesh<MeshType>();
        std::string basis = fmt::format( "Pch{}", MeshType::nOrder );
        //std::cout << "distance to range  create basis " << basis << std::endl;
        using distange_to_range_space_type = FunctionSpace<MeshType, bases<Lagrange<MeshType::nOrder,Scalar,Continuous,PointSetFekete>>>;
        auto Vh = M_mmeshCommon->template createFunctionSpace<distange_to_range_space_type>( basis );
        auto u = Vh->elementPtr();
        auto rangeFaces = dtrs.markers().empty()? boundaryfaces(themesh) : markedfaces( themesh, dtrs.markers() );
        *u = distanceToRange( _space=Vh, _range=rangeFaces, _max_distance = dtrs.maxDistance() );
        M_distanceToRanges[dtrs.name()] = u;

        for ( auto const& normalization : dtrs.normalizations() )
        {
            std::string const& normalizationType = normalization.type();
            if ( normalizationType == "min_max" || normalizationType == "mean" )
            {
                double uMax = u->max();
                double uMin = u->min();

                // min_max normalization : u_normalized \in [0,1] with u_normalized = (u-min(u))/(max(u)-min(u)
                if ( normalizationType == "min_max" )
                {
                    //std::cout << "uMax=" << uMax << " uMin=" << uMin << std::endl;
                    auto uNormalizedMinMax = Vh->elementPtr();
                    *uNormalizedMinMax = *u;
                    //*uNormalizedMinMax -= uMin;
                    uNormalizedMinMax->add( -uMin );
                    *uNormalizedMinMax *= 1./(uMax-uMin);
                    if ( true )
                    {
                        // u_normalized_ab = a + u_normalized*(b-a)
                        auto const& [aExpr,bExpr] = normalization.range();
                        double a = aExpr.template expr<1,1>().evaluate()(0,0);
                        double b = bExpr.template expr<1,1>().evaluate()(0,0);
                        if ( b < a )
                            throw std::runtime_error( fmt::format("b={} less than a={}",a,b) );
                        //double a = 0, b = 1;
                        *uNormalizedMinMax *= (b-a);
                        uNormalizedMinMax->add( a );
                    }
                    M_distanceToRanges[fmt::format("{}_normalized_{}",dtrs.name(),normalization.name())] = uNormalizedMinMax;
                }

                // mean normalization :  u_normalized = (u-average(u))/(max(u)-min(u))
                if ( normalizationType == "mean" )
                {
                    double average = mean(_range=elements(support(u->functionSpace())),_expr=idv(u))(0,0);
                    auto uNormalizedMean = Vh->elementPtr();
                    *uNormalizedMean = *u;
                    //*uNormalizedMean -= average;
                    uNormalizedMean->add( -average );
                    *uNormalizedMean *= 1./(uMax-uMin);
                    M_distanceToRanges[fmt::format("{}_normalized_{}",dtrs.name(),normalization.name())] = uNormalizedMean;
                }
            }
        }
    }
}

template <typename IndexType>
template <typename MeshType>
void
ModelMesh<IndexType>::applyRemesh( std::shared_ptr<MeshType> const& newMesh )
{
    auto oldmesh = this->mesh<MeshType>();
    if ( oldmesh == newMesh )
        return;

    this->setMesh( newMesh );

    M_mmeshCommon->clearFunctionSpaces();
    // TODO : common clear space + applyRemesh to pointmeasure

    M_distanceToRanges.clear();
    this->updateDistanceToRange<MeshType>();

    if constexpr ( MeshType::nDim > 1 )
    {
        if ( this->hasMeshMotion() )
        {
            auto mmt = this->meshMotionTool<MeshType>();

            std::vector<std::tuple<std::string,elements_reference_wrapper_t<MeshType>>> computationDomains;
            if ( M_meshMotionSetup )
            {
                auto const& cdm = M_meshMotionSetup->computationalDomainMarkers();
                if ( cdm.find( "@elements@" ) == cdm.end() )
                    computationDomains.push_back( std::make_tuple( M_name, markedelements(newMesh, cdm) ) );
            }
            mmt->applyRemesh( newMesh, computationDomains );
        }
    }
    // TODO : other fields
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

    // TODO here : other types
    using geoshape_list_type = boost::mp11::mp_list< boost::mp11::mp_identity_t<Simplex<2>>,
                                                     boost::mp11::mp_identity_t<Simplex<3>> >;
    boost::mp11::mp_for_each<geoshape_list_type>( [this,&p,&prefix_symbol]( auto I ){
                                                      using _mesh_type = typename Mesh<std::decay_t<decltype(I)>>::type;
                                                      if ( this->mesh<_mesh_type>() )
                                                          this->modelFields<_mesh_type>( "",prefix_symbol ).updateInformationObject( p["Fields"] );
                                                  } );
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


#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER1_LIST \
    BOOST_PP_TUPLE_TO_LIST(                                             \
        ( ( Simplex,2,1,2),                                             \
          ( Simplex,3,1,3),                                             \
          ( Simplex,1,1,2),                                             \
          ( Simplex,1,1,3) ) )                                          \
    /**/

#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER2_LIST \
    BOOST_PP_TUPLE_TO_LIST(                                             \
        ( ( Simplex,2,2,2),                                             \
          ( Simplex,3,2,3),                                             \
          ( Simplex,1,2,2),                                             \
          ( Simplex,1,2,3) ) )                                          \
    /**/
#else
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER2_LIST BOOST_PP_NIL
#endif

#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER3_LIST \
    BOOST_PP_TUPLE_TO_LIST(                                             \
        ( ( Simplex,2,3,2),                                             \
          ( Simplex,3,3,3) ) )                                          \
    /**/
#else
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER3_LIST BOOST_PP_NIL
#endif

#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER4_LIST \
    BOOST_PP_TUPLE_TO_LIST(                                             \
        ( ( Simplex,2,4,2),                                             \
          ( Simplex,3,4,3) ) )                                          \
    /**/
#else
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER4_LIST BOOST_PP_NIL
#endif


#define FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_CLASS_NAME(T)   BOOST_PP_TUPLE_ELEM(4, 0, T)
#define FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_DIM(T)   BOOST_PP_TUPLE_ELEM(4, 1, T)
#define FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_ORDER(T)   BOOST_PP_TUPLE_ELEM(4, 2, T)
#define FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_REALDIM(T)   BOOST_PP_TUPLE_ELEM(4, 3, T)

#define FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_CLASS(T)               \
    FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_CLASS_NAME(T)              \
    < FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_DIM(T),                  \
      FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_ORDER(T),                \
      FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_REALDIM(T) >             \
    /**/

#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_LIST     \
    BOOST_PP_LIST_APPEND( FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER1_LIST, \
                          BOOST_PP_LIST_APPEND( FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER2_LIST, \
                                                BOOST_PP_LIST_APPEND( FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER3_LIST, \
                                                                      FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_ORDER4_LIST ) ) ) \
    /**/

#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_INDEXTYPE_LIST    \
    BOOST_PP_ARRAY_TO_LIST((1, (uint32_type)))                          \
    /**/


/* Generates code for all binary operators and integral type pairs. */
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_METHODS_OP(_, IS) \
    FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_METHODS_OP_CODE IS    \
   /**/
#define FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_METHODS_OP_CODE(PP_I,PP_GS) \
    template void ModelMesh<PP_I>::updateForUse<Mesh<FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_CLASS(PP_GS)>>( ModelMeshes<PP_I> const& ); \
    template void ModelMesh<PP_I>::applyRemesh<Mesh<FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_CLASS(PP_GS)>>( std::shared_ptr<Mesh<FEELPP_TOOLBOXES_PP_MODELMESHES_GEOSHAPE_CLASS(PP_GS)>> const& );
    /**/

BOOST_PP_LIST_FOR_EACH_PRODUCT( FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_METHODS_OP, 2, (FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_INDEXTYPE_LIST,FEELPP_TOOLBOXES_PP_MODELMESHES_INSTANTIATION_GEOSHAPE_LIST) )
/**/

} // namespace FeelModels
} // namespace Feel
