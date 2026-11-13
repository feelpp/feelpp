/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2020-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#pragma once

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/analysis.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feelmodels/modelcore/modelmeshadaptation_impl.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename IndexType>
template <typename MeshType>
void
ModelMesh<IndexType>::updateForUse( ModelMeshes<IndexType> const& mMeshes, RuntimeOrder runtimeGeometryOrder )
{
    using mesh_type = MeshType;
    auto makeMesh = [&]( worldcomm_ptr_t const& wc ) {
        if constexpr ( mesh_type::is_order_dynamic )
            return new mesh_type( runtimeGeometryOrder, M_name, wc );
        else
            return new mesh_type( M_name, wc );
    };

    mMeshes.log("ModelMesh","updateForUse", "start");


    if ( !M_mmeshCommon->hasMesh() )
    {
        std::shared_ptr<mesh_type> meshLoaded;
        std::string meshFilename;
        std::string meshFilenameBase = fmt::format("{}.mesh",mMeshes.keyword());
        auto & importConfig = M_mmeshCommon->importConfig();
        importConfig.updateForUse( mMeshes );
        nl::json partConfig = {
            { "partitioner", this->partitioningSetup().json() }
        };


        auto wcPtr = ( importConfig.loadByMasterRankOnly() )? mMeshes.worldCommPtr()->subWorldCommSeqPtr() : mMeshes.worldCommPtr();
        if ( importConfig.importFromSubmesh() )
        {
            mMeshes.log("ModelMesh","updateForUse", fmt::format("load from submesh : nMarkers={}",importConfig.submeshMarkers().size()) );

            if ( !importConfig.submeshMarkers().empty() )
            {
                switch ( importConfig.submeshCoDim() )
                {
                case 0:
                {
                    auto inputMesh = importConfig.template submeshInputMesh<mesh_type>();
                    bool doRepartitioning = importConfig.generatePartitioning();//true;
                    size_type updateSubmeshCtx = importConfig.meshComponents();
                    if ( doRepartitioning )
                    {
                        updateSubmeshCtx = MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED;
                        if ( true ) // TODO: only required if partitioning has constraints
                            updateSubmeshCtx |= MESH_UPDATE_FACES_MINIMAL;
                    }

                    std::shared_ptr<mesh_type> submesh;
                    if ( inputMesh )
                        submesh = createSubmesh(_range=markedelements(inputMesh, importConfig.submeshMarkers() ),
                                                _context=size_type(EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT), // WARNING not use relation, put an option!
                                                _update=updateSubmeshCtx );

                    if ( doRepartitioning )
                    {
                        if ( submesh )
                            submesh->saveHDF5( "toto.json" );
                        meshLoaded = loadMesh(_mesh=makeMesh( wcPtr ),
                                              _filename="toto.json",
                                              _update=importConfig.meshComponents(),
                                              _partitioning=partConfig,
                                              //_partitions=importConfig.numberOfPartition(),
                                              _rebuild_partitions=true );
                    }
                    else
                        meshLoaded = submesh;

                    // TODO straighten
                    break;
                }
                default:
                    CHECK( false ) << "TODO";
                }
            }
            else
                CHECK( false ) << "missing data for mesh import with an inputMesh";
        }
        else if ( importConfig.hasMeshFilename() )
        {
            std::string const& inputMeshFilename = importConfig.meshFilename();
            mMeshes.log("ModelMesh","updateForUse", "load mesh file : " + inputMeshFilename);

            std::string rootpath = mMeshes.rootRepository();
            std::string meshPartitionedFilename = (fs::path( rootpath ) / (meshFilenameBase + ".json")).string();
            std::string meshFileExt = fs::path( inputMeshFilename ).extension().string();
            bool generatePartitioning = importConfig.generatePartitioning();
            // if ( generatePartitioning && meshFileExt != ".msh" )
            //     CHECK( false ) << "Can not rebuild at this time the mesh partitionining with other format than .msh : TODO";

            if ( !importConfig.loadByMasterRankOnly() || mMeshes.worldCommPtr()->isMasterRank() )
            {
                meshLoaded = loadMesh(_mesh=makeMesh( wcPtr ),
                                      _filename=inputMeshFilename,
                                      _prefix=mMeshes.prefix(),
                                      _vm=mMeshes.clovm(),
                                      _worldcomm=wcPtr/*mMeshes.worldCommPtr()*/,
                                      _straighten=importConfig.straightenMesh(),
                                      _partitioning=partConfig,
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
                meshLoaded = createGMSHMesh(_mesh=makeMesh( wcPtr ),
                                            _desc=geodesc,
                                            _prefix=mMeshes.prefix(),
                                            _vm=mMeshes.clovm(),
                                            _worldcomm=wcPtr/*mMeshes.worldCommPtr()*/,
                                            _h=importConfig.meshSize(),
                                            _straighten=importConfig.straightenMesh(),
                                            _partitions=importConfig.numberOfPartition(),
                                            _partitioning=partConfig,
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

    auto themesh = this->mesh<MeshType>();

    M_updateFieldInformations = [this]( nl::json& info, std::string const& prefix ) {
        this->template modelFields<MeshType>( "", prefix ).updateInformationObject( info );
    };

    // update data by mesh entity
    for ( auto & [name,data] : M_codbme )
    {
        data.setMesh( this->mesh<mesh_type>() );
        data.template updateForUse<mesh_type>();
    }


    // fields
    static constexpr auto tuple_t_basis = ModelMesh<IndexType>::template basisFieldTypeSupported<mesh_type>();
    for ( auto const& fs : M_fieldsSetup )
    {
        std::string const& basis = fs.basis();
        hana::for_each( tuple_t_basis, [this,&basis,&fs,&themesh]( auto const& b )
                        {
                            if ( basis == hana::at_c<0>( b ) )
                            {
                                using basis_type = typename std::decay_t<decltype(hana::at_c<1>( b ) )>::type;
                                using space_type = FunctionSpace<MeshType, bases<basis_type> >;

                                auto Vh = M_mmeshCommon->template createFunctionSpace<space_type>( basis );
                                auto u = Vh->elementPtr();

                                // TODO : manage sync of dof interprocess
                                for ( auto const& partSetup : fs.parts() )
                                {
                                    auto rangeMeshElt = partSetup.markers().empty() ? elements(themesh) : markedelements(themesh, (std::set<std::string>)partSetup.markers() );
                                    if ( partSetup.hasExpr() )
                                    {
                                        if ( !partSetup.mExpr().template hasExpr<space_type::nComponents1,space_type::nComponents2>() )
                                            throw std::runtime_error( fmt::format("ModelMesh<IndexType>::updateForUse fields basis and expr rank invalid, expr should be [{},{}]", space_type::nComponents1,space_type::nComponents2 ) );
                                        auto const& theexpr = partSetup.mExpr().template expr<space_type::nComponents1,space_type::nComponents2>();
                                        u->on(_range=rangeMeshElt, _expr=theexpr );
                                    }
                                    else
                                    {
                                        // TODO use rangeMeshElt
                                        CHECK( partSetup.markers().empty() ) << "TODO ModelMesh<IndexType>::updateForUse fields : requires the use of range elements with file loading";
                                        u->load(_path=partSetup.filename(),_type="default",_space_path=partSetup.spaceFilename());
                                    }
                                }
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
        themesh = this->mesh<MeshType>();
    }

#if defined(FEELPP_MODELS_HAS_MESHALE)
    if constexpr ( mesh_type::is_order_dynamic )
    {
        if ( M_meshMotionSetup )
            throw std::invalid_argument("Mesh motion is not yet supported with runtime geometry order");
    }
    else if constexpr ( mesh_type::nDim>1 )
    {
        if ( M_meshMotionSetup )
        {
            //auto themesh = this->mesh<MeshType>();
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

#else
    if ( M_meshMotionSetup )
        throw std::invalid_argument("Mesh motion support is disabled in this build");
#endif
}


template <typename IndexType>
template <typename MeshType>
void
ModelMesh<IndexType>::updateDistanceToRange()
{
    for ( auto const& dtrs : M_distanceToRangeSetup )
    {
        auto themesh = this->mesh<MeshType>();
        std::string basis = fmt::format( "Pch{}", themesh->order() );
        //std::cout << "distance to range  create basis " << basis << std::endl;
        using distange_to_range_space_type = FunctionSpace<MeshType, bases<Lagrange<MeshType::is_order_dynamic ? Dynamic : MeshType::nOrder,Scalar,Continuous,PointSetFekete>>>;
        auto Vh = M_mmeshCommon->template createFunctionSpace<distange_to_range_space_type>( basis, RuntimeOrder{ themesh->order() } );
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
                    double average = Feel::mean(_range=elements(support(u->functionSpace())),_expr=idv(u))(0,0);
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

#if defined(FEELPP_MODELS_HAS_MESHALE)
    if constexpr ( !MeshType::is_order_dynamic && MeshType::nDim > 1 )
    {
        if ( this->hasMeshMotion() )
        {
            auto mmt = this->meshMotionTool<MeshType>();

            std::vector<std::tuple<std::string,Range<MeshType,MESH_ELEMENTS>>> computationDomains;
            if ( M_meshMotionSetup )
            {
                auto const& cdm = M_meshMotionSetup->computationalDomainMarkers();
                if ( cdm.find( "@elements@" ) == cdm.end() )
                    computationDomains.push_back( std::make_tuple( M_name, markedelements(newMesh, cdm) ) );
            }
            mmt->applyRemesh( newMesh, computationDomains );
        }
    }
#endif
    // TODO : other fields
}


} // namespace FeelModels
} // namespace Feel
