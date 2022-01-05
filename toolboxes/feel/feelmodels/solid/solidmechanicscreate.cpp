/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::SolidMechanics( std::string const& prefix,
                                                    std::string const& keyword,
                                                    worldcomm_ptr_t const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<nDim>( "solid" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
{
    this->log("SolidMechanics","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("SolidMechanics","constructor", "finish" );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix,
                                          std::string const& keyword,
                                         worldcomm_ptr_t const& worldComm,
                                         std::string const& subPrefix,
                                         ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type>( prefix, keyword, worldComm, subPrefix, modelRep );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::expandStringFromSpec( std::string const& expr )
{
    std::string res = expr;
    boost::replace_all( res, "$solid_disp_order", (boost::format("%1%")%nOrderDisplacement).str() );
    boost::replace_all( res, "$solid_geo_order", (boost::format("%1%")%nOrderGeo).str() );
    std::string solidTag = (boost::format("P%1%G%2%")%nOrderDisplacement %nOrderGeo ).str();
    boost::replace_all( res, "$solid_tag", solidTag );
    return res;
}

// add members instatantiations need by static function expandStringFromSpec
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
const uint16_type SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::nOrderDisplacement;
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
const uint16_type SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::nOrderGeo;


//---------------------------------------------------------------------------------------------------//


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    this->log("SolidMechanics","loadParameterFromOptionsVm", "start" );
#if 0
    std::string formulation = soption(_name="formulation",_prefix=this->prefix());
    M_useDisplacementPressureFormulation = false;
    if ( formulation == "displacement-pressure" )
        M_useDisplacementPressureFormulation = true;
    M_mechanicalProperties->setUseDisplacementPressureFormulation(M_useDisplacementPressureFormulation);

    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        this->setModelName( soption(_name="model",_prefix=this->prefix()) );
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );
#endif
    this->setSolver( soption(_name="solver",_prefix=this->prefix()) );

    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "dirichlet-neumann";

    M_isHOVisu = nOrderGeo > 2;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());

    //time schema parameters
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix());
    M_timeSteppingUseMixedFormulation = false;
    M_genAlpha_alpha_m=1.0;
    M_genAlpha_alpha_f=1.0;
    if ( M_timeStepping == "Newmark" )
    {
        M_genAlpha_alpha_m=1.0;
        M_genAlpha_alpha_f=1.0;
    }
    else if ( M_timeStepping == "Generalized-Alpha" )
    {
#if 0
        M_genAlpha_rho=doption(_name="time-rho",_prefix=this->prefix());
        M_genAlpha_alpha_m=(2.- M_genAlpha_rho)/(1.+M_genAlpha_rho);
        M_genAlpha_alpha_f=1./(1.+M_genAlpha_rho);
#endif
    }
    else if ( M_timeStepping == "BDF" || M_timeStepping == "Theta" )
    {
        M_timeSteppingUseMixedFormulation = true;
        M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix());
    }
    else CHECK( false ) << "time stepping not supported : " << M_timeStepping << "\n";

    M_genAlpha_gamma=0.5+M_genAlpha_alpha_m-M_genAlpha_alpha_f;
    M_genAlpha_beta=0.25*(1+M_genAlpha_alpha_m-M_genAlpha_alpha_f)*(1+M_genAlpha_alpha_m-M_genAlpha_alpha_f);

    M_useMassMatrixLumped = false;

    this->log("SolidMechanics","loadParameterFromOptionsVm", "finish" );
}



//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_bcDirichletMarkerManagement.clearMarkerDirichletBC();
    M_bcNeumannMarkerManagement.clearMarkerNeumannBC();
    M_bcFSIMarkerManagement.clearMarkerFluidStructureInterfaceBC();
    M_bcRobinMarkerManagement.clearMarkerRobinBC();

    std::string dirichletbcType = "elimination";//soption(_name="dirichletbc.type",_prefix=this->prefix());
    this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        M_bcDirichletMarkerManagement.addMarkerDirichletBC( dirichletbcType, name(d), markers(d), ComponentType::NO_COMPONENT );
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        this->M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( (boost::format("displacement_%1%")%compTag).str(), "Dirichlet" );
        for( auto const& d : this->M_bcDirichletComponents.find(comp)->second )
            M_bcDirichletMarkerManagement.addMarkerDirichletBC( dirichletbcType, name(d), markers(d), comp );
    }

    this->M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "Neumann_scalar" );
    for( auto const& d : this->M_bcNeumannScalar )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d),markers(d));
    this->M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "Neumann_vectorial" );
    for( auto const& d : this->M_bcNeumannVectorial )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::VECTORIAL,name(d),markers(d));
    this->M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<nDim>( "displacement", "Neumann_tensor2" );
    for( auto const& d : this->M_bcNeumannTensor2 )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::TENSOR2,name(d),markers(d));

    this->M_bcInterfaceFSI = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "interface_fsi" );
    for( auto const& d : this->M_bcInterfaceFSI )
        M_bcFSIMarkerManagement.addMarkerFluidStructureInterfaceBC( markers(d) );

    this->M_bcRobin = this->modelProperties().boundaryConditions().template getVectorFieldsList<nDim>( "displacement", "robin" );
    for( auto const& d : this->M_bcRobin )
        M_bcRobinMarkerManagement.addMarkerRobinBC( name(d),markers(d) );

    this->M_bcNeumannEulerianFrameScalar = this->modelProperties().boundaryConditions().getScalarFields( { { "displacement", "Neumann_eulerian_scalar" },{ "displacement", "FollowerPressure" } } );
    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
        M_bcNeumannEulerianFrameMarkerManagement.addMarkerNeumannEulerianFrameBC(MarkerManagementNeumannEulerianFrameBC::NeumannEulerianFrameBCShape::SCALAR,name(d),markers(d));
    this->M_bcNeumannEulerianFrameVectorial = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "Neumann_eulerian_vectorial" );
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
        M_bcNeumannEulerianFrameMarkerManagement.addMarkerNeumannEulerianFrameBC(MarkerManagementNeumannEulerianFrameBC::NeumannEulerianFrameBCShape::VECTORIAL,name(d),markers(d));
    this->M_bcNeumannEulerianFrameTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<nDim>( "displacement", "Neumann_eulerian_tensor2" );
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
        M_bcNeumannEulerianFrameMarkerManagement.addMarkerNeumannEulerianFrameBC(MarkerManagementNeumannEulerianFrameBC::NeumannEulerianFrameBCShape::TENSOR2,name(d),markers(d));

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "VolumicForces" );



    auto XhDisp = this->functionSpaceDisplacement();
    auto XhCompDisp = XhDisp->compSpace();
    static const uint16_type nDofComponentsDisp = XhDisp->dof()->nDofComponents();
    auto mesh = this->mesh();

    std::set<std::string> dispMarkers;
    std::map<ComponentType,std::set<std::string> > compDispMarkers;

    // strong Dirichlet bc on displacement from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) );
        dispMarkers.insert( listMark.begin(), listMark.end() );
    }
    // strong Dirichlet bc on displacement component from expression
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d), comp );
            compDispMarkers[comp].insert( listMark.begin(), listMark.end() );
        }
    }



    //-------------------------------------//
    // distribute mesh markers by entity
    std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > meshMarkersDispByEntities;
    std::map<ComponentType, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > meshMarkersCompDispByEntities;
    meshMarkersDispByEntities = detail::distributeMarkerListOnSubEntity( mesh, dispMarkers );
    for ( auto const& compMarkerPair : compDispMarkers )
    {
        meshMarkersCompDispByEntities[compMarkerPair.first] = detail::distributeMarkerListOnSubEntity( mesh, compMarkerPair.second );
    }
    //-------------------------------------//
    // on topological faces
    auto const& listMarkedFacesDisp = std::get<0>( meshMarkersDispByEntities );
    if ( !listMarkedFacesDisp.empty() )
        this->updateDofEliminationIds( "displacement", XhDisp, markedfaces( mesh,listMarkedFacesDisp ) );
    // on marked edges (only 3d)
    if constexpr ( nDim == 3)
    {
        auto const& listMarkedEdgesDisp = std::get<1>( meshMarkersDispByEntities );
        if ( !listMarkedEdgesDisp.empty() )
            this->updateDofEliminationIds( "displacement", XhDisp, markededges( mesh,listMarkedEdgesDisp ) );
    }
    // on marked points
    auto const& listMarkedPointsDisp = std::get<2>( meshMarkersDispByEntities );
    if ( !listMarkedPointsDisp.empty() )
        this->updateDofEliminationIds( "displacement", XhDisp, markedpoints( mesh,listMarkedPointsDisp ) );
    //-------------------------------------//
    // on disp components
    for ( auto const& meshMarkersPair : meshMarkersCompDispByEntities )
    {
        ComponentType comp = meshMarkersPair.first;
        // topological faces
        auto const& listMarkedFacesCompDisp = std::get<0>( meshMarkersPair.second );
        if ( !listMarkedFacesCompDisp.empty() )
        {
            auto therange = markedfaces(mesh,listMarkedFacesCompDisp );
            auto dofsToAdd = XhDisp->dofs( therange, comp );
            XhDisp->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll("displacement",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = XhDisp->dofs( therange, comp, true );
            this->dofEliminationIdsMultiProcess("displacement",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
        // edges (only 3d)
        auto const& listMarkedEdgesCompDisp = std::get<1>( meshMarkersPair.second );
        if ( !listMarkedEdgesCompDisp.empty() )
        {
            auto therange = markededges(mesh,listMarkedEdgesCompDisp );
            auto dofsToAdd = XhDisp->dofs( therange, comp );
            XhDisp->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll("displacement",MESH_EDGES).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = XhDisp->dofs( therange, comp, true );
            this->dofEliminationIdsMultiProcess("displacement",MESH_EDGES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
        // points
        auto const& listMarkedPointsCompDisp = std::get<2>( meshMarkersPair.second );
        if ( !listMarkedPointsCompDisp.empty() )
        {
            auto therange = markedpoints(mesh,listMarkedPointsCompDisp );
            auto dofsToAdd = XhDisp->dofs( therange, comp );
            XhDisp->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll("displacement",MESH_POINTS).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = XhDisp->dofs( therange, comp, true );
            this->dofEliminationIdsMultiProcess("displacement",MESH_POINTS).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
    }


}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("SolidMechanics","initMesh", "start" );
    this->timerTool("Constructor").start();

    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    this->timerTool("Constructor").stop("createMesh");
    this->log("SolidMechanics","initMesh", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("SolidMechanics","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("SolidMechanics","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("SolidMechanics","initFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    auto mom = this->materialsProperties()->materialsOnMesh(this->mesh());
    //--------------------------------------------------------//
    // function space for displacement
    if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_XhDisplacement = space_displacement_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_XhDisplacement = space_displacement_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),
                                                         _range=M_rangeMeshElements );
    }
    // field displacement
    M_fieldDisplacement.reset( new element_displacement_type( M_XhDisplacement, "structure displacement" ));

    //--------------------------------------------------------//
    // init pressure space with displcaement-pressure formulation
    std::set<std::string> physicUseDisplacementPressureFormulation;
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        if ( physicSolidData->useDisplacementPressureFormulation() )
            physicUseDisplacementPressureFormulation.insert(physicName);
    }
    if ( !physicUseDisplacementPressureFormulation.empty() )
    {
        if ( mom->isDefinedOnWholeMesh( physicUseDisplacementPressureFormulation ) )
            M_XhPressure = space_pressure_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
        else
        {
            //auto matNamesWithDisplacementPressureFormulation = this->materialsProperties()->physicToMaterials( physicUseDisplacementPressureFormulation );
            auto rangePressure = markedelements(this->mesh(), mom->markers( physicUseDisplacementPressureFormulation ) );
            M_XhPressure = space_pressure_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),
                                                     _range=rangePressure );
        }
        M_fieldPressure.reset( new element_pressure_type( M_XhPressure, "pressure" ) );
    }
    if ( M_timeSteppingUseMixedFormulation )
        M_fieldVelocity.reset( new element_displacement_type( M_XhDisplacement, "velocity" ));

    //--------------------------------------------------------//
    // pre-stress ( not functional )
    if (false)
        U_displ_struct_prestress.reset(new element_vectorial_type( M_XhDisplacement, "structure displacement prestress" ));
    //--------------------------------------------------------//

    this->timerTool("Constructor").stop("createSpaces");
    this->log("SolidMechanics","initFunctionSpaces", "finish" );

}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesNormalStress()
{
    if ( !M_XhNormalStress )
        M_XhNormalStress = space_normal_stress_type::New( _mesh=this->mesh(), _worldscomm=this->localNonCompositeWorldsComm() );
    if ( !M_fieldNormalStressFromStruct )
        M_fieldNormalStressFromStruct.reset( new element_normal_stress_type( M_XhNormalStress ) );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createExporters()
{
    this->log("SolidMechanics","createExporters", "start" );
    this->timerTool("Constructor").start();

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";//change_coords_only, change, static
    if (!M_isHOVisu)
    {
        if constexpr ( nOrderGeo <= 2 /*&& doExport*/ )
        {
            M_exporter = exporter( _mesh=this->mesh(),
                                   //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                                   _name="Export",
                                   _geo=geoExportType,
                                   _worldcomm=M_XhDisplacement->worldComm(),
                                   _path=this->exporterPath() );
        }
    }
    else
    {
#if 1 //defined(FEELPP_HAS_VTK)
        std::shared_ptr<mesh_visu_ho_type> meshVisuHO;
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        if ( hovisuSpaceUsed == "displacement" )
        {
            //auto Xh_create_ho = space_create_ho_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
            auto Xh_create_ho = this->functionSpaceDisplacement()->compSpace();

            bool doLagP1parallel=false;
            auto opLagP1 = lagrangeP1(_space=Xh_create_ho,
                                      _backend=this->algebraicBackend(),
                                      //_worldscomm=this->localNonCompositeWorldsComm(),
                                      _path=this->rootRepository(),
                                      _prefix=this->prefix(),
                                      _rebuild=!this->doRestart(),
                                      _parallel=doLagP1parallel );
            meshVisuHO = opLagP1->mesh();
        }
        else if ( hovisuSpaceUsed == "pressure" )
        {
            CHECK( false ) << "not implement\n";
        }
        else if ( hovisuSpaceUsed == "p1" )
        {
            meshVisuHO = this->mesh()->createP1mesh();
        }
        else CHECK( false ) << "invalid hovisu.space-used " << hovisuSpaceUsed;

        M_exporter_ho = exporter( _mesh=meshVisuHO,
                                  //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"ExportHO")),
                                  _name="ExportHO",
                                  _geo=geoExportType,
                                  _worldcomm=M_XhDisplacement->worldComm(),
                                  _path=this->exporterPath() );

        M_XhVectorialVisuHO = space_vectorial_visu_ho_type::New( _mesh=meshVisuHO, _worldscomm=this->localNonCompositeWorldsComm());
        M_displacementVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"u_visuHO"));

        M_opIdisplacement = opInterpolation(_domainSpace=this->functionSpaceDisplacement(),
                                            _imageSpace=M_XhVectorialVisuHO,
                                            _range=elements(M_XhVectorialVisuHO->mesh()),
                                            _backend=this->algebraicBackend(),
                                            _type=InterpolationNonConforme(false) );

#if 0
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::NormalStress ) )
        {
            this->createAdditionalFunctionSpacesNormalStress();
            M_opInormalstress = opInterpolation(_domainSpace=M_XhNormalStress,
                                                _imageSpace=M_XhVectorialVisuHO,
                                                //_range=elements(M_XhVectorialVisuHO->mesh()),
                                                _range=boundaryfaces(M_XhVectorialVisuHO->mesh()),
                                                _backend=this->algebraicBackend(),
                                                _type=InterpolationNonConforme(false) );
        }
#endif
        if ( this->hasDisplacementPressureFormulation() )
        {
            //M_XhScalarVisuHO = space_scalar_visu_ho_type::New(_mesh=opLagP1->mesh(), _worldscomm=this->localNonCompositeWorldsComm());
            M_XhScalarVisuHO = M_XhVectorialVisuHO->compSpace();
            M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));

            M_opIpressure = opInterpolation(_domainSpace=M_XhPressure,
                                            _imageSpace=M_XhScalarVisuHO,
                                            _range=elements(M_XhScalarVisuHO->mesh()),
                                            _backend=this->algebraicBackend(),
                                            _type=InterpolationNonConforme(false) );
        }
#endif // FEELPP_HAS_VTK
    }
    this->timerTool("Constructor").stop("createExporters");
    this->log("SolidMechanics","createExporters", "finish" );
}


//---------------------------------------------------------------------------------------------------//

namespace detail
{
template <typename SpaceType>
NullSpace<double> getNullSpace( SpaceType const& space, mpl::int_<2> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( vec(Py(),-Px()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
    return userNullSpace;
}
template <typename SpaceType>
NullSpace<double> getNullSpace( SpaceType const& space, mpl::int_<3> /**/ )
{
    // auto mode1 = space->element( oneX() );
    // auto mode2 = space->element( oneY() );
    // auto mode3 = space->element( oneZ() );
    // auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    // auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    // auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    auto mode1 = space->element();
    mode1.on( _range=space->template rangeElements<0>(), _expr=oneX() );
    auto mode2 = space->element();
    mode2.on( _range=space->template rangeElements<0>(), _expr=oneY() );
    auto mode3 = space->element();
    mode3.on( _range=space->template rangeElements<0>(), _expr=oneZ() );
    auto mode4 = space->element();
    mode4.on( _range=space->template rangeElements<0>(), _expr=vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element();
    mode5.on( _range=space->template rangeElements<0>(), _expr=vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element();
    mode6.on( _range=space->template rangeElements<0>(), _expr=vec(cst(0.),Pz(),-Py()) );

    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    return userNullSpace;
}
template< typename TheBackendType >
NullSpace<double> extendNullSpace( NullSpace<double> const& ns,
                                   std::shared_ptr<TheBackendType> const& mybackend,
                                   std::shared_ptr<DataMap<>> const& dm )
{
    std::vector< typename NullSpace<double>::vector_ptrtype > myvecbasis(ns.size());
    for ( int k=0;k< ns.size();++k )
    {
        // TODO : use method in vectorblock.hpp : BlocksBaseVector<T>::setVector (can be static)
        myvecbasis[k] = mybackend->newVector(dm);
        auto const& subvec = ns.basisVector(k);
        auto const& basisGpToContainerGpSubVec = subvec.map().dofIdToContainerId( 0 ); //only one
        auto const& basisGpToContainerGpVec = dm->dofIdToContainerId( 0 ); // space index of disp
        for ( int i=0;i<basisGpToContainerGpSubVec.size();++i )
            myvecbasis[k]->set( basisGpToContainerGpVec[i], subvec( basisGpToContainerGpSubVec[i] ) );
#if 0
        for( int i = 0 ; i < ns.basisVector(k).map().nLocalDofWithGhost() ; ++i )
            myvecbasis[k]->set(i, ns.basisVector(k)(i) );
#endif
        myvecbasis[k]->close();
    }
    NullSpace<double> userNullSpace( myvecbasis, mybackend );
    return userNullSpace;
}
} // detail

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createsSolid1dReduced()
{
    M_solid1dReduced.reset( new solid_1dreduced_type( this->prefix(), this->keyword()+"_1d", this->worldCommPtr(), "" , this->repository() ) );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildAlgebraicFactory )
{
    if ( this->isUpdatedForUse() ) return;

    this->log("SolidMechanics","init", "start" );
    this->timerTool("Constructor").start();
#if 0
    if ( M_modelName.empty() )
    {
        std::string theSolidModel = this->modelProperties().models().model( this->keyword() ).equations();
        this->setModelName( theSolidModel );
    }
    if ( M_solverName.empty() )
    {
        if ( M_modelName == "Elasticity" || M_modelName == "Generalised-String" )
            this->setSolverName( "LinearSystem" );
        else
            this->setSolverName( "Newton" );
    }
#endif
    if ( this->physics().empty() )
        this->initPhysics( this->keyword(), this->modelProperties().models() );

    this->initMaterialProperties();

    // backend
    this->initAlgebraicBackend();

    bool hasSolid1dReduced = false;
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        if ( physicSolidData->equation() == "Generalised-String" )
            hasSolid1dReduced = true;
    }
    if ( hasSolid1dReduced )
    {
        if ( !M_solid1dReduced )
            this->createsSolid1dReduced();
        M_solid1dReduced->setPhysics( this->physicsFromCurrentType(), this->physicDefault() );
        M_solid1dReduced->setMaterialsProperties( this->materialsProperties() );

        M_solid1dReduced->setManageParameterValues( false );
        //if ( !M_solid1dReduced->modelPropertiesPtr() )   // TODO : not build modelPropertiesPtr
        {
            M_solid1dReduced->setModelProperties( this->modelPropertiesPtr() );
            M_solid1dReduced->setManageParameterValuesOfModelProperties( false );
        }

        M_solid1dReduced->init();
    }

    if ( this->hasSolidEquationStandard() )
    {
        if ( !this->mesh() )
            this->initMesh();

        this->materialsProperties()->addMesh( this->mesh() );

        this->initFunctionSpaces();

        this->initBoundaryConditions();

        // start or restart time step scheme
        if ( !this->isStationary() )
            this->initTimeStep();

        // update parameters values
        this->modelProperties().parameters().updateParameterValues();
        // init function defined in json
        this->initUserFunctions();
        // init post-processinig (exporter, measure at point, ...)
        this->initPostProcess();
    }
    else
    {
        // update missing info
        if ( !this->isStationary() )
        {
            if ( this->hasSolidEquation1dReduced() )
            {
                // up initial time
                this->setTimeInitial( M_solid1dReduced->timeInitial() );
                // up current time
                this->updateTime( M_solid1dReduced->currentTime() );
            }
        }
    }

    // update constant parameters
    this->updateParameterValues();


    if ( M_solverName == "automatic" )
    {
        bool isLinear = true;
        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            if ( !this->materialsProperties()->hasPhysic( physicName ) )
                continue;
            auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
            if ( physicSolidData->equation() =="Hyper-Elasticity" )
            {
                isLinear = false;
                break;
            }
        }
        if ( isLinear )
            M_solverName = "LinearSystem";
        else
            M_solverName = "Newton";
    }

    // define start dof index
    size_type currentStartIndex = 0;
    // prepare block vector
    int nBlock = this->nBlockMatrixGraph();
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    int cptBlock = 0;
    // update block vector (index + data struct)
    if ( this->hasSolidEquationStandard() )
    {
        this->setStartSubBlockSpaceIndex( "displacement", currentStartIndex++ );
        if ( this->hasDisplacementPressureFormulation() )
            this->setStartSubBlockSpaceIndex( "pressure", currentStartIndex++ );
        if ( M_timeSteppingUseMixedFormulation )
            this->setStartSubBlockSpaceIndex( "velocity", currentStartIndex++ );

        bvs->operator()(cptBlock++) = this->fieldDisplacementPtr();
        if ( this->hasDisplacementPressureFormulation() )
            bvs->operator()(cptBlock++) = M_fieldPressure;
        if ( M_timeSteppingUseMixedFormulation )
            bvs->operator()(cptBlock++) = M_fieldVelocity;
    }
    if ( this->hasSolidEquation1dReduced() )
    {
        this->setStartSubBlockSpaceIndex( "solid-1dreduced", currentStartIndex );
        auto blockVectorSolution1dReduced = *M_solid1dReduced->algebraicBlockVectorSolution();
        int nBlock1dReduced = blockVectorSolution1dReduced.size();
        int numberOfBlockSpace1dReduced = 0;
        for ( int k=0;k<nBlock1dReduced ;++k )
        {
            bvs->operator()(cptBlock+k) = blockVectorSolution1dReduced(k);
            numberOfBlockSpace1dReduced += blockVectorSolution1dReduced(k)->map().numberOfDofIdToContainerId();
        }
        cptBlock += nBlock1dReduced;
        currentStartIndex += numberOfBlockSpace1dReduced;
    }
    // init vector representing all blocs
    bvs->buildVector( this->backend() );

    // update algebraic model
    if (buildAlgebraicFactory)
    {
        auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
        this->setAlgebraicFactory( algebraicFactory );

        if ( this->hasSolidEquationStandard() )
        {
            if ( this->nBlockMatrixGraph() == 1 )
            {
                NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceDisplacement(), mpl::int_<nDim>() ) ;
                if ( boption(_name="use-null-space",_prefix=this->prefix() ) )
                    algebraicFactory->attachNullSpace( userNullSpace );
                if ( boption(_name="use-near-null-space",_prefix=this->prefix() ) )
                    algebraicFactory->attachNearNullSpace( userNullSpace );
            }
            else
            {
                NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceDisplacement(), mpl::int_<nDim>() ) ;
                NullSpace<double> userNullSpaceFull = detail::extendNullSpace( userNullSpace, algebraicFactory->backend(), algebraicFactory->sparsityMatrixGraph()->mapRowPtr() );
                if ( boption(_name="use-near-null-space",_prefix=this->prefix() ) )
                {
                    algebraicFactory->attachNearNullSpace( 0,userNullSpace ); // for block disp in fieldsplit
                    algebraicFactory->attachNearNullSpace( userNullSpaceFull ); // for multigrid on full system
                }
            }

#if 0
            if ( true )
            {
                auto massbf = form2( _trial=this->functionSpaceDisplacement(), _test=this->functionSpaceDisplacement());
                auto const& u = this->fieldDisplacement();
                auto const& rho = this->mechanicalProperties()->fieldRho();
                massbf += integrate( _range=elements( this->mesh() ),
                                     _expr=idv(rho)*inner( idt(u),id(u) ) );
                massbf.close();
                M_algebraicFactory->preconditionerTool()->attachAuxiliarySparseMatrix( "mass-matrix", massbf.matrixPtr() );
            }
#endif

            if ( M_timeStepping == "Theta" )
            {
                M_timeStepThetaSchemePreviousContrib = this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
                algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
                algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
            }
        }
    }

    this->setIsUpdatedForUse( true );

    this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("SolidMechanics","init", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("SolidMechanics","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    fs::path saveTsDir = fs::path(this->rootRepository())/fs::path( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ts")) );

    if ( this->hasSolidEquationStandard() )
    {
        if ( M_timeStepping == "Newmark" )
        {
            M_timeStepNewmark = newmark( _space=M_XhDisplacement,
                                         _name="displacement"+suffixName,
                                         _prefix=this->prefix(),
                                         _initial_time=ti, _final_time=tf, _time_step=dt,
                                         _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                         _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
            M_timeStepNewmark->setfileFormat( myFileFormat );
            M_timeStepNewmark->setPathSave( ( saveTsDir/"displacement" ).string() );
            M_fieldVelocity = M_timeStepNewmark->currentVelocityPtr();
            M_fieldAcceleration = M_timeStepNewmark->currentAccelerationPtr();
        }
        else if ( M_timeStepping == "BDF" || M_timeStepping == "Theta" )
        {
            int bdfOrder = 1;
            if ( M_timeStepping == "BDF" )
                bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order");
            int nConsecutiveSave = std::max( 2, bdfOrder ); // at least 2 is required by fsi when restart
            M_timeStepBdfDisplacement = this->createBdf( M_XhDisplacement,"displacement", bdfOrder, nConsecutiveSave, myFileFormat );
            M_timeStepBdfVelocity = this->createBdf( M_XhDisplacement,"velocity", bdfOrder, nConsecutiveSave, myFileFormat );
        }

        if ( this->hasDisplacementPressureFormulation() )
        {
            M_savetsPressure = bdf( _space=this->functionSpacePressure(),
                                    _name="pressure"+suffixName,
                                    _prefix=this->prefix(),
                                    _order=1,
                                    _initial_time=ti, _final_time=tf, _time_step=dt,
                                    _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                    _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
            M_savetsPressure->setfileFormat( myFileFormat );
            M_savetsPressure->setPathSave( ( saveTsDir/"pressure" ).string() );
        }
    }
    else// if ( this->is1dReducedModel() )
    {
    }

    this->log("SolidMechanics","initTimeStep", "create time stepping object done" );

    // update current time, initial solution and restart
    if ( this->hasSolidEquationStandard() )
    {
        if ( !this->doRestart() )
        {
            if ( Environment::vm().count(prefixvm(this->prefix(),"time-initial.displacement.files.directory").c_str()) )
            {
                std::string initialDispFilename = Environment::expand( soption( _name="time-initial.displacement.files.directory",_prefix=this->prefix() ) );
                std::string saveType = soption( _name="time-initial.displacement.files.format",_prefix=this->prefix() );
                M_fieldDisplacement->load(_path=initialDispFilename, _type=saveType );
            }
            if ( M_timeStepping == "Newmark" )
                this->updateTime( M_timeStepNewmark->timeInitial() );
            else
                this->updateTime( M_timeStepBdfDisplacement->timeInitial() );
        }
        else // do a restart
        {
            if ( M_timeStepping == "Newmark" )
            {
                // restart time step
                double tir = M_timeStepNewmark->restart();

                // load a previous solution as current solution
                *M_fieldDisplacement = M_timeStepNewmark->previousUnknown();
                // up initial time
                this->setTimeInitial( tir );
                // up current time
                this->updateTime( tir );
            }
            else
            {
                double tir = M_timeStepBdfDisplacement->restart();
                *M_fieldDisplacement = M_timeStepBdfDisplacement->unknown(0);
                M_timeStepBdfVelocity->restart();
                *M_fieldVelocity = M_timeStepBdfVelocity->unknown(0);
                this->setTimeInitial( tir );
                this->updateTime( tir );
            }
            if ( this->hasDisplacementPressureFormulation() )
            {
                M_savetsPressure->restart();
                *M_fieldPressure = M_savetsPressure->unknown(0);
            }
        }

    }
    else // if (this->is1dReducedModel())
    {
    }

    this->timerTool("Constructor").stop("initTimeStep");
    this->log("SolidMechanics","initTimeStep", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initUserFunctions()
{
    if ( this->modelProperties().functions().empty() )
        return;

    for ( auto const& modelfunc : this->modelProperties().functions() )
    {
        auto const& funcData = modelfunc.second;
        std::string funcName = funcData.name();

        if ( funcData.isScalar() )
        {
            if ( this->hasFieldUserScalar( funcName ) )
                continue;
            M_fieldsUserScalar[funcName] = this->functionSpaceDisplacement()->compSpace()->elementPtr();
        }
        else if ( funcData.isVectorial2() )
        {
            if ( nDim != 2 ) continue;
            if ( this->hasFieldUserVectorial( funcName ) )
                continue;
            M_fieldsUserVectorial[funcName] = this->functionSpaceDisplacement()->elementPtr();
        }
        else if ( funcData.isVectorial3() )
        {
            if ( nDim != 3 ) continue;
            if ( this->hasFieldUserVectorial( funcName ) )
                continue;
            M_fieldsUserVectorial[funcName] = this->functionSpaceDisplacement()->elementPtr();
        }
    }

    this->updateUserFunctions();
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateUserFunctions( bool onlyExprWithTimeSymbol )
{
    if ( this->modelProperties().functions().empty() )
        return;

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().functions().setParameterValues( paramValues );
    for ( auto const& modelfunc : this->modelProperties().functions() )
    {
        auto const& funcData = modelfunc.second;
        if ( onlyExprWithTimeSymbol && !funcData.hasSymbol("t") )
            continue;

        std::string funcName = funcData.name();
        if ( funcData.isScalar() )
        {
            CHECK( this->hasFieldUserScalar( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserScalar[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionScalar() );
        }
        else if ( funcData.isVectorial2() )
        {
            if constexpr( nDim == 2 )
            {
                CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
                M_fieldsUserVectorial[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionVectorial2() );
            }
            else CHECK( false ) << "TODO";
        }
        else if ( funcData.isVectorial3() )
        {
            if constexpr( nDim == 3 )
            {
                CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
                M_fieldsUserVectorial[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionVectorial3() );
            }
            else CHECK( false ) << "TODO";
        }
    }
}

//---------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix ) const
{
    std::set<std::string> res;
    for ( auto const& o : ifields )
    {
        if ( o == prefixvm(prefix,"displacement") || o == prefixvm(prefix,"all") )
            res.insert( "displacement" );
        if ( o == prefixvm(prefix,"velocity") || o == prefixvm(prefix,"all") )
            res.insert( "velocity" );
        if ( o == prefixvm(prefix,"acceleration") || o == prefixvm(prefix,"all") )
            res.insert( "acceleration" );
        if ( o == prefixvm(prefix,"normal-stress") || o == prefixvm(prefix,"all") )
            res.insert( "normal-stress" );
        if ( o == prefixvm(prefix,"pid") || o == prefixvm(prefix,"all") )
            res.insert( "pid" );
        if ( o == prefixvm(prefix,"pressure") || o == prefixvm(prefix,"all") )
            res.insert( "pressure" );
        if ( o == prefixvm(prefix,"material-properties") || o == prefixvm(prefix,"all") )
            res.insert( "material-properties" );
        if ( o == prefixvm(prefix,"Von-Mises") || o == prefixvm(prefix,"all") )
            res.insert( "Von-Mises" );
        if ( o == prefixvm(prefix,"Tresca") || o == prefixvm(prefix,"all") )
            res.insert( "Tresca" );
        if ( o == prefixvm(prefix,"principal-stresses") || o == prefixvm(prefix,"all") )
            res.insert( "principal-stresses" );

        // add user functions
        if ( this->hasFieldUserScalar( o ) || this->hasFieldUserVectorial( o ) )
            res.insert( o );
    }
    return res;
}
#endif
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("SolidMechanics","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::set<std::string> fieldsAvailable = { "displacement", "von-mises-criterion", "tresca-criterion", "principal-stresses" };
    if ( !this->isStationary() )
        fieldsAvailable.insert( "velocity" );
    if ( this->hasDisplacementPressureFormulation() )
        fieldsAvailable.insert( "pressure" );
    this->setPostProcessExportsAllFieldsAvailable( fieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );

    std::set<std::string> saveFieldsAvailable = { "displacement" };
    if ( !this->isStationary() )
        saveFieldsAvailable.insert( "velocity" );
    if ( this->hasDisplacementPressureFormulation() )
        saveFieldsAvailable.insert( "pressure" );
    this->setPostProcessSaveAllFieldsAvailable( saveFieldsAvailable );


    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
    {
        if (this->hasSolidEquationStandard())
            this->createExporters();
        // restart exporter
        if (this->doRestart())
            this->restartExporters( this->timeInitial() );
    }

    auto const& ptree = this->modelProperties().postProcess().pTree( this->keyword() );

    // volume variation
    std::string ppTypeMeasures = "Measures";
    std::string ppTypeMeasuresVolumeVariation = "VolumeVariation";
    for( auto const& ptreeLevel0 : ptree )
    {
        std::string ptreeLevel0Name = ptreeLevel0.first;
        if ( ptreeLevel0Name != ppTypeMeasures ) continue;
        for( auto const& ptreeLevel1 : ptreeLevel0.second )
        {
            std::string ptreeLevel1Name = ptreeLevel1.first;
            if ( ptreeLevel1Name == ppTypeMeasuresVolumeVariation )
            {
                // TODO : init name and markers from ptree
                std::string vvname = "volume_variation";
                std::set<std::string> markers; markers.insert( "" );
                if ( !markers.empty() )
                {
                    M_postProcessVolumeVariation[vvname] = markers;
                }
            }
        }
    }

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

#if 0
    std::set<std::string> fieldNameStressScalar = { "Von-Mises","Tresca","princial-stress-1","princial-stress-2","princial-stress-3",
                                                    "stress_xx","stress_xy","stress_xz","stress_yx","stress_yy","stress_yz","stress_zx","stress_zy","stress_zz" };
    // points evaluation
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
    {
        if (!this->isStandardModel()) break;// TODO

        auto const& ptPos = evalPoints.pointPosition();
        node_type ptCoord(3);
        for ( int c=0;c<3;++c )
            ptCoord[c]=ptPos.value()(c);

        auto const& fields = evalPoints.fields();
        for ( std::string const& field : fields )
        {
            if ( field == "displacement" || field == "velocity" || field == "acceleration" )
            {
                if ( !M_postProcessMeasuresContextDisplacement )
                    M_postProcessMeasuresContextDisplacement.reset( new context_displacement_type( this->functionSpaceDisplacement()->context() ) );
                int ctxId = M_postProcessMeasuresContextDisplacement->nPoints();
                M_postProcessMeasuresContextDisplacement->add( ptCoord );
                std::string ptNameExport = (boost::format("%1%_%2%")%field %ptPos.name()).str();
                this->postProcessMeasuresEvaluatorContext().add( field, ctxId, ptNameExport );
            }
            else if ( field == "pressure" )
            {
                if ( !M_useDisplacementPressureFormulation )
                    continue;
                if ( !M_postProcessMeasuresContextPressure )
                    M_postProcessMeasuresContextPressure.reset( new context_pressure_type( this->functionSpacePressure()->context() ) );
                int ctxId = M_postProcessMeasuresContextPressure->nPoints();
                M_postProcessMeasuresContextPressure->add( ptCoord );
                std::string ptNameExport = (boost::format("pressure_%1%")%ptPos.name()).str();
                this->postProcessMeasuresEvaluatorContext().add("pressure", ctxId, ptNameExport );
            }
            else if ( fieldNameStressScalar.find( field ) != fieldNameStressScalar.end() )
            {
                this->createAdditionalFunctionSpacesStressTensor();
                if (!M_postProcessMeasuresContextStressScalar )
                    M_postProcessMeasuresContextStressScalar.reset( new context_stress_scal_type( M_XhStressTensor->compSpace()->context() ) );
                int ctxId = M_postProcessMeasuresContextStressScalar->nPoints();
                M_postProcessMeasuresContextStressScalar->add( ptCoord );
                std::string ptNameExport = (boost::format("%1%_%2%")%field %ptPos.name()).str();
                this->postProcessMeasuresEvaluatorContext().add( field, ctxId, ptNameExport );
            }

        }
    }
#endif
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::restartExporters( double time )
{
    // restart exporter
    if (this->doRestart() && this->restartPath().empty())
    {
        if ( this->hasSolidEquationStandard() )
        {
            if (!M_isHOVisu)
            {
                if ( M_exporter && M_exporter->doExport() )
                    M_exporter->restart(this->timeInitial());
            }
            else
            {
#if 1 // defined(FEELPP_HAS_VTK)
                if ( M_exporter_ho && M_exporter_ho->doExport() )
                    M_exporter_ho->restart(this->timeInitial());
                #endif
            }
        }
        else
        {
        }
    }
}


} //FeelModels

} // Feel




