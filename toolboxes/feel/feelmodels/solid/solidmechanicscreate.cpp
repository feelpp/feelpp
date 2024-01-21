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
    M_boundaryConditions = std::make_shared<boundary_conditions_type>( this->shared_from_this() );
    if ( !this->modelProperties().boundaryConditions().hasSection( this->keyword() ) )
        return;
    M_boundaryConditions->setup( this->modelProperties().boundaryConditions().section( this->keyword() ) );

    for ( auto const& [bcId,bcData] : M_boundaryConditions->displacementImposed() )
        bcData->updateDofEliminationIds( *this, "displacement", this->functionSpaceDisplacement() );
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

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
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
    std::set<typename materialsproperties_type::physic_id_type> physicUseDisplacementPressureFormulation;
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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models )
{
    auto currentPhysic = std::dynamic_pointer_cast<ModelPhysicSolid<nDim>>( physicsTree.physic() );
    CHECK( currentPhysic ) << "wrong physic";
    if ( currentPhysic->equation() == "Generalised-String" )
    {
        if ( !M_solid1dReduced )
            this->createsSolid1dReduced();
        physicsTree.addChild( M_solid1dReduced, models );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildAlgebraicFactory )
{
    if ( this->isUpdatedForUse() ) return;

    this->log("SolidMechanics","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();
#if 0
    this->initPhysics( this->keyword(), this->modelProperties().models() );
#elif 0
    this->initPhysics( this->keyword(), this->modelProperties().models() );
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
        typename super_physics_type::PhysicsTree physicsTree( this->shared_from_this() );
        physicsTree.addLeaf( M_solid1dReduced );
        this->initPhysics( physicsTree, this->modelProperties().models() );
    }
#else
    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );
#endif

    this->initMaterialProperties();

    // backend
    this->initAlgebraicBackend();

#if 0
    bool hasSolid1dReduced = false;
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        if ( physicSolidData->equation() == "Generalised-String" )
            hasSolid1dReduced = true;
    }
#endif
    if ( M_solid1dReduced )//hasSolid1dReduced )
    {
        if ( !M_solid1dReduced )
            this->createsSolid1dReduced();
        //M_solid1dReduced->setPhysics( this->physicsFromCurrentType() );
        M_solid1dReduced->setMaterialsProperties( this->materialsProperties() );

        M_solid1dReduced->setManageParameterValues( false );
        //if ( !M_solid1dReduced->modelPropertiesPtr() )   // TODO : not build modelPropertiesPtr
        {
            M_solid1dReduced->setModelProperties( this->modelPropertiesPtr() );
            M_solid1dReduced->setManageParameterValuesOfModelProperties( false );
        }

        M_solid1dReduced->init();
    }
    else M_solid1dReduced.reset();// not very nice

    if ( this->hasSolidEquationStandard() )
    {
        //if ( !this->mesh() )
        this->initMesh();

        this->materialsProperties()->addMesh( this->mesh() );

        this->initFunctionSpaces();

        this->initBoundaryConditions();

        // start or restart time step scheme
        if ( !this->isStationary() )
            this->initTimeStep();

        // update parameters values
        this->modelProperties().parameters().updateParameterValues();
        
        
        // update initial conditions
        std::cout << "Update initial conditions" << std::endl;
        this->updateInitialConditions( this->symbolsExpr() );
        
        
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



    if ( this->modelProperties().postProcess().hasJsonProperties( this->keyword() ) )
    {
        auto const& j_pp = this->modelProperties().postProcess().jsonProperties( this->keyword() );
        std::string ppTypeMeasures = "Measures";
        if ( j_pp.contains( ppTypeMeasures ) )
        {
            auto j_pp_measures = j_pp.at( ppTypeMeasures );
            for ( auto const& [j_pp_measureskey,j_pp_measuresval] : j_pp_measures.items() )
            {
                if ( j_pp_measureskey == "VolumeVariation" )
                {
                    ModelMarkers _markers;
                    _markers.setup( j_pp_measuresval /*,indexes*/);
                    std::string vvname = "volume_variation";
                    if ( !_markers.empty() )
                        M_postProcessVolumeVariation[vvname] = _markers;
                }
            }
        }
    }

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

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




