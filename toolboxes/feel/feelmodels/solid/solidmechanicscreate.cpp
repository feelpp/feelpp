/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>


namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::SolidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    worldcomm_ptr_t const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldComm, subPrefix, modelRep ),
    M_hasBuildFromMesh( false ), M_hasBuildFromMesh1dReduced( false ), M_isUpdatedForUse( false ),
    M_mechanicalProperties( new mechanicalproperties_type( prefix ) )
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
    // load info from .json file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // set of worldComm for the function spaces
    this->createWorldsComm();
    //-----------------------------------------------------------------------------//
    // build  mesh, space,exporter,...
    if (buildMesh) this->build();
    //-----------------------------------------------------------------------------//
    this->log("SolidMechanics","constructor", "finish" );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix,
                                         bool buildMesh,
                                         worldcomm_ptr_t const& worldComm,
                                         std::string const& subPrefix,
                                         ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, modelRep );
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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::build()
{
    /**/ if ( M_pdeType == "Elasticity" )                   { buildStandardModel(); }
    else if ( M_pdeType == "Elasticity-Large-Deformation" ) { buildStandardModel(); }
    else if ( M_pdeType == "Hyper-Elasticity" )             { buildStandardModel(); }
    else if ( M_pdeType == "Generalised-String" )           { build1dReducedModel(); }
    else
        CHECK(false) << "invalid pdeType" << M_pdeType;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::build( mesh_ptrtype mesh )
{
    if ( M_pdeType == "Elasticity" || M_pdeType == "Elasticity-Large-Deformation" || M_pdeType == "Hyper-Elasticity" )
        this->buildStandardModel(mesh);
    else
        CHECK(false) << "invalid pdeType" << M_pdeType << "need standard model";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::build( mesh_1dreduced_ptrtype mesh )
{
    if ( M_pdeType == "Generalised-String" )
        this->build1dReducedModel(mesh);
    else
        CHECK(false) << "invalid pdeType" << M_pdeType << " (need 1d_reduced model)";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::buildStandardModel( mesh_ptrtype mesh )
{
    this->log("SolidMechanics","buildStandardModel", "start" );
    //-----------------------------------------------------------------------------//
    M_isStandardModel=true;
    M_is1dReducedModel=false;
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    if ( mesh )
        M_mesh = mesh;
    else
        this->createMesh();
    //-----------------------------------------------------------------------------//
    // functionSpaces and elements
    this->createFunctionSpaces();
    //-----------------------------------------------------------------------------//
    // time schema
    this->createTimeDiscretisation();
    //-----------------------------------------------------------------------------//
    // exporters
    this->createExporters();
    //-----------------------------------------------------------------------------//
    M_hasBuildFromMesh = true;
    this->log("SolidMechanics","buildStandardModel", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype mesh )
{
    if (this->doRestart())
        this->build();
    else
        this->build(mesh);
}
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadMesh( mesh_1dreduced_ptrtype mesh )
{
    // no restart case here :
    // we consider that the mesh given is identicaly (from createsubmesh of fluid mesh)
    this->build(mesh);
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::build1dReducedModel( mesh_1dreduced_ptrtype mesh )
{
    this->log("SolidMechanics","build1dReducedModel", "start" );
    //-----------------------------------------------------------------------------//
    M_isStandardModel=false;
    M_is1dReducedModel=true;
    //-----------------------------------------------------------------------------//
    if ( mesh )
        M_mesh_1dReduced = mesh;
    else
        this->createMesh1dReduced();
    //-----------------------------------------------------------------------------//
    this->createFunctionSpaces1dReduced();
    //-----------------------------------------------------------------------------//
    this->createTimeDiscretisation1dReduced();
    //-----------------------------------------------------------------------------//
    this->createExporters1dReduced();
    //-----------------------------------------------------------------------------//
    M_hasBuildFromMesh1dReduced = true;
    this->log("SolidMechanics","build1dReducedModel", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    this->log("SolidMechanics","loadParameterFromOptionsVm", "start" );

    std::string formulation = soption(_name="formulation",_prefix=this->prefix());
    M_useDisplacementPressureFormulation = false;
    if ( formulation == "displacement-pressure" )
        M_useDisplacementPressureFormulation = true;
    M_mechanicalProperties->setUseDisplacementPressureFormulation(M_useDisplacementPressureFormulation);

    std::string theSolidModel = this->modelProperties().models().model("fluid").equations();
    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        theSolidModel = soption(_name="model",_prefix=this->prefix());
    this->pdeType( theSolidModel );

    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        M_pdeSolver = soption(_name="solver",_prefix=this->prefix());

    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "dirichlet-neumann";

    M_isHOVisu = nOrderGeo > 1;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());

    // overwrite export field options in json if given in cfg
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_displacement").c_str()) )
        if ( boption(_name="do_export_displacement",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Displacement );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_velocity").c_str()) )
        if ( boption(_name="do_export_velocity",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Velocity );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_acceleration").c_str()) )
        if ( boption(_name="do_export_acceleration",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Acceleration );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_normalstress").c_str()) )
        if ( boption(_name="do_export_normalstress",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::NormalStress );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_pressure").c_str()) )
        if ( boption(_name="do_export_pressure",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pressure );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_material_properties").c_str()) )
        if ( boption(_name="do_export_material_properties",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::MaterialProperties );
    //if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_velocityinterfacefromfluid").c_str()) )
    //    M_doExportVelocityInterfaceFromFluid = boption(_name="do_export_velocityinterfacefromfluid",_prefix=this->prefix());
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_all").c_str()) )
        if ( boption(_name="do_export_all",_prefix=this->prefix()) )
        {
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Displacement );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Velocity );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Acceleration );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::NormalStress );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pressure );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::MaterialProperties );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::FSI );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pid );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::VonMises );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Tresca );
            this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::PrincipalStresses );
        }

    //time schema parameters
    std::string timeSchema = soption(_name="time-schema",_prefix=this->prefix());
    if (timeSchema == "Newmark")
    {
        M_genAlpha_alpha_m=1.0;
        M_genAlpha_alpha_f=1.0;
    }
    else if (timeSchema == "Generalized-Alpha")
    {
        M_genAlpha_rho=doption(_name="time-rho",_prefix=this->prefix());
        M_genAlpha_alpha_m=(2.- M_genAlpha_rho)/(1.+M_genAlpha_rho);
        M_genAlpha_alpha_f=1./(1.+M_genAlpha_rho);
    }
    else CHECK( false ) << "time scheme not supported : " << timeSchema << "\n";

    M_genAlpha_gamma=0.5+M_genAlpha_alpha_m-M_genAlpha_alpha_f;
    M_genAlpha_beta=0.25*(1+M_genAlpha_alpha_m-M_genAlpha_alpha_f)*(1+M_genAlpha_alpha_m-M_genAlpha_alpha_f);

    M_useMassMatrixLumped = false;

    // axi-sym
    M_thickness_1dReduced = doption(_name="1dreduced-thickness",_prefix=this->prefix());
    M_radius_1dReduced = doption(_name="1dreduced-radius",_prefix=this->prefix());

    this->log("SolidMechanics","loadParameterFromOptionsVm", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createWorldsComm()
{
    this->log("SolidMechanics","createWorldsComm", "start" );

    if (this->worldComm().localSize()==this->worldComm().globalSize())
    {
        auto vecWorldComm{ makeWorldsComm(space_displacement_type::nSpaces,this->worldCommPtr())};
        auto vecLocalWorldComm{ makeWorldsComm(1,this->worldCommPtr()) };
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
    }
    else
    {
#if 0
        // manage world comm : WARNING only true without the lagrange multiplier
        const int DisplacementWorld=0;
        const int PressureWorld=1;
        int CurrentWorld=0;
        if (this->worldComm().globalRank() < this->worldComm().globalSize()/2 )
            CurrentWorld=DisplacementWorld;
        else
            CurrentWorld=PressureWorld;

        std::vector<WorldComm> vecWorldComm(2);
        std::vector<WorldComm> vecLocalWorldComm(1);
        if (this->worldComm().globalSize()>1)
        {
            vecWorldComm[0]=this->worldComm().subWorldComm(DisplacementWorld);
            vecWorldComm[1]=this->worldComm().subWorldComm(PressureWorld);
            vecLocalWorldComm[0]=this->worldComm().subWorldComm(CurrentWorld);
        }
        else
        {
            vecWorldComm[0]=WorldComm();
            vecWorldComm[1]=WorldComm();
            vecLocalWorldComm[0]=WorldComm();
        }
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
#else
        // non composite case
        worldscomm_ptr_t vecWorldComm{ makeWorldsComm(1,this->worldCommPtr()) };
        worldscomm_ptr_t vecLocalWorldComm{ makeWorldsComm(1,this->worldCommPtr()) };
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
#endif

    }

    this->log("SolidMechanics","createWorldsComm", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerFluidStructureInterfaceBC();
    this->clearMarkerRobinBC();

    std::string dirichletbcType = "elimination";//soption(_name="dirichletbc.type",_prefix=this->prefix());
    this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::NO_COMPONENT );
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        this->M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( (boost::format("displacement_%1%")%compTag).str(), "Dirichlet" );
        for( auto const& d : this->M_bcDirichletComponents.find(comp)->second )
            this->addMarkerDirichletBC( dirichletbcType, marker(d), comp );
    }

    this->M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "Neumann_scalar" );
    for( auto const& d : this->M_bcNeumannScalar )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,marker(d));
    this->M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "Neumann_vectorial" );
    for( auto const& d : this->M_bcNeumannVectorial )
        this->addMarkerNeumannBC(NeumannBCShape::VECTORIAL,marker(d));
    this->M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<nDim>( "displacement", "Neumann_tensor2" );
    for( auto const& d : this->M_bcNeumannTensor2 )
        this->addMarkerNeumannBC(NeumannBCShape::TENSOR2,marker(d));

    this->M_bcInterfaceFSI = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "interface_fsi" );
    for( auto const& d : this->M_bcInterfaceFSI )
        this->addMarkerFluidStructureInterfaceBC( marker(d) );

    this->M_bcRobin = this->modelProperties().boundaryConditions().template getVectorFieldsList<nDim>( "displacement", "robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    this->M_bcNeumannEulerianFrameScalar = this->modelProperties().boundaryConditions().getScalarFields( { { "displacement", "Neumann_eulerian_scalar" },{ "displacement", "FollowerPressure" } } );
    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
        this->addMarkerNeumannEulerianFrameBC(NeumannEulerianFrameBCShape::SCALAR,marker(d));
    this->M_bcNeumannEulerianFrameVectorial = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "Neumann_eulerian_vectorial" );
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
        this->addMarkerNeumannEulerianFrameBC(NeumannEulerianFrameBCShape::VECTORIAL,marker(d));
    this->M_bcNeumannEulerianFrameTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<nDim>( "displacement", "Neumann_eulerian_tensor2" );
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
        this->addMarkerNeumannEulerianFrameBC(NeumannEulerianFrameBCShape::TENSOR2,marker(d));

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "displacement", "VolumicForces" );

}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("SolidMechanics","createMesh", "start" );
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    this->timerTool("Constructor").stop("createMesh");
    this->log("SolidMechanics","createMesh", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createMesh1dReduced()
{
    this->log("SolidMechanics","createMesh1dReduced", "start" );
    std::string prefix1dreduced = prefixvm(this->prefix(),"1dreduced");

    std::string modelMeshRestartFile = prefixvm(this->prefix(),"SolidMechanics1dreducedMesh.path");
    std::string smpath = (fs::path( this->rootRepository() ) / fs::path( modelMeshRestartFile)).string();

    if (this->doRestart())
    {
        this->log("SolidMechanics","createMesh1dReduced", "reload mesh (because restart)" );

        if ( !this->restartPath().empty() )
            smpath = (fs::path( this->restartPath() ) / fs::path( modelMeshRestartFile)).string();

#if defined(SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH)
        auto meshSM2dClass = reloadMesh<mesh_type>(smpath,this->worldCommPtr());
        SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH(meshSM2dClass);
        M_mesh_1dReduced=mesh;
#else
        M_mesh_1dReduced = reloadMesh<mesh_1dreduced_type>(smpath,this->worldCommPtr());
#endif
    }
    else
    {
        if (Environment::vm().count(prefixvm(this->prefix(),"1dreduced-geofile")))
        {
            this->log("SolidMechanics","createMesh1dReduced", "use 1dreduced-geofile" );
            std::string geofile=soption(_name="1dreduced-geofile",_prefix=this->prefix() );
            std::string path = this->rootRepository();
            std::string mshfile = path + "/" + prefix1dreduced + ".msh";
            this->setMeshFile(mshfile);

            fs::path curPath=fs::current_path();
            bool hasChangedRep=false;
            if ( curPath != fs::path(this->rootRepository()) )
            {
                this->log("createMeshModel","", "change repository (temporary) for build mesh from geo : "+ this->rootRepository() );
                bool hasChangedRep=true;
                Environment::changeRepository( _directory=boost::format(this->rootRepository()), _subdir=false );
            }

            gmsh_ptrtype geodesc = geo( _filename=geofile,
                                        _prefix=prefix1dreduced,
                                        _worldcomm=this->worldCommPtr() );
            // allow to have a geo and msh file with a filename equal to prefix
            geodesc->setPrefix(prefix1dreduced);
            M_mesh_1dReduced = createGMSHMesh(_mesh=new mesh_1dreduced_type,_desc=geodesc,
                                              _prefix=prefix1dreduced,_worldcomm=this->worldCommPtr(),
                                              _partitions=this->worldComm().localSize() );

            // go back to previous repository
            if ( hasChangedRep )
                Environment::changeRepository( _directory=boost::format(curPath.string()), _subdir=false );
        }
        else
        {
            this->loadConfigMeshFile1dReduced( prefix1dreduced );
        }
        this->saveMeshFile( smpath );
    }

    this->log("SolidMechanics","createMesh1dReduced", "finish" );
} // createMesh1dReduced


//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    this->log("SolidMechanics","createFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().materials().setParameterValues( paramValues );
    M_mechanicalProperties->updateForUse( this->mesh(), this->modelProperties().materials(),  this->localNonCompositeWorldsComm() );

    //--------------------------------------------------------//
    // function space for displacement
    if ( M_mechanicalProperties->isDefinedOnWholeMesh() )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_XhDisplacement = space_displacement_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), M_mechanicalProperties->markers());
        M_XhDisplacement = space_displacement_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),
                                                         _range=M_rangeMeshElements );
    }
    // field displacement
    M_fieldDisplacement.reset( new element_displacement_type( M_XhDisplacement, "structure displacement" ));
    //--------------------------------------------------------//
    if ( M_useDisplacementPressureFormulation )
    {
        if ( M_mechanicalProperties->isDefinedOnWholeMesh() )
            M_XhPressure = space_pressure_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
        else
            M_XhPressure = space_pressure_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),
                                                     _range=M_rangeMeshElements );
        M_fieldPressure.reset( new element_pressure_type( M_XhPressure, "pressure" ) );
    }
    //subfunctionspace vectorial
    //M_XhVectorial = M_Xh;

    //--------------------------------------------------------//
    // pre-stress ( not functional )
    if (false)
        U_displ_struct_prestress.reset(new element_vectorial_type( M_XhDisplacement, "structure displacement prestress" ));
    //--------------------------------------------------------//

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_XhDisplacement->worldCommPtr() );

    this->timerTool("Constructor").stop("createSpaces");
    this->log("SolidMechanics","createFunctionSpaces", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpaces1dReduced()
{
    this->log("SolidMechanics","createFunctionSpaces1dReduced", "start" );

    // function space and elements
    M_Xh_vect_1dReduced = space_vect_1dreduced_type::New(_mesh=M_mesh_1dReduced,
                                                         _worldscomm=makeWorldsComm(1,M_mesh_1dReduced->worldCommPtr()));
    M_Xh_1dReduced = M_Xh_vect_1dReduced->compSpace();
    // scalar field
    M_disp_1dReduced.reset( new element_1dreduced_type( M_Xh_1dReduced, "structure displacement" ));
    M_velocity_1dReduced.reset( new element_1dreduced_type( M_Xh_1dReduced, "structure velocity" ));
    M_acceleration_1dReduced.reset( new element_1dreduced_type( M_Xh_1dReduced, "structure acceleration" ));
    // vectorial field
    M_disp_vect_1dReduced.reset(new element_vect_1dreduced_type( M_Xh_vect_1dReduced, "structure vect 1d displacement" ));
    M_velocity_vect_1dReduced.reset(new element_vect_1dreduced_type( M_Xh_vect_1dReduced, "velocity vect 1d displacement" ));

    // backend : use worldComm of Xh_1dReduced
    M_backend_1dReduced = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh_1dReduced->worldCommPtr() );

    this->log("SolidMechanics","createFunctionSpaces1dReduced", "finish" );
}


//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesNormalStress()
{
    if ( !M_XhNormalStress )
        M_XhNormalStress = space_normal_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    if ( !M_fieldNormalStressFromStruct )
        M_fieldNormalStressFromStruct.reset( new element_normal_stress_type( M_XhNormalStress ) );
}
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesStressTensor()
{
    if ( !M_XhStressTensor )
        M_XhStressTensor = space_stress_tensor_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    if ( !M_fieldStressTensor )
        M_fieldStressTensor.reset( new element_stress_tensor_type( M_XhStressTensor ) );
    if ( !M_fieldVonMisesCriterions )
        M_fieldVonMisesCriterions.reset( new element_stress_scal_type( M_XhStressTensor->compSpace() ) );
    if ( !M_fieldTrescaCriterions )
        M_fieldTrescaCriterions.reset( new element_stress_scal_type( M_XhStressTensor->compSpace() ) );
    if ( M_fieldsPrincipalStresses.size() != 3 )
        M_fieldsPrincipalStresses.resize( 3 );
    for (int d=0;d<M_fieldsPrincipalStresses.size();++d)
        if( !M_fieldsPrincipalStresses[d] )
            M_fieldsPrincipalStresses[d].reset( new element_stress_scal_type( M_XhStressTensor->compSpace() ) );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesFSI()
{
    if ( this->isStandardModel() )
        this->createAdditionalFunctionSpacesFSIStandard();
    else
        this->createAdditionalFunctionSpacesFSI1dReduced();
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesFSIStandard()
{
    this->log("SolidMechanics","createAdditionalFunctionSpacesFSIStandard", "start" );

    //--------------------------------------------------------//
    // function space for normal stress
    if ( !M_XhNormalStress )
        this->createAdditionalFunctionSpacesNormalStress();
    if ( !M_fieldNormalStressFromFluid )
        M_fieldNormalStressFromFluid.reset(new element_normal_stress_type( M_XhNormalStress, "normalStressBoundaryFromFluid" ));

    //--------------------------------------------------------//
    if ( !M_fieldVelocityInterfaceFromFluid && ( this->couplingFSIcondition() == "robin-neumann" ||
                                                 this->couplingFSIcondition() == "robin-robin" ||
                                                 this->couplingFSIcondition() == "robin-robin-genuine" ||
                                                 this->couplingFSIcondition() == "nitsche" ) )
        M_fieldVelocityInterfaceFromFluid.reset( new element_vectorial_type( M_XhDisplacement, "velocityInterfaceFromFluid" ));


    if ( !M_XhSubMeshDispFSI && ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
                                  this->couplingFSIcondition() == "nitsche" ) )
    {
        auto subfsimesh = createSubmesh(this->mesh(),markedfaces(this->mesh(),this->markerNameFSI()) );
        M_XhSubMeshDispFSI = space_tracemesh_disp_type::New( _mesh=subfsimesh, _worldscomm=this->localNonCompositeWorldsComm() );
        M_fieldSubMeshDispFSI.reset( new element_tracemesh_disp_type(M_XhSubMeshDispFSI) );
    }

    this->log("SolidMechanics","createAdditionalFunctionSpacesFSIStandard", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesFSI1dReduced()
{
    if ( M_XhStressVect_1dReduced ) return;
    this->log("SolidMechanics","createAdditionalFunctionSpacesFSI1dReduced", "start" );

    // normal stress as source term
    M_XhStressVect_1dReduced = space_stress_vect_1dreduced_type::New(_mesh=M_mesh_1dReduced,
                                                                       _worldscomm=makeWorldsComm(1,M_mesh_1dReduced->worldCommPtr()));
    M_stress_1dReduced.reset( new element_stress_scal_1dreduced_type( M_XhStressVect_1dReduced->compSpace(), "structure stress" ));
    M_stress_vect_1dReduced.reset(new element_stress_vect_1dreduced_type( M_XhStressVect_1dReduced, "stress 1d vect displacement" ));

    this->log("SolidMechanics","createAdditionalFunctionSpacesFSI1dReduced", "finish" );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createTimeDiscretisation()
{
    this->log("SolidMechanics","createTimeDiscretisation", "start" );
    this->timerTool("Constructor").start();

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_timeStepNewmark = newmark( _vm=Environment::vm(), _space=M_XhDisplacement,
                                      _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"newmark"+suffixName)),
                                      _prefix=this->prefix(),
                                      _initial_time=ti, _final_time=tf, _time_step=dt,
                                      _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                      _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_timeStepNewmark->setfileFormat( myFileFormat );
    M_timeStepNewmark->setPathSave( (fs::path(this->rootRepository()) /
                                          fs::path( prefixvm(this->prefix(), (boost::format("newmark_dt_%1%")%dt).str() ) ) ).string() );

    if ( M_useDisplacementPressureFormulation )
    {
        M_savetsPressure = bdf( _vm=Environment::vm(), _space=this->functionSpacePressure(),
                                _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure"+suffixName)),
                                _prefix=this->prefix(),
                                _initial_time=ti, _final_time=tf, _time_step=dt,
                                _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
        M_savetsPressure->setfileFormat( myFileFormat );
        M_savetsPressure->setPathSave( (fs::path(this->rootRepository()) /
                                        fs::path( prefixvm(this->prefix(),"save-pressure" ) ) ).string() );
    }


    this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("SolidMechanics","createTimeDiscretisation", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createTimeDiscretisation1dReduced()
{
    this->log("SolidMechanics","createTimeDiscretisation1dReduced", "start" );

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_newmark_displ_1dReduced = newmark( _vm=Environment::vm(),
                                          _space=M_Xh_1dReduced,
                                          _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"structure-1dreduced."+suffixName)),
                                          _prefix=this->prefix(),
                                          _initial_time=ti, _final_time=tf, _time_step=dt,
                                          _restart=this->doRestart(),_restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                          _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_newmark_displ_1dReduced->setfileFormat( myFileFormat );
    M_newmark_displ_1dReduced->setPathSave( (fs::path(this->rootRepository()) /
                                             fs::path( prefixvm(this->prefix(), (boost::format("newmark_dt_%1%")%dt).str() ) ) ).string() );

    this->log("SolidMechanics","createTimeDiscretisation1dReduced", "finish" );
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
        M_exporter = exporter( _mesh=this->mesh(),
                               //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                               _name="Export",
                               _geo=geoExportType,
                               _worldcomm=M_XhDisplacement->worldComm(),
                               _path=this->exporterPath() );
    }
    else
    {
#if 1 //defined(FEELPP_HAS_VTK)
        std::shared_ptr<mesh_visu_ho_type> meshVisuHO;
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        bool doLagP1parallel=false;
        if ( hovisuSpaceUsed == "displacement" )
        {
            //auto Xh_create_ho = space_create_ho_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
            auto Xh_create_ho = this->functionSpaceDisplacement()->compSpace();

            auto opLagP1 = lagrangeP1(_space=Xh_create_ho,
                                      _backend=M_backend,
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
                                            _backend=M_backend,
                                            _type=InterpolationNonConforme(false) );

        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::NormalStress ) )
        {
            this->createAdditionalFunctionSpacesNormalStress();
            M_opInormalstress = opInterpolation(_domainSpace=M_XhNormalStress,
                                                _imageSpace=M_XhVectorialVisuHO,
                                                //_range=elements(M_XhVectorialVisuHO->mesh()),
                                                _range=boundaryfaces(M_XhVectorialVisuHO->mesh()),
                                                _backend=M_backend,
                                                _type=InterpolationNonConforme(false) );
        }

        if ( M_useDisplacementPressureFormulation )
        {
            //M_XhScalarVisuHO = space_scalar_visu_ho_type::New(_mesh=opLagP1->mesh(), _worldscomm=this->localNonCompositeWorldsComm());
            M_XhScalarVisuHO = M_XhVectorialVisuHO->compSpace();
            M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));

            M_opIpressure = opInterpolation(_domainSpace=M_XhPressure,
                                            _imageSpace=M_XhScalarVisuHO,
                                            _range=elements(M_XhScalarVisuHO->mesh()),
                                            _backend=M_backend,
                                            _type=InterpolationNonConforme(false) );
        }
#endif // FEELPP_HAS_VTK
    }
    this->timerTool("Constructor").stop("createExporters");
    this->log("SolidMechanics","createExporters", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::createExporters1dReduced()
{
    this->log("SolidMechanics","createExporters1dReduced", "start" );

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";

    if (!M_isHOVisu)
    {
        M_exporter_1dReduced = exporter( _mesh=this->mesh1dReduced(),
                                      //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export-1dReduced")),
                                      _name="Export-1dReduced",
                                      _geo=geoExportType,
                                      _worldcomm=M_Xh_1dReduced->worldComm(),
                                      _path=this->exporterPath() );
    }

    this->log("SolidMechanics","createExporters1dReduced", "finish" );
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
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( oneZ() );
    auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    return userNullSpace;
}
template< typename TheBackendType >
NullSpace<double> extendNullSpace( NullSpace<double> const& ns,
                                   std::shared_ptr<TheBackendType> const& mybackend,
                                   std::shared_ptr<DataMap> const& dm )
{
    std::vector< typename NullSpace<double>::vector_ptrtype > myvecbasis(ns.size());
    for ( int k=0;k< ns.size();++k )
    {
        myvecbasis[k] = mybackend->newVector(dm);
        for( int i = 0 ; i < ns.basisVector(k).map().nLocalDofWithGhost() ; ++i )
            myvecbasis[k]->set(i, ns.basisVector(k)(i) );
        myvecbasis[k]->close();
    }
    NullSpace<double> userNullSpace( myvecbasis, mybackend );
    return userNullSpace;
}
} // detail

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildAlgebraicFactory )
{
    if ( M_isUpdatedForUse ) return;

    this->log("SolidMechanics","init", "start" );
    this->timerTool("Constructor").start();

    if ( (this->isStandardModel() && !M_hasBuildFromMesh ) ||
         (this->is1dReducedModel() && !M_hasBuildFromMesh1dReduced ) )
        this->build();

    if ( this->markerNameFSI().size()>0 )
        this->createAdditionalFunctionSpacesFSI();

    //-------------------------------------------------//
    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // update parameters values
    this->modelProperties().parameters().updateParameterValues();
    // init function defined in json
    this->initUserFunctions();
    // init post-processinig (exporter, measure at point, ...)
    this->initPostProcess();

    this->updateBoundaryConditionsForUse();

    // update block vector (index + data struct)
    if (this->isStandardModel())
    {
        // define start dof index
        size_type currentStartIndex = 0;
        this->setStartSubBlockSpaceIndex( "displacement", currentStartIndex++ );
        if ( M_useDisplacementPressureFormulation )
            this->setStartSubBlockSpaceIndex( "pressure", currentStartIndex++ );

        // prepare block vector
        int nBlock = this->nBlockMatrixGraph();
        M_blockVectorSolution.resize( nBlock );
        M_blockVectorSolution(0) = this->fieldDisplacementPtr();
        int cptBlock=1;
        if ( M_useDisplacementPressureFormulation )
        {
            M_blockVectorSolution(cptBlock) = M_fieldPressure;
            ++cptBlock;
        }
        // init vector associated to the block
        M_blockVectorSolution.buildVector( this->backend() );
    }

    // update algebraic model
    if (buildAlgebraicFactory)
    {
        if (this->isStandardModel())
        {
            M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );

            if ( this->nBlockMatrixGraph() == 1 )
            {
                NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceDisplacement(), mpl::int_<nDim>() ) ;
                if ( boption(_name="use-null-space",_prefix=this->prefix() ) )
                    M_algebraicFactory->attachNullSpace( userNullSpace );
                if ( boption(_name="use-near-null-space",_prefix=this->prefix() ) )
                    M_algebraicFactory->attachNearNullSpace( userNullSpace );
            }
            else
            {
                NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceDisplacement(), mpl::int_<nDim>() ) ;
                NullSpace<double> userNullSpaceFull = detail::extendNullSpace( userNullSpace, M_algebraicFactory->backend(), M_algebraicFactory->sparsityMatrixGraph()->mapRowPtr() );
                if ( boption(_name="use-near-null-space",_prefix=this->prefix() ) )
                {
                    M_algebraicFactory->attachNearNullSpace( 0,userNullSpace ); // for block disp in fieldsplit
                    M_algebraicFactory->attachNearNullSpace( userNullSpaceFull ); // for multigrid on full system
                }
            }

        }
        else if (this->is1dReducedModel())
        {
            M_algebraicFactory_1dReduced.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend1dReduced()) );
            M_algebraicFactory_1dReduced->initFromFunctionSpace( this->functionSpace1dReduced() );
        }
    }

    M_isUpdatedForUse = true;

    this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("SolidMechanics","init", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    // update timediscr and exporters
    if (!this->doRestart())
    {
        if (this->isStandardModel())
        {
            if ( Environment::vm().count(prefixvm(this->prefix(),"time-initial.displacement.files.directory").c_str()) )
            {
                std::string initialDispFilename = Environment::expand( soption( _name="time-initial.displacement.files.directory",_prefix=this->prefix() ) );
                std::string saveType = soption( _name="time-initial.displacement.files.format",_prefix=this->prefix() );
                M_fieldDisplacement->load(_path=initialDispFilename, _type=saveType );
            }

            // start time step
            M_timeStepNewmark->start(*M_fieldDisplacement);
            if ( M_useDisplacementPressureFormulation ) M_savetsPressure->start( M_XhPressure );
            // up current time
            this->updateTime( M_timeStepNewmark->time() );
        }
        else if (this->is1dReducedModel())
        {
            // start time step
            M_newmark_displ_1dReduced->start(*M_disp_1dReduced);
            // up current time
            this->updateTime( M_newmark_displ_1dReduced->time() );
        }
    }
    else // do a restart
    {
        if (this->isStandardModel())
        {
            // restart time step
            M_timeStepNewmark->restart();
            if ( M_useDisplacementPressureFormulation ) M_savetsPressure->restart();
            // load a previous solution as current solution
            *M_fieldDisplacement = M_timeStepNewmark->previousUnknown();
            if ( M_useDisplacementPressureFormulation ) *M_fieldPressure = M_savetsPressure->unknown(0);
            // up initial time
            this->setTimeInitial( M_timeStepNewmark->timeInitial() );
            // up current time
            this->updateTime( M_timeStepNewmark->time() );

        }
        else  if (this->is1dReducedModel())
        {
            // restart time step
            M_newmark_displ_1dReduced->restart();
            // load a previous solution as current solution
            *M_disp_1dReduced = M_newmark_displ_1dReduced->previousUnknown();
            this->updateInterfaceDispFrom1dDisp();
            // up initial time
            this->setTimeInitial(M_newmark_displ_1dReduced->timeInitial());
            // up current time
            this->updateTime( M_newmark_displ_1dReduced->time() );
        }
    }

    this->log("SolidMechanics","init", "start/restart timeStep scheme done" );


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
            M_fieldsUserScalar[funcName]->on(_range=elements(this->mesh()),_expr=funcData.expressionScalar() );
        }
        else if ( funcData.isVectorial2() )
        {
            if ( nDim != 2 ) continue;
            CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserVectorial[funcName]->on(_range=elements(this->mesh()),_expr=funcData.expressionVectorial2() );
        }
        else if ( funcData.isVectorial3() )
        {
            if ( nDim != 3 ) continue;
            CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserVectorial[funcName]->on(_range=elements(this->mesh()),_expr=funcData.expressionVectorial3() );
        }
    }
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    std::string modelName = "solid";

    // update post-process expression
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    for ( auto const& o : this->modelProperties().postProcess().exports( modelName ).fields() )
    {
        if ( o == "displacement" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Displacement );
        if ( o == "velocity" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Velocity );
        if ( o == "acceleration" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Acceleration );
        if ( o == "stress" || o == "normal-stress" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::NormalStress );
        if ( o == "pressure" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pressure );
        if ( o == "material-properties" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::MaterialProperties );
        if ( o == "pid" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pid );
        if ( o == "fsi" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::FSI );
        if ( o == "Von-Mises" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::VonMises );
        if ( o == "Tresca" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Tresca );
        if ( o == "principal-stresses" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::PrincipalStresses );

        // add user functions
        if ( this->hasFieldUserScalar( o ) || this->hasFieldUserVectorial( o ) )
            M_postProcessUserFieldExported.insert( o );
    }
    // clean doExport with fields not available
    if ( !M_useDisplacementPressureFormulation )
        M_postProcessFieldExported.erase( SolidMechanicsPostProcessFieldExported::Pressure );
    if ( this->is1dReducedModel() )
        M_postProcessFieldExported.erase( SolidMechanicsPostProcessFieldExported::NormalStress );
#if 0
    if ( !this->fieldVelocityInterfaceFromFluidPtr() )
        M_postProcessFieldExported.erase( SolidMechanicsPostProcessFieldExported::FSI );
#endif

    if (this->isStandardModel())
        this->createExporters();
    else  if (this->is1dReducedModel())
        this->createExporters1dReduced();

    // restart exporter
    if (this->doRestart())
        this->restartExporters( this->timeInitial() );


    auto const& ptree = this->modelProperties().postProcess().pTree( modelName );

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

    std::set<std::string> fieldNameStressScalar = { "Von-Mises","Tresca","princial-stress-1","princial-stress-2","princial-stress-3",
                                                    "stress_xx","stress_xy","stress_xz","stress_yx","stress_yy","stress_yz","stress_zx","stress_zy","stress_zz" };
    // points evaluation
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( modelName ) )
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

    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just for have time in first column
    }

}

} //FeelModels

} // Feel




