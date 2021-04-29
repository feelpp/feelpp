/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

//#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

namespace Feel
{
namespace FeelModels
{

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::CoefficientFormPDE( typename super_type::super2_type::infos_type const& infosPDE,
                                                            std::string const& prefix,
                                                            std::string const& keyword,
                                                            worldcomm_ptr_t const& worldComm,
                                                            std::string const& subPrefix,
                                                            ModelBaseRepository const& modelRep )
    :
    super_type( infosPDE, prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep, ModelBaseCommandLineOptions( coefficientformpde_options( prefix ) ) )
{}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("CoefficientFormPDE","init", "start" );
    this->timerTool("Constructor").start();

    this->initMaterialProperties();

    if ( !this->mesh() )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // stabilization gls
    if ( this->M_applyStabilization && !this->M_stabilizationGLSParameter )
    {
        typedef StabilizationGLSParameter<mesh_type, nOrderUnknown> stab_gls_parameter_impl_type;
        this->M_stabilizationGLSParameter.reset( new stab_gls_parameter_impl_type( this->mesh(), prefixvm(this->prefix(),"stabilization.gls.parameter"), this->clovm() ) );
        this->M_stabilizationGLSParameter->init();
    }

    // update constant parameters into
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // post-process
    this->initPostProcess();


    // backend
    this->M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr(), this->clovm() );

    // subspaces index
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( this->unknownName(), currentStartIndex++ );

     // vector solution
    this->M_blockVectorSolution.resize( 1 );
    this->M_blockVectorSolution(0) = this->fieldUnknownPtr();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("CoefficientFormPDE","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}




COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("CoefficientFormPDE","initFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    // functionspace
    if ( mom->isDefinedOnWholeMesh( this->physic() ) )
    {
        this->M_rangeMeshElements = elements(this->mesh());
        M_Xh = space_unknown_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        this->M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physic() ));
        M_Xh = space_unknown_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=this->M_rangeMeshElements );
    }

    M_fieldUnknown.reset( new element_unknown_type( M_Xh,this->unknownName() ) );

    double tElpased = this->timerTool("Constructor").stop("initFunctionSpaces");
    this->log("CoefficientFormPDE","initFunctionSpaces",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_bcDirichletMarkerManagement.clearMarkerDirichletBC();
    M_bcNeumannMarkerManagement.clearMarkerNeumannBC();
    M_bcRobinMarkerManagement.clearMarkerRobinBC();

    if constexpr ( unknown_is_scalar )
        M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( { { this->physic(), std::string("Dirichlet") } } );
    else
        M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( { { this->physic(), std::string("Dirichlet") } } );

    for( auto const& d : this->M_bcDirichlet )
        M_bcDirichletMarkerManagement.addMarkerDirichletBC("elimination", name(d), markers(d) );

    if constexpr ( unknown_is_vectorial )
    {
        for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
        {
            std::string compTag = ( comp == ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
            this->M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( { { this->physic(), (boost::format("Dirichlet_%1%")%compTag).str() } } );
            for( auto const& d : this->M_bcDirichletComponents.find(comp)->second )
                M_bcDirichletMarkerManagement.addMarkerDirichletBC( "elimination", name(d), markers(d), comp );
        }
    }

    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( { { this->physic(),  std::string("Neumann") }  } );
    for( auto const& d : this->M_bcNeumann )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d),markers(d));

    std::string tmp = this->physic();
    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( std::move(tmp), "Robin" );
    for( auto const& d : this->M_bcRobin )
        M_bcRobinMarkerManagement.addMarkerRobinBC( name(d),markers(d) );

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();



    // TODO : improve lisibility/usablity by adding a struct for the marker distribution
    M_meshMarkersDofEliminationUnknown.clear();
    std::map<ComponentType, std::decay_t<decltype(detail::distributeMarkerListOnSubEntity( mesh, std::set<std::string>{} ) )> > meshMarkersUnknownByEntitiesAndComp;
    // strong Dirichlet bc
    std::set<std::string> unknownMarkers;
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) );
        if ( listMark.empty() )
            continue;
        unknownMarkers.insert( listMark.begin(), listMark.end() );
        auto markersDistrib = detail::distributeMarkerListOnSubEntity( mesh, std::set<std::string>( listMark.begin(), listMark.end() ) );
        M_meshMarkersDofEliminationUnknown[ComponentType::NO_COMPONENT].emplace( std::make_pair( name(d), std::move( markersDistrib ) ) );
    }
    if ( !unknownMarkers.empty() )
        meshMarkersUnknownByEntitiesAndComp[ComponentType::NO_COMPONENT] = detail::distributeMarkerListOnSubEntity( mesh, unknownMarkers );

    // strong Dirichlet bc on unknown component
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        unknownMarkers.clear();
        for( auto const& d : bcDirComp.second )
        {
            auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d), comp );
            if ( listMark.empty() )
                continue;
            unknownMarkers.insert( listMark.begin(), listMark.end() );
            auto markersDistrib = detail::distributeMarkerListOnSubEntity( mesh, std::set<std::string>( listMark.begin(), listMark.end() ) );
            M_meshMarkersDofEliminationUnknown[comp].emplace( std::make_pair( name(d), std::move( markersDistrib ) ) );
        }
        if ( !unknownMarkers.empty() )
            meshMarkersUnknownByEntitiesAndComp[comp] = detail::distributeMarkerListOnSubEntity( mesh, unknownMarkers );
    }

    // update dof eliminition id
    for ( auto const& [comp,meshMarkersUnknownByEntities] : meshMarkersUnknownByEntitiesAndComp )
    {
        // on element
        auto const& listMarkedElementUnknown = std::get<3>( meshMarkersUnknownByEntities );
        if ( !listMarkedElementUnknown.empty() )
            this->updateDofEliminationIds( this->unknownName(), Xh, markedelements( mesh,listMarkedElementUnknown ), comp );
        // on topological faces
        auto const& listMarkedFacesUnknown = std::get<0>( meshMarkersUnknownByEntities );
        if ( !listMarkedFacesUnknown.empty() )
            this->updateDofEliminationIds( this->unknownName(), Xh, markedfaces( mesh,listMarkedFacesUnknown ), comp );
        // on marked edges (only 3d)
        if constexpr ( nDim == 3)
        {
            auto const& listMarkedEdgesUnknown = std::get<1>( meshMarkersUnknownByEntities );
            if ( !listMarkedEdgesUnknown.empty() )
                this->updateDofEliminationIds( this->unknownName(), Xh, markededges( mesh,listMarkedEdgesUnknown ), comp );
        }
        // on marked points
        auto const& listMarkedPointsUnknown = std::get<2>( meshMarkersUnknownByEntities );
        if ( !listMarkedPointsUnknown.empty() )
            this->updateDofEliminationIds( this->unknownName(), Xh, markedpoints( mesh,listMarkedPointsUnknown ), comp );
    }
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("CoefficientFormPDE","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    CHECK( this->timeStepping() == "BDF" || this->timeStepping() == "Theta" ) << "invalid time-stepping";

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix

    int bdfOrder = 1;
    if ( this->timeStepping() == "BDF" )
        bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order",_vm=this->clovm());
    int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme

    M_bdfUnknown = this->createBdf( this->spaceUnknown(),this->unknownName(), bdfOrder, nConsecutiveSave, myFileFormat );

    if (!this->doRestart())
    {
        // up current time
        this->updateTime( M_bdfUnknown->timeInitial() );
    }
    else
    {
        // start time step
        double tir = M_bdfUnknown->restart();
        // load a previous solution as current solution
        *this->fieldUnknownPtr() = M_bdfUnknown->unknown(0);
        // up initial time
        this->setTimeInitial( tir );
        // up current time
        this->updateTime( tir );
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("CoefficientFormPDE","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("CoefficientFormPDE","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->initBasePostProcess();

    // point measures
    auto fieldNamesWithSpaceUnknown = std::make_pair( std::set<std::string>({this->unknownName()}), this->spaceUnknown() );
    auto fieldNamesWithSpaces = hana::make_tuple( fieldNamesWithSpaceUnknown );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
        M_measurePointsEvaluation->init( evalPoints );

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("CoefficientFormPDE","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    // init petsc vector associated to the block
    this->M_blockVectorSolution.buildVector( this->backend() );

    this->M_algebraicFactory.reset( new typename super_type::model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = 1;//this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceUnknown(),
                                _trial=this->spaceUnknown() )->graph();
    return myblockGraph;
}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );


    // FunctionSpace
    nl::json subPt;
    subPt[this->unknownName()] = M_Xh->journalSection().to_string();
    p.emplace( "Function Spaces",  subPt );

    this->modelFields().updateInformationObject( p["Fields"] );

#if 0
    if ( M_algebraicFactory )
        M_algebraicFactory->updateInformationObject( p["Algebraic Solver"] );
#endif
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );
        nl::json::json_pointer jsonPointerSpaceUnknown( jsonInfoFunctionSpaces.at( this->unknownName() ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceUnknown ) )
            tabInfoFunctionSpaces->add( this->unknownName(), TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceUnknown ), tabInfoProp ) );
        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );

    return tabInfo;
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------Info : CoefficientFormPDE----------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- name      : " << this->physic()
           << "\n     -- time mode : " << std::string( (this->isStationary())?"Stationary":"Transient");
    *_ostr << "\n   Boundary conditions"
           << M_bcDirichletMarkerManagement.getInfoDirichletBC()
           << M_bcNeumannMarkerManagement.getInfoNeumannBC()
           << M_bcRobinMarkerManagement.getInfoRobinBC();
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();
#if 0
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << this->mesh()->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
#endif
    *_ostr << "\n   Space Unknown Discretization"
           << "\n     -- name of unknown         : " << this->unknownName()
           << "\n     -- symbol of unknown : " << this->unknownSymbol()
           << "\n     -- basis : " << this->unknownBasis()
           << "\n     -- number of dof : " << this->spaceUnknown()->nDof() << " (" << this->spaceUnknown()->nLocalDof() << ")";
    if ( !this->isStationary() )
    {
        *_ostr << "\n   Time Discretization"
               << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
               << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
               << "\n     -- time step    : " << this->timeStepBase()->timeStep()
               << "\n     -- type : " << this->timeStepping();
        if ( this->timeStepping() == "BDF" )
            *_ostr << " ( order=" << this->timeStepBdfUnknown()->timeOrder() << " )";
    }
    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";
#endif
    return _ostr;
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("CoefficientFormPDE","startTimeStep", "start");

    // start time step
    if (!this->doRestart())
        M_bdfUnknown->start( M_bdfUnknown->unknowns() );
     // up current time
    this->updateTime( M_bdfUnknown->time() );

    // update all expressions in bc or in house prec
    this->updateParameterValues();

    this->log("CoefficientFormPDE","startTimeStep", "finish");
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("CoefficientFormPDE","updateTimeStep", "start");
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    bool rebuildCstAssembly = false;
    if ( this->timeStepping() == "BDF" )
    {
        int previousTimeOrder = this->timeStepBdfUnknown()->timeOrder();
        M_bdfUnknown->next( this->fieldUnknown() );
        int currentTimeOrder = this->timeStepBdfUnknown()->timeOrder();
        rebuildCstAssembly = previousTimeOrder != currentTimeOrder && this->timeStepBase()->strategy() == TS_STRATEGY_DT_CONSTANT;
        this->updateTime( this->timeStepBdfUnknown()->time() );
    }
    else if ( this->timeStepping() == "Theta" )
    {
        M_bdfUnknown->next( this->fieldUnknown() );
        this->updateTime( this->timeStepBdfUnknown()->time() );
    }

    if ( rebuildCstAssembly )
        this->setNeedToRebuildCstPart(true);

    this->updateParameterValues();

    this->timerTool("TimeStepping").stop("updateTimeStep");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("CoefficientFormPDE","updateTimeStep", "finish");
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    this->log("CoefficientFormPDE","setParameterValues", "start");

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : this->M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );

    this->log("CoefficientFormPDE","setParameterValues", "finish");
}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::solve()
{
    // TODO
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
#if 0
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
#endif
}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( ModelAlgebraic::DataUpdateJacobian & data ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("CoefficientFormPDE","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( this->unknownName(), data );

    this->log("CoefficientFormPDE","updateJacobianDofElimination","finish" );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( ModelAlgebraic::DataUpdateResidual & data ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("CoefficientFormPDE","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( this->unknownName(), data );

    this->log("CoefficientFormPDE","updateResidualDofElimination","finish" );
}


} // namespace Feel
} // namespace FeelModels
