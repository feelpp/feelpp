//!

#include <feel/feelmodels/electromagnetic/electromagnetic.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/modelcore/utils.hpp>

namespace Feel::FeelModels
{

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::Electromagnetic( std::string const& prefix,
                                                      std::string const& keyword,
                                                      worldcomm_ptr_t const& worldComm,
                                                      ModelBaseRepository const& modelRep,
                                                      ModelBaseCommandLineOptions const& modelOptions)
    :
    super_type( prefix, keyword, worldComm, "", modelRep, modelOptions ),
    ModelPhysics<mesh_type::nDim>( "electromagnetic" ),
    ModelBase( prefix, keyword, worldComm, "", modelRep, modelOptions )
{
    this->log("Electromagnetic","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectromagneticConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectromagneticSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectromagneticPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectromagneticTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("Electromagnetic","constructor", "finish");
}


ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_solverName = soption(_prefix=this->prefix(),_name="solver");
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("Electromagnetic","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("Electromagnetic","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()


ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlockElectric = M_electricModel->startSubBlockSpaceIndices().size();//nBlockMatrixGraph();
    int nBlockMagnetic = M_magneticModel->startSubBlockSpaceIndices().size();//nBlockMatrixGraph();
    int nBlock = nBlockElectric + nBlockMagnetic;
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);

    int indexBlock=0;

    auto blockMatElectric = M_electricModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockElectric ;++tk1 )
        for (int tk2=0;tk2<nBlockElectric ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = blockMatElectric(tk1,tk2);
#if 0 // TODO coupling
    BlocksStencilPattern patCoupling1(1,nBlockHeat,size_type(Pattern::ZERO));
    patCoupling1(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(indexBlock,indexBlock+nBlockHeat) = stencil(_test=M_heatModel->spaceTemperature(),
                                                                     _trial=M_electricModel->spaceElectricPotential(),
                                                                     _pattern_block=patCoupling1,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    if ( true )
    {
        BlocksStencilPattern patCoupling2(nBlockHeat,1,size_type(Pattern::ZERO));
        patCoupling2(0,0) = size_type(Pattern::COUPLED);
        myblockGraph(indexBlock+nBlockHeat,indexBlock) = stencil(_test=M_electricModel->spaceElectricPotential(),
                                                                         _trial=M_heatModel->spaceTemperature(),
                                                                         _pattern_block=patCoupling2,
                                                                         _diag_is_nonzero=false,_close=false)->graph();
    }
#endif
    indexBlock += nBlockElectric;

    auto blockMatMagnetic = M_magneticModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockMagnetic ;++tk1 )
        for (int tk2=0;tk2<nBlockMagnetic ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = blockMatMagnetic(tk1,tk2);

    myblockGraph.close();

    return myblockGraph;
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models )
{
    if ( !M_electricModel )
    {
        M_electricModel = electric_model_type::New( _prefix=prefixvm(this->prefix(),"electric"),
                                                    _keyword="electric",
                                                    _worldcomm=this->worldCommPtr(),
                                                    _repository=this->repository(),
                                                    _vm=this->clovm() );
    }

    if ( !M_magneticModel )
    {
        M_magneticModel = magnetic_model_type::New( _prefix=prefixvm(this->prefix(),"magnetic"),
                                                    _keyword="magnetic",
                                                    _worldcomm=this->worldCommPtr(),
                                                    _repository=this->repository(),
                                                    _vm=this->clovm() );
    }

    physicsTree.addChild( M_electricModel, models );
    physicsTree.addChild( M_magneticModel, models );

    physicsTree.updateMaterialSupportFromChildren( "intersect" );
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("Electromagnetic","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();

    // physics
    this->initPhysics(
        this->shared_from_this(),
        [this]( typename super_physics_type::PhysicsTree & physicsTree ) {
            physicsTree.updatePhysics( this->shared_from_this(), this->modelProperties().models() );
            CHECK( M_electricModel && M_magneticModel ) << "missing initialization of electric and magnetic models";
            physicsTree.updatePhysics( M_electricModel, this->modelProperties().models() );
            physicsTree.updatePhysics( M_magneticModel, this->modelProperties().models() );
        } );

    // physical properties
    if ( !M_materialsProperties )
    {
        //auto paramValues = this->modelProperties().parameters().toParameterValues();
        //this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

        // init electric toolbox
    M_electricModel->setManageParameterValues( false );
    if ( !M_electricModel->modelPropertiesPtr() )
    {
        M_electricModel->setModelProperties( this->modelPropertiesPtr() );
        M_electricModel->setManageParameterValuesOfModelProperties( false );
    }
    M_electricModel->setModelMeshAsShared( this->modelMesh() );
    M_electricModel->setMaterialsProperties( M_materialsProperties );
    M_electricModel->init( false );

    // init magnetic toolbox
    M_magneticModel->setManageParameterValues( false );
    if ( !M_magneticModel->modelPropertiesPtr() )
    {
        M_magneticModel->setModelProperties( this->modelPropertiesPtr() );
        M_magneticModel->setManageParameterValuesOfModelProperties( false );
    }
    M_magneticModel->setModelMeshAsShared( this->modelMesh() );
    M_magneticModel->setMaterialsProperties( M_materialsProperties );
    M_magneticModel->init( false );


#if 0 // TODO
    M_modelName = "Electromagnetic";
    if ( M_solverName == "automatic" )
    {
        if ( this->materialsProperties()->hasElectricConductivityDependingOnSymbol( "heat_T" ) )
            M_solverName = "Newton";
        else
            M_solverName = "Linear";
    }
    M_modelUseJouleEffect = true;

    if ( M_solverName == "Linear" || M_solverNewtonInitialGuessUseLinearHeat )
    {
        M_heatModel->initAlgebraicFactory();
        M_heatModel->algebraicFactory()->setFunctionLinearAssembly( [this]( auto & data ) {
                                                                        return this->updateLinear_Heat( data );
                                                                    } );
        M_heatModel->algebraicFactory()->setFunctionResidualAssembly( [this]( auto & data ) {
                                                                        return this->updateResidual_Heat( data );
                                                                    } );
    }
    if ( M_solverName == "Linear" || M_solverNewtonInitialGuessUseLinearElectric )
    {
        M_electricModel->initAlgebraicFactory();
        M_electricModel->algebraicFactory()->setFunctionLinearAssembly( [this]( auto & data ) {
                                                                            return this->updateLinear_Electric( data );
                                                                        } );                                    
    }
#endif

#if 0
    M_rangeMeshElements = ( M_heatModel->thermalProperties()->isDefinedOnWholeMesh() && M_electricModel->electricProperties()->isDefinedOnWholeMesh() )?
        elements(this->mesh() ) :
        intersect( M_heatModel->rangeMeshElements(), M_electricModel->rangeMeshElements() );
#endif

    // post-process
    this->initPostProcess();

    // update constant parameters into
    this->updateParameterValues();

    // backend
    this->initAlgebraicBackend();

    // block vector solution
    auto const& blockVectorSolutionElectric = *M_electricModel->algebraicBlockVectorSolution();
    auto const& blockVectorSolutionMagnetic = *M_magneticModel->algebraicBlockVectorSolution();
    int nBlockElectric = blockVectorSolutionElectric.size();
    int nBlockMagnetic = blockVectorSolutionMagnetic.size();
    int nBlock = nBlockMagnetic + nBlockElectric;
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    int indexBlock=0;
    int numberOfBlockSpaceElectric = 0;
    for ( int k=0;k<nBlockElectric ;++k )
    {
        bvs->operator()(indexBlock+k) = blockVectorSolutionElectric(k);
        numberOfBlockSpaceElectric += blockVectorSolutionElectric(k)->map().numberOfDofIdToContainerId();
    }
    indexBlock += nBlockElectric;
    for ( int k=0;k<nBlockMagnetic ;++k )
        bvs->operator()(indexBlock+k) = blockVectorSolutionMagnetic(k);
    indexBlock += nBlockMagnetic;
    // init monolithic vector associated to the block vector
    bvs->buildVector( this->backend() );

    size_type currentStartBlockSpaceIndex = 0;
    this->setStartSubBlockSpaceIndex( "electric", currentStartBlockSpaceIndex );
    currentStartBlockSpaceIndex += numberOfBlockSpaceElectric;
    this->setStartSubBlockSpaceIndex( "magnetic", currentStartBlockSpaceIndex );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( M_solverName == "Newton" || M_solverName == "Picard" )
        {
            auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
            this->setAlgebraicFactory( algebraicFactory );
        }
    }

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Electromagnetic","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Electromagnetic","initPostProcess", "start");
    this->timerTool("Constructor").start();


    // need to not include export fields of material of subphysics
    std::set<std::string> ppExportsAllFieldsAvailableMagnetic = Feel::FeelModels::detail::set_difference( this->magneticModel()->postProcessExportsAllFieldsAvailable(),
                                                                                                      this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->magneticModel()->physicsAvailable() ) );
    std::set<std::string> ppExportsAllFieldsAvailableElectric = Feel::FeelModels::detail::set_difference( this->electricModel()->postProcessExportsAllFieldsAvailable(),
                                                                                                          this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->electricModel()->physicsAvailable() ) );
    std::set<std::string> ppExportsAllFieldsAvailable;
    for ( auto const& s : ppExportsAllFieldsAvailableMagnetic )
        ppExportsAllFieldsAvailable.insert( prefixvm( this->magneticModel()->keyword(), s) );
    for ( auto const& s : ppExportsAllFieldsAvailableElectric )
        ppExportsAllFieldsAvailable.insert( prefixvm( this->electricModel()->keyword(), s) );

    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        if ( this->doRestart() && this->restartPath().empty() )
        {
            if ( M_exporter->doExport() )
                M_exporter->restart(this->timeInitial());
        }
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("Electromagnetic","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    // p.put( "toolbox-magnetic", M_magneticModel->journalSectionName() );
    // p.put( "toolbox-electric", M_electricModel->journalSectionName() );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    // Numerical Solver
    nl::json subPt;
    subPt.emplace( "solver", M_solverName );
    p["Numerical Solver"] = subPt;

    // Exporter
#if 0
    if ( M_exporter )
    {
        subPt.clear();
        subPt.put( "type",M_exporter->type() );
        subPt.put( "freq save",M_exporter->freq() );
        pt::ptree subPt2;
        for ( std::string const& fieldName : this->postProcessExportsFields() )
            subPt2.push_back( std::make_pair("", pt::ptree( fieldName ) ) );
        subPt.put_child( "fields", subPt2 );
        p.put_child( "Exporter", subPt );
    }
#endif

    // Algebraic Solver
    if ( this->algebraicFactory() )
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );

    p["Toolbox Magnetic"] = M_magneticModel->journalSection().to_string();
    p["Toolbox Electric"] = M_electricModel->journalSection().to_string();
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    // Numerical Solver
    if ( jsonInfo.contains( "Numerical Solver" ) )
    {
        Feel::Table tabInfoNumSolver;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoNumSolver, jsonInfo.at("Numerical Solver"), tabInfoProp );
        tabInfo->add( "Numerical Solver",  TabulateInformations::New( tabInfoNumSolver, tabInfoProp ) );
    }

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    // generate sub toolboxes info
    if ( M_magneticModel && jsonInfo.contains( "Toolbox Magnetic" ) )
    {
        nl::json::json_pointer jsonPointerMagnetic( jsonInfo.at( "Toolbox Magnetic" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerMagnetic ) )
        {
            auto tabInfos_magnetic = M_magneticModel->tabulateInformations( JournalManager::journalData().at( jsonPointerMagnetic ), tabInfoProp );
            TabulateInformationsSections::cast( tabInfos_magnetic )->erase( "Materials Properties" );
            tabInfo->add( "Toolbox Magnetic", tabInfos_magnetic );
        }
    }
    if ( M_electricModel && jsonInfo.contains( "Toolbox Electric" ) )
    {
        nl::json::json_pointer jsonPointerElectric( jsonInfo.at( "Toolbox Electric" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerElectric ) )
        {
            auto tabInfos_electric = M_electricModel->tabulateInformations( JournalManager::journalData().at( jsonPointerElectric ), tabInfoProp );
            TabulateInformationsSections::cast( tabInfos_electric )->erase( "Materials Properties" );
            tabInfo->add( "Toolbox Electric", tabInfos_electric );
        }
    }

    return tabInfo;
}
ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::startTimeStep()
{
#if 0 // TODO
    this->heatModel()->startTimeStep();
    this->updateTime( this->heatModel()->time() );
    this->updateParameterValues();
#endif
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
#if 0 // TODO
    this->heatModel()->updateTimeStep();
    this->updateTime( this->heatModel()->time() );
    this->updateParameterValues();
#endif
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
#if 0 // TODO
    this->log("Electromagnetic","exportResults", "start");
    this->timerTool("PostProcessing").start();

    auto mfields = this->modelFields();
    auto symbolExpr = this->symbolsExpr( mfields );
    //std::cout << "holalla \n "<< symbolExpr.names() << std::endl;
    M_heatModel->exportResults( time, symbolExpr );
    M_electricModel->exportResults( time, symbolExpr );

    auto exprExport =  hana::concat( M_materialsProperties->exprPostProcessExports( this->mesh(),this->physicsAvailable(),symbolExpr ),
                                     hana::concat( M_heatModel->exprPostProcessExportsToolbox( symbolExpr,M_heatModel->keyword() ),
                                                   M_electricModel->exprPostProcessExportsToolbox( symbolExpr,M_electricModel->keyword() ) ) );
    this->executePostProcessExports( M_exporter, time, mfields, symbolExpr, exprExport );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Electromagnetic","exportResults", "finish");
#endif
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }
    M_electricModel->setParameterValues( paramValues );
    M_magneticModel->setParameterValues( paramValues );
}


ELECTROMAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTROMAGNETIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Electromagnetic","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    if ( M_solverName == "Linear" )
    {
        M_electricModel->solve();
        M_magneticModel->solve();
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
    }
    else if ( M_solverName == "Newton" || M_solverName == "Picard" )
    {
#if 0 // TODO
        // initial guess
        if ( M_solverNewtonInitialGuessUseLinearElectric )
            M_electricModel->solve();
        if ( M_solverNewtonInitialGuessUseLinearHeat )
            M_heatModel->solve();

        // solve non linear monolithic system
        M_electricModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("electric") );
        M_magenticModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("magnetic") );
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );
        this->algebraicBlockVectorSolution()->localize();
#endif
    }

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("Electromagnetic","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

} // end namespace Feel::FeelModels
