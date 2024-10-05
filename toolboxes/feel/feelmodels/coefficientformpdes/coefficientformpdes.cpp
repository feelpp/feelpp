/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>

//#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/modelcore/utils.hpp>

#include <feel/feelvf/operators.hpp>

namespace Feel
{
namespace FeelModels
{

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::CoefficientFormPDEs( std::string const& prefix,
                                                              std::string const& keyword,
                                                              worldcomm_ptr_t const& worldComm,
                                                              std::string const& subPrefix,
                                                              ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    super_physics_type( "GenericPDEs" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep/*, ModelBaseCommandLineOptions( coefficientformpdes_options( prefix ) )*/ )
{
    M_solverName = soption(_prefix=this->prefix(),_name="solver",_vm=this->clovm());
}


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models )
{
    auto currentPhysic = std::dynamic_pointer_cast<ModelPhysicCoefficientFormPDEs<nDim>>( physicsTree.physic() );
    CHECK( currentPhysic ) << "wrong physic";

    using model_physic_cfpde_type = ModelPhysicCoefficientFormPDE<nDim>;
    for ( auto & eqData : currentPhysic->pdes() )
    {
        std::string const& eqName = std::get<0>( eqData );
        auto & eqInfos = std::get<1>( eqData );
        //std::string const& eqName = eqInfos->equationName();
        if ( models.hasType( eqName ) )
        {
            for ( auto const& [_nameFromModels,_model] : models.models( eqName ) )
            {
                if ( !eqInfos )
                    eqInfos = std::make_shared<typename model_physic_cfpde_type::infos_type>( eqName, _model.setup() );
                break;
            }
        }

        CHECK( eqInfos ) << "something is missing";
        std::string const& eqBasisTag = eqInfos->unknownBasis();
        std::shared_ptr<coefficient_form_pde_base_type> newCoefficientFormPDE;
        hana::for_each( tuple_type_unknown_basis, [this,&eqInfos,&eqName,&eqBasisTag,&newCoefficientFormPDE]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) == eqBasisTag )
                            {
                                using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                                newCoefficientFormPDE.reset( new coefficient_form_pde_type( eqInfos, prefixvm( this->prefix(),eqName ), eqName/*this->keyword()*/,
                                                                                            this->worldCommPtr(), this->subPrefix(), this->repository() ) );
                            }
                        });

        //std::get<1>( eq ) = newCoefficientFormPDE;
        M_coefficientFormPDEs.push_back( newCoefficientFormPDE );
        physicsTree.addChild( newCoefficientFormPDE, models );
    }

#if 0
    physicsTree.updateMaterialSupportFromChildren( "intersect" );
#endif
}


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("CoefficientFormPDEs","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();
    DCHECK( this->hasModelProperties() ) << "no model properties";

    this->initPhysics(
        this->shared_from_this(),
        [this]( typename super_physics_type::PhysicsTree & physicsTree ) {
            physicsTree.updatePhysics( this->shared_from_this(), this->modelProperties().models() );
#if 0
            CHECK( M_heatModel && M_electricModel ) << "aiai";
            physicsTree.updatePhysics( M_heatModel, this->modelProperties().models() );
            physicsTree.updatePhysics( M_electricModel, this->modelProperties().models() );
#endif
        } );


    this->initMaterialProperties();

    this->initMesh();

    for ( auto & cfpdeBase : M_coefficientFormPDEs )
    {
        hana::for_each( tuple_type_unknown_basis, [this,&cfpdeBase]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                return;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                            if ( !cfpde ) CHECK( false ) << "failure in dynamic_pointer_cast";
                            cfpde->setManageParameterValues( false );
                            if ( !cfpde->hasModelProperties() )
                            {
                                cfpde->setModelProperties( this->modelPropertiesPtr() );
                                cfpde->setManageParameterValuesOfModelProperties( false );
                            }
                            cfpde->setMaterialsProperties( M_materialsProperties );
                            cfpde->setModelMeshAsShared( this->modelMesh() );

                            // TODO check if the same space has already built
                            cfpde->init( false );
                        });
    }

    if ( !this->isStationary() )
    {
        CHECK( !M_coefficientFormPDEs.empty() ) << "no equation";
        // up initial time
        this->setTimeInitial( M_coefficientFormPDEs.front()->timeInitial() );
        // up current time
        this->updateTime( M_coefficientFormPDEs.front()->currentTime() );
    }

    // update constant parameters into
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // post-process
    this->initPostProcess();


    // backend
    this->initAlgebraicBackend();

    int nBlock = 0;
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
        nBlock += cfpdeBase->algebraicBlockVectorSolution()->size();
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    int indexBlock = 0, startBlockSpace = 0;
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
    {
        this->setStartSubBlockSpaceIndex( cfpdeBase->equationName(), startBlockSpace );
        auto const& blockVectorSolutionPDE = *(cfpdeBase->algebraicBlockVectorSolution());
        int nBlockPDE = blockVectorSolutionPDE.size();
        for ( int k=0;k<nBlockPDE ;++k )
        {
            bvs->operator()(indexBlock+k) = blockVectorSolutionPDE(k);
            startBlockSpace += blockVectorSolutionPDE(k)->map().numberOfDofIdToContainerId();
        }
        indexBlock += nBlockPDE;
    }
    // init petsc vector associated to the block
    bvs->buildVector( this->backend() );

    if ( M_solverName == "automatic" )
        this->updateAutomaticSolverSelection();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();


}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("CoefficientFormPDEs","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
         super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );
    CHECK( this->mesh() ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("CoefficientFormPDEs","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("CoefficientFormPDEs","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("CoefficientFormPDEs","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("CoefficientFormPDEs","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::set<std::string> ppExportsAllFieldsAvailable;
    for (auto const& cfpde : M_coefficientFormPDEs )
    {
        std::set<std::string> ppExportsAllFieldsAvailableInCFPDE = Feel::FeelModels::detail::set_difference( cfpde->postProcessExportsAllFieldsAvailable(),
                                                                                                             this->materialsProperties()->postProcessExportsAllFieldsAvailable( cfpde->mesh(), cfpde->physicsAvailable() ) );
        for ( auto const& s : ppExportsAllFieldsAvailableInCFPDE )
            ppExportsAllFieldsAvailable.insert( prefixvm( cfpde->keyword(), s) );
    }

    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() || this->hasPostProcessExportsExpr() )
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

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("CoefficientFormPDEs","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );


#if 0

    int nEq = M_coefficientFormPDEs.size();
    int nBlock = nEq;//this->nBlockMatrixGraph();
    BlocksBaseSparseMatrix<double> myblockMat( nBlock,nBlock );
    //BlocksBaseGraphCSR myblockGraph(nEq,nEq);
    for ( int k=0;k<nEq;++k )
    {
        auto const& cfpdeBase = M_coefficientFormPDEs[k];
        hana::for_each( tuple_type_unknown_basis, [this,&nEq,&k,&myblockMat,&cfpdeBase]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                return;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                            if ( !cfpde ) CHECK( false ) << "failure in dynamic_pointer_cast";

                            int rowId = this->startSubBlockSpaceIndex( cfpde->equationName() );

                            myblockMat(k,k) = this->backend()->newIdentityMatrix( cfpde->spaceUnknown()->dof(), cfpde->spaceUnknown()->dof() );
                            if ( k == 0 )
                            {
                                auto dofsBody = cfpde->spaceUnknown()->dofs( markedfaces(this->mesh(),"Interface") );
                                auto matFI_Id = myblockMat(k,k);
                                for ( auto dofid : dofsBody )
                                {
                                    matFI_Id->set( dofid,dofid, 0.);
                                }
                            }
                            std::cout << "TOTO k="<<k<< " : " << cfpdeBase->unknownBasis() << std::endl;
                            //auto Matp = this->backend()->newIdentityMatrix( algebraicFactory->matrix()->mapColPtr(), algebraicFactory->matrix()->mapRowPtr() );
// #if 0
//                             myblockGraph(rowId,rowId) = stencil(_test=cfpde->spaceUnknown(),
//                                                                 _trial=cfpde->spaceUnknown() )->graph();
// #else
//                             auto blockGraph = cfpde->buildBlockMatrixGraph();
//                             for (int bg1 = 0 ; bg1< blockGraph.nRow() ; ++bg1)
//                                 for (int bg2 = 0 ; bg2< blockGraph.nCol() ; ++bg2)
//                                     myblockGraph(rowId+bg1,rowId+bg2) = blockGraph(bg1,bg2);
// #endif
                            if ( k != 0 )
                                return;

                            // maybe coupling with other equation in the row
                            for ( int k2=0;k2<nEq;++k2 )
                            {
                                if ( k2 != 1 )
                                    continue;
                                // if ( k == k2 )
                                //     continue;
                                // if ( k != 0 )
                                //     continue;
                                auto const& cfpdeBase2 = M_coefficientFormPDEs[k2];
                                hana::for_each( tuple_type_unknown_basis, [this,&rowId,&myblockMat,&cfpde,&cfpdeBase2,&k,&k2]( auto const& e2 )
                                                {
                                                    if ( this->unknowBasisTag( e2 ) != cfpdeBase2->unknownBasis() )
                                                        return;

                                                    using coefficient_form_pde_2_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e2)>;
                                                    auto cfpde2 = std::dynamic_pointer_cast<coefficient_form_pde_2_type>( cfpdeBase2 );
                                                    if ( !cfpde2 ) CHECK( false ) << "failure in dynamic_pointer_cast";
                                                    std::cout << "TITI k2="<<k2<< " : " << cfpdeBase2->unknownBasis() << std::endl;


#if 1
                                                    sparse_matrix_ptrtype mat;
                                                    bool buildNewMatrix = true;
                                                    OperatorInterpolationMatrixSetup matSetup( mat, buildNewMatrix? Feel::DIFFERENT_NONZERO_PATTERN : Feel::SAME_NONZERO_PATTERN,
                                                                                               //startBlockIndexVelocity, startBlockIndexAngularVelocity );
                                                                                               0,0/*rowId, startBlockIndexAngularVelocity*/ );
                                                    auto w = cfpde2->spaceUnknown()->element();

                                                    if constexpr( std::decay_t<decltype(*cfpde2)>::unknown_is_scalar )
                                                    {
                                                        auto opI_AngularVelocity = opInterpolation( _domainSpace=cfpde2->spaceUnknown(),_imageSpace=cfpde->spaceUnknown(),
                                                                                                    _range= markedfaces(this->mesh(),"Bound_super"/*"Interface"*/),
                                                                                                    //_type= makeExprInterpolation( id(w)*vec(-Py()+massCenter(1,0),Px()-massCenter(0,0) ), nonconforming_t() ),
                                                                                                    //_type= makeExprInterpolation( id(w), nonconforming_t() ),
                                                                                                    //_type= makeExprInterpolation( -cross( N(), trans(grad(w)) ), nonconforming_t() ),
                                                                                                    //_type= makeExprInterpolation( -cross( trans(grad(w)),N() ), nonconforming_t() ),
                                                                                                    _type= makeExprInterpolation( -trans(grad(w)), nonconforming_t() ),

                                                                                                    _matrix=matSetup );

                                                        if ( buildNewMatrix )
                                                            mat = opI_AngularVelocity->matPtr();
                                                        myblockMat(k,k2) = opI_AngularVelocity->matPtr();
                                                    }
                                                    // M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
#endif

                                                    // auto tse = this->trialSymbolsExpr( mfields, cfpde2->trialSelectorModelFields() );
                                                    // auto trialSymbolNames = tse.names();

                                                    // if ( cfpde->hasSymbolDependencyInCoefficients( trialSymbolNames, se ) )
                                                    // {
                                                    //     int colId = this->startSubBlockSpaceIndex( cfpde2->equationName() );
                                                    //     myblockGraph(rowId,colId) = stencil(_test=cfpde->spaceUnknown(),
                                                    //                                  _trial=cfpde2->spaceUnknown() )->graph();
                                                    // }
                                                });
                            }

                        });
    }




    auto matP = this->backend()->newBlockMatrix(_block=myblockMat, _copy_values=true);

    //auto matP = this->backend()->newIdentityMatrix( algebraicFactory->matrix()->mapColPtr(), algebraicFactory->matrix()->mapRowPtr() );
    auto matQ = this->backend()->newIdentityMatrix( matP->mapColPtr(), matP->mapRowPtr() );
    // for ( auto & [bpname,bbc] : *this )
    // {
    //     auto dofsBody = fluidToolbox.functionSpaceVelocity()->dofs( bbc.rangeMarkedFacesOnFluid() );
    //     auto const& basisToContainerGpVelocity = matQ->mapCol().dofIdToContainerId( startBlockIndexVelocity );
    //     for ( auto dofid : dofsBody )
    //     {
    //         matQ->set( basisToContainerGpVelocity[dofid],basisToContainerGpVelocity[dofid], 0.);
    //     }
    // }
    // matQ->close();

    auto applyQ = [this,matQ](vector_ptrtype const& Ud, vector_ptrtype & Ui){
                      matQ->multVector( *Ud, *Ui );
                      //this->solverPtAP_applyQ(Ud,Ui);
                  };
    algebraicFactory->initSolverPtAP( matP, applyQ );

#endif

    bool hasTimeSteppingTheta = false;
    for (auto & cfpde : M_coefficientFormPDEs )
    {
        if ( cfpde->timeStepping() == "Theta" )
        {
            hasTimeSteppingTheta = true;
            break;
        }
    }

    if ( hasTimeSteppingTheta )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
        algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );

        bool hasStabilizationGLS = false;
        for (auto const& cfpde : M_coefficientFormPDEs )
        {
            if ( cfpde->applyStabilization() )
            {
                hasStabilizationGLS = true;
                break;
            }
        }
        if ( hasStabilizationGLS )
            algebraicFactory->dataInfos().addVectorInfo( "time-stepping.previous-solution", this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() ) );
    }

}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("CoefficientFormPDEs","buildBlockMatrixGraph", "start" );
    int nEq = M_coefficientFormPDEs.size();
    BlocksBaseGraphCSR myblockGraph(nEq,nEq);

    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    for ( int k=0;k<nEq;++k )
    {
        auto const& cfpdeBase = M_coefficientFormPDEs[k];
        hana::for_each( tuple_type_unknown_basis, [this,&nEq,&myblockGraph,&k,&cfpdeBase,&mfields,&se]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                return;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                            if ( !cfpde ) CHECK( false ) << "failure in dynamic_pointer_cast";

                            int rowId = this->startSubBlockSpaceIndex( cfpde->equationName() );
#if 0
                            myblockGraph(rowId,rowId) = stencil(_test=cfpde->spaceUnknown(),
                                                                _trial=cfpde->spaceUnknown() )->graph();
#else
                            auto blockGraph = cfpde->buildBlockMatrixGraph();
                            for (int bg1 = 0 ; bg1< blockGraph.nRow() ; ++bg1)
                                for (int bg2 = 0 ; bg2< blockGraph.nCol() ; ++bg2)
                                    myblockGraph(rowId+bg1,rowId+bg2) = blockGraph(bg1,bg2);
#endif

                            // maybe coupling with other equation in the row
                            for ( int k2=0;k2<nEq;++k2 )
                            {
                                if ( k == k2 )
                                    continue;
                                auto const& cfpdeBase2 = M_coefficientFormPDEs[k2];
                                hana::for_each( tuple_type_unknown_basis, [this,&myblockGraph,&rowId,&cfpde,&cfpdeBase2,&mfields,&se]( auto const& e2 )
                                                {
                                                    if ( this->unknowBasisTag( e2 ) != cfpdeBase2->unknownBasis() )
                                                        return;

                                                    using coefficient_form_pde_2_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e2)>;
                                                    auto cfpde2 = std::dynamic_pointer_cast<coefficient_form_pde_2_type>( cfpdeBase2 );
                                                    if ( !cfpde2 ) CHECK( false ) << "failure in dynamic_pointer_cast";

                                                    auto tse = this->trialSymbolsExpr( mfields, cfpde2->trialSelectorModelFields() );
                                                    auto trialSymbolNames = tse.names();

                                                    if ( cfpde->hasSymbolDependencyInCoefficients( trialSymbolNames, se ) )
                                                    {
                                                        int colId = this->startSubBlockSpaceIndex( cfpde2->equationName() );
                                                        myblockGraph(rowId,colId) = stencil(_test=cfpde->spaceUnknown(),
                                                                                     _trial=cfpde2->spaceUnknown() )->graph();
                                                    }
                                                });
                            }

                        });
    }
    myblockGraph.close();
    this->log("CoefficientFormPDEs","buildBlockMatrixGraph", "finish" );
    return myblockGraph;

}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateAutomaticSolverSelection()
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( 0/*rowStartInVector*/ ) );
    auto trialSymbolNames = tse.names();
    int nEq = M_coefficientFormPDEs.size();
    bool isLinear = true;
    for ( int k=0;k<nEq;++k )
    {
        auto const& cfpdeBase = M_coefficientFormPDEs[k];

        hana::for_each( tuple_type_unknown_basis, [this,&cfpdeBase,&isLinear,&se,&trialSymbolNames]( auto const& e )
        {
            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                return;

            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
            if ( !cfpde ) CHECK( false ) << "failure in dynamic_pointer_cast";

            if ( cfpde->hasSymbolDependencyInCoefficients( trialSymbolNames, se ) ||
                 cfpde->hasSymbolDependencyInBoundaryConditions( trialSymbolNames, se ) ||
                 ( cfpde->applyStabilization() && cfpde->stabilizationGLS_applyShockCapturing() )
                 )
            {
                isLinear = false;
                return;
            }
        });
        if ( !isLinear )
            break;
    }

    M_solverName = ( isLinear )? "Linear" : "Newton";
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
std::string const&
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::unknowBasisTag( variant_unknown_basis_type const& vb )
{
    return S_unknownBasisTags[vb.index()];
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
const std::vector<std::string> COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::S_unknownBasisTags = { COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TAG };


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    this->modelFields().updateInformationObject( p["Fields"] );

    if ( this->algebraicFactory() )
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );

    if ( this->hasModelProperties() )
        this->modelProperties().parameters().updateInformationObject( p["Parameters"] );

    nl::json subPt;
    for (auto & cfpde : M_coefficientFormPDEs )
        subPt[cfpde->keyword()] = cfpde->journalSection().to_string();
    p.emplace( "Coefficient Form PDE", subPt );

    //this->symbolsExpr().updateInformationObject( p["Symbols Expression"] );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp.newByIncreasingVerboseLevel() ) );

    if ( jsonInfo.contains("Parameters") )
        tabInfo->add( "Parameters", TabulateInformationTools::FromJSON::tabulateInformationsSymbolsExpr( jsonInfo.at("Parameters"), tabInfoProp, true ) );

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    //if ( jsonInfo.contains( "Symbols Expression" ) )
    //tabInfo->add( "Symbols Expression", TabulateInformationTools::FromJSON::tabulateInformationsSymbolsExpr( jsonInfo.at("Symbols Expression"), tabInfoProp/*.newByIncreasingVerboseLevel()*/ ) );

    if ( jsonInfo.contains("Coefficient Form PDE") )
    {
        auto const& jsonInfo_cfpde = jsonInfo.at("Coefficient Form PDE");
        for (auto & cfpde : M_coefficientFormPDEs )
        {
            if ( jsonInfo_cfpde.contains(cfpde->keyword()) )
            {
                nl::json::json_pointer jsonPointerCFPDE( jsonInfo_cfpde.at( cfpde->keyword() ).template get<std::string>() );
                if ( JournalManager::journalData().contains( jsonPointerCFPDE ) )
                {
                    auto tabInfos_cfpde = cfpde->tabulateInformations(  JournalManager::journalData().at( jsonPointerCFPDE ), tabInfoProp );
                    tabInfo->add( (boost::format("Toolbox Coefficient Form PDE : %1%")%cfpde->keyword()).str(), tabInfos_cfpde );
                }
            }
        }
    }
    return tabInfo;
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<TSBase>
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::timeStepBase() const
{
    CHECK( !M_coefficientFormPDEs.empty() ) << "no equation";
    return M_coefficientFormPDEs.front()->timeStepBase();
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    for (auto & cfpde : M_coefficientFormPDEs )
        cfpde->startTimeStep();

    // up current time
    this->updateTime( this->timeStepBase()->time() );

    this->updateParameterValues();
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    for (auto & cfpde : M_coefficientFormPDEs )
        cfpde->updateTimeStep();

    // up current time
    this->updateTime( this->timeStepBase()->time() );

    this->updateParameterValues();
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateTimeStepCurrentResidual()
{
    if ( this->isStationary() )
        return;
    auto algebraicFactory = this->algebraicFactory();
    if ( !algebraicFactory )
        return;

    bool hasTimeSteppingTheta = false;
    for (auto const& cfpde : M_coefficientFormPDEs )
    {
        if ( cfpde->timeStepping() == "Theta" )
        {
            hasTimeSteppingTheta = true;
            break;
        }
    }

    if ( hasTimeSteppingTheta )
    {
        M_timeStepThetaSchemePreviousContrib->zero();
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        ModelAlgebraic::DataUpdateResidual dataResidual( this->algebraicBlockVectorSolution()->vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
        dataResidual.addInfo( "time-stepping.evaluate-residual-without-time-derivative" );
        algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        algebraicFactory->evaluateResidual( dataResidual );
        algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        bool hasStabilizationGLS = false;
        for (auto const& cfpde : M_coefficientFormPDEs )
        {
            if ( cfpde->applyStabilization() )
            {
                hasStabilizationGLS = true;
                break;
            }
        }
        if ( hasStabilizationGLS )
        {
            auto & dataInfos = algebraicFactory->dataInfos();
            *dataInfos.vectorInfo( "time-stepping.previous-solution" ) = *this->algebraicBlockVectorSolution()->vectorMonolithic();
            dataInfos.addParameterValuesInfo( "time-stepping.previous-parameter-values", M_currentParameterValues );
        }
    }
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->materialsProperties()->exprPostProcessExports( this->mesh(),this->physicsAvailable(),se ) );
}


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    for ( auto [physicName,physicData] : this->physics/*FromCurrentType*/() )
        physicData->updateParameterValues( paramValues );


    this->updateParameterValues_postProcess( paramValues, prefixvm("postprocess",this->keyword(),"_" ) );
    for (auto const& cfpdeBase : M_coefficientFormPDEs )
        cfpdeBase->updateParameterValues_postProcess( paramValues, prefixvm("postprocess",cfpdeBase->keyword(),"_" ) );

    this->setParameterValues( paramValues );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto const& [param,val] : paramValues )
        M_currentParameterValues[param] = val;

    if ( this->manageParameterValuesOfModelProperties() )
    {
        //std::cout << "JJJ paramValues : " << paramValues << std::endl;
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for (auto & cfpdeBase : M_coefficientFormPDEs )
        cfpdeBase->setParameterValues( paramValues );
}


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("CoefficientFormPDEs","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();

    for (auto & cfpdeBase : M_coefficientFormPDEs )
        cfpdeBase->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex( cfpdeBase->equationName() ) );

    this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );

    this->algebraicBlockVectorSolution()->localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("CoefficientFormPDEs","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
bool
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::checkResults() const
{
    auto se = this->symbolsExpr();
    // several calls (not do in on line) to be sure that all check have been run
    bool checkValue = super_type::checkResults( se );
    std::vector<bool> checkCFPDE;
    checkCFPDE.reserve(M_coefficientFormPDEs.size());
    for (auto & cfpdeBase : M_coefficientFormPDEs )
        checkCFPDE.push_back( cfpdeBase->checkResults( se ) );
    checkValue = checkValue && (std::find(std::begin(checkCFPDE), std::end(checkCFPDE), false) == std::end(checkCFPDE));
    return checkValue;
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol, this->rowStartInVector() );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateLinearPDE_spec<FilterBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateLinearPDEDofElimination_spec<FilterBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    vector_ptrtype& U = data.initialGuess();
    auto mctx = this->modelContext( U, this->rowStartInVector() );
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateNewtonInitialGuess_spec<FilterBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateJacobian_spec<FilterBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
        cfpdeBase->updateJacobianDofElimination( data );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol, this->rowStartInVector() );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateResidual_spec<FilterBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
        cfpdeBase->updateResidualDofElimination( data );
}

} // namespace Feel
} // namespace FeelModels
