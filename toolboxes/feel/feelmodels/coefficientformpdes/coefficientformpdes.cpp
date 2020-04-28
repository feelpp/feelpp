/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/modelcore/utils.hpp>

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
    super_type( prefix, keyword, worldComm, subPrefix, modelRep )
{
    M_solverName = soption(_prefix=this->prefix(),_name="solver");
}


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("CoefficientFormPDEs","init", "start" );
    this->timerTool("Constructor").start();

    if ( this->physics().empty() )
        this->setupGenericPDEs( this->keyword(), this->modelProperties().models().model( this->keyword() ).ptree() );

    if ( !this->M_mesh )
        this->initMesh();

    CHECK( this->hasModelProperties() ) << "no model properties";

    this->initMaterialProperties();

    for ( auto const& eq : this->pdes() )
    {
        std::string const eqBasisTag = eq.unknownBasis();
        std::shared_ptr<coefficient_form_pde_base_type> newCoefficientFormPDE;
        hana::for_each( tuple_type_unknown_basis, [this,&eq,&eqBasisTag,&newCoefficientFormPDE]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) == eqBasisTag )
                            {
                                using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                                std::shared_ptr<coefficient_form_pde_type> _newCoefficientFormPDE( new coefficient_form_pde_type( eq, prefixvm( this->prefix(),eq.physicDefault())/*this->prefix()*/, eq.physicDefault()/*this->keyword()*/,
                                                                                                                                  this->worldCommPtr(), this->subPrefix(), this->repository() ) );
                                _newCoefficientFormPDE->setManageParameterValues( false );
                                if ( !_newCoefficientFormPDE->hasModelProperties() )
                                {
                                    _newCoefficientFormPDE->setModelProperties( this->modelPropertiesPtr() );
                                    _newCoefficientFormPDE->setManageParameterValuesOfModelProperties( false );
                                }
                                _newCoefficientFormPDE->setMaterialsProperties( M_materialsProperties );
                                _newCoefficientFormPDE->setMesh( this->mesh() );

                                // TODO check if the same space has already built
                                _newCoefficientFormPDE->init( false );
                                newCoefficientFormPDE = _newCoefficientFormPDE;
                            }
                        });
        M_coefficientFormPDEs.push_back( newCoefficientFormPDE );
    }

    // post-process
    this->initPostProcess();

    // update constant parameters into
    this->updateParameterValues();

    // backend
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

    int nBlock = 0;
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
        nBlock += cfpdeBase->blockVectorSolution().size();
    M_blockVectorSolution.resize( nBlock );
    int indexBlock = 0, startBlockSpace = 0;
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
    {
        this->setStartSubBlockSpaceIndex( cfpdeBase->physicDefault(), startBlockSpace );
        auto const& blockVectorSolutionPDE = cfpdeBase->blockVectorSolution();
        int nBlockPDE = blockVectorSolutionPDE.size();
        for ( int k=0;k<nBlockPDE ;++k )
        {
            M_blockVectorSolution(indexBlock+k) = blockVectorSolutionPDE(k);
            startBlockSpace += blockVectorSolutionPDE(k)->map().numberOfDofIdToContainerId();
        }
        indexBlock += nBlockPDE;
    }

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

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

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
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->prefix(), this->repository().expr() ) );
        M_materialsProperties->updateForUse( M_mesh, this->modelProperties().materials(), *this );
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
                                                                                                             this->materialsProperties()->postProcessExportsAllFieldsAvailable( cfpde->physicsAvailable() ) );
        for ( auto const& s : ppExportsAllFieldsAvailableInCFPDE )
            ppExportsAllFieldsAvailable.insert( prefixvm( cfpde->keyword(), s) );
    }

    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->physicsAvailable() ) );
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
    this->log("CoefficientFormPDEs","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );

    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("CoefficientFormPDEs","buildBlockMatrixGraph", "start" );
    int nEq = M_coefficientFormPDEs.size();
    BlocksBaseGraphCSR myblockGraph(nEq,nEq);
    for ( int k=0;k<nEq;++k )
    {
        auto const& cfpdeBase = M_coefficientFormPDEs[k];
        hana::for_each( tuple_type_unknown_basis, [this,&nEq,&myblockGraph,&k,&cfpdeBase]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                return;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                            if ( !cfpde ) CHECK( false ) << "failure in dynamic_pointer_cast";

                            int rowId = this->startSubBlockSpaceIndex( cfpde->physicDefault() );
                            myblockGraph(rowId,rowId) = stencil(_test=cfpde->spaceUnknown(),
                                                                _trial=cfpde->spaceUnknown() )->graph();

                            // maybe coupling with other equation in the row
                            for ( int k2=0;k2<nEq;++k2 )
                            {
                                if ( k == k2 )
                                    continue;
                                auto const& cfpdeBase2 = M_coefficientFormPDEs[k2];
                                hana::for_each( tuple_type_unknown_basis, [this,&myblockGraph,&rowId,&cfpde,&cfpdeBase2]( auto const& e2 )
                                                {
                                                    if ( this->unknowBasisTag( e2 ) != cfpdeBase2->unknownBasis() )
                                                        return;

                                                    using coefficient_form_pde_2_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e2)>;
                                                    auto cfpde2 = std::dynamic_pointer_cast<coefficient_form_pde_2_type>( cfpdeBase2 );
                                                    if ( !cfpde2 ) CHECK( false ) << "failure in dynamic_pointer_cast";

                                                    auto tse = this->trialSymbolsExpr( this->modelFields(), cfpde2->trialSelectorModelFields() );
                                                    auto trialSymbolNames = tse.names();

                                                    bool coeffDependOnUnknown = false;
                                                    for ( std::string const& matName : cfpde->materialsProperties()->physicToMaterials( cfpde->physicDefault() ) )
                                                    {
                                                        if ( ( cfpde->materialsProperties()->hasProperty( matName, cfpde->convectionCoefficientName() ) &&
                                                               cfpde->materialsProperties()->materialProperty( matName, cfpde->convectionCoefficientName() ).hasSymbolDependency( trialSymbolNames ) ) ||
                                                             ( cfpde->materialsProperties()->hasProperty( matName, cfpde->diffusionCoefficientName() ) &&
                                                               cfpde->materialsProperties()->materialProperty( matName, cfpde->diffusionCoefficientName() ).hasSymbolDependency( trialSymbolNames ) ) ||
                                                             ( cfpde->materialsProperties()->hasProperty( matName, cfpde->reactionCoefficientName() ) &&
                                                               cfpde->materialsProperties()->materialProperty( matName, cfpde->reactionCoefficientName() ).hasSymbolDependency( trialSymbolNames ) ) ||
                                                             ( cfpde->materialsProperties()->hasProperty( matName, cfpde->sourceCoefficientName() ) &&
                                                               cfpde->materialsProperties()->materialProperty( matName, cfpde->sourceCoefficientName() ).hasSymbolDependency( trialSymbolNames ) ) ||
                                                             ( cfpde->materialsProperties()->hasProperty( matName, cfpde->firstTimeDerivativeCoefficientName() ) &&
                                                               cfpde->materialsProperties()->materialProperty( matName, cfpde->firstTimeDerivativeCoefficientName() ).hasSymbolDependency( trialSymbolNames ) )
                                                             )
                                                        {
                                                            coeffDependOnUnknown = true;
                                                            break;
                                                        }
                                                    }

                                                    if ( coeffDependOnUnknown )
                                                    {
                                                        int colId = this->startSubBlockSpaceIndex( cfpde2->physicDefault() );
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
std::string const&
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::unknowBasisTag( variant_unknown_basis_type const& vb )
{
    return S_unknownBasisTags[vb.index()];
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
const std::vector<std::string> COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::S_unknownBasisTags = { COEFFICIENTFORMPDES_LIST_UNKNOWN_BASIS_TAG };


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    for (auto const& cfpde  : M_coefficientFormPDEs )
        *_ostr << cfpde->getInfo()->str();
    return _ostr;
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateInformationObject( pt::ptree & p )
{
    // TODO
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->materialsProperties()->exprPostProcessExports( this->physicsAvailable(),se ) );
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

    this->setParameterValues( paramValues );
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    if ( this->manageParameterValuesOfModelProperties() )
    {
        //std::cout << "JJJ paramValues : " << paramValues << std::endl;
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
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

    M_blockVectorSolution.updateVectorFromSubVectors();

    for (auto & cfpdeBase : M_coefficientFormPDEs )
        cfpdeBase->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex( cfpdeBase->physicDefault() ) );

    // if ( this->materialsProperties()->hasThermalConductivityDependingOnSymbol( "heat_T" ) )
    //     M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );
    // else
    if ( M_solverName == "automatic" ) // TODO define automatic solver from pde/option
    {
        M_algebraicFactory->solve( "Linear", M_blockVectorSolution.vectorMonolithic() );
        // M_blockVectorSolution.localize();
        // M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );
    }
    else
        M_algebraicFactory->solve( M_solverName, M_blockVectorSolution.vectorMonolithic() );

    M_blockVectorSolution.localize();

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
    // several calls (not do in on line) to be sure that all check have been run
    bool checkValue = super_type::checkResults();
    std::vector<bool> checkCFPDE;
    checkCFPDE.reserve(M_coefficientFormPDEs.size());
    for (auto & cfpdeBase : M_coefficientFormPDEs )
        checkCFPDE.push_back( cfpdeBase->checkResults() );
    checkValue = checkValue && (std::find(std::begin(checkCFPDE), std::end(checkCFPDE), false) == std::end(checkCFPDE));
    return checkValue;
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    auto mctx = this->modelContext();
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
    auto mctx = this->modelContext();
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
