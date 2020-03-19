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
{}


COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("CoefficientFormPDEs","init", "start" );
    this->timerTool("Constructor").start();

    if ( this->physic().empty() )
        this->setupGenericPDEs( this->keyword(), this->modelProperties().models().model( this->keyword() ).ptree() );

    if ( !this->M_mesh )
        this->initMesh();

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
                                std::shared_ptr<coefficient_form_pde_type> _newCoefficientFormPDE( new coefficient_form_pde_type( eq, this->prefix(), eq.physic()/*this->keyword()*/, this->worldCommPtr(), this->subPrefix(), this->repository() ) );
                                //newCoefficientFormPDE.reset( new typename self_type::traits::template coefficient_form_pde_t<decltype(e)>( this->prefix(), this->keyword(), this->worldCommPtr(), this->subPrefix(), this->repository() ) );
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

    // TODO
    // update constant parameters into
    //this->updateParameterValues();

    // backend
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

    int nBlock = 0;
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
        nBlock += cfpdeBase->blockVectorSolution().size();
    M_blockVectorSolution.resize( nBlock );
    int indexBlock = 0, startBlockSpace = 0;
    for ( auto const& cfpdeBase : M_coefficientFormPDEs )
    {
        this->setStartSubBlockSpaceIndex( cfpdeBase->physic(), startBlockSpace );
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
                                                                                                             this->materialsProperties()->postProcessExportsAllFieldsAvailable( cfpde->physics() ) );
        for ( auto const& s : ppExportsAllFieldsAvailableInCFPDE )
            ppExportsAllFieldsAvailable.insert( prefixvm( cfpde->keyword(), s) );
    }

    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->physics() ) );
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
    int nBlock = M_coefficientFormPDEs.size();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    for (int k=0;k<nBlock;++k )
    {
        auto const& cfpdeBase = M_coefficientFormPDEs[k];
        hana::for_each( tuple_type_unknown_basis, [this,&myblockGraph,&k,&cfpdeBase]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) == cfpdeBase->unknownBasis() )
                            {
                                using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                                auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                                if ( !cfpde ) CHECK( false ) << "ai";
                                myblockGraph(k,k) = stencil(_test=cfpde->spaceUnknown(),
                                                            _trial=cfpde->spaceUnknown() )->graph();
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
    auto se = this->symbolsExpr();
    this->exportResults( time, se, this->materialsProperties()->exprPostProcessExports( this->physics(),se ) );
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
        cfpdeBase->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex( cfpdeBase->physic() ) );

    // if ( this->materialsProperties()->hasThermalConductivityDependingOnSymbol( "heat_T" ) )
    //     M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );
    // else
    M_algebraicFactory->solve( "LinearSystem", M_blockVectorSolution.vectorMonolithic() );

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

} // namespace Feel
} // namespace FeelModels
