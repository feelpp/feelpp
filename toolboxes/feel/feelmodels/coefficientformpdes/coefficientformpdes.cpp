/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

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
        hana::for_each( tuple_type_unknown_basis, [this,&eqBasisTag,&newCoefficientFormPDE]( auto const& e )
                        {
                            if ( this->unknowBasisTag( e ) == eqBasisTag )
                            {
                                newCoefficientFormPDE.reset( new typename self_type::traits::template coefficient_form_pde_t<decltype(e)>( this->prefix(), this->keyword(), this->worldCommPtr(), this->subPrefix(), this->repository() ) );
                                //std::cout << "unknow basis name " << this->unknowBasisTag( e ) << std::endl;
                            }
                        });
    }


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
    int nBlock = 1;//this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
#if 0
    myblockGraph(0,0) = stencil(_test=this->spaceTemperature(),
                                _trial=this->spaceTemperature() )->graph();
#endif
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

} // namespace Feel
} // namespace FeelModels
