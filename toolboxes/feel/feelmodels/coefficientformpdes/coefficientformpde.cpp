/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

//#include <feel/feelmodels/modelmesh/createmesh.hpp>
//#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

namespace Feel
{
namespace FeelModels
{

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::CoefficientFormPDE( typename super_type::super2_type const& genericPDE,
                                                            std::string const& prefix,
                                                            std::string const& keyword,
                                                            worldcomm_ptr_t const& worldComm,
                                                            std::string const& subPrefix,
                                                            ModelBaseRepository const& modelRep )
    :
    super_type( genericPDE, prefix, keyword, worldComm, subPrefix, modelRep )
{}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("CoefficientFormPDE","init", "start" );
    this->timerTool("Constructor").start();

    if ( !this->M_mesh )
        this->initMesh();

    this->initMaterialProperties();

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    this->initInitialConditions();
#if 0
    // stabilization gls
    if ( M_stabilizationGLS )
    {
        typedef StabilizationGLSParameter<mesh_type, nOrderTemperature> stab_gls_parameter_impl_type;
        M_stabilizationGLSParameter.reset( new stab_gls_parameter_impl_type( this->mesh(),prefixvm(this->prefix(),"stabilization-gls.parameter") ) );
        M_stabilizationGLSParameter->init();
    }
#endif
    // post-process
    this->initPostProcess();

    // update constant parameters into
    this->updateParameterValues();

    // backend
    this->M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

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

    // functionspace
    if ( this->materialsProperties()->isDefinedOnWholeMesh( this->physic() ) )
    {
        this->M_rangeMeshElements = elements(this->mesh());
        M_Xh = space_unknown_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        this->M_rangeMeshElements = markedelements(this->mesh(), this->materialsProperties()->markers( this->physic() ));
        M_Xh = space_unknown_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=this->M_rangeMeshElements );
    }

    M_fieldUnknown.reset( new element_unknown_type(M_Xh,"temperature"));

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

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( { { this->physic(), std::string("Dirichlet") } } );
    for( auto const& d : this->M_bcDirichlet )
        M_bcDirichletMarkerManagement.addMarkerDirichletBC("elimination", name(d), markers(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( { { this->physic(),  std::string("Neumann") }  } );
    for( auto const& d : this->M_bcNeumann )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d),markers(d));

    std::string tmp = this->physic();
    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( std::move(tmp), "Robin" );
    for( auto const& d : this->M_bcRobin )
        M_bcRobinMarkerManagement.addMarkerRobinBC( name(d),markers(d) );

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    std::set<std::string> unknownMarkers;

    // strong Dirichlet bc on temperature from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) );
        unknownMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersTemperatureByEntities = detail::distributeMarkerListOnSubEntity( mesh, unknownMarkers );

    // on topological faces
    auto const& listMarkedFacesTemperature = std::get<0>( meshMarkersTemperatureByEntities );
    if ( !listMarkedFacesTemperature.empty() )
        this->updateDofEliminationIds( "temperature", Xh, markedfaces( mesh,listMarkedFacesTemperature ) );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initInitialConditions()
{
    // TODO
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    // TODO
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->initBasePostProcess();
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
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateInformationObject( pt::ptree & p )
{
    // TODO
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------------Info : Heat------------------||"
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
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << this->mesh()->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space Unknown Discretization"
           << "\n     -- name of unknown         : " << this->unknownName()
           << "\n     -- symbol of unknown : " << this->unknownSymbol()
           << "\n     -- basis : " << this->unknownBasis()
           << "\n     -- number of dof : " << this->spaceUnknown()->nDof() << " (" << this->spaceUnknown()->nLocalDof() << ")";
    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
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
        this->materialsProperties()->setParameterValues( paramValues );
    }
#if 0
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
#endif
    this->log("CoefficientFormPDE","setParameterValues", "finish");
}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::solve()
{
    // TODO
}


} // namespace Feel
} // namespace FeelModels
