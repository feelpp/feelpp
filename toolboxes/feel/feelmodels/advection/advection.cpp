/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/advection/advection.hpp>

namespace Feel {
namespace FeelModels {

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
ADVECTION_CLASS_TEMPLATE_TYPE::Advection( 
        std::string const& prefix,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
: super_type( prefix, worldComm, subPrefix, modelRep )
{
    this->log("Advection","constructor", "start" );

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"Advection.info") );
    //-----------------------------------------------------------------------------//
    // load info from .bc file
    this->loadConfigBCFile();
    // get periodicity from options (if needed)
    this->loadPeriodicityFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh, space, exporter,...
    //if ( buildMesh ) this->build();
    //-----------------------------------------------------------------------------//
    this->log("Advection","constructor", "finish");
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
typename ADVECTION_CLASS_TEMPLATE_TYPE::self_ptrtype 
ADVECTION_CLASS_TEMPLATE_TYPE::New( 
        std::string const& prefix,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type>( prefix, worldComm, subPrefix, modelRep );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    // Init super_type
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
    // Set initial value
    this->setInitialValue();
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::setInitialValue()
{
    this->modelProperties().parameters().updateParameterValues();
    if( !this->M_icValue.empty() )
    {
        auto const& phi = this->fieldSolutionPtr();
        this->M_icValue.setParameterValues( this->modelProperties().parameters().toParameterValues() );
        for( auto const& iv : this->M_icValue )
        {
            if( markers(iv).empty() )
            {
                phi->on( _range=elements(this->mesh()),
                         _expr=expression(iv),
                         _geomap=this->geomap() );
            }
            else
            {
                phi->on( _range=markedelements(this->mesh(), markers(iv)),
                         _expr=expression(iv),
                         _geomap=this->geomap() );
            }
        }
    }
}

namespace detail {

template<int Dim, bool isVectorial, 
        typename = typename std::enable_if<!isVectorial>::type
        >
map_scalar_field<2> getBCFields( 
        BoundaryConditions const& bc, 
        std::string const& field, std::string const& type
        )
{
    return bc.getScalarFields( std::string(field), std::string(type) );
}

template<int Dim, bool isVectorial, 
        typename = typename std::enable_if<isVectorial>::type
        >
map_vector_field<Dim, 1, 2> getBCFields( 
        BoundaryConditions const& bc, 
        std::string const& field, std::string const& type
        )
{
    return bc.getVectorFields<Dim>( std::string(field), std::string(type) );
}

} // namespace detail

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::loadConfigICFile()
{
    this->M_icValue = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().initialConditions(), this->prefix(), "InitialValue" );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->M_bcInflowMarkers.clear();

    //this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Dirichlet" );
    this->M_bcDirichlet = detail::getBCFields<nDim, is_vectorial>( 
            this->modelProperties().boundaryConditions(), this->prefix(), "Dirichlet");
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", name(d), markers(d) );

    //this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Neumann" );
    this->M_bcNeumann = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().boundaryConditions(), this->prefix(), "Neumann"
            );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,name(d),markers(d));

    //this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "advection", "Robin" );
    //for( auto const& d : this->M_bcRobin )
        //this->addMarkerRobinBC( marker(d) );

    for( std::string const& bcMarker: this->modelProperties().boundaryConditions().markers( this->prefix(), "inflow" ) )
        if( std::find(this->M_bcInflowMarkers.begin(), this->M_bcInflowMarkers.end(), bcMarker) == this->M_bcInflowMarkers.end() )
            this->M_bcInflowMarkers.push_back( bcMarker );

    //M_sources = this->modelProperties().boundaryConditions().template getScalarFields( "advection", "Sources" );
    this->M_sources = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().boundaryConditions(), this->prefix(), "Sources"
            );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::solve()
{
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();

    this->M_bcDirichlet.setParameterValues( paramValues );
    this->M_bcNeumann.setParameterValues( paramValues );
    //M_bcRobin.setParameterValues( paramValues );
    this->M_sources.setParameterValues( paramValues );

    super_type::solve(); 
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{
// TODO
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    if ( this->M_bcDirichlet.empty() && this->M_bcInflowMarkers.empty() ) return;

    this->log("Advection","updateBCStrongDirichletLinearPDE","start" );

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto const& u = this->fieldSolution();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    // Dirichlet bc
    for( auto const& d : this->M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=expression(d) );
    }

    // Inflow bc
    if( !this->isStationary() )
    {
        for( auto const& bcMarker: this->M_bcInflowMarkers )
        {
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(mesh, bcMarker),
                        _element=u,
                        _rhs=F,
                        _expr=(
                            // Transient part
                            idv(this->timeStepBDF()->polyDeriv())
                            // Advection part
                            - gradv(u)*idv(this->fieldAdvectionVelocity())
                            // Diffusion part
                            + (this->hasDiffusion()) * idv(this->diffusionReactionModel()->fieldDiffusionCoeff())*laplacianv(u)
                            // Reaction part
                            - (this->hasReaction()) * idv(this->diffusionReactionModel()->fieldReactionCoeff())*idv(u)
                            // Source part
                            + (this->hasSourceAdded() || this->hasSourceTerm()) * idv(*this->M_fieldSource)
                            )/(this->timeStepBDF()->polyDerivCoefficient(0))
                  );
        }
    }

    this->log("Advection","updateBCStrongDirichletLinearPDE","finish" );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( ModelAlgebraic::DataUpdateLinear & data ) const
{
    if( this->M_sources.empty() ) return;

    bool BuildCstPart = data.buildCstPart();
    if ( BuildCstPart )
        return;
    vector_ptrtype& F = data.rhs();

    auto mesh = this->mesh();
    auto space = this->functionSpace();
    auto const& v = this->fieldSolution();
    auto linearForm = form1( _test=space, _vector=F );
    for( auto const& d : this->M_sources )
    {
        auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
        linearForm += integrate( _range=rangeUsed,
                                 _expr= inner(expression(d),id(v)),
                                 _geomap=this->geomap() );

    }
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTION_CLASS_TEMPLATE_TYPE::hasSourceTerm() const
{
    return !this->M_sources.empty(); 
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::addMarkerInflowBC( std::string const& markerName )
{
    if( std::find(this->M_bcInflowMarkers.begin(), this->M_bcInflowMarkers.end(), markerName) == this->M_bcInflowMarkers.end() )
        this->M_bcInflowMarkers.push_back( markerName );
}

namespace detail {

template<typename AdvT> 
void
loadPeriodicityFromOptionsVm(AdvT &, mpl::false_)
{}

template<typename AdvT>
void
loadPeriodicityFromOptionsVm(AdvT & adv, mpl::true_)
{
    node_type translat( AdvT::nDim );
    translat[0] = doption(_name="periodicity.translate-x",_prefix=adv.prefix());
    if ( AdvT::nDim >=2 )
        translat[1] = doption(_name="periodicity.translate-y",_prefix=adv.prefix());
    if ( AdvT::nDim == 3 )
        translat[2]= doption(_name="periodicity.translate-z",_prefix=adv.prefix());
    std::string marker1 = soption(_name="periodicity.marker1",_prefix=adv.prefix());
    std::string marker2 = soption(_name="periodicity.marker2",_prefix=adv.prefix());
    auto theperiodicity = Periodic<>( adv.mesh()->markerName(marker1),adv.mesh()->markerName(marker2), translat );

    adv.setPeriodicity( theperiodicity );
}

}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::loadPeriodicityFromOptionsVm()
{
    detail::loadPeriodicityFromOptionsVm( *this, mpl::bool_<PeriodicityType::is_periodic>() );
}

} // namespace FeelModels
} // namespace Feel
