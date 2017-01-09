/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/advection/advection.hpp>

namespace Feel {
namespace FeelModels {

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
ADVECTION_CLASS_TEMPLATE_TYPE::Advection( 
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
: super_type( prefix, worldComm, subPrefix, rootRepository)
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
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
{
    return boost::make_shared<self_type>( prefix, worldComm, subPrefix, rootRepository );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    //-----------------------------------------------------------------------------//
    // Set model from options
    std::string advection_model = this->modelProperties().model();
    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        advection_model = soption(_name="model",_prefix=this->prefix());
    if( !advection_model.empty() )
        this->setModelName( advection_model );
    // Init super_type
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
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
ADVECTION_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();

    //this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Dirichlet" );
    this->M_bcDirichlet = detail::getBCFields<nDim, is_vectorial>( 
            this->modelProperties().boundaryConditions(), this->prefix(), "Dirichlet");
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", marker(d) );

    //this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Neumann" );
    this->M_bcNeumann = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().boundaryConditions(), this->prefix(), "Neumann"
            );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));

    //this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "advection", "Robin" );
    //for( auto const& d : this->M_bcRobin )
        //this->addMarkerRobinBC( marker(d) );

    //M_sources = this->modelProperties().boundaryConditions().template getScalarFields( "advection", "Sources" );
    M_sources = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().boundaryConditions(), this->prefix(), "Sources"
            );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::solve()
{
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();

    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    //M_bcRobin.setParameterValues( paramValues );
    M_sources.setParameterValues( paramValues );

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
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("Advection","updateBCStrongDirichletLinearPDE","start" );

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto const& u = this->fieldSolution();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                _element=u,_rhs=F,_expr=expression(d) );
    }

    this->log("Advection","updateBCStrongDirichletLinearPDE","finish" );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE(element_advection_ptrtype& fieldSource, bool BuildCstPart) const
{
    if( this->M_sources.empty() ) return;

    bool BuildSourceTerm = !BuildCstPart;

    if ( BuildSourceTerm )
    {
        //auto linearForm = form1( 
                //_test=this->functionSpace(), 
                //_vector=F,
                //_rowstart=this->rowStartInVector() 
                //);
        //auto const& v = this->fieldSolution();
        
        auto fieldSourceAux = this->functionSpace()->element();

        for( auto const& d : M_sources )
        {
            if ( marker(d).empty() )
            {
                //linearForm +=
                    //integrate( _range=elements(this->mesh()),
                               //_expr= inner( expression(d),id(v) ),
                               //_geomap=this->geomap() );
                fieldSourceAux.on( 
                        _range=elements(this->mesh()), 
                        _expr=expression(d), 
                        _geomap=this->geomap()
                        );
            }
            else
            {
                //linearForm +=
                    //integrate( _range=markedelements(this->mesh(),marker(d)),
                               //_expr= inner( expression(d),id(v) ),
                               //_geomap=this->geomap() );
                fieldSourceAux.on( 
                        _range=markedelements(this->mesh(),marker(d)), 
                        _expr=expression(d), 
                        _geomap=this->geomap()
                        );
            }
            
            *fieldSource += fieldSourceAux;
        }

    }
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTION_CLASS_TEMPLATE_TYPE::hasSourceTerm() const
{
    return !this->M_sources.empty(); 
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
