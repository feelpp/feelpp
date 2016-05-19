
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
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", marker(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));

    //this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "advection", "Robin" );
    //for( auto const& d : this->M_bcRobin )
        //this->addMarkerRobinBC( marker(d) );

    M_sources = this->modelProperties().boundaryConditions().template getScalarFields( "advection", "Sources" );
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{

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

} // namespace FeelModels
} // namespace Feel
