/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/heat/heat.hpp>

#include <feel/feelvf/vf.hpp>

//#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>

namespace Feel
{
namespace FeelModels
{


HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data, this->symbolsExpr() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Heat","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    auto u = this->spaceTemperature()->element( U, this->rowStartInVector() );

    for( auto const& d : M_bcDirichlet )
    {
        auto theExpr = expression(d,this->symbolsExpr());
        u.on(_range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=theExpr );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateNewtonInitialGuess","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto const t = this->spaceTemperature()->element(XVec, this->rowStartInVector());
    this->updateJacobian( data, this->symbolsExpr(t) );
}
HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto const t = this->spaceTemperature()->element(XVec, this->rowStartInVector());
    this->updateResidual( data, this->symbolsExpr(t) );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateResidualDofElimination","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualWeakBC( DataUpdateResidual & data, element_temperature_external_storage_type const& u ) const
{
    if ( this->M_bcRobin.empty() && this->M_bcNeumann.empty() ) return;

    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& v = this->fieldTemperature();

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() );

    for( auto const& d : this->M_bcNeumann )
    {
        if ( buildCstPart )
        {
            auto theExpr = expression(d,this->symbolsExpr());
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

    for( auto const& d : this->M_bcRobin )
    {
        auto theExpr1 = expression1( d,this->symbolsExpr() );
        if ( buildNonCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idv(u)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            auto theExpr2 = expression2( d,this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= -timeSteppingScaling*theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualSourceTerm( DataUpdateResidual & data  ) const
{
    if ( this->M_volumicForcesProperties.empty() ) return;

    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    if ( buildCstPart )
    {
        auto myLinearForm = form1( _test=this->spaceTemperature(), _vector=R,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldTemperature();

        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto theExpr = expression(d,this->symbolsExpr());
            auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= -timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateJacobianDofElimination","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobianRobinBC(  DataUpdateJacobian & data ) const
{
    if ( this->M_bcRobin.empty() ) return;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );


    if ( buildNonCstPart )
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = this->fieldTemperature();

        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );
        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expression1( d,this->symbolsExpr() );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const 
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcDirichlet )
    {
        auto theExpr = expression(d,this->symbolsExpr());
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=theExpr );
    }

    this->log("Heat","updateLinearPDEDofElimination","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDESourceTerm( DataUpdateLinear & data ) const
{
    if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        this->M_overwritemethod_updateSourceTermLinearPDE(data);
        return;
    }

    if ( this->M_volumicForcesProperties.empty() ) return;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    if ( buildNonCstPart )
    {
        auto myLinearForm = form1( _test=this->spaceTemperature(), _vector=F,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldTemperature();

        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto theExpr = expression( d,this->symbolsExpr() );
            auto rangeBodyForceUsed = ( markers(d).empty() )? this->rangeMeshElements() : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEWeakBC( DataUpdateLinear & data ) const
{
    if ( this->M_bcNeumann.empty() && this->M_bcRobin.empty() ) return;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    if ( buildNonCstPart )
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = this->fieldTemperature();

        auto myLinearForm = form1( _test=Xh, _vector=F,
                                   _rowstart=this->rowStartInVector() );
        for( auto const& d : this->M_bcNeumann )
        {
            auto theExpr = expression( d,this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }

        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );
        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expression1( d,this->symbolsExpr() );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
            auto theExpr2 = expression2( d,this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }

    }
}




} // end namespace FeelModels
} // end namespace Feel
