/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>

//#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>


namespace Feel
{
namespace FeelModels
{

#if 0
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualWeakBC( DataUpdateResidual & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const
{
    this->log("FluidMechanics","updateResidualModel", "start" );

    boost::mpi::timer thetimer,thetimer2;
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !BuildCstPart;

    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );
    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto const& v = u;
    auto const& q = p;

    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormV = form1( _test=XhV, _vector=R,
                              _pattern=size_type(Pattern::COUPLED),
                              _rowstart=this->rowStartInVector() );


    auto Id = eye<nDim,nDim>();

    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateResidualModel","finish in "+(boost::format("%1% s") % timeElapsed).str() );


} // updateResidual
#endif
} // namespace FeelModels
} // namespace Feel
