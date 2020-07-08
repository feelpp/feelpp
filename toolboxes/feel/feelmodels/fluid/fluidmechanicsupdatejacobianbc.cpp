/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>



namespace Feel
{
namespace FeelModels
{

#if 0
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianWeakBC( DataUpdateJacobian & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const
{
    using namespace Feel::vf;

    this->log("FluidMechanics","updateJacobianWeakBC", "start" );
    boost::mpi::timer t1;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool BuildCstPart = data.buildCstPart();

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();
    auto const& v = u;
    auto const& q = p;

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    auto bilinearFormVV = form2( _test=XhV,_trial=XhV,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix,
                                 _colstart=colStartInMatrix );
    // identity Matrix
    auto const Id = eye<nDim,nDim>();

    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    double timeElapsed=t1.elapsed();
    this->log("FluidMechanics","updateJacobianWeakBC","finish in "+(boost::format("%1% s") % timeElapsed).str() );

} // updateJacobian

#endif




} // namespace FeelModels

} // namespace Feel
