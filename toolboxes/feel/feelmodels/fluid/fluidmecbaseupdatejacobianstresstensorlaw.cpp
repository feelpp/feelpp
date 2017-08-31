
#include <feel/feelmodels/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>


namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateJacobianModel( DataUpdateJacobian & data, element_fluid_external_storage_type const& U ) const
{
    using namespace Feel::vf;

    this->log("FluidMechanics","updateJacobianModel", "start" );
    boost::mpi::timer t1;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool _BuildCstPart = data.buildCstPart();

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    if ( this->densityViscosityModel()->dynamicViscosityLaw() == "newtonian")
    {
        //--------------------------------------------------------------------------------------------------//
        // identity matrix
        auto const Id = eye<nDim,nDim>();
        //--------------------------------------------------------------------------------------------------//
        // shear tensor
        //bool useSymTensor = option(_name="strain_tensor.use-sym-tensor",_prefix=this->prefix()).as<bool>();
        //auto const defv = Feel::vf::FSI::fluidStrainTensor(gradv(u),useSymTensor); // D(u) = 0.5*(a+at) or D(u)=a
        //auto const deft = Feel::vf::FSI::fluidStrainTensor(gradt(u),useSymTensor); // D(u) = 0.5*(a+at) or D(u)=a
        auto const deft = sym(gradt(u));
        auto const defv = sym(gradv(u));
        //--------------------------------------------------------------------------------------------------//
        // newtonian law
        auto const& mu = this->densityViscosityModel()->fieldMu();
        auto const sigma_newtonian_viscous = idv(mu)*deft;
        auto const Sigmat_newtonian = -idt(p)*Id + 2*idv(mu)*deft;
        //--------------------------------------------------------------------------------------------------//
        if (_BuildCstPart)
        {
#if 1
            bilinearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= inner(Sigmat_newtonian,grad(v)),
                           _geomap=this->geomap() );
#else
            //auto StressTensorExprJac = Feel::vf::FSI::fluidMecNewtonianStressTensorJacobian(u,p,viscosityModel,false/*true*/);
            bilinearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= 2*idv(mu)*inner(deft,grad(v)),
                           //_expr= inner( StressTensorExprJac, grad(v) ),
                           _geomap=this->geomap() );
            bilinearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -idt(p)*div(v),
                           _geomap=this->geomap() );
#endif

        }
    }
    else
    {
        if ( _BuildCstPart )
            bilinearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -idt(p)*div(v),
                           _geomap=this->geomap() );

        if (!_BuildCstPart)
        {
            auto StressTensorExprJac = Feel::vf::FeelModels::fluidMecNewtonianStressTensorJacobian<2*nOrderVelocity>(u,p,*this->densityViscosityModel(),false/*true*/);
            bilinearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           //_expr= inner( 2*sigma_powerlaw_viscous/*Sigmat_powerlaw*/,grad(v) ),
                           _expr= inner( StressTensorExprJac,grad(v) ),
                           _geomap=this->geomap() );
        }
    } // non newtonian
    //--------------------------------------------------------------------------------------------------//

    double timeElapsed=t1.elapsed();
    this->log("FluidMechanics","updateJacobianModel","finish in "+(boost::format("%1% s") % timeElapsed).str() );

} // updateJacobian

} // namespace FeelModels

} // namespace Feel
