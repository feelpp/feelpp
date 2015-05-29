
#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels2/modelvf/fluidmecstresstensor.hpp>

#define FSI_FLUID_USE_OPT_EXPR 1


namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateJacobianModel( element_fluid_type const& U,
                                                             sparse_matrix_ptrtype& J, vector_ptrtype& R,
                                                             bool _BuildCstPart ) const
{
#if defined(FEELMODELS_FLUID_BUILD_JACOBIAN_CODE)
    using namespace Feel::vf;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateJacobianModel", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer t1;

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    auto const rowStartInMatrix = this->rowStartInMatrix();
    auto const colStartInMatrix = this->colStartInMatrix();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

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
    auto const gammapoint2t = 2.0*inner(deft,defv);
    auto const gammapoint2v = 2.0*inner(defv,defv,mpl::int_<InnerProperties::IS_SAME>());
    //--------------------------------------------------------------------------------------------------//
    // newtonian law
    auto const sigma_newtonian_viscous = idv(*M_P0Mu)*deft;
    auto const Sigmat_newtonian = -idt(p)*Id + 2*idv(*M_P0Mu)*deft;
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    // boundaries conditions
    //auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    //--------------------------------------------------------------------------------------------------//

    if (this->stressTensorLawType() == "newtonian")
    {
        if (_BuildCstPart)
        {
#if 1
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= inner(Sigmat_newtonian,grad(v)),
                           _geomap=this->geomap() );
#else
            //auto StressTensorExprJac = Feel::vf::FSI::fluidMecNewtonianStressTensorJacobian(u,p,viscosityModel,false/*true*/);
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= 2*idv(*M_P0Mu)*inner(deft,grad(v)),
                           //_expr= inner( StressTensorExprJac, grad(v) ),
                           _geomap=this->geomap() );
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= -idt(p)*div(v),
                           _geomap=this->geomap() );
#endif

            //pressure bc condition
            if ( !this->markerPressureBC().empty() )
            {
                //ForEachBC( bcDef,cl::pressure,
                bilinearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,this->markerPressureBC()),
                               _expr= -inner( 2*sigma_newtonian_viscous*N(),id(v) ),
                               _geomap=this->geomap() );
            }
        }
    }
    else
    {
        if ( _BuildCstPart )
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= -idt(p)*div(v),
                           _geomap=this->geomap() );
#if FSI_FLUID_USE_OPT_EXPR
        if (!_BuildCstPart)
        {
#if 0
            auto viscosityModel = viscosityModelDesc( this->stressTensorLawType(),
                                                      *M_P0Mu,
                                                      this->powerLaw_n(),this->powerLaw_k(),
                                                      this->nonNewtonianMu0(),this->nonNewtonianMuInf(),
                                                      this->carreauLaw_lambda(), this->carreauLaw_n(),
                                                      this->carreauYasudaLaw_lambda(), this->carreauYasudaLaw_n(), this->carreauYasudaLaw_a(),
                                                      this->walburnSchneckLaw_C1(),this->walburnSchneckLaw_C2(),this->walburnSchneckLaw_C3(),this->walburnSchneckLaw_C4(),
                                                      this->nonNewtonianHematocrit(), this->nonNewtonianTPMA() );
#endif
            auto StressTensorExprJac = Feel::vf::FeelModels::fluidMecNewtonianStressTensorJacobian<2*nOrderVelocity>(u,p,*this->viscosityModel()/*viscosityModel*/,false/*true*/);
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           //_expr= inner( 2*sigma_powerlaw_viscous/*Sigmat_powerlaw*/,grad(v) ),
                           _expr= inner( StressTensorExprJac,grad(v) ),
                           _geomap=this->geomap() );
            //pressure bc condition
            if ( !this->markerPressureBC().empty() )
            {
                //ForEachBC( bcDef,cl::pressure,
                bilinearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,this->markerPressureBC()),
                               _expr= -inner(StressTensorExprJac*N(),id(v)),
                               _geomap=this->geomap() );
            }
        }

#else
        if (!_BuildCstPart)
        {
            // parameters
            double power_n=this->powerLaw_n(), power_k=this->powerLaw_k();
            double power_n_generic=(power_n-1.)/2., power_k_generic=power_k;
            double mu_0=this->nonNewtonianMu0();
            double mu_inf=this->nonNewtonianMuInf();
            double carreau_lambda=this->carreauLaw_lambda();
            double carreau_n=this->carreauLaw_n();
            double carreauYasuda_lambda=this->carreauYasudaLaw_lambda();
            double carreauYasuda_n=this->carreauYasudaLaw_n();
            double carreauYasuda_a=this->carreauYasudaLaw_a();
            if (this->stressTensorLawType() == "walburn-schneck_law")
            {
                double hematocrit = this->nonNewtonianHematocrit();
                power_k_generic=this->walburnSchneckLaw_C1()*math::exp(hematocrit*this->walburnSchneckLaw_C2())*math::exp( this->walburnSchneckLaw_C4()*this->nonNewtonianTPMA()/math::pow(hematocrit,2) );
                power_n_generic=-this->walburnSchneckLaw_C3()*hematocrit;
            }
            // power law
            //auto chiSup = chi(trace(defv*trans(defv)) >1e-9);
            //auto chiInv = chi(trace(defv*trans(defv)) <=1e-9);
            auto const muv_powerlaw = power_k_generic*pow( gammapoint2v /*+chiInv*/ ,  power_n_generic ) /**chiSup*/;
            auto const mut_powerlaw = gammapoint2t*power_k_generic*power_n_generic*pow( gammapoint2v , power_n_generic-1.0 );
            auto const sigma_powerlaw_viscous = val(muv_powerlaw)*deft  + mut_powerlaw*defv;
            //--------------------------------------------------------------------------------------------------//
            // carreau
            const double carreau_lambda2 = math::pow(carreau_lambda,2);
            auto const muv_carreauLaw = cst(mu_inf) + (cst(mu_0) - cst(mu_inf))*pow( cst(1.) + carreau_lambda2*gammapoint2v, (carreau_n-1)/2.0 );
            auto const mut_carreauLaw = gammapoint2t*( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2*pow( cst(1.) + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
            auto const sigma_carreauLaw_viscous = val(muv_carreauLaw)*deft  + mut_carreauLaw*defv;
            //--------------------------------------------------------------------------------------------------//
            // carreau-yasuda
            const double carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            auto const muv_carreauYasudaLaw = cst(mu_inf) + (cst(mu_0) - cst(mu_inf))*pow( cst(1.) + carreauYasuda_lambdaA*pow( gammapoint2v, carreauYasuda_a/2.) , (carreauYasuda_n-1)/carreauYasuda_a );
            auto const part1_carreauYasudaLaw = cst( ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf) );
            auto const part2_carreauYasudaLaw = cst(1.) + carreauYasuda_lambdaA*pow( gammapoint2v, carreauYasuda_a/2.);
            auto const part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
            auto const mut_carreauYasudaLaw = gammapoint2t*part1_carreauYasudaLaw*pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
            auto const sigma_carreauYasudaLaw_viscous = val(muv_carreauYasudaLaw)*deft  + mut_carreauYasudaLaw*defv;
            //--------------------------------------------------------------------------------------------------//
            // strain tensor
            auto const Sigmat_powerlaw = -idt(p)*Id + 2*sigma_powerlaw_viscous;
            auto const Sigmat_carreauLaw = -idt(p)*Id + 2*sigma_carreauLaw_viscous;
            auto const Sigmat_carreauYasudaLaw = -idt(p)*Id + 2*sigma_carreauYasudaLaw_viscous;

            if (this->stressTensorLawType() == "power_law" || this->stressTensorLawType() == "walburn-schneck_law")
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               //_expr= inner( Sigmat_powerlaw,grad(v) ),
                               _expr= inner( 2*sigma_powerlaw_viscous,grad(v) ),
                               _quad=_Q<2*nOrderVelocity+nOrderVelocity-1>(),
                               _geomap=this->geomap() );
                if ( M_pressureBCType.size() > 0 )
                {
                    //ForEachBC( bcDef,cl::pressure,
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,M_pressureBCType),
                                   _expr= -inner(2*sigma_powerlaw_viscous*N(),id(v)),
                                   _geomap=this->geomap() );
                }
            }
            else if (this->stressTensorLawType() == "carreau_law")
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               //_expr= inner( Sigmat_carreauLaw,grad(v) ),
                               _expr= inner( 2*sigma_carreauLaw_viscous,grad(v) ),
                               _geomap=this->geomap() );
                if ( M_pressureBCType.size() > 0 )
                {
                    //ForEachBC( bcDef,cl::pressure,
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,M_pressureBCType),
                                   _expr= -inner( 2*sigma_carreauLaw_viscous*N(),id(v) ),
                                   _quad=_Q<2*nOrderVelocity+nOrderVelocity-1>(),
                                   _geomap=this->geomap() );
                }
            }
            else if (this->stressTensorLawType() == "carreau-yasuda_law")
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               //_expr= inner( Sigmat_carreauYasudaLaw,grad(v) ),
                               _expr= inner( 2*sigma_carreauYasudaLaw_viscous,grad(v) ),
                               _quad=_Q<2*nOrderVelocity+nOrderVelocity-1>(),
                               _geomap=this->geomap() );
                if ( M_pressureBCType.size() > 0 )
                {
                    //ForEachBC( bcDef,cl::pressure,
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,M_pressureBCType),
                                   _expr= -inner( 2*sigma_carreauYasudaLaw_viscous*N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
        } // if (!_BuildCstPart)
#endif
    } // non newtonian
    //--------------------------------------------------------------------------------------------------//

    double timeElapsed=t1.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateJacobianModel",
                                               "finish in "+(boost::format("%1% s") % timeElapsed).str(),
                                               this->worldComm(),this->verboseAllProc());

#endif // defined(FEELMODELS_FLUID_BUILD_JACOBIAN_CODE)
} // updateJacobian

} // namespace FeelModels

} // namespace Feel
