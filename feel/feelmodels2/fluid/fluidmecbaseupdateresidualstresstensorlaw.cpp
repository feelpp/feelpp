/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

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
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidualModel( element_fluid_type const& U,
                                                             vector_ptrtype& R,
                                                             bool BuildCstPart,
                                                             bool UseJacobianLinearTerms ) const
{
#if defined(FEELMODELS_FLUID_BUILD_RESIDUAL_CODE)
    using namespace Feel::vf;

    //if (this->stressTensorLawType() == "newtonian" && UseJacobianLinearTerms ) return;

    //!this->velocityDivIsEqualToZero() && BuildCstPart

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidualModel", "start for " + this->stressTensorLawType(),
                                               this->worldComm(),this->verboseAllProc());

    boost::mpi::timer thetimer,thetimer2;
    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    //--------------------------------------------------------------------------------------------------//

    auto const rowStartInMatrix = this->rowStartInMatrix();
    auto const colStartInMatrix = this->colStartInMatrix();
    auto const rowStartInVector = this->rowStartInVector();
    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=rowStartInVector );

    //--------------------------------------------------------------------------------------------------//

    //--------------------------------------------------------------------------------------------------//
    // identity matrix
    auto const Id = eye<nDim,nDim>();
    //--------------------------------------------------------------------------------------------------//
    // strain tensor
    //bool useSymTensor = option(_name="strain_tensor.use-sym-tensor",_prefix=this->prefix()).as<bool>();
    //auto const defv = Feel::vf::FSI::fluidStrainTensor(gradv(u),useSymTensor); // D(u) = 0.5*(a+at) or D(u)=a
    auto const defv = sym(gradv(u));
    //--------------------------------------------------------------------------------------------------//
    // newtonian law
    auto const mu_newtonian = idv(*M_P0Mu);
    auto const Sigmav_newtonian = -idv(p)*Id + 2*mu_newtonian*defv;
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    // boundaries conditions
    //auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    //--------------------------------------------------------------------------------------------------//

    double timeElapsedBis = thetimer2.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidualModel",
                                               "finish init in "+(boost::format("%1% s") % timeElapsedBis).str(),
                                               this->worldComm(),this->verboseAllProc());

    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    if (this->stressTensorLawType() == "newtonian")
    {
        // sigma : grad(v) on Omega
        if (!BuildCstPart && !UseJacobianLinearTerms )
        {
            Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidualModel",
                                  "build stress stensor law", this->worldComm(),this->verboseAllProc());
#if 1
            linearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           //_expr= inner( StressTensorExpr,grad(v) ),
                           _expr= inner( val(Sigmav_newtonian),grad(v) ),
                           _geomap=this->geomap() );
#else
            form1( Xh, R ) +=
                integrate( _range=elements(mesh),
                           _expr= 2*idv(*M_P0Mu)*trace(trans(defv)*grad(v)),
                           _geomap=this->geomap() );
            form1( Xh, R ) +=
                integrate( _range=elements(mesh),
                           _expr= -idv(p)*div(v),
                           _geomap=this->geomap() );
#endif

            if ( !this->markerPressureBC().empty() )
            {
                linearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,this->markerPressureBC()),
                               //_expr= -inner( ViscousStressTensorExpr*N(),id(v) ),
                               _expr= -inner( 2*mu_newtonian*defv*N(),id(v) ),
                               _geomap=this->geomap() );
            }
        }
    }
    else
    {
        if (!BuildCstPart && !UseJacobianLinearTerms )
        {
            linearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= -idv(p)*div(v),
                           _geomap=this->geomap() );
        }
#if FSI_FLUID_USE_OPT_EXPR
        if (!BuildCstPart)
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
            auto const StressTensorExpr = Feel::vf::FeelModels::fluidMecNewtonianStressTensor<2*nOrderVelocity>(u,p,*this->viscosityModel()/*viscosityModel*/,false/*true*/);
            // sigma : grad(v) on Omega
            linearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= inner( StressTensorExpr,grad(v) ),
                           _geomap=this->geomap() );
            //pressure condition
            if ( !this->markerPressureBC().empty() )
            {
                linearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,this->markerPressureBC()),
                               _expr= -inner( StressTensorExpr*N(),id(v) ),
                               _geomap=this->geomap() );
            }
        }
#else
        if (!BuildCstPart)
        {
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
            auto const gammapoint2v = 2.0*inner(defv,defv,mpl::int_<InnerProperties::IS_SAME>());
            // power law
            //auto chiSup = chi(trace(defv*trans(defv)) >1e-9);
            //auto chiInv = chi(trace(defv*trans(defv)) <=1e-9);
            auto const mu_powerlaw = power_k_generic*pow( gammapoint2v /*+chiInv*/ , cst( power_n_generic ) )/**chiSup*/;
            //--------------------------------------------------------------------------------------------------//
            // carreau
            auto const mu_carreauLaw = cst(mu_inf) + (cst(mu_0) - cst(mu_inf))*pow( cst(1.) + math::pow(carreau_lambda,2.)*gammapoint2v, (carreau_n-1)/2.0 );
            //--------------------------------------------------------------------------------------------------//
            // carreau-yasuda
            auto const mu_carreauYasudaLaw = cst(mu_inf) + (cst(mu_0) - cst(mu_inf))*
                pow( cst(1.) + math::pow(carreauYasuda_lambda,carreauYasuda_a)*pow( gammapoint2v, carreauYasuda_a/2.) , (carreauYasuda_n-1)/carreauYasuda_a );
            //--------------------------------------------------------------------------------------------------//
            // strain tensor
            auto const Sigmav_powerlaw = -idv(p)*Id + 2*mu_powerlaw*defv;
            auto const Sigmav_carreauLaw = -idv(p)*Id + 2*mu_carreauLaw*defv;
            auto const Sigmav_carreauYasudaLaw = -idv(p)*Id + 2*mu_carreauYasudaLaw*defv;

            if (this->stressTensorLawType() == "power_law" || this->stressTensorLawType() == "walburn-schneck_law")
            {
                linearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               //_expr= inner( val(Sigmav_powerlaw),grad(v) ),
                               _expr= inner( val(2*mu_powerlaw*defv),grad(v) ),
                               _quad=_Q<2*nOrderVelocity+nOrderVelocity-1>(),
                               _geomap=this->geomap() );
                if ( M_pressureBCType.size() > 0 )
                {
                    //ForEachBC( bcDef,cl::pressure,
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,M_pressureBCType),
                                   _expr= -inner( 2*mu_powerlaw*defv*N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
            else if (this->stressTensorLawType() == "carreau_law")
            {
                linearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               //_expr= inner( val(Sigmav_carreauLaw),grad(v) ),
                               _expr= inner( val(2*mu_carreauLaw*defv),grad(v) ),
                               _quad=_Q<2*nOrderVelocity+nOrderVelocity-1>(),
                               _geomap=this->geomap() );
                if ( M_pressureBCType.size() > 0 )
                {
                    //ForEachBC( bcDef,cl::pressure,
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,M_pressureBCType),
                                   _expr= -inner( 2*mu_carreauLaw*defv*N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
            else if (this->stressTensorLawType() == "carreau-yasuda_law")
            {
                linearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               //_expr= inner( val(Sigmav_carreauYasudaLaw),grad(v) ),
                               _expr= inner( val(2*mu_carreauYasudaLaw*defv),grad(v) ),
                               _quad=_Q<2*nOrderVelocity+nOrderVelocity-1>(),
                               _geomap=this->geomap() );
                if ( M_pressureBCType.size() > 0 )
                {
                    //ForEachBC( bcDef,cl::pressure,
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,M_pressureBCType),
                                   _expr= -inner( 2*mu_carreauYasudaLaw*defv*N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
        } // if (!BuildCstPart)
#endif
    } // non newtonian
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    // take into account that div u != 0
    if (!this->velocityDivIsEqualToZero() && BuildCstPart)
    {
        linearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= -idv(this->velocityDiv())*id(q),
                       _geomap=this->geomap() );

        auto coeffDiv = (2./3.)*idv(*M_P0Mu); //(eps-2mu/3)
        linearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= val(-coeffDiv*gradv(this->velocityDiv()))*id(v),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidualModel",
                                               "finish in "+(boost::format("%1% s") % timeElapsed).str(),
                                               this->worldComm(),this->verboseAllProc());

#endif // defined(FEELMODELS_FLUID_BUILD_RESIDUAL_CODE)

} // updateResidual

} // namespace FeelModels
} // namespace Feel
