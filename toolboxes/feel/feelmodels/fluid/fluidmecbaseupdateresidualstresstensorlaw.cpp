/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>


namespace Feel
{
namespace FeelModels
{


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidualModel( DataUpdateResidual & data, element_fluid_external_storage_type const& U ) const
{
    using namespace Feel::vf;

    this->log("FluidMechanics","updateResidualModel", "start for " + this->densityViscosityModel()->dynamicViscosityLaw() );

    boost::mpi::timer thetimer,thetimer2;
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//
    double timeElapsedBis = thetimer2.elapsed();
    this->log("FluidMechanics","updateResidualModel",
              "finish init in "+(boost::format("%1% s") % timeElapsedBis).str() );

    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    if ( this->densityViscosityModel()->dynamicViscosityLaw() == "newtonian")
    {
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
        auto const& mu = this->densityViscosityModel()->fieldMu();
        auto const mu_newtonian = idv(mu);
        auto const Sigmav_newtonian = -idv(p)*Id + 2*mu_newtonian*defv;
        //--------------------------------------------------------------------------------------------------//
        // sigma : grad(v) on Omega
        if (!BuildCstPart && !UseJacobianLinearTerms )
        {
            this->log("FluidMechanics","updateResidualModel","assembly with newtonian viscosity" );
#if 1
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           //_expr= inner( StressTensorExpr,grad(v) ),
                           _expr= inner( val(Sigmav_newtonian),grad(v) ),
                           _geomap=this->geomap() );
#else
            form1( Xh, R ) +=
                integrate( _range=M_rangeMeshElements,
                           _expr= 2*idv(*M_P0Mu)*trace(trans(defv)*grad(v)),
                           _geomap=this->geomap() );
            form1( Xh, R ) +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -idv(p)*div(v),
                           _geomap=this->geomap() );
#endif

        }
    }
    else
    {
        if (!BuildCstPart && !UseJacobianLinearTerms )
        {
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -idv(p)*div(v),
                           _geomap=this->geomap() );
        }

        if (!BuildCstPart)
        {
            auto const StressTensorExpr = Feel::vf::FeelModels::fluidMecNewtonianStressTensor<2*nOrderVelocity>(u,p,*this->densityViscosityModel(),false/*true*/);
            // sigma : grad(v) on Omega
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= inner( StressTensorExpr,grad(v) ),
                           _geomap=this->geomap() );
        }
    } // non newtonian
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    // take into account that div u != 0
    if (!this->velocityDivIsEqualToZero() && BuildCstPart)
    {
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= idv(this->velocityDiv())*id(q),
                       _geomap=this->geomap() );

        auto coeffDiv = (2./3.)*idv(this->densityViscosityModel()->fieldMu()); //(eps-2mu/3)
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= val(-coeffDiv*gradv(this->velocityDiv()))*id(v),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateResidualModel","finish in "+(boost::format("%1% s") % timeElapsed).str() );


} // updateResidual

} // namespace FeelModels
} // namespace Feel
