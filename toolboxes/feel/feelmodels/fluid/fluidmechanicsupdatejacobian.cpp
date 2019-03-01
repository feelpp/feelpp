/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelvf/fluidmecconvection.hpp>
#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool _BuildCstPart = data.buildCstPart();

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("FluidMechanics","updateJacobian",(boost::format("start %1%") %sc).str() );
    boost::mpi::timer thetimer;

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    //bool BuildNonCstPart_robinFSI = BuildNonCstPart;
    //if (this->useFSISemiImplicitScheme()) BuildNonCstPart_robinFSI=BuildCstPart;

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    size_type rowStartInVector = this->rowStartInVector();
    auto bilinearForm_PatternDefault = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::DEFAULT),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

    auto U = Xh->element(XVec, rowStartInVector);
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    //--------------------------------------------------------------------------------------------------//

    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // strain tensor
    auto const deft = sym(gradt(u));
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();
    auto const& rho = this->materialProperties()->fieldRho();
    // stress tensor
    auto const Sigmat = -idt(p)*Id + 2*idv(mu)*deft;
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    boost::mpi::timer timerAssemble;

    //--------------------------------------------------------------------------------------------------//
    // convection terms
    if ( BuildNonCstPart )
    {
        if (this->doStabConvectionEnergy())
        {
            // convection term + stabilisation energy of convection with neumann bc (until outflow bc) ( see Nobile thesis)
            // auto const convecTerm = (trans(val(gradv(u)*idv(*M_P0Rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(*M_P0Rho)) ) )*id(v);
            // stabTerm = trans(divt(u)*val(0.5*idv(*M_P0Rho)*idv(u))+val(0.5*idv(*M_P0Rho)*divv(u))*idt(u))*id(v)

            auto const convecTerm = Feel::FeelModels::fluidMecConvectionJacobianWithEnergyStab(u,rho);
            bilinearForm_PatternCoupled +=
                //bilinearForm_PatternDefault +=
                integrate ( _range=M_rangeMeshElements,
                            _expr=convecTerm,
                            _geomap=this->geomap() );
        }
        else
        {
#if 0
            auto const convecTerm = (trans(val(gradv(u)*idv(rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(rho)) ) )*id(v);
#else
            auto const convecTerm = Feel::FeelModels::fluidMecConvectionJacobian(u,rho);
#endif
            bilinearForm_PatternCoupled +=
                //bilinearForm_PatternDefault +=
                integrate ( _range=M_rangeMeshElements,
                            _expr=convecTerm,
                            _geomap=this->geomap() );
        }
    }

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain() && BuildCstPart )
    {
        bilinearForm_PatternCoupled +=
            //bilinearForm_PatternDefault +=
            integrate (_range=M_rangeMeshElements,
                       _expr= -trans(gradt(u)*idv(rho)*idv( this->meshVelocity() ))*id(v),
                       _geomap=this->geomap() );
    }
#endif

    double timeElapsed = timerAssemble.elapsed();
    this->log("FluidMechanics","updateJacobian",(boost::format("assemble convection term in %1% s") %timeElapsed).str() );

    //--------------------------------------------------------------------------------------------------//
    // sigma : grad(v) on Omega
    for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& dynamicViscosity = this->materialProperties()->dynamicViscosity(matName);
        if ( dynamicViscosity.isNewtonianLaw() )
        {
            //auto const deft = sym(gradt(u));
            //--------------------------------------------------------------------------------------------------//
            // newtonian law
            auto const& mu = this->materialProperties()->fieldMu();
            auto const sigma_newtonian_viscous = idv(mu)*deft;
            auto const Sigmat_newtonian = -idt(p)*Id + 2*idv(mu)*deft;
            //--------------------------------------------------------------------------------------------------//
            if ( BuildCstPart )
            {
#if 1
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= inner(Sigmat_newtonian,grad(v)),
                               _geomap=this->geomap() );
#else
                //auto StressTensorExprJac = Feel::vf::FSI::fluidMecNewtonianStressTensorJacobian(u,p,viscosityModel,false/*true*/);
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= 2*idv(mu)*inner(deft,grad(v)),
                               //_expr= inner( StressTensorExprJac, grad(v) ),
                               _geomap=this->geomap() );
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= -idt(p)*div(v),
                               _geomap=this->geomap() );
#endif

            }
        }
        else
        {
            if ( BuildCstPart )
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= -idt(p)*div(v),
                               _geomap=this->geomap() );

            if ( BuildNonCstPart )
            {
                auto StressTensorExprJac = Feel::FeelModels::fluidMecNewtonianStressTensorJacobian<2*nOrderVelocity>(u,p,*this->materialProperties(),matName);
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               //_expr= inner( 2*sigma_powerlaw_viscous/*Sigmat_powerlaw*/,grad(v) ),
                               _expr= inner( StressTensorExprJac,grad(v) ),
                               _geomap=this->geomap() );
            }
        } // non newtonian
    }

    //--------------------------------------------------------------------------------------------------//
    // incompressibility term
    if (BuildCstPart)
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -divt(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    //transients terms
    bool Build_TransientTerm = !BuildCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=BuildCstPart;

    if (!this->isStationaryModel() && Build_TransientTerm/*BuildCstPart*/)
    {
        bilinearForm_PatternDefault +=
            integrate( _range=M_rangeMeshElements,
                       _expr= idv(rho)*trans(idt(u))*id(v)*M_bdf_fluid->polyDerivCoefficient(0),
                       _geomap=this->geomap() );
    }

    // peusdo transient continuation
    if ( BuildNonCstPart && data.hasInfo( "use-pseudo-transient-continuation" ) )
    {
        double pseudoTimeStepDelta = data.doubleInfo("pseudo-transient-continuation.delta");
        auto norm2_uu = this->materialProperties()->fieldRho().functionSpace()->element(); // TODO : improve this (maybe create an expression instead)
        //norm2_uu.on(_range=M_rangeMeshElements,_expr=norm2(idv(u))/h());
        auto fieldNormu = u.functionSpace()->compSpace()->element( norm2(idv(u)) );
        auto maxu = fieldNormu.max( this->materialProperties()->fieldRho().functionSpace() );
        //auto maxux = u[ComponentType::X].max( this->materialProperties()->fieldRho().functionSpace() );
        //auto maxuy = u[ComponentType::Y].max( this->materialProperties()->fieldRho().functionSpace() );
        //norm2_uu.on(_range=M_rangeMeshElements,_expr=norm2(vec(idv(maxux),idv(maxux)))/h());
        norm2_uu.on(_range=M_rangeMeshElements,_expr=idv(maxu)/h());
        
        bilinearForm_PatternDefault +=
            integrate(_range=M_rangeMeshElements,
                      _expr=(1./pseudoTimeStepDelta)*idv(norm2_uu)*inner(idt(u),id(u)),
                      //_expr=(1./pseudoTimeStepDelta)*(norm2(idv(u))/h())*inner(idt(u),id(u)),
                      //_expr=(1./pseudoTimeStepDelta)*inner(idt(u),id(u)),
                      _geomap=this->geomap() );
    }
    //--------------------------------------------------------------------------------------------------//
    // user-defined additional terms
    this->updateJacobianAdditional( J, _BuildCstPart );

    //--------------------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && BuildCstPart )
        {
            double beta = this->definePressureCstPenalisationBeta();
            for ( auto const& rangeElt : M_definePressureCstMeshRanges )
                bilinearForm_PatternCoupled +=
                    integrate( _range=rangeElt,
                               _expr=beta*idt(p)*id(q),
                               _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" && BuildCstPart )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("define-pressure-cst-lm") ) << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startBlockIndexDefinePressureCstLM = this->startSubBlockSpaceIndex("define-pressure-cst-lm");
            for ( int k=0;k<M_XhMeanPressureLM.size();++k )
            {
                auto lambda = M_XhMeanPressureLM[k]->element();
                form2( _test=Xh, _trial=M_XhMeanPressureLM[k], _matrix=J,
                       _rowstart=this->rowStartInMatrix(),
                       _colstart=this->colStartInMatrix()+startBlockIndexDefinePressureCstLM+k ) +=
                    integrate( _range=M_definePressureCstMeshRanges[k],
                               _expr= id(p)*idt(lambda) /*+ idt(p)*id(lambda)*/,
                               _geomap=this->geomap() );

                form2( _test=M_XhMeanPressureLM[k], _trial=Xh, _matrix=J,
                       _rowstart=this->rowStartInMatrix()+startBlockIndexDefinePressureCstLM+k,
                       _colstart=this->colStartInMatrix() ) +=
                    integrate( _range=M_definePressureCstMeshRanges[k],
                               _expr= + idt(p)*id(lambda),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->updateJacobianStabilisation( data, U );

    //--------------------------------------------------------------------------------------------------//

    this->updateJacobianWeakBC( data, U );

    //--------------------------------------------------------------------------------------------------//
    /*double*/ timeElapsed=thetimer.elapsed();
    this->log("FluidMechanics","updateJacobian",(boost::format("finish %1% in %2% s") %sc %timeElapsed).str() );

} // updateJacobian

} // namespace FeelModels

} // namespace Feel


