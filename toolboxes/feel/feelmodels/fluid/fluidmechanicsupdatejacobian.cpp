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

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");
    size_type startBlockIndexPressure = this->startSubBlockSpaceIndex("pressure");
    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    size_type rowStartInVector = this->rowStartInVector();
    auto bilinearFormVV_PatternDefault = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::DEFAULT),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto bilinearFormVV_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix+0,
                                 _colstart=colStartInMatrix+1 );
    auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix+1,
                                 _colstart=colStartInMatrix+0 );

    auto u = XhV->element(XVec, rowStartInVector+0);
    auto const& v = u;
    auto p = XhP->element(XVec, rowStartInVector+1);
    auto const& q = p;

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
    CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);

        // convection terms
        if ( physicFluidData->equation() == "Navier-Stokes" )
        {
            if ( BuildNonCstPart )
            {
                if (this->doStabConvectionEnergy())
                {
                    // convection term + stabilisation energy of convection with neumann bc (until outflow bc) ( see Nobile thesis)
                    // auto const convecTerm = (trans(val(gradv(u)*idv(*M_P0Rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(*M_P0Rho)) ) )*id(v);
                    // stabTerm = trans(divt(u)*val(0.5*idv(*M_P0Rho)*idv(u))+val(0.5*idv(*M_P0Rho)*divv(u))*idt(u))*id(v)

                    auto const convecTerm = Feel::FeelModels::fluidMecConvectionJacobianWithEnergyStab(u,rho);
                    bilinearFormVV_PatternCoupled +=
                        //bilinearForm_PatternDefault +=
                        integrate ( _range=M_rangeMeshElements,
                                    _expr=timeSteppingScaling*convecTerm,
                                    _geomap=this->geomap() );
                }
                else
                {
#if 0
                    auto const convecTerm = (trans(val(gradv(u)*idv(rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(rho)) ) )*id(v);
#else
                    auto const convecTerm = Feel::FeelModels::fluidMecConvectionJacobian(u,rho);
#endif
                    bilinearFormVV_PatternCoupled +=
                        //bilinearForm_PatternDefault +=
                        integrate ( _range=M_rangeMeshElements,
                                    _expr=timeSteppingScaling*convecTerm,
                                    _geomap=this->geomap() );
                }
            }
            if ( data.hasVectorInfo( "explicit-part-of-solution" ) && BuildCstPart )
            {
                auto uExplicitPartOfSolution = XhV->element( data.vectorInfo( "explicit-part-of-solution" ), rowStartInVector+startBlockIndexVelocity );
                bilinearFormVV_PatternCoupled +=
                    integrate ( _range=M_rangeMeshElements,
                                _expr= timeSteppingScaling*idv(rho)*trans( gradt(u)*idv(uExplicitPartOfSolution) + gradv(uExplicitPartOfSolution )*idt(u) )*id(v),
                                _geomap=this->geomap() );
            }
        }
    } // foreach physic

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain() && BuildCstPart )
    {
        bilinearFormVV_PatternCoupled +=
            //bilinearForm_PatternDefault +=
            integrate (_range=M_rangeMeshElements,
                       _expr= -timeSteppingScaling*trans(gradt(u)*idv(rho)*idv( this->meshVelocity() ))*id(v),
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
        if ( BuildCstPart )
            bilinearFormVP +=
                integrate( _range=range,
                           _expr= -idt(p)*div(v),
                           _geomap=this->geomap() );

        if ( ( dynamicViscosity.isNewtonianLaw() && BuildCstPart ) ||
             ( !dynamicViscosity.isNewtonianLaw() && BuildNonCstPart ) )
        {
            auto StressTensorExprJac = Feel::FeelModels::fluidMecNewtonianViscousStressTensorJacobian(gradv(u),u,*this->materialProperties(),matName);
            bilinearFormVV_PatternCoupled +=
                integrate( _range=range,
                           _expr= timeSteppingScaling*inner( StressTensorExprJac,grad(v) ),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // incompressibility term
    if (BuildCstPart)
    {
        bilinearFormPV +=
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
        bilinearFormVV_PatternDefault +=
            integrate( _range=M_rangeMeshElements,
                       _expr= idv(rho)*trans(idt(u))*id(v)*M_bdfVelocity->polyDerivCoefficient(0),
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

        bilinearFormVV_PatternDefault +=
            integrate(_range=M_rangeMeshElements,
                      _expr=(1./pseudoTimeStepDelta)*idv(norm2_uu)*inner(idt(u),id(u)),
                      //_expr=(1./pseudoTimeStepDelta)*(norm2(idv(u))/h())*inner(idt(u),id(u)),
                      //_expr=(1./pseudoTimeStepDelta)*inner(idt(u),id(u)),
                      _geomap=this->geomap() );
    }
    //--------------------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && BuildCstPart )
        {
            auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=J,
                                         _pattern=size_type(Pattern::COUPLED),
                                         _rowstart=rowStartInMatrix+1,
                                         _colstart=colStartInMatrix+1 );
            double beta = this->definePressureCstPenalisationBeta();
            for ( auto const& rangeElt : M_definePressureCstMeshRanges )
                bilinearFormPP +=
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
                form2( _test=XhP, _trial=M_XhMeanPressureLM[k], _matrix=J,
                       _rowstart=this->rowStartInMatrix()+1,
                       _colstart=this->colStartInMatrix()+startBlockIndexDefinePressureCstLM+k ) +=
                    integrate( _range=M_definePressureCstMeshRanges[k],
                               _expr= id(p)*idt(lambda) /*+ idt(p)*id(lambda)*/,
                               _geomap=this->geomap() );

                form2( _test=M_XhMeanPressureLM[k], _trial=XhP, _matrix=J,
                       _rowstart=this->rowStartInMatrix()+startBlockIndexDefinePressureCstLM+k,
                       _colstart=this->colStartInMatrix()+1 ) +=
                    integrate( _range=M_definePressureCstMeshRanges[k],
                               _expr= + idt(p)*id(lambda),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->updateJacobianStabilisation( data, u,p );

    //--------------------------------------------------------------------------------------------------//

    this->updateJacobianWeakBC( data, u,p );

    //--------------------------------------------------------------------------------------------------//
    /*double*/ timeElapsed=thetimer.elapsed();
    this->log("FluidMechanics","updateJacobian",(boost::format("finish %1% in %2% s") %sc %timeElapsed).str() );

} // updateJacobian

} // namespace FeelModels

} // namespace Feel


