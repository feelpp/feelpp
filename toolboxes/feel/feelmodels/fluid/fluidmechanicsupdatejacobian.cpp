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

    bool Build_TransientTerm = BuildNonCstPart;
    if ( !this->isStationaryModel() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        Build_TransientTerm=BuildCstPart;


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

    auto se = this->symbolsExpr();
    //--------------------------------------------------------------------------------------------------//

    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // strain tensor
    auto const deft = sym(gradt(u));
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    boost::mpi::timer timerAssemble;

    //--------------------------------------------------------------------------------------------------//
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );

            // stress tensor sigma : grad(v)
            if ( BuildCstPart )
                bilinearFormVP +=
                    integrate( _range=range,
                               _expr= -idt(p)*div(v),
                               _geomap=this->geomap() );

            if ( ( physicFluidData->dynamicViscosity().isNewtonianLaw() && BuildCstPart ) ||
                 ( !physicFluidData->dynamicViscosity().isNewtonianLaw() && BuildNonCstPart ) )
            {
                auto StressTensorExprJac = Feel::FeelModels::fluidMecNewtonianViscousStressTensorJacobian(gradv(u),u,*physicFluidData,matProps);
                bilinearFormVV_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( StressTensorExprJac,grad(v) ),
                               _geomap=this->geomap() );
            }

            // convection terms
            if ( physicFluidData->equation() == "Navier-Stokes" )
            {
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                if ( BuildNonCstPart )
                {
                    if (this->doStabConvectionEnergy())
                    {
                        // convection term + stabilisation energy of convection with neumann bc (until outflow bc) ( see Nobile thesis)
                        // auto const convecTerm = (trans(val(gradv(u)*idv(*M_P0Rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(*M_P0Rho)) ) )*id(v);
                        // stabTerm = trans(divt(u)*val(0.5*idv(*M_P0Rho)*idv(u))+val(0.5*idv(*M_P0Rho)*divv(u))*idt(u))*id(v)

                        auto const convecTerm = densityExpr*Feel::FeelModels::fluidMecConvectionJacobianWithEnergyStab(u);
                        bilinearFormVV_PatternCoupled +=
                            //bilinearForm_PatternDefault +=
                            integrate ( _range=range,
                                        _expr=timeSteppingScaling*convecTerm,
                                        _geomap=this->geomap() );
                    }
                    else
                    {
                        //auto const convecTerm = (trans(val(gradv(u)*idv(rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(rho)) ) )*id(v);
                        auto const convecTerm = densityExpr*Feel::FeelModels::fluidMecConvectionJacobian(u);

                        bilinearFormVV_PatternCoupled +=
                            //bilinearForm_PatternDefault +=
                            integrate ( _range=range,
                                        _expr=timeSteppingScaling*convecTerm,
                                        _geomap=this->geomap() );
                    }
                }
#if defined( FEELPP_MODELS_HAS_MESHALE )
                if (this->isMoveDomain() && BuildCstPart )
                {
                    bilinearFormVV_PatternCoupled +=
                        //bilinearForm_PatternDefault +=
                        integrate (_range=range,
                                   _expr= -timeSteppingScaling*trans(gradt(u)*densityExpr*idv( this->meshVelocity() ))*id(v),
                                   _geomap=this->geomap() );
                }
#endif

                if ( data.hasVectorInfo( "explicit-part-of-solution" ) && BuildCstPart )
                {
                    auto uExplicitPartOfSolution = XhV->element( data.vectorInfo( "explicit-part-of-solution" ), rowStartInVector+startBlockIndexVelocity );
                    bilinearFormVV_PatternCoupled +=
                        integrate ( _range=range,
                                    _expr= timeSteppingScaling*densityExpr*trans( gradt(u)*idv(uExplicitPartOfSolution) + gradv(uExplicitPartOfSolution )*idt(u) )*id(v),
                                    _geomap=this->geomap() );
                }
            }

            //transients terms
            if ( !this->isStationaryModel() && Build_TransientTerm )
            {
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                bilinearFormVV_PatternDefault +=
                    integrate( _range=range,
                               _expr= densityExpr*inner(idt(u),id(v))*M_bdfVelocity->polyDerivCoefficient(0),
                               _geomap=this->geomap() );
            }

        } // foreach mat
    } // foreach physic


    double timeElapsed = timerAssemble.elapsed();
    this->log("FluidMechanics","updateJacobian",(boost::format("assemble convection term in %1% s") %timeElapsed).str() );

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
    // peusdo transient continuation
    if ( BuildNonCstPart && data.hasInfo( "use-pseudo-transient-continuation" ) )
    {
#if 0
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
#else
        CHECK( false ) << "TODO VINCENT";
#endif
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


