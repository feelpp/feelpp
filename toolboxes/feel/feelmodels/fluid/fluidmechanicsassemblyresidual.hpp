#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_ASSEMBLY_RESIDUAL_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_ASSEMBLY_RESIDUAL_HPP 1

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !BuildCstPart;

    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("FluidMechanics","updateResidual", "start"+sc );
    this->timerTool("Solve").start();

    //--------------------------------------------------------------------------------------------------//

    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );
    // if ( !doAssemblyRhs )
    //     std::cout << "hola \n";

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationaryModel() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( M_timeStepping == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - M_timeStepThetaValue;
            else
                timeSteppingScaling = M_timeStepThetaValue;
        }
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");
    size_type startBlockIndexPressure = this->startSubBlockSpaceIndex("pressure");
    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormV_PatternDefault = form1( _test=XhV, _vector=R,
                                            _pattern=size_type(Pattern::DEFAULT),
                                            _rowstart=rowStartInVector+startBlockIndexVelocity );
    auto linearFormV = form1( _test=XhV, _vector=R,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=rowStartInVector+startBlockIndexVelocity );
    auto linearFormP = form1( _test=XhP, _vector=R,
                              _pattern=size_type(Pattern::COUPLED),
                              _rowstart=rowStartInVector+startBlockIndexPressure );

    auto const& u = mctx.field( FieldTag::velocity(this), "velocity" );
    auto const& p = mctx.field( FieldTag::pressure(this), "pressure" );
    auto const& v = u;
    auto const& q = p;
    auto const& se = mctx.symbolsExpr();

    //--------------------------------------------------------------------------------------------------//

    // identity Matrix
    auto Id = eye<nDim,nDim>();

    //--------------------------------------------------------------------------------------------------//

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto dynamicViscosityLawPtr = std::static_pointer_cast<dynamic_viscosity_law_type>( matProps.law( "dynamic-viscosity" ) );

            // stress tensor sigma : grad(v)
            if ( !timeSteppingEvaluateResidualWithoutTimeDerivative && BuildNonCstPart && !UseJacobianLinearTerms )
            {
                linearFormV +=
                    integrate( _range=range,
                               _expr= -idv(p)*div(v),
                               _geomap=this->geomap() );
            }

            bool doAssemblyStressTensor = ( dynamicViscosityLawPtr->isConstant() && !physicFluidData->turbulence().isEnabled() )? BuildNonCstPart && !UseJacobianLinearTerms : BuildNonCstPart;
            if ( doAssemblyStressTensor )
            {
                auto const StressTensorExpr = Feel::FeelModels::fluidMecStressTensor(u/*gradv(u)*/,idv(p),*physicFluidData,matProps,false/*true*/,se);
                // sigma : grad(v) on Omega
                linearFormV +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( StressTensorExpr,grad(v) ),
                               _geomap=this->geomap() );
            }
            
#if 0
            if ( physicFluidData->turbulence().isEnabled() && BuildNonCstPart )
            {
                auto lmix = min( 0.41*idv(M_fieldDist2Wall), cst(0.09)*cst(0.0635/2.) );
                auto mut = pow(lmix,2)*abs(gradv(u)(0,1));
                linearFormV +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( 2*mut*sym(gradv(u)),grad(v) ),
                               _geomap=this->geomap() );
            }
#endif
            


            // convection
            if ( BuildNonCstPart && physicFluidData->equation() == "Navier-Stokes" )
            {
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                if ( this->useSemiImplicitTimeScheme() )
                {
                    auto const& beta_u = this->useSemiImplicitTimeScheme()? mctx.field( FieldTag::velocity_extrapolated(this), "velocity_extrapolated" ) : u;
                    linearFormV +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*densityExpr*trans( gradv(u)*idv(beta_u) )*id(v),
                                   _geomap=this->geomap() );
                    if ( this->doStabConvectionEnergy() )
                        CHECK( false ) << "TODO";
                }
                else
                {
                    if ( this->doStabConvectionEnergy() )
                    {
                        linearFormV +=
                            integrate( _range=range,
                                       //_expr= /*idv(*M_P0Rho)**/inner( Feel::vf::FSI::fluidMecConvection(u,*M_P0Rho) + idv(*M_P0Rho)*0.5*divv(u)*idv(u), id(v) ),
                                       _expr=timeSteppingScaling*densityExpr*inner( Feel::FeelModels::fluidMecConvectionWithEnergyStab(u), id(v) ),
                                       _geomap=this->geomap() );
                    }
                    else
                    {
                        // convection term
                        // auto convecTerm = val( idv(rho)*trans( gradv(u)*idv(u) ))*id(v);
                        auto convecTerm = densityExpr*inner( Feel::FeelModels::fluidMecConvection(u),id(v) );
                        linearFormV +=
                            integrate( _range=range,
                                       _expr=timeSteppingScaling*convecTerm,
                                       _geomap=this->geomap() );
                    }
                }
            }

#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( this->hasMeshMotion() && !BuildCstPart && !UseJacobianLinearTerms && ( physicFluidData->equation() == "Navier-Stokes" ||  physicFluidData->equation() == "StokesTransient" ) )
            {
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                // mesh velocity (convection) term
                linearFormV +=
                    integrate( _range=range,
                               _expr= -timeSteppingScaling*val(densityExpr*trans( gradv(u)*( idv( this->meshMotionTool()->velocity() ))))*id(v),
                               _geomap=this->geomap() );
            }
#endif
            if ( !BuildCstPart && !UseJacobianLinearTerms && physicFluidData->equation() == "Navier-Stokes" && data.hasVectorInfo( "explicit-part-of-solution" ) )
            {
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                auto uExplicitPartOfSolution = XhV->element( data.vectorInfo( "explicit-part-of-solution" ), rowStartInVector+startBlockIndexVelocity );
                linearFormV +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= timeSteppingScaling*val( densityExpr*trans( gradv(u)*idv(uExplicitPartOfSolution) + gradv(uExplicitPartOfSolution )*idv(u)  ))*id(v),
                               _geomap=this->geomap() );
            }


            //------------------------------------------------------------------------------------//
            //transients terms
            if ( !this->isStationaryModel() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
            {
                bool Build_TransientTerm = !BuildCstPart;
                if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=!BuildCstPart && !UseJacobianLinearTerms;

                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                if (Build_TransientTerm) //  !BuildCstPart && !UseJacobianLinearTerms )
                {
                    linearFormV_PatternDefault +=
                        integrate( _range=range,
                                   _expr= densityExpr*inner(idv(u),id(v))*M_bdfVelocity->polyDerivCoefficient(0),
                                   _geomap=this->geomap() );
                }

                if (BuildCstPart && doAssemblyRhs)
                {
                    auto buzz = M_bdfVelocity->polyDeriv();
                    linearFormV_PatternDefault +=
                        integrate( _range=range,
                                   _expr= -densityExpr*inner(idv(buzz),id(v)),
                                   _geomap=this->geomap() );
                }
            }

            //------------------------------------------------------------------------------------//
            //! gravity force
            if ( doAssemblyRhs && physicFluidData->gravityForceEnabled() && false )
            {
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                auto const& gravityForce = physicFluidData->gravityForceExpr();
                bool assembleGravityTerm = gravityForce.expression().isNumericExpression()? BuildCstPart : BuildNonCstPart;
                if ( assembleGravityTerm )
                {
                    auto const& gravityForceExpr = gravityForce; // TODO apply symbols expr
                    linearFormV +=
                        integrate( _range=M_rangeMeshElements,
                                   _expr= -timeSteppingScaling*densityExpr*inner(gravityForceExpr,id(u)),
                                   _geomap=this->geomap() );
                }
            }

            //--------------------------------------------------------------------------------------------------//
            // take into account that div u != 0
            if (!this->velocityDivIsEqualToZero() && BuildCstPart && doAssemblyRhs )
            {
                linearFormP +=
                    integrate( _range=range,
                               _expr= idv(this->velocityDiv())*id(q),
                               _geomap=this->geomap() );

                auto muExpr = expr( matProps.property("dynamic-viscosity").template expr<1,1>(), se );
                auto coeffDiv = (2./3.)*muExpr;
                linearFormV +=
                    integrate( _range=range,
                               _expr= val(-timeSteppingScaling*coeffDiv*gradv(this->velocityDiv()))*id(v),
                               _geomap=this->geomap() );
            }

            //--------------------------------------------------------------------------------------------------//
            // apply turbulent wall function : wall shear stress bc
            if ( physicFluidData->turbulence().isEnabled() && (physicFluidData->turbulence().model() == "k-epsilon" )  && BuildNonCstPart )
            {
                auto frictionVelocityExpr = physicFluidData->turbulence().frictionVelocityWallFunctionExpr( matName, se );
                auto u_plusExpr = cst(11.06); // TODO
                auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
                //auto expr_turbulence_wst = expr( expr( "physics_fluid_fluid_fluid_Omega_turbulence_k_epsilon_u_tauBC_V2/11.06:physics_fluid_fluid_fluid_Omega_turbulence_k_epsilon_u_tauBC_V2" ) , se );
                for ( auto const& [bcName,bcWall] : M_turbulenceModelBoundaryConditions.wall() )
                {
                    linearFormV +=
                        integrate( _range=markedfaces(this->mesh(), bcWall.markers() /* {"Gamma1","Gamma3"}*/ ),
                                   //_expr= timeSteppingScaling*expr_turbulence_wst*inner( idt(u),id(v) ),
                                   _expr=timeSteppingScaling*densityExpr*(frictionVelocityExpr/u_plusExpr)*inner( idv(u),id(v) ),
                                   _geomap=this->geomap() );
                }
            }

            // if ( true && BuildNonCstPart )
            // {
            //     auto expr_turbulence_wst = expr( expr( "physics_fluid_fluid_fluid_Omega_turbulence_k_epsilon_u_tauBC_V2/11.06:physics_fluid_fluid_fluid_Omega_turbulence_k_epsilon_u_tauBC_V2" ) , se );
            //     linearFormV +=
            //         integrate( _range=markedfaces(this->mesh(),{"Gamma1","Gamma3"} ),
            //                    _expr= timeSteppingScaling*expr_turbulence_wst*inner( idv(u),id(v) ),
            //                    _geomap=this->geomap() );
            // }



            // stabilization gls
            if ( M_stabilizationGLS && M_stabilizationGLSDoAssembly )
            {
                this->updateResidualStabilizationGLS( data, mctx, *physicFluidData, matProps, range );
            }

        } // foreach mat
    } // foreach physic

    //--------------------------------------------------------------------------------------------------//

    // incompressibility term
    if (!BuildCstPart && !UseJacobianLinearTerms )
    {
        linearFormP +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -divv(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    // body forces
#if 0 // VINCENT
    if (BuildCstPart && doAssemblyRhs)
    {
        if ( this->M_overwritemethod_updateSourceTermResidual != NULL )
        {
            this->M_overwritemethod_updateSourceTermResidual(R);
        }
        else
        {
            for( auto const& d : this->M_volumicForcesProperties )
            {
                auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
                linearFormV +=
                    integrate( _range=rangeBodyForceUsed,
                               _expr= -timeSteppingScaling*inner( expression(d,se),id(v) ),
                               _geomap=this->geomap() );
            }
        }

        if (M_haveSourceAdded)
        {
            linearFormV +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -timeSteppingScaling*trans(idv(*M_SourceAdded))*id(v),
                           _geomap=this->geomap() );
        }
    }
#endif

    //------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && !BuildCstPart && !UseJacobianLinearTerms )
        {
            double beta = this->definePressureCstPenalisationBeta();
            for ( auto const& rangeElt : M_definePressureCstMeshRanges )
                linearFormP +=
                    integrate( _range=rangeElt,
                               _expr=beta*idv(p)*id(q),
                               _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("define-pressure-cst-lm") ) << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startBlockIndexDefinePressureCstLM = this->startSubBlockSpaceIndex("define-pressure-cst-lm");

            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                for ( int k=0;k<M_XhMeanPressureLM.size();++k )
                {
                    auto lambda = M_XhMeanPressureLM[k]->element(XVec,rowStartInVector+startBlockIndexDefinePressureCstLM+k);
                    //M_blockVectorSolution.setSubVector( lambda, *XVec, rowStartInVector+startBlockIndexDefinePressureCstLM+k );
                    //for ( size_type k=0;k<M_XhMeanPressureLM->nLocalDofWithGhost();++k )
                    //    lambda( k ) = XVec->operator()( startDofIndexDefinePressureCstLM + k);

                    form1( _test=M_XhMeanPressureLM[k],_vector=R,
                           _rowstart=rowStartInVector+startBlockIndexDefinePressureCstLM+k ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= id(p)*idv(lambda) + idv(p)*id(lambda),
                                   _geomap=this->geomap() );
                }
            }
#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE)
            if ( BuildCstPart && doAssemblyRhs )
            {
                for ( int k=0;k<M_XhMeanPressureLM.size();++k )
                {
                    auto lambda = M_XhMeanPressureLM[k]->element();
                    form1( _test=M_XhMeanPressureLM[k],_vector=R,
                           _rowstart=rowStartInVector+startDofIndexDefinePressureCstLM+k ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= -(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE(this->shared_from_this()))*id(lambda),
                                   _geomap=this->geomap() );
                }
            }
#endif
        }
    }



    //------------------------------------------------------------------------------------//
#if 0
    if ( UsePeriodicity && !BuildCstPart )
    {
        std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
        double pressureJump = doption(_name="periodicity.pressure-jump",_prefix=this->prefix());
        linearForm_PatternCoupled +=
            integrate( _range=markedfaces( this->mesh(),this->mesh()->markerName(marker1) ),
                       _expr=-inner(pressureJump*N(),id(v) ) );
    }
#endif
    //------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------//


    //--------------------------------------------------------------------------------------------------//
    // Neumann boundary condition
    if ( BuildCstPart && doAssemblyRhs )
    {
        for ( auto const& [bcId,bcData] : M_boundaryConditions->normalStress() )
        {
            if ( bcData->isScalarExpr() )
            {
                auto normalStessExpr = bcData->exprScalar( se );
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*normalStessExpr*inner( N(),id(v) ),
                               _geomap=this->geomap() );
            }
            else if ( bcData->isVectorialExpr() )
            {
                auto normalStessExpr = bcData->exprVectorial( se );
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*inner( normalStessExpr,id(v) ),
                               _geomap=this->geomap() );
            }
            else if ( bcData->isMatrixExpr() )
            {
                auto normalStessExpr = bcData->exprMatrix( se );
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*inner( normalStessExpr*N(),id(v) ),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
    int indexOutletWindkessel = 0;
    for ( auto const& [bcName,bcData] : M_boundaryConditions->outletWindkessel() )
    {
        if ( bcData->useImplicitCoupling() )
        {
            if ( BuildCstPart || (BuildNonCstPart && !UseJacobianLinearTerms) )
            {
                size_type startBlockIndexWindkessel = this->startSubBlockSpaceIndex("windkessel");

                bool hasWindkesselActiveDof = M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0;
                int blockStartWindkesselVec = rowStartInVector + startBlockIndexWindkessel + 2*indexOutletWindkessel;
                auto const& basisToContainerGpPressureDistalVec = R->map().dofIdToContainerId( blockStartWindkesselVec );
                auto const& basisToContainerGpPressureProximalVec = R->map().dofIdToContainerId( blockStartWindkesselVec+1 );
                if ( hasWindkesselActiveDof )
                    CHECK( !basisToContainerGpPressureDistalVec.empty() && !basisToContainerGpPressureProximalVec.empty() ) << "incomplete datamap info";
                const size_type gpPressureDistalVec = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalVec[0] : 0;
                const size_type gpPressureProximalVec = (hasWindkesselActiveDof)? basisToContainerGpPressureProximalVec[0] : 0;

                // windkessel parameters
                double Rd = bcData->expr_Rd( se ).evaluate()(0,0);
                double Rp = bcData->expr_Rp( se ).evaluate()(0,0);
                double Cd = bcData->expr_Cd( se ).evaluate()(0,0);

                auto const& windkesselData = M_fluidOutletWindkesselData.at(bcName);
                auto bdfDistal = std::get<1>( windkesselData ).at(0);

                R->setIsClosed( false );
                if ( BuildCstPart && hasWindkesselActiveDof && doAssemblyRhs )
                {
                    auto rhsTimeDerivDistal = bdfDistal->polyDeriv();
                    R->add( gpPressureDistalVec, -Cd*rhsTimeDerivDistal(0) );
                }

                if ( BuildNonCstPart && !UseJacobianLinearTerms )
                {
                    auto rangeFaceFluidOutlet = markedfaces(mesh,bcData->markers());

                    // TODO : not correct, use view here
                    auto pDistal = std::get<0>( windkesselData ).at(0);
                    auto pProximal = std::get<0>( windkesselData ).at(1);

                    // first equation
                    if ( hasWindkesselActiveDof )
                    {
                        R->add( gpPressureDistalVec,
                                unwrap_ptr(pDistal)(0)* Cd*bdfDistal->polyDerivCoefficient(0)+ 1./Rd );
                    }
                    form1( _test=M_fluidOutletWindkesselSpace,_vector=R,
                           _rowstart=blockStartWindkesselVec ) +=
                        integrate( _range=rangeFaceFluidOutlet,
                                   _expr=-inner(idv(u),N())*id(pDistal),
                                   _geomap=this->geomap() );

                    // second equation
                    if ( hasWindkesselActiveDof )
                    {
                        R->add( gpPressureProximalVec, unwrap_ptr(pProximal)(0) - unwrap_ptr(pDistal)(0) );
                    }

                    // really correct?
                    form1( _test=M_fluidOutletWindkesselSpace,_vector=R,
                           _rowstart=blockStartWindkesselVec+1 )+=
                        integrate( _range=rangeFaceFluidOutlet,
                                   _expr=-Rp*inner(idv(u),N())*id(pProximal),
                                   _geomap=this->geomap() );

                    // coupling with fluid model
                    form1( _test=XhV, _vector=R,
                           _rowstart=rowStartInVector ) +=
                        integrate( _range=rangeFaceFluidOutlet,
                                   _expr= timeSteppingScaling*idv(pProximal)*inner(N(),id(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
        else
        {
            if ( BuildCstPart && doAssemblyRhs )
            {
                auto const& beta = M_bdfVelocity->poly();
                // windkessel parameters
                double Rd = bcData->expr_Rd( se ).evaluate()(0,0);
                double Rp = bcData->expr_Rp( se ).evaluate()(0,0);
                double Cd = bcData->expr_Cd( se ).evaluate()(0,0);

                auto const& windkesselData = M_fluidOutletWindkesselData.at(bcName);
                auto pDistal = std::get<0>( windkesselData ).at(0);
                auto pProximal = std::get<0>( windkesselData ).at(1);
                auto bdfDistal = std::get<1>( windkesselData ).at(0);

                auto rangeFaceFluidOutlet = markedfaces(mesh,bcData->markers());

                auto outletQ = integrate(_range=rangeFaceFluidOutlet,
                                         _expr=trans(idv(beta))*N() ).evaluate()(0,0);


                auto rhsTimeDerivDistal = bdfDistal->polyDeriv();

                double denum = Rd*Cd*bdfDistal->polyDerivCoefficient(0) + 1;
                pDistal->setConstant( Rd*outletQ/denum );
                pDistal->add( Rd*Cd/denum, rhsTimeDerivDistal );
                pProximal->setConstant( Rp*outletQ );
                pProximal->add( 1., *pDistal );

                linearFormV +=
                    integrate( _range=rangeFaceFluidOutlet,
                               _expr= timeSteppingScaling*idv(pProximal)*inner(N(),id(v)),
                               _geomap=this->geomap() );
            }
        }
        ++indexOutletWindkessel;
    }

#if 0 // VINCENT
    //--------------------------------------------------------------------------------------------------//

    if (!BuildCstPart && !UseJacobianLinearTerms && !this->markerSlipBC().empty() )
    {
        auto muExpr = this->dynamicViscosityExpr( u,se );
        auto densityExpr = this->materialsProperties()->template materialPropertyExpr<1,1>( "density", se );
        auto P = Id-N()*trans(N());
        double gammaN = doption(_name="bc-slip-gammaN",_prefix=this->prefix());
        double gammaTau = doption(_name="bc-slip-gammaTau",_prefix=this->prefix());
        auto betaExtrapolated = M_bdfVelocity->poly();
        auto beta = vf::project( _space=XhV,
                                 _range=boundaryfaces(mesh),
                                 _expr=densityExpr*idv(betaExtrapolated) );
        auto Cn = gammaN*max(abs(trans(idv(beta))*N()),muExpr/vf::h());
        auto Ctau = gammaTau*muExpr/vf::h() + max( -trans(idv(beta))*N(),cst(0.) );

        linearFormV +=
            integrate( _range=markedfaces(mesh,this->markerSlipBC()),
                       _expr=
                       val(timeSteppingScaling*Cn*(trans(idv(u))*N()))*(trans(id(v))*N())+
                       val(timeSteppingScaling*Ctau*trans(idv(u)))*id(v),
                       //+ trans(idv(p)*Id*N())*id(v)
                       //- trans(id(v))*N()* trans(2*idv(mu)*defv*N())*N()
                       _geomap=this->geomap()
                       );
    }
#endif // VINCENT
    //--------------------------------------------------------------------------------------------------//

    // weak formulation of the boundaries conditions
    if ( M_boundaryConditions->hasVelocityImposedNitsche() )
    {
        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            std::set<std::string> allmarkers;
            for ( auto const& [bcId,bcData] : M_boundaryConditions->velocityImposedNitsche() )
                allmarkers.insert( bcData->markers().begin(), bcData->markers().end() );

            auto Sigmav = this->stressTensorExpr(u,p,se);
            linearFormV +=
                integrate( _range=markedfaces(mesh,allmarkers),
                           _expr= -timeSteppingScaling*trans(Sigmav*N())*id(v)
                           /**/   + timeSteppingScaling*this->dirichletBCnitscheGamma()*inner( idv(u),id(v) )/hFace(),
                           _geomap=this->geomap() );
        }
        if ( BuildCstPart && doAssemblyRhs )
        {
            for ( auto const& [bcId,bcData] : M_boundaryConditions->velocityImposedNitsche() )
            {
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*this->dirichletBCnitscheGamma()*inner( bcData->expr(se),id(v) )/hFace(),
                               _geomap=this->geomap() );
            }
        }
    }

    //------------------------------------------------------------------------------------//
    // Dirichlet with Lagrange-mulitplier
    if ( M_boundaryConditions->hasVelocityImposedLagrangeMultiplier() )
    {
        CHECK( this->hasStartSubBlockSpaceIndex("dirichletlm") ) << " start dof index for dirichletlm is not present\n";
        size_type startBlockIndexDirichletLM = this->startSubBlockSpaceIndex("dirichletlm");

        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            std::set<std::string> allmarkers;
            for ( auto const& [bcId,bcData] : M_boundaryConditions->velocityImposedLagrangeMultiplier() )
                allmarkers.insert( bcData->markers().begin(), bcData->markers().end() );

            auto lambdaBC = this->XhDirichletLM()->element( XVec, rowStartInVector+startBlockIndexDirichletLM );
            //int dataBaseIdLM = XVec->map().basisIndexFromGp( rowStartInVector+startBlockIndexDirichletLM );
            //this->algebraicBlockVectorSolution()->setSubVector( lambdaBC, *XVec, rowStartInVector+startBlockIndexDirichletLM );

            linearFormV +=
                integrate( _range=markedfaces(mesh,allmarkers),
                           _expr= inner( idv(lambdaBC),id(u) ) );

            form1( _test=this->XhDirichletLM(),_vector=R,
                   _rowstart=rowStartInVector+startBlockIndexDirichletLM ) +=
                integrate( //_range=elements(this->meshDirichletLM()),
                    _range=markedfaces( this->mesh(),allmarkers),
                           _expr= inner(idv(u),id(lambdaBC) ) );
        }
#if 1
        if ( BuildCstPart && doAssemblyRhs )
        {
            auto lambdaBC = this->XhDirichletLM()->element();
            for ( auto const& [bcId,bcData] : M_boundaryConditions->velocityImposedLagrangeMultiplier() )
            {
                form1( _test=this->XhDirichletLM(),_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexDirichletLM ) +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               //_range=markedelements(this->meshDirichletLM(),PhysicalName),
                               _expr= -inner( bcData->expr(se),id(lambdaBC) ),
                               _geomap=this->geomap() );
            }
        }
#endif

    }

    //------------------------------------------------------------------------------------//

    if ( !M_boundaryConditions->pressureImposed().empty() )
    {

        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            std::set<std::string> allmarkers;
            for ( auto const& [bcName,bcData] : M_boundaryConditions->pressureImposed() )
                allmarkers.insert( bcData->markers().begin(), bcData->markers().end() );
            auto rangeFacesPressureBC = markedfaces( this->mesh(),allmarkers );

            CHECK( this->hasStartSubBlockSpaceIndex("pressurelm1") ) << " start dof index for pressurelm1 is not present\n";
            size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");

            auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( XVec, rowStartInVector+startBlockIndexPressureLM1 );

            if ( nDim==2 )
            {
                linearFormV +=
                    integrate( _range=rangeFacesPressureBC,
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,0)*idv(lambdaPressure1),
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM1 ) +=
                    integrate( _range=rangeFacesPressureBC,
                               _expr=-trans(cross(idv(u),N()))(0,0)*id(lambdaPressure1),
                               _geomap=this->geomap() );
            }
            else if ( nDim==3 )
            {
                auto alpha = 1./sqrt(1-Nz()*Nz());
                linearFormV +=
                    integrate( _range=rangeFacesPressureBC,
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,2)*idv(lambdaPressure1)*alpha,
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM1 ) +=
                    integrate( _range=rangeFacesPressureBC,
                               _expr=-trans(cross(idv(u),N()))(0,2)*id(lambdaPressure1)*alpha,
                               _geomap=this->geomap() );

                CHECK( this->hasStartSubBlockSpaceIndex("pressurelm2") ) << " start dof index for pressurelm2 is not present\n";
                size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");

                auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( XVec, rowStartInVector+startBlockIndexPressureLM2 );

                linearFormV +=
                    integrate( _range=rangeFacesPressureBC,
                               _expr= -timeSteppingScaling*trans(cross(id(u),N()))(0,0)*alpha*idv(lambdaPressure2)*Ny()
                               +timeSteppingScaling*trans(cross(id(u),N()))(0,1)*alpha*idv(lambdaPressure2)*Nx(),
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM2 ) +=
                    integrate( _range=rangeFacesPressureBC,
                               _expr= -trans(cross(idv(u),N()))(0,0)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +trans(cross(idv(u),N()))(0,1)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );
            }
        }
        if ( BuildCstPart && doAssemblyRhs )
        {
            for ( auto const& [bcName,bcData] : M_boundaryConditions->pressureImposed() )
            {
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= timeSteppingScaling*bcData->expr(se)*trans(N())*id(v),
                               _geomap=this->geomap() );
            }
        }

    }

    //--------------------------------------------------------------------------------------------------//

    if ( !M_bodySetBC.empty() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
    {
        this->log("FluidMechanics","updateJacobianWeakBC","assembly of body bc");

        for ( auto const& [bpname,bpbc] : M_bodySetBC )
        {
            size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
            bool hasActiveDofTranslationalVelocity = bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0;
            double massBody = bpbc.body().mass();
            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                R->setIsClosed( false );
                if ( hasActiveDofTranslationalVelocity)
                {
                    auto uTranslationalVelocity = bpbc.spaceTranslationalVelocity()->element( XVec, rowStartInVector+startBlockIndexTranslationalVelocity );
                    auto const& basisToContainerGpTranslationalVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    for (int d=0;d<nDim;++d)
                    {
                        R->add( basisToContainerGpTranslationalVelocityVector[d],
                                uTranslationalVelocity(d)*bpbc.bdfTranslationalVelocity()->polyDerivCoefficient(0)*massBody );
                    }
                }
            }
            if ( BuildCstPart && doAssemblyRhs )
            {
                R->setIsClosed( false );
                if ( hasActiveDofTranslationalVelocity )
                {
                    auto const& basisToContainerGpTranslationalVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    auto translationalVelocityPolyDeriv = bpbc.bdfTranslationalVelocity()->polyDeriv();
                    for (int d=0;d<nDim;++d)
                    {
                        R->add( basisToContainerGpTranslationalVelocityVector[d],
                                -massBody*translationalVelocityPolyDeriv(d) );

                        if ( bpbc.gravityForceEnabled() )
                        {
                            R->add( basisToContainerGpTranslationalVelocityVector[d],
                                    -bpbc.gravityForceWithMass()(d) );
                        }
                    }
                }
            }

            if ( !bpbc.isInNBodyArticulated() || ( bpbc.getNBodyArticulated().masterBodyBC().name() == bpbc.name() ) )
            {
                size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
                int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();
                bool hasActiveDofAngularVelocity = nLocalDofAngularVelocity > 0;
                R->setIsClosed( false );

                if ( hasActiveDofAngularVelocity )
                {
                    auto uAngularVelocity = bpbc.spaceAngularVelocity()->element( XVec, rowStartInVector+startBlockIndexAngularVelocity );

                    auto w_eval = idv(uAngularVelocity).evaluate(false);
                    auto const& basisToContainerGpAngularVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexAngularVelocity );
                    auto const& momentOfInertia = bpbc.momentOfInertia_inertialFrame();
                    //auto timeDerivativeOfMomentOfInertia = bpbc.timeDerivativeOfMomentOfInertia_bodyFrame(this->timeStep());
                    if ( !BuildCstPart && !UseJacobianLinearTerms )
                    {
                        typename Body::moment_of_inertia_type termWithTimeDerivativeOfMomentOfInertia;
                        if constexpr ( nDim == 2 )
                            termWithTimeDerivativeOfMomentOfInertia = bpbc.timeDerivativeOfMomentOfInertia_bodyFrame(this->timeStep());
                        else
                        {
                            auto rotationMat = bpbc.rigidRotationMatrix();
                            termWithTimeDerivativeOfMomentOfInertia = rotationMat*bpbc.timeDerivativeOfMomentOfInertia_bodyFrame(this->timeStep())*(rotationMat.transpose());
                        }

                        //auto contribLhsAngularVelocity = (bpbc.bdfAngularVelocity()->polyDerivCoefficient(0)*((momentOfInertiaExpr*idv(uAngularVelocity)).evaluate(false))).eval(); // not works if not call eval(), probably aliasing issue
                        auto contribLhsAngularVelocity = (bpbc.bdfAngularVelocity()->polyDerivCoefficient(0)*momentOfInertia + termWithTimeDerivativeOfMomentOfInertia)*w_eval;
                        for (int i=0;i<nLocalDofAngularVelocity;++i)
                        {
                            double val = contribLhsAngularVelocity(i,0);
                            R->add( basisToContainerGpAngularVelocityVector[i], val );
                        }
                    }
                    if ( BuildNonCstPart )
                    {
                        if constexpr ( nDim == 3 )
                        {
                            // term : w x Iw
                            auto Iw_eval = momentOfInertia*w_eval;
                            auto wCrossIw = w_eval.cross( Iw_eval );
                            for (int i=0;i<nLocalDofAngularVelocity;++i)
                                R->add( basisToContainerGpAngularVelocityVector[i], wCrossIw(i) );
                        }
                    }
                    if ( BuildCstPart && doAssemblyRhs )
                    {
                        //auto const& basisToContainerGpAngularVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexAngularVelocity );
                        auto angularVelocityPolyDeriv = bpbc.bdfAngularVelocity()->polyDeriv();
                        //auto momentOfInertiaExpr = bpbc.momentOfInertiaExpr();
                        auto contribRhsAngularVelocity = momentOfInertia*(idv(angularVelocityPolyDeriv).evaluate(false));
                        for (int i=0;i<nLocalDofAngularVelocity;++i)
                        {
                            R->add( basisToContainerGpAngularVelocityVector[i],
                                    -contribRhsAngularVelocity(i,0)
                                    //momentOfInertia*angularVelocityPolyDeriv(d)
                                    );
                        }
                    }
                } // hasActiveDofAngularVelocity
            }

        } // for ( auto const& [bpname,bpbc] : M_bodySetBC )

        for ( auto const& nba : M_bodySetBC.nbodyArticulated() )
        {
            if ( nba.articulationMethod() != "lm" )
                continue;

            for ( auto const& ba : nba.articulations() )
            {
                auto const& bbc1 = ba.body1();
                auto const& bbc2 = ba.body2();
                size_type startBlockIndexTranslationalVelocityBody1 = this->startSubBlockSpaceIndex("body-bc."+bbc1.name()+".translational-velocity");
                size_type startBlockIndexTranslationalVelocityBody2 = this->startSubBlockSpaceIndex("body-bc."+bbc2.name()+".translational-velocity");
                bool hasActiveDofTranslationalVelocityBody1 = bbc1.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0;
                bool hasActiveDofTranslationalVelocityBody2 = bbc2.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0;
                size_type startBlockIndexArticulationLMTranslationalVelocity = this->startSubBlockSpaceIndex( "body-bc.articulation-lm."+ba.name()+".translational-velocity");
                if ( !BuildCstPart && !UseJacobianLinearTerms )
                {
                    R->setIsClosed( false );
                    if ( hasActiveDofTranslationalVelocityBody1 )
                    {
                        auto const& basisToContainerGpTranslationalVelocityBody1 = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocityBody1 );
                        auto const& basisToContainerGpArticulationLMTranslationalVelocity = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexArticulationLMTranslationalVelocity );
                        for (int d=0;d<nDim;++d)
                        {
                            R->add( basisToContainerGpTranslationalVelocityBody1[d], XVec->operator()( basisToContainerGpArticulationLMTranslationalVelocity[d] ) );
                            R->add( basisToContainerGpArticulationLMTranslationalVelocity[d], XVec->operator()( basisToContainerGpTranslationalVelocityBody1[d] ) );
                        }
                    }
                    if ( hasActiveDofTranslationalVelocityBody2 )
                    {
                        auto const& basisToContainerGpTranslationalVelocityBody2 = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocityBody2 );
                        auto const& basisToContainerGpArticulationLMTranslationalVelocity = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexArticulationLMTranslationalVelocity );
                        for (int d=0;d<nDim;++d)
                        {
                            R->add( basisToContainerGpArticulationLMTranslationalVelocity[d], -XVec->operator()( basisToContainerGpTranslationalVelocityBody2[d] ) );
                            R->add( basisToContainerGpTranslationalVelocityBody2[d], -XVec->operator()( basisToContainerGpArticulationLMTranslationalVelocity[d] ) );
                        }
                    }
                }
                if ( BuildCstPart && doAssemblyRhs )
                {
                    R->setIsClosed( false );
                    if ( hasActiveDofTranslationalVelocityBody1 )
                    {
                        auto const& basisToContainerGpArticulationLMTranslationalVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexArticulationLMTranslationalVelocity );
                        auto articulationTranslationalVelocityExpr = ba.translationalVelocityExpr( se ).evaluate(false);
                        for (int d=0;d<nDim;++d)
                        {
                            R->add( basisToContainerGpArticulationLMTranslationalVelocityVector[d],
                                    -articulationTranslationalVelocityExpr(d) );
                        }
                    }
                }
            } // ba
        } // nba

    }

    //------------------------------------------------------------------------------------//

    this->updateResidualStabilisation( data, unwrap_ptr(u), unwrap_ptr(p) );

    //------------------------------------------------------------------------------------//

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("FluidMechanics","updateResidual","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

} // namespace Feel
} // namespace FeelModels

#endif
