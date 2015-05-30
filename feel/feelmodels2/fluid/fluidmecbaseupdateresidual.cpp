/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */


#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels2/modelvf/fluidmecconvection.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidual( const vector_ptrtype& XVec, vector_ptrtype& R, bool BuildCstPart,
                                                        bool UseJacobianLinearTerms,
                                                        bool _doClose, bool _doBCStrongDirichlet) const
{
#if defined(FEELMODELS_FLUID_BUILD_RESIDUAL_CODE)
    using namespace Feel::vf;

    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    if (this->verbose()) Feel::FeelModels::Log("--------------------------------------------------\n",
                                               this->prefix()+".FluidMechanics","updateResidual","start"+sc,
                                               this->worldComm(),this->verboseAllProc());

    boost::mpi::timer thetimer, thetimerBis;

    bool BuildNonCstPart = !BuildCstPart;

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    //auto const rowStartInMatrix = this->rowStartInMatrix();
    //auto const colStartInMatrix = this->colStartInMatrix();
    size_type rowStartInVector = this->rowStartInVector();
    auto linearForm_PatternDefault = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::DEFAULT),
                                            _rowstart=rowStartInVector );
    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=rowStartInVector );

    auto U = Xh->element("u");//U = *XVec;
    // copy vector values in fluid element
    for ( size_type k=0;k<Xh->nLocalDofWithGhost();++k )
        U(k) = XVec->operator()(/*rowStartInVector+*/k);

    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    //--------------------------------------------------------------------------------------------------//

    // strain tensor (trial)
    auto defv = sym(gradv(u));
    // identity Matrix
    auto Id = eye<nDim,nDim>();
    // dynamic viscosity
    auto const& mu = this->densityViscosityModel()->fieldMu();
    auto const& rho = this->densityViscosityModel()->fieldRho();
    // stress tensor (eval)
    auto Sigmav = -idv(p)*Id + 2*idv(mu)*defv;

    double timeElapsedBis=thetimerBis.elapsed();
    this->log("FluidMechanics","updateResidual","init done in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//

    thetimerBis.restart();
    if (!BuildCstPart)
    {
        if ( this->doStabConvectionEnergy() )
        {
            linearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           //_expr= /*idv(*M_P0Rho)**/inner( Feel::vf::FSI::fluidMecConvection(u,*M_P0Rho) + idv(*M_P0Rho)*0.5*divv(u)*idv(u), id(v) ),
                           _expr=inner( Feel::vf::FeelModels::fluidMecConvectionWithEnergyStab(u,rho), id(v) ),
                           _geomap=this->geomap() );

            /*if (this->isMoveDomain()  && !BuildCstPart && !UseJacobianLinearTerms)
             {
             linearForm_PatternCoupled +=
             integrate( _range=elements(mesh),
             _expr= -0.5*idv(M_P0Rho)*divv(this->meshVelocity())*trans(idv(u))*id(v),
             _geomap=this->geomap() );
             }*/
        }
        else
        {
            // convection term
#if 0
            auto convecTerm = val( idv(rho)*trans( gradv(u)*idv(u) ))*id(v);
#else
            auto convecTerm = inner( Feel::vf::FeelModels::fluidMecConvection(u,rho),id(v) );
#endif
            linearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr=convecTerm,
                           _geomap=this->geomap() );
        }

    }

    timeElapsedBis=thetimerBis.elapsed();thetimerBis.restart();
    this->log("FluidMechanics","updateResidual","build convective--1-- term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (M_isMoveDomain && !BuildCstPart && !UseJacobianLinearTerms )
    {
        // mesh velocity (convection) term
        linearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= -val(idv(rho)*trans( gradv(u)*( idv( this->meshVelocity() ))))*id(v),
                       _geomap=this->geomap() );
    }
#endif
    timeElapsedBis=thetimerBis.elapsed();
    this->log("FluidMechanics","updateResidual","build convective--2-- term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//

    this->updateResidualModel(U/*XVec*/, R, BuildCstPart, UseJacobianLinearTerms);

    //--------------------------------------------------------------------------------------------------//

    // incompressibility term
    if (!BuildCstPart && !UseJacobianLinearTerms )
    {
        linearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= divv(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( BuildCstPart )
    {
        //neumann boundary condition
        this->updateBCNeumannResidual( R );
        //pressure boundary condition
        this->updateBCPressureResidual( R );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" )
    {
        if ( this->fluidOutletWindkesselCoupling() == "explicit" )
        {
            if ( BuildCstPart )
            {
                auto const Beta = M_bdf_fluid->poly();

                for (int k=0;k<this->nFluidOutlet();++k)
                {
                    // Windkessel model
                    double Rd=this->fluidOutletWindkesselCoeffRd()[k];// 6.2e3;
                    double Rp=this->fluidOutletWindkesselCoeffRp()[k];//400;
                    double Cd=this->fluidOutletWindkesselCoeffCd()[k];//2.72e-4;
                    double Deltat = this->timeStepBDF()->timeStep();

                    double xiBF = Rd*Cd*this->timeStepBDF()->polyDerivCoefficient(0)*Deltat+Deltat;
                    double alphaBF = Rd*Cd/(xiBF);
                    double gammaBF = Rd*Deltat/xiBF;
                    double kappaBF = (Rp*xiBF+ Rd*Deltat)/xiBF;

                    std::string markerOutlet=this->fluidOutletMarkerName(k);
                    auto outletQ = integrate(_range=markedfaces(mesh,markerOutlet),
                                             _expr=trans(idv(Beta.template element<0>()))*N() ).evaluate()(0,0);

                    double pressureDistalOld  = 0;
                    for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                        pressureDistalOld += Deltat*this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld()[k][i];

                    M_fluidOutletWindkesselPressureDistal[k] = alphaBF*pressureDistalOld + gammaBF*outletQ;
                    M_fluidOutletWindkesselPressureProximal[k] = kappaBF*outletQ + alphaBF*pressureDistalOld;


                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= M_fluidOutletWindkesselPressureProximal[k]*trans(N())*id(v),
                                   _geomap=this->geomap() );
                }
            }
        } // explicit
        else if ( this->fluidOutletWindkesselCoupling() == "implicit" )
        {
            CHECK( this->startDofIndexFieldsInMatrix().find("windkessel") != this->startDofIndexFieldsInMatrix().end() )
                << " start dof index for windkessel is not present\n";
            size_type startDofIndexWindkessel = this->startDofIndexFieldsInMatrix().find("windkessel")->second;

            if ( BuildCstPart || (!BuildCstPart && !UseJacobianLinearTerms) )
            {
                auto presDistalProximal = M_fluidOutletWindkesselSpace->element();
                auto presDistal = presDistalProximal.template element<0>();
                auto presProximal = presDistalProximal.template element<1>();

                for (int k=0;k<this->nFluidOutlet();++k)
                {
                    // Windkessel model
                    double Rd=this->fluidOutletWindkesselCoeffRd()[k];// 6.2e3;
                    double Rp=this->fluidOutletWindkesselCoeffRp()[k];//400;
                    double Cd=this->fluidOutletWindkesselCoeffCd()[k];//2.72e-4;
                    double Deltat = this->timeStepBDF()->timeStep();
                    std::string markerOutlet = this->fluidOutletMarkerName(k);

                    const size_type rowStartWindkessel = startDofIndexWindkessel + 2*k;
                    const size_type rowStartInVectorWindkessel = rowStartInVector + rowStartWindkessel;
                    //----------------------------------------------------//
                    if ( BuildCstPart  && M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                    {
                        double pressureDistalOld  = 0;
                        for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                            pressureDistalOld += this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld()[k][i];
                        // add in vector
                        R->add( rowStartInVectorWindkessel, -Cd*pressureDistalOld);
                    }
                    //----------------------------------------------------//
                    if ( !BuildCstPart && !UseJacobianLinearTerms )
                    {
                        for ( size_type kk=0;kk<M_fluidOutletWindkesselSpace->nLocalDofWithGhost();++kk )
                            presDistalProximal( kk ) = XVec->operator()( rowStartWindkessel + kk);
                        //----------------------------------------------------//
                        // 1ere ligne
                        if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                        {
                            const double value = presDistalProximal(0)*(Cd*this->timeStepBDF()->polyDerivCoefficient(0)+1./Rd);
                            R->add( rowStartInVectorWindkessel,value );
                        }
                        form1( _test=M_fluidOutletWindkesselSpace,_vector=R,
                               _rowstart=rowStartInVectorWindkessel ) +=
                            integrate( _range=markedfaces(mesh,markerOutlet),
                                       _expr=-(trans(idv(u))*N())*id(presDistal),
                                       _geomap=this->geomap() );
                        //----------------------------------------------------//
                        // 2eme ligne
                        if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                        {
                            R->add( rowStartInVectorWindkessel+1,  presDistalProximal(1)-presDistalProximal(0) );
                        }
                        form1( _test=M_fluidOutletWindkesselSpace,_vector=R,
                               _rowstart=rowStartInVectorWindkessel )+=
                            integrate( _range=markedfaces(mesh,markerOutlet),
                                       _expr=-Rp*(trans(idv(u))*N())*id(presProximal),
                                       _geomap=this->geomap() );
                        //----------------------------------------------------//
                        // coupling with fluid model
                        form1( _test=Xh, _vector=R,
                               _rowstart=rowStartInVector ) +=
                            integrate( _range=markedfaces(mesh,markerOutlet),
                                       _expr= idv(presProximal)*trans(N())*id(v),
                                       _geomap=this->geomap() );
                    }
                }
            }

        } // implicit

    }

    //--------------------------------------------------------------------------------------------------//

    if (!BuildCstPart && !UseJacobianLinearTerms && !this->markerSlipBC().empty() )
    {
        // slip condition :
        auto P = Id-N()*trans(N());
        double gammaN = doption(_name="bc-slip-gammaN",_prefix=this->prefix());
        double gammaTau = doption(_name="bc-slip-gammaTau",_prefix=this->prefix());
        auto Beta = M_bdf_fluid->poly();
        auto beta = vf::project( _space=Beta.template element<0>().functionSpace(),
                                 _range=boundaryfaces(Beta.template element<0>().mesh()),
                                 _expr=idv(rho)*idv(Beta.template element<0>()) );
        auto Cn = gammaN*max(abs(trans(idv(beta))*N()),idv(mu)/vf::h());
        auto Ctau = gammaTau*idv(mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) );

        linearForm_PatternCoupled +=
            integrate( _range=markedfaces(mesh,this->markerSlipBC()),
                       _expr= val(Cn*(trans(idv(u))*N()))*(trans(id(v))*N())+
                       val(Ctau*trans(idv(u)))*id(v),
                       //+ trans(idv(p)*Id*N())*id(v)
                       //- trans(id(v))*N()* trans(2*idv(mu)*defv*N())*N()
                       _geomap=this->geomap()
                       );
    }

    //--------------------------------------------------------------------------------------------------//

    // weak formulation of the boundaries conditions
    if ( this->hasMarkerDirichletBCnitsche() )
    {
        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            linearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerDirichletBCnitsche() ),
                           _expr= -trans(Sigmav*N())*id(v)
                           /**/   + this->dirichletBCnitscheGamma()*inner( idv(u),id(v) )/hFace(),
                           _geomap=this->geomap() );
        }

        if ( BuildCstPart )
        {
            this->updateBCDirichletNitscheResidual( R );

#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet" )
            {
                // compute integrate range (intersection with nitsche and moving marker)
                std::list<std::string> movingBCmarkers;
                for ( std::string const& marker : this->markersNameMovingBoundary() )
                {
                    auto it = std::find_if( this->markerDirichletBCnitsche().begin(), this->markerDirichletBCnitsche().end(),
                                            [&marker]( std::string const& m ) { return m == marker; } );
                    if ( it != this->markerDirichletBCnitsche().end() )
                        movingBCmarkers.push_back( marker );
                }

                linearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,movingBCmarkers),
                               _expr= - this->dirichletBCnitscheGamma()*inner( idv(this->meshVelocity2()),id(v) )/hFace(),
                               _geomap=this->geomap() );
            }
#endif
        }
    }

    //--------------------------------------------------------------------------------------------------//

    // volume force
    if (BuildCstPart)
    {
        this->updateSourceTermResidual( R );
        /*#if defined(FLUIDMECHANICS_VOLUME_FORCE)
         linearForm_PatternCoupled +=
         integrate( _range=elements(mesh),
         _expr= -trans(f)*id(v),
         _geomap=this->geomap() );
         #endif*/
        if (M_haveSourceAdded)
        {
            linearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= -trans(idv(*M_SourceAdded))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //------------------------------------------------------------------------------------//

    //transients terms
    if (!this->isStationary())
    {
        bool Build_TransientTerm = !BuildCstPart;
        if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=!BuildCstPart && !UseJacobianLinearTerms;

        if (Build_TransientTerm) //  !BuildCstPart && !UseJacobianLinearTerms )
        {
            linearForm_PatternDefault +=
                integrate( _range=elements(mesh),
                           _expr= val(idv(rho)*trans(idv(u))*M_bdf_fluid->polyDerivCoefficient(0))*id(v),
                           _geomap=this->geomap() );
        }

        if (BuildCstPart)
        {
            auto Buzz = M_bdf_fluid->polyDeriv();
            auto buzz = Buzz.template element<0>();
            linearForm_PatternDefault +=
                integrate( _range=elements(mesh),
                           _expr= val(-idv(rho)*trans(idv(buzz)))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && !BuildCstPart && !UseJacobianLinearTerms )
        {
            double beta = this->definePressureCstPenalisationBeta();
            linearForm_PatternCoupled +=
                integrate( _range=elements(Xh->mesh()),
                           _expr=beta*idv(p)*id(q),
                           _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" )
        {
            CHECK( this->startDofIndexFieldsInMatrix().find("define-pressure-cst-lm") != this->startDofIndexFieldsInMatrix().end() )
                << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startDofIndexDefinePressureCstLM = this->startDofIndexFieldsInMatrix().find("define-pressure-cst-lm")->second;

#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_ONLY_ON_BOUNDARY)
            auto therange = boundaryfaces(mesh);
#else
            auto therange = elements(mesh);
#endif
            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                auto lambda = M_XhMeanPressureLM->element();
                for ( size_type k=0;k<M_XhMeanPressureLM->nLocalDofWithGhost();++k )
                    lambda( k ) = XVec->operator()( startDofIndexDefinePressureCstLM + k);

                form1( _test=M_XhMeanPressureLM,_vector=R,
                       _rowstart=rowStartInVector+startDofIndexDefinePressureCstLM ) +=
                    integrate( _range=therange,//elements(mesh),
                               _expr= id(p)*idv(lambda) + idv(p)*id(lambda),
                               _geomap=this->geomap() );
            }
#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE)
            if (BuildCstPart)
            {
                auto lambda = M_XhMeanPressureLM->element();
                form1( _test=M_XhMeanPressureLM,_vector=R,
                       _rowstart=rowStartInVector+startDofIndexDefinePressureCstLM ) +=
                    integrate( _range=elements(mesh),
                               _expr= -(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE(this->shared_from_this()))*id(lambda),
                               _geomap=this->geomap() );
            }
#endif
        }
    }


    //------------------------------------------------------------------------------------//

    this->updateResidualStabilisation(U/*XVec*/, R, BuildCstPart, UseJacobianLinearTerms);

    //------------------------------------------------------------------------------------//
    if (this->hasMarkerDirichletBClm())
    {
        CHECK( this->startDofIndexFieldsInMatrix().find("dirichletlm") != this->startDofIndexFieldsInMatrix().end() )
            << " start dof index for dirichletlm is not present\n";
        size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;

        auto lambdaBC = this->XhDirichletLM()->element();
        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            //size_type rowStartDirichletLM = startDofIndexDirichletLM;
            for ( size_type kk=0;kk<this->XhDirichletLM()->nLocalDofWithGhost();++kk )
                lambdaBC( kk ) = XVec->operator()( startDofIndexDirichletLM + kk);

            form1( _test=Xh,_vector=R,
                   _rowstart=rowStartInVector )+=
                integrate( _range=markedfaces(mesh,this->markerDirichletBClm() ), //elements(this->meshDirichletLM()),
                           _expr= inner( idv(lambdaBC),id(u) ) );

            form1( _test=this->XhDirichletLM(),_vector=R,
                   _rowstart=rowStartInVector+startDofIndexDirichletLM ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner(idv(u),id(lambdaBC) ) );
        }

        if ( BuildCstPart )
        {
            this->updateBCDirichletLagMultResidual( R );

#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet" )
            {
                std::list<std::string> movingBCmarkers = Feel::FSI::detail::intersectionList( this->markersNameMovingBoundary(),
                                                                                              this->markerDirichletBClm() );
                form1( _test=this->XhDirichletLM(),_vector=R,
                       _rowstart=rowStartInVector+startDofIndexDirichletLM ) +=
                    integrate( _range=markedfaces(mesh,movingBCmarkers), //markedelements(this->meshDirichletLM(),movingBCmarkers),
                               _expr= -inner( idv(this->meshVelocity2()),id(lambdaBC) ) );
            }
#endif
        }

    }

    //------------------------------------------------------------------------------------//
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() && this->couplingFSIcondition() == "robin" )
    {
        double gammaRobinFSI = this->gammaNitschFSI();//2500;//10;
        double muFluid = this->mu();//0.03;

        if ( (!BuildCstPart && !UseJacobianLinearTerms) )
        {
            linearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= gammaRobinFSI*muFluid*inner(idv(u),id(v))/hFace(),
                           _geomap=this->geomap() );
        }

        if ( !BuildCstPart )
        {
            linearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= -gammaRobinFSI*muFluid*inner(idv(this->meshVelocity2()),id(u))/hFace(),
                           _geomap=this->geomap() );

            //Deformations tensor (trial)
            auto uEval = this->getSolution()->template element<0>();
            auto pEval = this->getSolution()->template element<1>();
            auto defv = sym(gradv(uEval));
            // Strain tensor (trial)
            auto Sigmav = -idv(pEval)*Id + 2*idv(mu)*defv;
            linearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= -inner( Sigmav*N(),id(u)),
                           _geomap=this->geomap() );
        }
    }
#endif
    //------------------------------------------------------------------------------------//

#if FLUIDMECHANICS_USE_PERIODICITY
    if ( !BuildCstPart )
    {
        std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
        double pressureJump = doption(_name="periodicity.pressure-jump",_prefix=this->prefix());
        linearForm_PatternCoupled +=
            integrate( _range=markedfaces( this->mesh(),this->mesh()->markerName(marker1) ),
                       _expr=-inner(pressureJump*N(),id(v) ) );
    }
#endif

    //------------------------------------------------------------------------------------//


    //if ( _doClose ) R->close();

    if (this->hasMarkerDirichletBCelimination() && !BuildCstPart && _doBCStrongDirichlet)
    {
        this->updateBCStrongDirichletResidual(R);
    }

    //if ( _doClose ) R->close();

    //------------------------------------------------------------------------------------//

    double timeElapsed=thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidual",
                                               "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str()+
                                               "\n--------------------------------------------------",
                                               this->worldComm(),this->verboseAllProc());

#endif // defined(FEELMODELS_FLUID_BUILD_RESIDUAL_CODE)

} // updateResidual

} // namespace FeelModels
} // namespace Feel


