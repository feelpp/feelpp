/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */


#include <feel/feelmodels/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelvf/fluidmecconvection.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    using namespace Feel::vf;
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

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

#if 0
    auto U = Xh->element("u");//U = *XVec;
    // copy vector values in fluid element
#if 0
    for ( size_type k=0;k<Xh->nLocalDofWithGhost();++k )
        U(k) = XVec->operator()(/*rowStartInVector+*/k);
#else
    M_blockVectorSolution.setSubVector( U, *XVec, rowStartInVector );
#endif
#else
    auto U = Xh->element(XVec, rowStartInVector);
#endif

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
                integrate( _range=M_rangeMeshElements,
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
                integrate( _range=M_rangeMeshElements,
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
            integrate( _range=M_rangeMeshElements,
                       _expr= -val(idv(rho)*trans( gradv(u)*( idv( this->meshVelocity() ))))*id(v),
                       _geomap=this->geomap() );
        timeElapsedBis=thetimerBis.elapsed();
        this->log("FluidMechanics","updateResidual","build convective--2-- term in "+(boost::format("%1% s") % timeElapsedBis ).str() );
    }
#endif

    //--------------------------------------------------------------------------------------------------//

    this->updateResidualModel( data, U );

    //--------------------------------------------------------------------------------------------------//

    // incompressibility term
    if (!BuildCstPart && !UseJacobianLinearTerms )
    {
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -idv(rho)*divv(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( BuildCstPart )
    {
        //neumann boundary condition
        this->updateBCNeumannResidual( R );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasFluidOutletWindkessel() )
    {
        if ( this->hasFluidOutletWindkesselExplicit() )
        {
            if ( BuildCstPart )
            {
                auto const Beta = M_bdf_fluid->poly();

                for (int k=0;k<this->nFluidOutlet();++k)
                {
                    if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" || std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) != "explicit" )
                        continue;

                    // Windkessel model
                    std::string markerOutlet = std::get<0>( M_fluidOutletsBCType[k] );
                    auto const& windkesselParam = std::get<2>( M_fluidOutletsBCType[k] );
                    double Rd=std::get<1>(windkesselParam);
                    double Rp=std::get<2>(windkesselParam);
                    double Cd=std::get<3>(windkesselParam);
                    double Deltat = this->timeStepBDF()->timeStep();

                    double xiBF = Rd*Cd*this->timeStepBDF()->polyDerivCoefficient(0)*Deltat+Deltat;
                    double alphaBF = Rd*Cd/(xiBF);
                    double gammaBF = Rd*Deltat/xiBF;
                    double kappaBF = (Rp*xiBF+ Rd*Deltat)/xiBF;

                    auto outletQ = integrate(_range=markedfaces(mesh,markerOutlet),
                                             _expr=trans(idv(Beta.template element<0>()))*N() ).evaluate()(0,0);

                    double pressureDistalOld  = 0;
                    for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                        pressureDistalOld += Deltat*this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld().find(k)->second[i];

                    M_fluidOutletWindkesselPressureDistal[k] = alphaBF*pressureDistalOld + gammaBF*outletQ;
                    M_fluidOutletWindkesselPressureProximal[k] = kappaBF*outletQ + alphaBF*pressureDistalOld;

                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= M_fluidOutletWindkesselPressureProximal[k]*trans(N())*id(v),
                                   _geomap=this->geomap() );
                }
            }
        } // explicit
        if ( this->hasFluidOutletWindkesselImplicit() )
        {
            CHECK( this->startBlockIndexFieldsInMatrix().find("windkessel") != this->startBlockIndexFieldsInMatrix().end() )
                << " start dof index for windkessel is not present\n";
            size_type startBlockIndexWindkessel = this->startBlockIndexFieldsInMatrix().find("windkessel")->second;

            if ( BuildCstPart || (!BuildCstPart && !UseJacobianLinearTerms) )
            {
                //auto presDistalProximal = M_fluidOutletWindkesselSpace->element();
                //auto presDistal = presDistalProximal.template element<0>();
                //auto presProximal = presDistalProximal.template element<1>();

                int cptOutletUsed = 0;
                for (int k=0;k<this->nFluidOutlet();++k)
                {
                    if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" || std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) != "implicit" )
                        continue;

                    // Windkessel model
                    std::string markerOutlet = std::get<0>( M_fluidOutletsBCType[k] );
                    auto const& windkesselParam = std::get<2>( M_fluidOutletsBCType[k] );
                    double Rd=std::get<1>(windkesselParam);
                    double Rp=std::get<2>(windkesselParam);
                    double Cd=std::get<3>(windkesselParam);
                    double Deltat = this->timeStepBDF()->timeStep();

                    bool hasWindkesselActiveDof = M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0;
                    int blockStartWindkesselVec = rowStartInVector + startBlockIndexWindkessel + 2*cptOutletUsed;
                    auto const& basisToContainerGpPressureDistalVec = R->map().dofIdToContainerId( blockStartWindkesselVec );
                    auto const& basisToContainerGpPressureProximalVec = R->map().dofIdToContainerId( blockStartWindkesselVec+1 );
                    if ( hasWindkesselActiveDof )
                        CHECK( !basisToContainerGpPressureDistalVec.empty() && !basisToContainerGpPressureProximalVec.empty() ) << "incomplete datamap info";
                    const size_type gpPressureDistalVec = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalVec[0] : 0;
                    const size_type gpPressureProximalVec = (hasWindkesselActiveDof)? basisToContainerGpPressureProximalVec[0] : 0;

                    //const size_type rowStartWindkessel = startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                    //const size_type rowStartInVectorWindkessel = rowStartInVector + rowStartWindkessel;
                    ++cptOutletUsed;
                    //----------------------------------------------------//
                    if ( BuildCstPart  && hasWindkesselActiveDof )
                    {
                        double pressureDistalOld  = 0;
                        for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                            pressureDistalOld += this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld().find(k)->second[i];
                        // add in vector
                        R->add( gpPressureDistalVec/*rowStartInVectorWindkessel*/, -Cd*pressureDistalOld);
                    }
                    //----------------------------------------------------//
                    if ( !BuildCstPart && !UseJacobianLinearTerms )
                    {
                        auto presDistalProximal = M_fluidOutletWindkesselSpace->element(XVec,blockStartWindkesselVec);
                        auto presDistal = presDistalProximal.template element<0>();
                        auto presProximal = presDistalProximal.template element<1>();
#if 0
                        for ( size_type kk=0;kk<M_fluidOutletWindkesselSpace->nLocalDofWithGhost();++kk )
                            presDistalProximal( kk ) = XVec->operator()( rowStartWindkessel + kk);
#endif
                        //----------------------------------------------------//
                        // 1ere ligne
                        if ( hasWindkesselActiveDof )
                        {
                            const double value = presDistalProximal(0)*(Cd*this->timeStepBDF()->polyDerivCoefficient(0)+1./Rd);
                            R->add( gpPressureDistalVec/*rowStartInVectorWindkessel*/,value );
                        }
                        form1( _test=M_fluidOutletWindkesselSpace,_vector=R,
                               _rowstart=blockStartWindkesselVec/*rowStartInVectorWindkessel*/ ) +=
                            integrate( _range=markedfaces(mesh,markerOutlet),
                                       _expr=-(trans(idv(u))*N())*id(presDistal),
                                       _geomap=this->geomap() );
                        //----------------------------------------------------//
                        // 2eme ligne
                        if ( hasWindkesselActiveDof )
                        {
                            R->add( gpPressureProximalVec/*rowStartInVectorWindkessel+1*/,  presDistalProximal(1)-presDistalProximal(0) );
                        }
                        form1( _test=M_fluidOutletWindkesselSpace,_vector=R,
                               _rowstart=blockStartWindkesselVec/*rowStartInVectorWindkessel*/ )+=
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
            if ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann" && false )
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
                if ( !movingBCmarkers.empty() )
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

        if (M_haveSourceAdded)
        {
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -trans(idv(*M_SourceAdded))*id(v),
                           _geomap=this->geomap() );
        }
        if ( M_useGravityForce )
        {
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -idv(rho)*inner(M_gravityForce,id(u)),
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
                integrate( _range=M_rangeMeshElements,
                           _expr= val(idv(rho)*trans(idv(u))*M_bdf_fluid->polyDerivCoefficient(0))*id(v),
                           _geomap=this->geomap() );
        }

        if (BuildCstPart)
        {
            auto Buzz = M_bdf_fluid->polyDeriv();
            auto buzz = Buzz.template element<0>();
            linearForm_PatternDefault +=
                integrate( _range=M_rangeMeshElements,
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
                integrate( _range=M_rangeMeshElements,
                           _expr=beta*idv(p)*id(q),
                           _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" )
        {
            CHECK( this->startBlockIndexFieldsInMatrix().find("define-pressure-cst-lm") != this->startBlockIndexFieldsInMatrix().end() )
                << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startBlockIndexDefinePressureCstLM = this->startBlockIndexFieldsInMatrix().find("define-pressure-cst-lm")->second;

            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                auto lambda = M_XhMeanPressureLM->element();
                M_blockVectorSolution.setSubVector( lambda, *XVec, rowStartInVector+startBlockIndexDefinePressureCstLM );

                //for ( size_type k=0;k<M_XhMeanPressureLM->nLocalDofWithGhost();++k )
                //    lambda( k ) = XVec->operator()( startDofIndexDefinePressureCstLM + k);

                form1( _test=M_XhMeanPressureLM,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexDefinePressureCstLM ) +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= id(p)*idv(lambda) + idv(p)*id(lambda),
                               _geomap=this->geomap() );
            }
#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE)
            if (BuildCstPart)
            {
                auto lambda = M_XhMeanPressureLM->element();
                form1( _test=M_XhMeanPressureLM,_vector=R,
                       _rowstart=rowStartInVector+startDofIndexDefinePressureCstLM ) +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= -(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE(this->shared_from_this()))*id(lambda),
                               _geomap=this->geomap() );
            }
#endif
        }
    }


    //------------------------------------------------------------------------------------//

    this->updateResidualStabilisation( data, U );

    //------------------------------------------------------------------------------------//
    if (this->hasMarkerDirichletBClm())
    {
        CHECK( this->startBlockIndexFieldsInMatrix().find("dirichletlm") != this->startBlockIndexFieldsInMatrix().end() )
            << " start dof index for dirichletlm is not present\n";
        size_type startBlockIndexDirichletLM = this->startBlockIndexFieldsInMatrix().find("dirichletlm")->second;

        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            auto lambdaBC = this->XhDirichletLM()->element();
            //int dataBaseIdLM = XVec->map().basisIndexFromGp( rowStartInVector+startBlockIndexDirichletLM );
            M_blockVectorSolution.setSubVector( lambdaBC, *XVec, rowStartInVector+startBlockIndexDirichletLM );

            form1( _test=Xh,_vector=R,
                   _rowstart=rowStartInVector )+=
                integrate( _range=markedfaces(mesh,this->markerDirichletBClm() ),
                           _expr= inner( idv(lambdaBC),id(u) ) );

            form1( _test=this->XhDirichletLM(),_vector=R,
                   _rowstart=rowStartInVector+startBlockIndexDirichletLM ) +=
                integrate( //_range=elements(this->meshDirichletLM()),
                    _range=markedfaces( this->mesh(),this->markerDirichletBClm() ),
                           _expr= inner(idv(u),id(lambdaBC) ) );
        }

        if ( BuildCstPart )
        {
            this->updateBCDirichletLagMultResidual( R );

#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann" && false )
            {
                auto lambdaBC = this->XhDirichletLM()->element();
                form1( _test=this->XhDirichletLM(),_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexDirichletLM ) +=
                    integrate( _range=markedfaces(mesh,this->markersNameMovingBoundary() ),
                               _expr= -inner( idv(this->meshVelocity2()),id(lambdaBC) ) );
            }
#endif
        }

    }
    //------------------------------------------------------------------------------------//

    if ( this->hasMarkerPressureBC() )
    {

        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            CHECK( this->startBlockIndexFieldsInMatrix().find("pressurelm1") != this->startBlockIndexFieldsInMatrix().end() )
                << " start dof index for pressurelm1 is not present\n";
            size_type startBlockIndexPressureLM1 = this->startBlockIndexFieldsInMatrix().find("pressurelm1")->second;

            auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( XVec, rowStartInVector+startBlockIndexPressureLM1 );

            if ( nDim==2 )
            {
                form1( _test=Xh,_vector=R,
                       _rowstart=rowStartInVector ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(id(u),N()))(0,0)*idv(lambdaPressure1),
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idv(u),N()))(0,0)*id(lambdaPressure1),
                               _geomap=this->geomap() );
            }
            else if ( nDim==3 )
            {
                auto alpha = 1./sqrt(1-Nz()*Nz());
                form1( _test=Xh,_vector=R,
                       _rowstart=rowStartInVector ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(id(u),N()))(0,2)*idv(lambdaPressure1)*alpha,
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idv(u),N()))(0,2)*id(lambdaPressure1)*alpha,
                               _geomap=this->geomap() );

                CHECK( this->startBlockIndexFieldsInMatrix().find("pressurelm2") != this->startBlockIndexFieldsInMatrix().end() )
                    << " start dof index for pressurelm2 is not present\n";
                size_type startBlockIndexPressureLM2 = this->startBlockIndexFieldsInMatrix().find("pressurelm2")->second;

                auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( XVec, rowStartInVector+startBlockIndexPressureLM2 );

                form1( _test=Xh,_vector=R,
                       _rowstart=rowStartInVector ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -trans(cross(id(u),N()))(0,0)*alpha*idv(lambdaPressure2)*Ny()
                               +trans(cross(id(u),N()))(0,1)*alpha*idv(lambdaPressure2)*Nx(),
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM2 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -trans(cross(idv(u),N()))(0,0)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +trans(cross(idv(u),N()))(0,1)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );
            }
        }
        if ( BuildCstPart )
        {
            this->updateBCPressureResidual( R );
        }

    }


    //------------------------------------------------------------------------------------//
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() && ( this->couplingFSIcondition() == "robin-neumann" || this->couplingFSIcondition() == "robin-neumann-genuine" ||
                                   this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
                                   this->couplingFSIcondition() == "robin-neumann-generalized" || this->couplingFSIcondition() == "nitsche" ) )
    {
        auto const& uEval = (true)? this->fieldVelocity() : this->timeStepBDF()->unknown(0).template element<0>();
        auto const& pEval = (true)? this->fieldPressure() : this->timeStepBDF()->unknown(0).template element<1>();

        if ( !BuildCstPart )
        {
            auto defv = sym(gradv(uEval));
            // stress tensor
            auto Sigmav = -idv(pEval)*Id + 2*idv(mu)*defv;
            linearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= -inner( Sigmav*N(),id(v)),
                           _geomap=this->geomap() );
        }

        //---------------------------------------------------------------------------//
        if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
             this->couplingFSIcondition() == "robin-neumann" || this->couplingFSIcondition() == "robin-neumann-genuine" ||
             this->couplingFSIcondition() == "nitsche" )
        //---------------------------------------------------------------------------//
        {
            double gammaRobinFSI = this->couplingFSI_Nitsche_gamma();

            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                linearForm_PatternCoupled +=
                    integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                               _expr= (gammaRobinFSI*idv(mu)/hFace())*inner(idv(u),id(v)),
                               _geomap=this->geomap() );
            }

            if ( !BuildCstPart )
            {
                linearForm_PatternCoupled +=
                    integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                               _expr= -(gammaRobinFSI*idv(mu)/hFace())*inner(idv(this->meshVelocity2()),id(u)),
                               _geomap=this->geomap() );
            }
            if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" ||
                 this->couplingFSIcondition() == "nitsche" )
            {
                double alpha = this->couplingFSI_Nitsche_alpha();
                double gamma0RobinFSI = this->couplingFSI_Nitsche_gamma0();
                auto mysigma = id(q)*Id+2*alpha*idv(mu)*sym(grad(u));
                if ( !BuildCstPart && !UseJacobianLinearTerms )
                {
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                   _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*idv(mu)) )*idv(p)*id(q),
                                   _geomap=this->geomap() );
                    if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" )
                    {
                        linearForm_PatternCoupled +=
                            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                       _expr= -inner( idv(u), vf::N() )*id(q),
                                       _geomap=this->geomap() );
                    }
                    else if ( this->couplingFSIcondition() == "nitsche" )
                    {
                        linearForm_PatternCoupled +=
                            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                       _expr= -inner( idv(u), mysigma*vf::N() ),
                                       _geomap=this->geomap() );
                    }
                }
                if ( !BuildCstPart )
                {
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                   _expr= -( gamma0RobinFSI*hFace()/(gammaRobinFSI*idv(mu)) )*idv(pEval)*id(q),
                                   _geomap=this->geomap() );

                    if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" )
                    {
                        linearForm_PatternCoupled +=
                            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                       _expr= inner(idv(this->meshVelocity2()),vf::N())*id(q),
                                       _geomap=this->geomap() );
                    }
                    else if ( this->couplingFSIcondition() == "nitsche" )
                    {
                        linearForm_PatternCoupled +=
                            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                       _expr= inner(idv(this->meshVelocity2()),mysigma*vf::N()),
                                       _geomap=this->geomap() );
                    }
                }
            }
        }
        //---------------------------------------------------------------------------//
        else if ( this->couplingFSIcondition() == "robin-neumann-generalized" )
        //---------------------------------------------------------------------------//
        {
            bool useInterfaceOperator = this->couplingFSI_RNG_useInterfaceOperator() && !this->couplingFSI_solidIs1dReduced();

            if ( !useInterfaceOperator )
            {
                if ( !BuildCstPart && !UseJacobianLinearTerms )
                {
                    this->meshALE()->revertReferenceMesh();
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                   _expr=this->couplingFSI_RNG_coeffForm2()*inner(idv(u),id(v)),
                                   _geomap=this->geomap() );
                }

                // todo cst in implicit non cst in semi-implicit
                if ( BuildNonCstPart )
                {
                    this->meshALE()->revertReferenceMesh();
                    linearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                   _expr= inner(idv(this->couplingFSI_RNG_evalForm1()),id(v)),
                                   _geomap=this->geomap() );
                }
                this->meshALE()->revertMovingMesh();
            }
            else // useInterfaceOperator
            {
                CHECK(false) << "TODO";
            }
        }

    }
#endif
    //------------------------------------------------------------------------------------//

    if ( UsePeriodicity && !BuildCstPart )
    {
        std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
        double pressureJump = doption(_name="periodicity.pressure-jump",_prefix=this->prefix());
        linearForm_PatternCoupled +=
            integrate( _range=markedfaces( this->mesh(),this->mesh()->markerName(marker1) ),
                       _expr=-inner(pressureJump*N(),id(v) ) );
    }

    //------------------------------------------------------------------------------------//
    if ( M_useThermodynModel && M_useGravityForce )
    {
        DataUpdateResidual dataThermo( data );
        dataThermo.setDoBCStrongDirichlet( false );
        M_thermodynModel->updateResidual( dataThermo );

        if ( !BuildCstPart )
        {
            //auto const& t = M_thermodynModel->fieldTemperature();
            auto XhT = M_thermodynModel->spaceTemperature();
            auto t = XhT->element(XVec, M_thermodynModel->rowStartInVector() );
            auto const& thermalProperties = M_thermodynModel->thermalProperties();

            auto thecoeff = idv(thermalProperties->fieldRho())*idv(thermalProperties->fieldHeatCapacity());
            form1( _test=M_thermodynModel->spaceTemperature(), _vector=R,
                   _pattern=size_type(Pattern::COUPLED),
                   _rowstart=M_thermodynModel->rowStartInVector() ) +=
                integrate( _range=M_rangeMeshElementsAeroThermal,
                           _expr= thecoeff*(gradv(t)*idv(u))*id(t),
                           _geomap=this->geomap() );

            double T0 = M_BoussinesqRefTemperature;
            auto betaFluid = idv(thermalProperties->fieldThermalExpansion());
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElementsAeroThermal,
                           _expr= idv(thermalProperties->fieldRho())*(betaFluid*(idv(t)-T0))*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );
        }
    }
    //------------------------------------------------------------------------------------//

    bool hasStrongDirichletBC = this->hasMarkerDirichletBCelimination() || this->hasFluidInlet();
#if defined( FEELPP_MODELS_HAS_MESHALE )
    hasStrongDirichletBC = hasStrongDirichletBC || ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann" );
#endif
    if ( M_useThermodynModel && M_useGravityForce )
        hasStrongDirichletBC = hasStrongDirichletBC || M_thermodynModel->hasMarkerDirichletBCelimination();

    if (!BuildCstPart && _doBCStrongDirichlet && hasStrongDirichletBC)
    {
        R->close();

        auto resFeView = Xh->element(R,rowStartInVector);
        auto resFeViewVelocity = resFeView.template element<0>();

        auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("velocity");
        auto const& dofsWithValueImposedVelocity = ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposed->second : std::set<size_type>();
        for ( size_type thedof : dofsWithValueImposedVelocity )
            resFeViewVelocity.set( thedof,0. );
        sync( resFeViewVelocity, "=", dofsWithValueImposedVelocity );


        if ( this->hasMarkerPressureBC() )
        {
#if 0
            size_type startBlockIndexPressureLM1 = this->startBlockIndexFieldsInMatrix().find("pressurelm1")->second;
            auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM1 );
            lambdaPressure1.on(_range=boundaryfaces(M_meshLagrangeMultiplierPressureBC),
                               _expr=vf::zero<1,1>() );
            if ( nDim == 3 )
            {
                size_type startBlockIndexPressureLM2 = this->startBlockIndexFieldsInMatrix().find("pressurelm2")->second;
                auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM2 );
                lambdaPressure2.on(_range=boundaryfaces(M_meshLagrangeMultiplierPressureBC),
                                   _expr=vf::zero<1,1>() );
            }
#else
            auto itFindDofsWithValueImposedPressureBC = M_dofsWithValueImposed.find("pressurebc-lm");
            auto const& dofsWithValueImposedPressureBC = ( itFindDofsWithValueImposedPressureBC != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposedPressureBC->second : std::set<size_type>();
            size_type startBlockIndexPressureLM1 = this->startBlockIndexFieldsInMatrix().find("pressurelm1")->second;
            auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM1 );
            for ( size_type thedof : dofsWithValueImposedPressureBC )
                lambdaPressure1.set( thedof,0. );
            sync( lambdaPressure1, "=", dofsWithValueImposedPressureBC );
            if ( nDim == 3 )
            {
                size_type startBlockIndexPressureLM2 = this->startBlockIndexFieldsInMatrix().find("pressurelm2")->second;
                auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM2 );
                for ( size_type thedof : dofsWithValueImposedPressureBC )
                    lambdaPressure2.set( thedof,0. );
                sync( lambdaPressure2, "=", dofsWithValueImposedPressureBC );
            }

#endif
        }

        if ( M_useThermodynModel && M_useGravityForce )
            M_thermodynModel->updateBCDirichletStrongResidual( R );
    }

    //------------------------------------------------------------------------------------//

    double timeElapsed=thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidual",
                                               "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str()+
                                               "\n--------------------------------------------------",
                                               this->worldComm(),this->verboseAllProc());


} // updateResidual

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess(vector_ptrtype& U) const
{
    this->log("FluidMechanics","updateNewtonInitialGuess","start");

    size_type rowStartInVector = this->rowStartInVector();
    auto Xh = this->functionSpace();
    auto up = Xh->element( U, rowStartInVector );
    auto u = up.template element<0>();
    auto mesh = this->mesh();

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann")
    {
        this->log("FluidMechanics","updateNewtonInitialGuess","update moving boundary with strong Dirichlet");
        u.on(_range=markedfaces(mesh, this->markersNameMovingBoundary()),
             _expr=idv(this->meshVelocity2()) );
    }
#endif

    if ( this->hasMarkerDirichletBCelimination() )
    {
        // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
        std::map<std::string, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> > > mapMarkerBCToEntitiesMeshMarker;
        for( auto const& d : M_bcDirichlet )
        {
            mapMarkerBCToEntitiesMeshMarker[marker(d)] =
                detail::distributeMarkerListOnSubEntity( mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        }
        std::map<std::pair<std::string,ComponentType>, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> > > mapCompMarkerBCToEntitiesMeshMarker;
        for ( auto const& bcDirComp : M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                mapCompMarkerBCToEntitiesMeshMarker[std::make_pair(marker(d),comp)] =
                    detail::distributeMarkerListOnSubEntity(mesh, this->markerDirichletBCByNameId( "elimination",marker(d), comp )   );
            }
        }

        // strong Dirichlet bc with vectorial velocity
        for( auto const& d : M_bcDirichlet )
        {
            auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( marker(d) );
            if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                continue;
            auto const& listMarkerFaces = std::get<0>( itFindMarker->second );
            if ( !listMarkerFaces.empty() )
                u.on(_range=markedfaces(mesh,listMarkerFaces ),
                     _expr=expression(d) );
            auto const& listMarkerEdges = std::get<1>( itFindMarker->second );
            if ( !listMarkerEdges.empty() )
                u.on(_range=markededges(mesh,listMarkerEdges ),
                     _expr=expression(d) );
            auto const& listMarkerPoints = std::get<2>( itFindMarker->second );
            if ( !listMarkerPoints.empty() )
                u.on(_range=markedpoints(mesh,listMarkerPoints ),
                     _expr=expression(d) );
        }
        // strong Dirichlet bc with velocity components
        for ( auto const& bcDirComp : M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                auto itFindMarker = mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) );
                if ( itFindMarker == mapCompMarkerBCToEntitiesMeshMarker.end() )
                    continue;
                auto const& listMarkerFaces = std::get<0>( itFindMarker->second );
                if ( !listMarkerFaces.empty() )
                    u[comp].on(_range=markedfaces(mesh,listMarkerFaces ),
                               _expr=expression(d) );
                auto const& listMarkerEdges = std::get<1>( itFindMarker->second );
                if ( !listMarkerEdges.empty() )
                    u[comp].on(_range=markededges(mesh,listMarkerEdges ),
                               _expr=expression(d) );
                auto const& listMarkerPoints = std::get<2>( itFindMarker->second );
                if ( !listMarkerPoints.empty() )
                    u[comp].on(_range=markedpoints(mesh,listMarkerPoints ),
                               _expr=expression(d) );
            }
        }
    }

    for ( auto const& inletbc : M_fluidInletDesc )
    {
        std::string const& marker = std::get<0>( inletbc );
        auto itFindMark = M_fluidInletVelocityInterpolated.find(marker);
        if ( itFindMark == M_fluidInletVelocityInterpolated.end() )
            continue;
        auto const& inletVel = std::get<0>( itFindMark->second );
        u.on(_range=markedfaces(mesh, marker),
             _expr=-idv(inletVel)*N() );
    }
    // synchronize velocity dof on interprocess
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("velocity");
    if ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )
        sync( u, "=", itFindDofsWithValueImposed->second );

    if ( this->definePressureCst() && this->definePressureCstMethod() == "algebraic" )
    {
        auto upSol = this->functionSpace()->element( U, this->rowStartInVector() );
        auto pSol = upSol.template element<1>();
        CHECK( M_definePressureCstAlgebraicOperatorMeanPressure ) << "mean pressure operator does not init";
        double meanPressureCurrent = inner_product( *M_definePressureCstAlgebraicOperatorMeanPressure, pSol );
        double meanPressureImposed = 0;
        pSol.add( meanPressureImposed - meanPressureCurrent );
    }

    if ( M_useThermodynModel && M_useGravityForce )
    {
        M_thermodynModel->updateNewtonInitialGuess( U );
    }

    this->log("FluidMechanics","updateNewtonInitialGuess","finish");
}




} // namespace FeelModels
} // namespace Feel


