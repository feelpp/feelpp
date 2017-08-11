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
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool _BuildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    if (this->verbose()) Feel::FeelModels::Log("--------------------------------------------------\n",
                                               this->prefix()+".FluidMechanics","updateJacobian", "start"+sc,
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer thetimer;

    bool BuildNonCstPart = !_BuildCstPart;
    //bool BuildCstPart = ( _BuildCstPart  && this->application()->vm()[prefixvm(this->prefix(),"jacobian-linear-update")].as<bool>() ) || _isFirstExperience ;
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
#if 0
    auto U = Xh->element("u"); //U = *XVec;
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

    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // strain tensor
    auto const deft = sym(gradt(u));
    // dynamic viscosity
    auto const& mu = this->densityViscosityModel()->fieldMu();
    auto const& rho = this->densityViscosityModel()->fieldRho();
    // stress tensor
    auto const Sigmat = -idt(p)*Id + 2*idv(mu)*deft;
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    boost::mpi::timer timerAssemble;

    //--------------------------------------------------------------------------------------------------//
    // convection terms
    if (BuildNonCstPart)
    {
        if (this->doStabConvectionEnergy())
        {
            // convection term + stabilisation energy of convection with neumann bc (until outflow bc) ( see Nobile thesis)
            // auto const convecTerm = (trans(val(gradv(u)*idv(*M_P0Rho))*idt(u)) + trans(gradt(u)*val(idv(u)*idv(*M_P0Rho)) ) )*id(v);
            // stabTerm = trans(divt(u)*val(0.5*idv(*M_P0Rho)*idv(u))+val(0.5*idv(*M_P0Rho)*divv(u))*idt(u))*id(v)

            auto const convecTerm = Feel::vf::FeelModels::fluidMecConvectionJacobianWithEnergyStab(u,rho);
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
            auto const convecTerm = Feel::vf::FeelModels::fluidMecConvectionJacobian(u,rho);
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

#if 0
    //--------------------------------------------------------------------------------------------------//
    // stabilisation energy of convection  with neumann bc (until outflow bc) ( see Nobile thesis)
    if (this->doStabConvectionEnergy())
    {
        if (BuildNonCstPart)
        {
            bilinearForm_PatternCoupled +=
                integrate (_range=M_rangeMeshElements,
                           _expr= trans(divt(u)*val(0.5*idv(rho)*idv(u))+val(0.5*idv(rho)*divv(u))*idt(u))*id(v),
                           _geomap=this->geomap() );
        }

        /*if ( this->isMoveDomain() && BuildCstPart)
         {
         bilinearForm_PatternDefault +=
         integrate (_range=elements(mesh),
         _expr= -0.5*idv(M_P0Rho)*divv(this->meshVelocity())*trans(idt(u))*id(v),
         _geomap=this->geomap() );
         }*/
    }
#endif
    double timeElapsed=timerAssemble.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateJacobian",
                                               "assemble convection term in "+(boost::format("%1% s") % timeElapsed).str(),
                                               this->worldComm(),this->verboseAllProc());

    //--------------------------------------------------------------------------------------------------//
    // sigma : grad(v) on Omega
    this->updateJacobianModel( data, U );

    //--------------------------------------------------------------------------------------------------//
    // incompressibility term
    if (BuildCstPart)
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -idv(rho)*divt(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    //transients terms
    bool Build_TransientTerm = !BuildCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=BuildCstPart;

    if (!this->isStationary() && Build_TransientTerm/*BuildCstPart*/)
    {
        bilinearForm_PatternDefault +=
            integrate( _range=M_rangeMeshElements,
                       _expr= idv(rho)*trans(idt(u))*id(v)*M_bdf_fluid->polyDerivCoefficient(0),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && BuildCstPart )
        {
            double beta = this->definePressureCstPenalisationBeta();
            bilinearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr=beta*idt(p)*id(q),
                           _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" && BuildCstPart )
        {
            CHECK( this->startBlockIndexFieldsInMatrix().find("define-pressure-cst-lm") != this->startBlockIndexFieldsInMatrix().end() )
                << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startBlockIndexDefinePressureCstLM = this->startBlockIndexFieldsInMatrix().find("define-pressure-cst-lm")->second;

            auto lambda = M_XhMeanPressureLM->element();
            form2( _test=Xh, _trial=M_XhMeanPressureLM, _matrix=J,
                   _rowstart=this->rowStartInMatrix(),
                   _colstart=this->colStartInMatrix()+startBlockIndexDefinePressureCstLM ) +=
                integrate( _range=M_rangeMeshElements,
                           _expr= id(p)*idt(lambda) /*+ idt(p)*id(lambda)*/,
                           _geomap=this->geomap() );

            form2( _test=M_XhMeanPressureLM, _trial=Xh, _matrix=J,
                   _rowstart=this->rowStartInMatrix()+startBlockIndexDefinePressureCstLM,
                   _colstart=this->colStartInMatrix() ) +=
                integrate( _range=M_rangeMeshElements,
                           _expr= + idt(p)*id(lambda),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->updateJacobianStabilisation( data, U );

    //--------------------------------------------------------------------------------------------------//

    if (BuildCstPart && !this->markerSlipBC().empty() )
    {
        // slip condition :
        auto P = Id-N()*trans(N());
        double gammaN = doption(_name="bc-slip-gammaN",_prefix=this->prefix());
        double gammaTau = doption(_name="bc-slip-gammaTau",_prefix=this->prefix());
        auto Beta = M_bdf_fluid->poly();
        //auto beta = Beta.element<0>();
        auto beta = vf::project( _space=Beta.template element<0>().functionSpace(),
                                 _range=boundaryfaces(Beta.template element<0>().mesh()),
                                 _expr=idv(rho)*idv(Beta.template element<0>()) );
        auto Cn = val(gammaN*max(abs(trans(idv(beta))*N()),idv(mu)/vf::h()));
        auto Ctau = val(gammaTau*idv(mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) ));

        bilinearForm_PatternCoupled +=
            integrate( _range= markedfaces(mesh,this->markerSlipBC()),
                       _expr= Cn*(trans(idt(u))*N())*(trans(id(v))*N())+
                       Ctau*trans(idt(u))*id(v),
                       //+ trans(idt(p)*Id*N())*id(v)
                       //- trans(id(v))*N()* trans(2*idv(mu)*deft*N())*N()
                       _geomap=this->geomap()
                       );
    }

    //--------------------------------------------------------------------------------------------------//
    // weak formulation of the boundaries conditions
    if ( this->hasMarkerDirichletBCnitsche() && BuildCstPart)
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=markedfaces(mesh, this->markerDirichletBCnitsche()),
                       _expr= -trans(Sigmat*N())*id(v)
                       /**/   + this->dirichletBCnitscheGamma()*trans(idt(u))*id(v)/hFace(),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    // windkessel implicit
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        CHECK( this->startBlockIndexFieldsInMatrix().find("windkessel") != this->startBlockIndexFieldsInMatrix().end() )
            << " start dof index for windkessel is not present\n";
        size_type startBlockIndexWindkessel = this->startBlockIndexFieldsInMatrix().find("windkessel")->second;

        if (BuildCstPart)
        {
            auto presDistalProximal = M_fluidOutletWindkesselSpace->element();
            auto presDistal = presDistalProximal.template element<0>();
            auto presProximal = presDistalProximal.template element<1>();

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
                int blockStartWindkesselRow = rowStartInMatrix + startBlockIndexWindkessel + 2*cptOutletUsed;
                int blockStartWindkesselCol = colStartInMatrix + startBlockIndexWindkessel + 2*cptOutletUsed;
                auto const& basisToContainerGpPressureDistalRow = J->mapRow().dofIdToContainerId( blockStartWindkesselRow );
                auto const& basisToContainerGpPressureDistalCol = J->mapCol().dofIdToContainerId( blockStartWindkesselCol );
                auto const& basisToContainerGpPressureProximalRow = J->mapRow().dofIdToContainerId( blockStartWindkesselRow+1 );
                auto const& basisToContainerGpPressureProximalCol = J->mapCol().dofIdToContainerId( blockStartWindkesselCol+1 );
                if ( hasWindkesselActiveDof )
                    CHECK( !basisToContainerGpPressureDistalRow.empty() && !basisToContainerGpPressureDistalCol.empty() &&
                           !basisToContainerGpPressureProximalRow.empty() && !basisToContainerGpPressureProximalCol.empty() ) << "incomplete datamap info";
                const size_type gpPressureDistalRow = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalRow[0] : 0;
                const size_type gpPressureDistalCol = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalCol[0] : 0;
                const size_type gpPressureProximalRow = (hasWindkesselActiveDof)? basisToContainerGpPressureProximalRow[0] : 0;
                const size_type gpPressureProximalCol = (hasWindkesselActiveDof)? basisToContainerGpPressureProximalCol[0] : 0;
                //const size_type rowStartInMatrixWindkessel = rowStartInMatrix + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                //const size_type colStartInMatrixWindkessel = colStartInMatrix + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                ++cptOutletUsed;
                //--------------------//
                // first line
                if ( hasWindkesselActiveDof )
                {
                    J->add( gpPressureDistalRow/*rowStartInMatrixWindkessel*/, gpPressureDistalCol/*colStartInMatrixWindkessel*/,
                            Cd*this->timeStepBDF()->polyDerivCoefficient(0)+1./Rd );
                }

                form2( _test=M_fluidOutletWindkesselSpace,_trial=Xh,_matrix=J,
                       _rowstart=blockStartWindkesselRow/*rowStartInMatrixWindkessel*/,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr=-(trans(idt(u))*N())*id(presDistal),
                               _geomap=this->geomap() );

                //--------------------//
                // second line
                if ( hasWindkesselActiveDof )
                {
                    J->add( gpPressureProximalRow/*rowStartInMatrixWindkessel+1*/, gpPressureProximalCol/*colStartInMatrixWindkessel+1*/,  1.);

                    J->add( gpPressureProximalRow/*rowStartInMatrixWindkessel+1*/, gpPressureDistalCol/*colStartInMatrixWindkessel*/  , -1.);
                }

                form2( _test=M_fluidOutletWindkesselSpace,_trial=Xh,_matrix=J,
                       _rowstart=blockStartWindkesselRow/*rowStartInMatrixWindkessel*/,
                       _colstart=colStartInMatrix )+=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr=-Rp*(trans(idt(u))*N())*id(presProximal),
                               _geomap=this->geomap() );
                //--------------------//
                // coupling with fluid model
                form2( _test=Xh, _trial=M_fluidOutletWindkesselSpace, _matrix=J,
                       _rowstart=rowStartInMatrix,
                       _colstart=blockStartWindkesselCol/*colStartInMatrixWindkessel*/ ) +=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr= idt(presProximal)*trans(N())*id(v),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//

    if (this->hasMarkerDirichletBClm())
    {
        if (BuildCstPart)
        {
            CHECK( this->startBlockIndexFieldsInMatrix().find("dirichletlm") != this->startBlockIndexFieldsInMatrix().end() )
                << " start dof index for dirichletlm is not present\n";
            size_type startBlockIndexDirichletLM = this->startBlockIndexFieldsInMatrix().find("dirichletlm")->second;

            auto lambdaBC = this->XhDirichletLM()->element();
            form2( _test=Xh,_trial=this->XhDirichletLM(),_matrix=J,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix,
                   _colstart=colStartInMatrix+startBlockIndexDirichletLM )+=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(lambdaBC),id(u) ) );

            form2( _test=this->XhDirichletLM(),_trial=Xh,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix+startBlockIndexDirichletLM,
                   _colstart=colStartInMatrix ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(u),id(lambdaBC) ) );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasMarkerPressureBC() )
    {
        CHECK( this->startBlockIndexFieldsInMatrix().find("pressurelm1") != this->startBlockIndexFieldsInMatrix().end() )
            << " start dof index for pressurelm1 is not present\n";

        size_type startBlockIndexPressureLM1 = this->startBlockIndexFieldsInMatrix().find("pressurelm1")->second;
        if (BuildCstPart)
        {
            if ( nDim==2 )
            {
                form2( _test=Xh,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(id(u),N()))(0,0)*idt(M_fieldLagrangeMultiplierPressureBC1),
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=Xh,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idt(u),N()))(0,0)*id(M_fieldLagrangeMultiplierPressureBC1),
                               _geomap=this->geomap() );
            }
            else if ( nDim==3 )
            {
                auto alpha = 1./sqrt(1-Nz()*Nz());
                form2( _test=Xh,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(id(u),N()))(0,2)*idt(M_fieldLagrangeMultiplierPressureBC1)*alpha,
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=Xh,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idt(u),N()))(0,2)*id(M_fieldLagrangeMultiplierPressureBC1)*alpha,
                               _geomap=this->geomap() );

                CHECK( this->startBlockIndexFieldsInMatrix().find("pressurelm2") != this->startBlockIndexFieldsInMatrix().end() )
                    << " start dof index for pressurelm2 is not present\n";
                size_type startBlockIndexPressureLM2 = this->startBlockIndexFieldsInMatrix().find("pressurelm2")->second;

                form2( _test=Xh,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM2 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -trans(cross(id(u),N()))(0,0)*alpha*idt(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +trans(cross(id(u),N()))(0,1)*alpha*idt(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=Xh,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM2,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -trans(cross(idt(u),N()))(0,0)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +trans(cross(idt(u),N()))(0,1)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );
            }
        }

    }


    //--------------------------------------------------------------------------------------------------//

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() && ( this->couplingFSIcondition() == "robin-neumann" || this->couplingFSIcondition() == "robin-neumann-genuine" ||
                                   this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
                                   this->couplingFSIcondition() == "robin-neumann-generalized" || this->couplingFSIcondition() == "nitsche" ) )
    {
        //---------------------------------------------------------------------------//
        if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
             this->couplingFSIcondition() == "robin-neumann" || this->couplingFSIcondition() == "robin-neumann-genuine" ||
             this->couplingFSIcondition() == "nitsche" )
        //---------------------------------------------------------------------------//
        {
            double gammaRobinFSI = this->couplingFSI_Nitsche_gamma();
            if ( BuildCstPart )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                               _expr= ( gammaRobinFSI*idv(mu)/hFace() )*inner(idt(u),id(u)),
                               _geomap=this->geomap() );

                if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" ||
                     this->couplingFSIcondition() == "nitsche" )
                {
                    double alpha = this->couplingFSI_Nitsche_alpha();
                    double gamma0RobinFSI = this->couplingFSI_Nitsche_gamma0();
                    auto mysigma = id(q)*Id+2*alpha*idv(mu)*sym(grad(u));

                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                   _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*idv(mu)) )*idt(p)*id(q),
                                   _geomap=this->geomap() );

                    if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" )
                    {
                        bilinearForm_PatternCoupled +=
                            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                       _expr= -inner( idt(u), vf::N() )*id(q),
                                       _geomap=this->geomap() );
                    }
                    else if ( this->couplingFSIcondition() == "nitsche" )
                    {
                        bilinearForm_PatternCoupled +=
                            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                       _expr= -inner( idt(u), mysigma*vf::N() ),
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
                if ( BuildCstPart )
                {
                    this->meshALE()->revertReferenceMesh();
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                   _expr=this->couplingFSI_RNG_coeffForm2()*inner(idt(u),id(u)),
                                   _geomap=this->geomap() );
                    this->meshALE()->revertMovingMesh();
                }
            }
            else // useInterfaceOperator
            {
                CHECK( false ) << "TODO";
            }
        }
    }
#endif // FEELPP_MODELS_HAS_MESHALE

    //--------------------------------------------------------------------------------------------------//
    if ( M_useThermodynModel && M_useGravityForce )
    {
        DataUpdateJacobian dataThermo( data );
        dataThermo.setDoBCStrongDirichlet( false );
        M_thermodynModel->updateJacobian( dataThermo );

        if ( BuildNonCstPart )
        {
            auto XhT = M_thermodynModel->spaceTemperature();
            auto t = XhT->element(XVec, M_thermodynModel->rowStartInVector() );
            auto const& thermalProperties = M_thermodynModel->thermalProperties();

            auto thecoeff = idv(thermalProperties->fieldRho())*idv(thermalProperties->fieldHeatCapacity());
            form2( _test=XhT,_trial=XhT,_matrix=J,
                   _rowstart=M_thermodynModel->rowStartInMatrix(),
                   _colstart=M_thermodynModel->colStartInMatrix() ) +=
                integrate( _range=M_rangeMeshElementsAeroThermal,
                           _expr= thecoeff*(gradt(t)*idv(u))*id(t),
                       _geomap=this->geomap() );
            form2( _test=XhT,_trial=Xh,_matrix=J,
                   _rowstart=M_thermodynModel->rowStartInMatrix(),
                   _colstart=this->colStartInMatrix() ) +=
                integrate( _range=M_rangeMeshElementsAeroThermal,
                           _expr= thecoeff*(gradv(t)*idt(u))*id(t),
                           _geomap=this->geomap() );

            auto betaFluid = idv(thermalProperties->fieldThermalExpansion() );
            form2( _test=Xh,_trial=XhT,_matrix=J,
                   _rowstart=this->rowStartInMatrix(),
                   _colstart=M_thermodynModel->colStartInMatrix() ) +=
                integrate( _range=M_rangeMeshElementsAeroThermal,
                           _expr= idv(thermalProperties->fieldRho())*betaFluid*(idt(t))*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );
        }


    }
    //--------------------------------------------------------------------------------------------------//

    if ( BuildNonCstPart && _doBCStrongDirichlet )
    {
        if ( this->hasMarkerDirichletBCelimination() )
            this->updateBCStrongDirichletJacobian(J,RBis);

        std::list<std::string> markerDirichletEliminationOthers;
#if defined( FEELPP_MODELS_HAS_MESHALE )
        if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann")
        {
            for (std::string const& marker : this->markersNameMovingBoundary() )
                markerDirichletEliminationOthers.push_back( marker );
        }
#endif
        for ( auto const& inletbc : M_fluidInletDesc )
        {
            std::string const& marker = std::get<0>( inletbc );
            markerDirichletEliminationOthers.push_back( marker );
        }

        if ( !markerDirichletEliminationOthers.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(mesh, markerDirichletEliminationOthers ),
                    _element=u,_rhs=RBis,
                    _expr= vf::zero<nDim,1>() );


        if ( this->hasMarkerPressureBC() )
        {
            size_type startBlockIndexPressureLM1 = this->startBlockIndexFieldsInMatrix().find("pressurelm1")->second;
            form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,
                   _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                   _colstart=rowStartInMatrix+startBlockIndexPressureLM1 ) +=
                on( _range=boundaryfaces(M_meshLagrangeMultiplierPressureBC), _rhs=RBis,
                    _element=*M_fieldLagrangeMultiplierPressureBC1, _expr=cst(0.));
            if ( nDim == 3 )
            {
                size_type startBlockIndexPressureLM2 = this->startBlockIndexFieldsInMatrix().find("pressurelm2")->second;
                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM2,
                       _colstart=rowStartInMatrix+startBlockIndexPressureLM2 ) +=
                    on( _range=boundaryfaces(M_meshLagrangeMultiplierPressureBC), _rhs=RBis,
                        _element=*M_fieldLagrangeMultiplierPressureBC2, _expr=cst(0.));
            }
        }


        if ( M_useThermodynModel && M_useGravityForce )
        {
            M_thermodynModel->updateBCStrongDirichletJacobian( J,RBis );
        }

    }

    //--------------------------------------------------------------------------------------------------//

    /*double*/ timeElapsed=thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateJacobian",
                                               "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str()+
                                               "\n--------------------------------------------------",
                                               this->worldComm(),this->verboseAllProc());

} // updateJacobian

} // namespace FeelModels

} // namespace Feel


