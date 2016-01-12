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

    auto U = Xh->element("u"); //U = *XVec;
    // copy vector values in fluid element
    for ( size_type k=0;k<Xh->nLocalDofWithGhost();++k )
        U(k) = XVec->operator()(/*rowStartInVector+*/k);

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
                integrate ( _range=elements(mesh),
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
                integrate ( _range=elements(mesh),
                            _expr=convecTerm,
                            _geomap=this->geomap() );
        }
    }

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain() && BuildCstPart )
    {
        bilinearForm_PatternCoupled +=
            //bilinearForm_PatternDefault +=
            integrate (_range=elements(mesh),
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
                integrate (_range=elements(mesh),
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
    this->updateJacobianModel(U, J, RBis, _BuildCstPart);

    //--------------------------------------------------------------------------------------------------//
    // incompressibility term
    if (BuildCstPart)
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= -divt(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    //transients terms
    bool Build_TransientTerm = !BuildCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=BuildCstPart;

    if (!this->isStationary() && Build_TransientTerm/*BuildCstPart*/)
    {
        bilinearForm_PatternDefault +=
            integrate( _range=elements(mesh),
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
                integrate( _range=elements(Xh->mesh()),
                           _expr=beta*idt(p)*id(q),
                           _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" && BuildCstPart )
        {
            CHECK( this->startDofIndexFieldsInMatrix().find("define-pressure-cst-lm") != this->startDofIndexFieldsInMatrix().end() )
                << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startDofIndexDefinePressureCstLM = this->startDofIndexFieldsInMatrix().find("define-pressure-cst-lm")->second;

#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_ONLY_ON_BOUNDARY)
            auto therange=boundaryfaces(mesh);
#else
            auto therange=elements(mesh);
#endif
            auto lambda = M_XhMeanPressureLM->element();
            form2( _test=Xh, _trial=M_XhMeanPressureLM, _matrix=J,
                   _rowstart=this->rowStartInMatrix(),
                   _colstart=this->colStartInMatrix()+startDofIndexDefinePressureCstLM ) +=
                integrate( _range=therange,//elements(mesh),
                           _expr= id(p)*idt(lambda) /*+ idt(p)*id(lambda)*/,
                           _geomap=this->geomap() );

            form2( _test=M_XhMeanPressureLM, _trial=Xh, _matrix=J,
                   _rowstart=this->rowStartInMatrix()+startDofIndexDefinePressureCstLM,
                   _colstart=this->colStartInMatrix() ) +=
                integrate( _range=therange,//elements(mesh),
                           _expr= + idt(p)*id(lambda),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->updateJacobianStabilisation(U, J, RBis, _BuildCstPart);

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
        CHECK( this->startDofIndexFieldsInMatrix().find("windkessel") != this->startDofIndexFieldsInMatrix().end() )
            << " start dof index for windkessel is not present\n";
        size_type startDofIndexWindkessel = this->startDofIndexFieldsInMatrix().find("windkessel")->second;

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

                const size_type rowStartInMatrixWindkessel = rowStartInMatrix + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                const size_type colStartInMatrixWindkessel = colStartInMatrix + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                ++cptOutletUsed;
                //--------------------//
                // 1ere ligne
                if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                {
                    J->add( rowStartInMatrixWindkessel,colStartInMatrixWindkessel,
                            Cd*this->timeStepBDF()->polyDerivCoefficient(0)+1./Rd );
                }

                form2( _test=M_fluidOutletWindkesselSpace,_trial=Xh,_matrix=J,
                       _rowstart=rowStartInMatrixWindkessel,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr=-(trans(idt(u))*N())*id(presDistal),
                               _geomap=this->geomap() );

                //--------------------//
                // 2eme ligne
                if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                {
                    J->add( rowStartInMatrixWindkessel+1, colStartInMatrixWindkessel+1,  1.);

                    J->add( rowStartInMatrixWindkessel+1, colStartInMatrixWindkessel  , -1.);
                }

                form2( _test=M_fluidOutletWindkesselSpace,_trial=Xh,_matrix=J,
                       _rowstart=rowStartInMatrixWindkessel,
                       _colstart=colStartInMatrix )+=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr=-Rp*(trans(idt(u))*N())*id(presProximal),
                               _geomap=this->geomap() );
                //--------------------//
                // coupling with fluid model
                form2( _test=Xh, _trial=M_fluidOutletWindkesselSpace, _matrix=J,
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrixWindkessel ) +=
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
            CHECK( this->startDofIndexFieldsInMatrix().find("dirichletlm") != this->startDofIndexFieldsInMatrix().end() )
                << " start dof index for dirichletlm is not present\n";
            size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;

            auto lambdaBC = this->XhDirichletLM()->element();
            form2( _test=Xh,_trial=this->XhDirichletLM(),_matrix=J,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix,
                   _colstart=colStartInMatrix+startDofIndexDirichletLM )+=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(lambdaBC),id(u) ) );

            form2( _test=this->XhDirichletLM(),_trial=Xh,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix+startDofIndexDirichletLM,
                   _colstart=colStartInMatrix ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(u),id(lambdaBC) ) );
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


