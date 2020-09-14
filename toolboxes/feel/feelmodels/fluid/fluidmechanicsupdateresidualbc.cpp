/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>

//#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>


namespace Feel
{
namespace FeelModels
{


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualWeakBC( DataUpdateResidual & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const
{
    this->log("FluidMechanics","updateResidualModel", "start" );

    boost::mpi::timer thetimer,thetimer2;
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !BuildCstPart;

    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );
    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto const& v = u;
    auto const& q = p;

    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormV = form1( _test=XhV, _vector=R,
                              _pattern=size_type(Pattern::COUPLED),
                              _rowstart=this->rowStartInVector() );


    auto Id = eye<nDim,nDim>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();
    // density
    auto const& rho = this->materialProperties()->fieldRho();

    //--------------------------------------------------------------------------------------------------//
    // Neumann boundary condition
    if ( BuildCstPart && doAssemblyRhs )
    {
        for( auto const& d : this->M_bcNeumannScalar )
            linearFormV +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -timeSteppingScaling*expression(d,this->symbolsExpr())*inner( N(),id(v) ),
                           _geomap=this->geomap() );
        for( auto const& d : this->M_bcNeumannVectorial )
            linearFormV +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d)) ),
                           _expr= -timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                           _geomap=this->geomap() );
        for( auto const& d : this->M_bcNeumannTensor2 )
            linearFormV +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d)) ),
                           _expr= -timeSteppingScaling*inner( expression(d,this->symbolsExpr())*N(),id(v) ),
                           _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasFluidOutletWindkessel() )
    {
        if ( this->hasFluidOutletWindkesselExplicit() )
        {
            if ( BuildCstPart && doAssemblyRhs )
            {
                auto const beta = M_bdfVelocity->poly();

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
                                             _expr=trans(idv(beta))*N() ).evaluate()(0,0);

                    double pressureDistalOld  = 0;
                    for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                        pressureDistalOld += Deltat*this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld().find(k)->second[i];

                    M_fluidOutletWindkesselPressureDistal[k] = alphaBF*pressureDistalOld + gammaBF*outletQ;
                    M_fluidOutletWindkesselPressureProximal[k] = kappaBF*outletQ + alphaBF*pressureDistalOld;

                    linearFormV +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= timeSteppingScaling*M_fluidOutletWindkesselPressureProximal[k]*trans(N())*id(v),
                                   _geomap=this->geomap() );
                }
            }
        } // explicit
        if ( this->hasFluidOutletWindkesselImplicit() )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("windkessel") ) << " start dof index for windkessel is not present\n";
            size_type startBlockIndexWindkessel = this->startSubBlockSpaceIndex("windkessel");

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
                    if ( BuildCstPart  && hasWindkesselActiveDof && doAssemblyRhs )
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
                        form1( _test=XhV, _vector=R,
                               _rowstart=rowStartInVector ) +=
                            integrate( _range=markedfaces(mesh,markerOutlet),
                                       _expr= timeSteppingScaling*idv(presProximal)*trans(N())*id(v),
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
        auto betaExtrapolated = M_bdfVelocity->poly();
        auto beta = vf::project( _space=XhV,
                                 _range=boundaryfaces(mesh),
                                 _expr=idv(rho)*idv(betaExtrapolated) );
        auto Cn = gammaN*max(abs(trans(idv(beta))*N()),idv(mu)/vf::h());
        auto Ctau = gammaTau*idv(mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) );

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

    //--------------------------------------------------------------------------------------------------//

    // weak formulation of the boundaries conditions
    if ( this->hasMarkerDirichletBCnitsche() )
    {
        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            // deformation and stress tensor (eval)
            auto defv = sym(gradv(u));
            auto Sigmav = -idv(p)*Id + 2*idv(mu)*defv;
            linearFormV +=
                integrate( _range=markedfaces(mesh,this->markerDirichletBCnitsche() ),
                           _expr= -timeSteppingScaling*trans(Sigmav*N())*id(v)
                           /**/   + timeSteppingScaling*this->dirichletBCnitscheGamma()*inner( idv(u),id(v) )/hFace(),
                           _geomap=this->geomap() );
        }
        if ( BuildCstPart && doAssemblyRhs )
        {
            for( auto const& d : this->M_bcDirichlet )
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "nitsche",name(d) ) ),
                               _expr= -timeSteppingScaling*this->dirichletBCnitscheGamma()*inner( expression(d,this->symbolsExpr()),id(v) )/hFace(),
                               _geomap=this->geomap() );
        }
    }

    //------------------------------------------------------------------------------------//
    // Dirichlet with Lagrange-mulitplier
    if ( this->hasMarkerDirichletBClm() )
    {
        CHECK( this->hasStartSubBlockSpaceIndex("dirichletlm") ) << " start dof index for dirichletlm is not present\n";
        size_type startBlockIndexDirichletLM = this->startSubBlockSpaceIndex("dirichletlm");

        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            auto lambdaBC = this->XhDirichletLM()->element();
            //int dataBaseIdLM = XVec->map().basisIndexFromGp( rowStartInVector+startBlockIndexDirichletLM );
            M_blockVectorSolution.setSubVector( lambdaBC, *XVec, rowStartInVector+startBlockIndexDirichletLM );

            linearFormV +=
                integrate( _range=markedfaces(mesh,this->markerDirichletBClm() ),
                           _expr= inner( idv(lambdaBC),id(u) ) );

            form1( _test=this->XhDirichletLM(),_vector=R,
                   _rowstart=rowStartInVector+startBlockIndexDirichletLM ) +=
                integrate( //_range=elements(this->meshDirichletLM()),
                    _range=markedfaces( this->mesh(),this->markerDirichletBClm() ),
                           _expr= inner(idv(u),id(lambdaBC) ) );
        }
#if 1
        if ( BuildCstPart && doAssemblyRhs )
        {
            auto lambdaBC = this->XhDirichletLM()->element();
            for( auto const& d : this->M_bcDirichlet )
                form1( _test=this->XhDirichletLM(),_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexDirichletLM ) +=
                    integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "lm",name(d) ) ),
                               //_range=markedelements(this->meshDirichletLM(),PhysicalName),
                               _expr= -inner( expression(d,this->symbolsExpr()),id(lambdaBC) ),
                               _geomap=this->geomap() );
        }
#endif

    }

        //------------------------------------------------------------------------------------//

    if ( this->hasMarkerPressureBC() )
    {

        if ( !BuildCstPart && !UseJacobianLinearTerms )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("pressurelm1") ) << " start dof index for pressurelm1 is not present\n";
            size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");

            auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( XVec, rowStartInVector+startBlockIndexPressureLM1 );

            if ( nDim==2 )
            {
                linearFormV +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,0)*idv(lambdaPressure1),
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
                linearFormV +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,2)*idv(lambdaPressure1)*alpha,
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idv(u),N()))(0,2)*id(lambdaPressure1)*alpha,
                               _geomap=this->geomap() );

                CHECK( this->hasStartSubBlockSpaceIndex("pressurelm2") ) << " start dof index for pressurelm2 is not present\n";
                size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");

                auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( XVec, rowStartInVector+startBlockIndexPressureLM2 );

                linearFormV +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -timeSteppingScaling*trans(cross(id(u),N()))(0,0)*alpha*idv(lambdaPressure2)*Ny()
                               +timeSteppingScaling*trans(cross(id(u),N()))(0,1)*alpha*idv(lambdaPressure2)*Nx(),
                               _geomap=this->geomap() );

                form1( _test=M_spaceLagrangeMultiplierPressureBC,_vector=R,
                       _rowstart=rowStartInVector+startBlockIndexPressureLM2 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -trans(cross(idv(u),N()))(0,0)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +trans(cross(idv(u),N()))(0,1)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );
            }
        }
        if ( BuildCstPart && doAssemblyRhs )
        {
            for( auto const& d : this->M_bcPressure )
            {
                linearFormV +=
                    integrate( _range=markedfaces(this->mesh(),this->markerPressureBC(name(d)) ),
                               _expr= timeSteppingScaling*expression(d,this->symbolsExpr())*trans(N())*id(v),
                               _geomap=this->geomap() );
            }
        }

    }

    //--------------------------------------------------------------------------------------------------//

    if ( !M_bodySetBC.empty() )
    {
        this->log("FluidMechanics","updateJacobianWeakBC","assembly of body bc");

        for ( auto const& [bpname,bpbc] : M_bodySetBC )
        {
            size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
            size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");

            double massBody = bpbc.massExpr().evaluate()(0,0);
            auto momentOfInertiaExpr = bpbc.momentOfInertiaExpr();
            auto const& momentOfInertia = bpbc.body().momentOfInertia();
            bool hasActiveDofTranslationalVelocity = bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0;
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();
            bool hasActiveDofAngularVelocity = nLocalDofAngularVelocity > 0;

            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                if ( hasActiveDofTranslationalVelocity )
                {
                    auto uTranslationalVelocity = bpbc.spaceTranslationalVelocity()->element( XVec, rowStartInVector+startBlockIndexTranslationalVelocity );
                    auto const& basisToContainerGpTranslationalVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    for (int d=0;d<nDim;++d)
                    {
                        R->add( basisToContainerGpTranslationalVelocityVector[d],
                                uTranslationalVelocity(d)*bpbc.bdfTranslationalVelocity()->polyDerivCoefficient(0)*massBody );
                    }
                }
                if ( hasActiveDofAngularVelocity )
                {
                    auto uAngularVelocity = bpbc.spaceAngularVelocity()->element( XVec, rowStartInVector+startBlockIndexAngularVelocity );
                    auto const& basisToContainerGpAngularVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexAngularVelocity );
                    auto contribLhsAngularVelocity = bpbc.bdfAngularVelocity()->polyDerivCoefficient(0)*((momentOfInertiaExpr*idv(uAngularVelocity)).evaluate(false));
                    for (int i=0;i<nLocalDofAngularVelocity;++i)
                    {
                        R->add( basisToContainerGpAngularVelocityVector[i],
                                //uAngularVelocity(d)*bpbc.bdfAngularVelocity()->polyDerivCoefficient(0)*momentOfInertia
                                contribLhsAngularVelocity(i,0)
                                );
                    }
                }
            }
            if ( BuildCstPart && doAssemblyRhs )
            {
                if ( hasActiveDofTranslationalVelocity )
                {
                    auto const& basisToContainerGpTranslationalVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    //std::vector<double> _gravity = { 0., 2. };
                    //double massTilde = 10;
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
                if ( hasActiveDofAngularVelocity )
                {
                    auto const& basisToContainerGpAngularVelocityVector = R->map().dofIdToContainerId( rowStartInVector+startBlockIndexAngularVelocity );
                    auto angularVelocityPolyDeriv = bpbc.bdfAngularVelocity()->polyDeriv();
                    auto contribRhsAngularVelocity = (momentOfInertiaExpr*idv(angularVelocityPolyDeriv)).evaluate(false);
                    for (int i=0;i<nLocalDofAngularVelocity;++i)
                    {
                        R->add( basisToContainerGpAngularVelocityVector[i],
                                -contribRhsAngularVelocity(i,0)
                                //momentOfInertia*angularVelocityPolyDeriv(d)
                                );
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateResidualModel","finish in "+(boost::format("%1% s") % timeElapsed).str() );


} // updateResidual

} // namespace FeelModels
} // namespace Feel
