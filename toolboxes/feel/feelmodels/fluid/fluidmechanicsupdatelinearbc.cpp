/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEWeakBC( DataUpdateLinear & data ) const
{
    this->log("FluidMechanics","updateLinearPDEWeakBC", "start" );

    boost::timer thetimer;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;

    bool build_BoundaryNeumannTerm = BuildNonCstPart;
    bool BuildNonCstPart_robinFSI = BuildNonCstPart;
    if ( this->useFSISemiImplicitScheme() )
    {
        BuildNonCstPart_robinFSI = BuildCstPart;
        build_BoundaryNeumannTerm = BuildCstPart;
    }

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto const& u = this->fieldVelocity();
    auto const& v = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    auto const& q = this->fieldPressure();

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearFormVV = form2( _test=XhV,_trial=XhV,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix+0,
                                 _colstart=colStartInMatrix+0 );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix+0,
                                 _colstart=colStartInMatrix+1 );
    auto myLinearFormV = form1( _test=XhV, _vector=F,
                                _rowstart=this->rowStartInVector()+0 );

    //Deformations tensor (trial)
    auto deft = sym(gradt(u));
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();
    auto const& rho = this->materialProperties()->fieldRho();

    //--------------------------------------------------------------------------------------------------//

    if (BuildNonCstPart && !this->markerSlipBC().empty() )
    {
        auto P = Id-N()*trans(N());
        double gammaN = doption(_name="bc-slip-gammaN",_prefix=this->prefix());
        double gammaTau = doption(_name="bc-slip-gammaTau",_prefix=this->prefix());
        auto const& betaExtrapolate = M_bdfVelocity->poly();
        auto beta = vf::project( _space=XhV,
                                 _range=boundaryfaces(mesh),
                                 _expr=idv(rho)*idv(betaExtrapolate) );
        //auto beta = Beta.element<0>();
        auto Cn = gammaN*max(abs(trans(idv(beta))*N()),idv(mu)/vf::h());
        auto Ctau = gammaTau*idv(mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) );

        bilinearFormVV +=
            integrate( _range= markedfaces(mesh,this->markerSlipBC()),
                       _expr= timeSteppingScaling*val(Cn)*(trans(idt(u))*N())*(trans(id(v))*N())+
                       timeSteppingScaling*val(Ctau)*trans(idt(u))*id(v),
                       //+ trans(idt(p)*Id*N())*id(v)
                       //- trans(id(v))*N()* trans(2*idv(*M_P0Mu)*deft*N())*N()
                       _geomap=this->geomap()
                       );

    }

    //--------------------------------------------------------------------------------------------------//

    if (this->hasMarkerDirichletBClm())
    {
        CHECK( this->hasStartSubBlockSpaceIndex("dirichletlm") ) << " start dof index for dirichletlm is not present\n";
        size_type startBlockIndexDirichletLM = this->startSubBlockSpaceIndex("dirichletlm");
        auto lambdaBC = this->XhDirichletLM()->element();
        if (BuildCstPart)
        {
            form2( _test=XhV,_trial=this->XhDirichletLM(),_matrix=A,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix,
                   _colstart=colStartInMatrix+startBlockIndexDirichletLM ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(lambdaBC),id(u) ) );

            form2( _test=this->XhDirichletLM(),_trial=XhV,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix+startBlockIndexDirichletLM,
                   _colstart=colStartInMatrix ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(u),id(lambdaBC) ) );
        }
        if ( BuildNonCstPart)
        {
            for( auto const& d : this->M_bcDirichlet )
            {
                form1( _test=this->XhDirichletLM(),_vector=F,
                       _rowstart=this->rowStartInVector()+startBlockIndexDirichletLM ) +=
                    integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "lm",name(d) ) ),
                               _expr= inner( expression(d,this->symbolsExpr()),id(lambdaBC) ),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // weak formulation of the boundaries conditions
    if ( this->hasMarkerDirichletBCnitsche() )
    {
        if ( BuildCstPart)
        {
            auto viscousStressTensor = 2*idv(mu)*deft;
            bilinearFormVV +=
                integrate( _range=markedfaces(mesh,this->markerDirichletBCnitsche() ),
                           _expr= -timeSteppingScaling*inner(viscousStressTensor*N(),id(v) )
                           /**/   + timeSteppingScaling*this->dirichletBCnitscheGamma()*inner(idt(u),id(v))/hFace(),
                           _geomap=this->geomap() );
            bilinearFormVP +=
                integrate( _range=markedfaces(mesh,this->markerDirichletBCnitsche() ),
                           _expr= timeSteppingScaling*inner( idt(p)*N(), id(v) ),
                           _geomap=this->geomap() );
        }
        if ( BuildNonCstPart)
        {
            for( auto const& d : this->M_bcDirichlet )
            {
                myLinearFormV +=
                    integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "nitsche",name(d) ) ),
                               _expr= timeSteppingScaling*this->dirichletBCnitscheGamma()*inner( expression(d,this->symbolsExpr()),id(v) )/hFace(),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // Neumann bc
    if ( build_BoundaryNeumannTerm )
    {
        for( auto const& d : this->M_bcNeumannScalar )
            myLinearFormV +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= timeSteppingScaling*expression(d,this->symbolsExpr())*inner( N(),id(v) ),
                           _geomap=this->geomap() );
        for( auto const& d : this->M_bcNeumannVectorial )
            myLinearFormV +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d)) ),
                           _expr= timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                           _geomap=this->geomap() );
        for( auto const& d : this->M_bcNeumannTensor2 )
            myLinearFormV +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d)) ),
                           _expr= timeSteppingScaling*inner( expression(d,this->symbolsExpr())*N(),id(v) ),
                           _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasMarkerPressureBC() )
    {
        CHECK( this->hasStartSubBlockSpaceIndex("pressurelm1") ) << " start dof index for pressurelm1 is not present\n";
        size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");
        if (BuildCstPart)
        {
            if ( nDim == 2 )
            {
                form2( _test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,0)*idt(M_fieldLagrangeMultiplierPressureBC1),
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idt(u),N()))(0,0)*id(M_fieldLagrangeMultiplierPressureBC1),
                               _geomap=this->geomap() );
            }
            else if ( nDim == 3 )
            {
                auto alpha = 1./sqrt(1-Nz()*Nz());
                form2( _test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,2)*idt(M_fieldLagrangeMultiplierPressureBC1)*alpha,
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idt(u),N()))(0,2)*id(M_fieldLagrangeMultiplierPressureBC1)*alpha,
                               _geomap=this->geomap() );

                CHECK( this->hasStartSubBlockSpaceIndex("pressurelm2") ) << " start dof index for pressurelm2 is not present\n";
                size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");

                form2( _test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM2 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -timeSteppingScaling*trans(cross(id(u),N()))(0,0)*alpha*idt(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +timeSteppingScaling*trans(cross(id(u),N()))(0,1)*alpha*idt(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM2,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -trans(cross(idt(u),N()))(0,0)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +trans(cross(idt(u),N()))(0,1)*alpha*id(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );
            }
        }
        if ( BuildNonCstPart )
        {
            for( auto const& d : this->M_bcPressure )
            {
                myLinearFormV +=
                    integrate( _range=markedfaces(this->mesh(),this->markerPressureBC(name(d)) ),
                               _expr= -timeSteppingScaling*expression(d,this->symbolsExpr())*trans(N())*id(v),
                               _geomap=this->geomap() );
            }
        }

    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasFluidOutletWindkessel() )
    {
        this->timerTool("Solve").start();
        if ( this->hasFluidOutletWindkesselExplicit() )
        {
            if (BuildNonCstPart)
            {
                auto const& beta = M_bdfVelocity->poly();

                for (int k=0;k<this->nFluidOutlet();++k)
                {
                    if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" || std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) != "explicit" )
                        continue;

                    // Windkessel model
                    std::string markerOutlet = std::get<0>( M_fluidOutletsBCType[k] );
                    auto const& windkesselParam = std::get<2>( M_fluidOutletsBCType[k] );
                    double Rd=std::get<1>(windkesselParam);// 6.2e3;
                    double Rp=std::get<2>(windkesselParam);//400;
                    double Cd=std::get<3>(windkesselParam);//2.72e-4;

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

                    myLinearFormV +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= -timeSteppingScaling*M_fluidOutletWindkesselPressureProximal[k]*trans(N())*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }
        if ( this->hasFluidOutletWindkesselImplicit() )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("windkessel") ) << " start dof index for windkessel is not present\n";
            size_type startBlockIndexWindkessel = this->startSubBlockSpaceIndex("windkessel");

            if (BuildNonCstPart)
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
                    int blockStartWindkesselVec = rowStartInVector + startBlockIndexWindkessel + 2*cptOutletUsed;
                    auto const& basisToContainerGpPressureDistalRow = A->mapRow().dofIdToContainerId( blockStartWindkesselRow );
                    auto const& basisToContainerGpPressureDistalCol = A->mapCol().dofIdToContainerId( blockStartWindkesselCol );
                    auto const& basisToContainerGpPressureDistalVec = F->map().dofIdToContainerId( blockStartWindkesselVec );
                    auto const& basisToContainerGpPressureProximalRow = A->mapRow().dofIdToContainerId( blockStartWindkesselRow+1 );
                    auto const& basisToContainerGpPressureProximalCol = A->mapCol().dofIdToContainerId( blockStartWindkesselCol+1 );
                    auto const& basisToContainerGpPressureProximalVec = F->map().dofIdToContainerId( blockStartWindkesselVec+1 );
                    if ( hasWindkesselActiveDof )
                        CHECK( !basisToContainerGpPressureDistalRow.empty() && !basisToContainerGpPressureDistalCol.empty() &&
                               !basisToContainerGpPressureProximalRow.empty() && !basisToContainerGpPressureProximalCol.empty() &&
                               !basisToContainerGpPressureDistalVec.empty() && !basisToContainerGpPressureProximalVec.empty() ) << "incomplete datamap info";
                    const size_type gpPressureDistalRow = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalRow[0] : 0;
                    const size_type gpPressureDistalCol = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalCol[0] : 0;
                    const size_type gpPressureDistalVec = (hasWindkesselActiveDof)? basisToContainerGpPressureDistalVec[0] : 0;
                    const size_type gpPressureProximalRow = (hasWindkesselActiveDof)? basisToContainerGpPressureProximalRow[0] : 0;
                    const size_type gpPressureProximalCol = (hasWindkesselActiveDof)? basisToContainerGpPressureProximalCol[0] : 0;

                    //const size_type rowStartInMatrixWindkessel = rowStartInMatrix + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                    //const size_type colStartInMatrixWindkessel = colStartInMatrix + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                    //const size_type rowStartInVectorWindkessel = rowStartInVector + startDofIndexWindkessel + 2*cptOutletUsed/*k*/;
                    ++cptOutletUsed;
                    //--------------------//
                    // 1ere ligne
                    if ( hasWindkesselActiveDof )
                    {
                        A->add( gpPressureDistalRow/*rowStartInMatrixWindkessel*/,gpPressureDistalCol/*colStartInMatrixWindkessel*/,
                                Cd*this->timeStepBDF()->polyDerivCoefficient(0)+1./Rd );
                    }

                    form2( _test=M_fluidOutletWindkesselSpace,_trial=XhV,_matrix=A,
                           _rowstart=blockStartWindkesselRow/*rowStartInMatrixWindkessel*/,
                           _colstart=colStartInMatrix ) +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr=-(trans(idt(u))*N())*id(presDistal) );

                    if ( hasWindkesselActiveDof )
                    {
                        double pressureDistalOld  = 0;
                        for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                            pressureDistalOld += this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld().find(k)->second[i];
                        // add in vector
                        F->add( gpPressureDistalVec/*rowStartInVectorWindkessel*/, Cd*pressureDistalOld);
                    }
                    //--------------------//
                    // 2eme ligne
                    if ( hasWindkesselActiveDof )
                    {
                        A->add( gpPressureProximalRow/*rowStartInMatrixWindkessel+1*/, gpPressureProximalCol/*colStartInMatrixWindkessel+1*/,  1.);

                        A->add( gpPressureProximalRow/*rowStartInMatrixWindkessel+1*/, gpPressureDistalCol/*colStartInMatrixWindkessel*/  , -1.);
                    }

                    form2( _test=M_fluidOutletWindkesselSpace,_trial=XhV,_matrix=A,
                           _rowstart=blockStartWindkesselRow/*rowStartInMatrixWindkessel*/,
                           _colstart=colStartInMatrix )+=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr=-Rp*(trans(idt(u))*N())*id(presProximal) );
                    //--------------------//
                    // coupling with fluid model
                    form2( _test=XhV, _trial=M_fluidOutletWindkesselSpace, _matrix=A,
                           _rowstart=rowStartInMatrix,
                           _colstart=blockStartWindkesselCol/*colStartInMatrixWindkessel*/ ) +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= timeSteppingScaling*idt(presProximal)*trans(N())*id(v),
                                   _geomap=this->geomap() );



                }
            }

        }
        double timeElapsedBC = this->timerTool("Solve").stop();
        this->log("FluidMechanics","updateLinearPDE","assembly windkessel bc in "+(boost::format("%1% s") %timeElapsedBC).str() );
    }

    //--------------------------------------------------------------------------------------------------//

    if ( !M_bodyParticlesBC.empty() )
    {
        this->log("FluidMechanics","updateLinearPDE","assembly of body particles");

        for ( auto const& [bpname,bpbc] : M_bodyParticlesBC )
        {

            //CHECK( this->hasStartSubBlockSpaceIndex("particles-bc.translational-velocity") ) << " start dof index for particles-bc.translational-velocity is not present\n";
            //CHECK( this->hasStartSubBlockSpaceIndex("particles-bc.angular-velocity") ) << " start dof index for particles-bc.angular-velocity is not present\n";
            size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex("particles-bc."+bpbc.name()+".translational-velocity");
            size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("particles-bc."+bpbc.name()+".angular-velocity");

            // equation on translational velocity
#if 0
            // this is not good, not use an integrate but currently all row are eliminated
            form2( _test=M_XhNoSlipRigidParticlesTranslationalVelocity,_trial=M_XhNoSlipRigidParticlesTranslationalVelocity,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix+startBlockIndexTranslationalVelocity,
                   _colstart=colStartInMatrix+startBlockIndexTranslationalVelocity ) +=
                integrate( _range=elements(M_meshNoSlipRigidParticles),
                           _expr= M_bdfNoSlipRigidParticlesTranslationalVelocity->polyDerivCoefficient(0)*massParticle*inner( idt(M_fieldNoSlipRigidParticlesTranslationalVelocity), id(M_fieldNoSlipRigidParticlesTranslationalVelocity) ),
                           _geomap=this->geomap() );
            form1( _test=M_XhNoSlipRigidParticlesTranslationalVelocity, _vector=F,
                   _rowstart=rowStartInVector+startBlockIndexTranslationalVelocity ) +=
                integrate( _range=elements(M_meshNoSlipRigidParticles),
                           _expr= massParticle*inner( idv(M_bdfNoSlipRigidParticlesTranslationalVelocity->polyDeriv()),id(M_fieldNoSlipRigidParticlesTranslationalVelocity) ),
                           _geomap=this->geomap() );
#endif
            double massParticle = bpbc.massExpr().evaluate()(0,0);
            double momentOfInertia = bpbc.momentOfInertiaExpr().evaluate()(0,0);// NOT TRUE in 3D
            bool hasActiveDofTranslationalVelocity = bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0;
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();
            bool hasActiveDofAngularVelocity = nLocalDofAngularVelocity > 0;
            if ( BuildCstPart)
            {
                if ( hasActiveDofTranslationalVelocity )
                {
                    auto const& basisToContainerGpTranslationalVelocityRow = A->mapRow().dofIdToContainerId( rowStartInMatrix+startBlockIndexTranslationalVelocity );
                    auto const& basisToContainerGpTranslationalVelocityCol = A->mapCol().dofIdToContainerId( colStartInMatrix+startBlockIndexTranslationalVelocity );
                    for (int d=0;d<nDim;++d)
                    {
                        A->add( basisToContainerGpTranslationalVelocityRow[d], basisToContainerGpTranslationalVelocityCol[d],
                                bpbc.bdfTranslationalVelocity()->polyDerivCoefficient(0)*massParticle );
                    }
                }
                if ( hasActiveDofAngularVelocity )
                {
                    auto const& basisToContainerGpAngularVelocityRow = A->mapRow().dofIdToContainerId( rowStartInMatrix+startBlockIndexAngularVelocity );
                    auto const& basisToContainerGpAngularVelocityCol = A->mapCol().dofIdToContainerId( colStartInMatrix+startBlockIndexAngularVelocity );
                    for (int d=0;d<nLocalDofAngularVelocity;++d)
                    {
                        A->add( basisToContainerGpAngularVelocityRow[d], basisToContainerGpAngularVelocityCol[d],
                                bpbc.bdfAngularVelocity()->polyDerivCoefficient(0)*momentOfInertia );
                    }
                }
            }

            if ( BuildNonCstPart )
            {
                if ( hasActiveDofTranslationalVelocity )
                {
                    auto const& basisToContainerGpTranslationalVelocityVector = F->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    //std::vector<double> _gravity = { 0., 2. };
                    //double massTilde = 10;
                    for (int d=0;d<nDim;++d)
                    {
                        auto translationalVelocityPolyDeriv = bpbc.bdfTranslationalVelocity()->polyDeriv();
                        F->add( basisToContainerGpTranslationalVelocityVector[d],
                                massParticle*translationalVelocityPolyDeriv(d) );
#if 0
                        F->add( basisToContainerGpTranslationalVelocityVector[d],
                                massTilde*_gravity[d] );
#endif
                    }
                }
                if ( hasActiveDofAngularVelocity )
                {
                    auto const& basisToContainerGpAngularVelocityVector = F->map().dofIdToContainerId( rowStartInVector+startBlockIndexAngularVelocity );
                    for (int d=0;d<nLocalDofAngularVelocity;++d)
                    {
                        auto angularVelocityPolyDeriv = bpbc.bdfAngularVelocity()->polyDeriv();
                        F->add( basisToContainerGpAngularVelocityVector[d],
                                momentOfInertia*angularVelocityPolyDeriv(d) );
                    }
                }
            }
        } //  for ( auto const& [bpname,bpbc] : M_bodyParticlesBC )

    }
    //--------------------------------------------------------------------------------------------------//

#if 0
    if ( UsePeriodicity && BuildNonCstPart )
    {
        std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
        double pressureJump = doption(_name="periodicity.pressure-jump",_prefix=this->prefix());
        form1( _test=Xh, _vector=F,
               _rowstart=rowStartInVector ) +=
            integrate( _range=markedfaces( this->mesh(),this->mesh()->markerName(marker1) ),
                       _expr=inner(pressureJump*N(),id(v) ) );
    }
#endif

    //--------------------------------------------------------------------------------------------------//

    std::ostringstream ostr;ostr<<thetimer.elapsed()<<"s";
    this->log("FluidMechanics","updateLinearPDEWeakBC", "finish in "+ostr.str() );


} // updateLinearPDEWeakBC

} // end namespace FeelModels
} // end namespace Feel
