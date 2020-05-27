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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianWeakBC( DataUpdateJacobian & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const
{
    using namespace Feel::vf;

    this->log("FluidMechanics","updateJacobianWeakBC", "start" );
    boost::mpi::timer t1;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool BuildCstPart = data.buildCstPart();

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();
    auto const& v = u;
    auto const& q = p;

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    auto bilinearFormVV = form2( _test=XhV,_trial=XhV,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix,
                                 _colstart=colStartInMatrix );
    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // density
    auto const& rho = this->materialProperties()->fieldRho();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();

    //--------------------------------------------------------------------------------------------------//
    // Dirichlet bc by using Nitsche formulation
    if ( this->hasMarkerDirichletBCnitsche() && BuildCstPart )
    {
        auto deft = sym(gradt(u));
        auto viscousStressTensor = 2*idv(mu)*deft;
         bilinearFormVV +=
             integrate( _range=markedfaces(mesh, this->markerDirichletBCnitsche()),
                        _expr= -timeSteppingScaling*inner(viscousStressTensor*N(),id(v) )
                        /**/   + timeSteppingScaling*this->dirichletBCnitscheGamma()*inner(idt(u),id(v))/hFace(),
                        _geomap=this->geomap() );
        auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=J,
                                     _pattern=size_type(Pattern::COUPLED),
                                     _rowstart=rowStartInMatrix,
                                     _colstart=colStartInMatrix+1 );
         bilinearFormVP +=
             integrate( _range=markedfaces(mesh, this->markerDirichletBCnitsche()),
                        _expr= timeSteppingScaling*inner( idt(p)*N(),id(v) ),
                        _geomap=this->geomap() );
    }
    //--------------------------------------------------------------------------------------------------//
    // Dirichlet bc by using Lagrange-multiplier
    if ( this->hasMarkerDirichletBClm() )
    {
        if ( BuildCstPart )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("dirichletlm") ) << " start dof index for dirichletlm is not present\n";
            size_type startBlockIndexDirichletLM = this->startSubBlockSpaceIndex("dirichletlm");

            auto lambdaBC = this->XhDirichletLM()->element();
            std::cout << "----lm-1-----"<<std::endl;
            form2( _test=XhV,_trial=this->XhDirichletLM(),_matrix=J,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix,
                   _colstart=colStartInMatrix+startBlockIndexDirichletLM )+=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(lambdaBC),id(u) ) );

            std::cout << "----lm-2-----"<<std::endl;
            form2( _test=this->XhDirichletLM(),_trial=XhV,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix+startBlockIndexDirichletLM,
                   _colstart=colStartInMatrix ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(u),id(lambdaBC) ) );
            std::cout << "----lm-3-----"<<std::endl;
        }
    }
    //--------------------------------------------------------------------------------------------------//
    // pressure bc
    if ( this->hasMarkerPressureBC() )
    {
        CHECK( this->hasStartSubBlockSpaceIndex("pressurelm1") ) << " start dof index for pressurelm1 is not present\n";
        size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");
        if (BuildCstPart)
        {
            if ( nDim==2 )
            {
                form2( _test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,0)*idt(M_fieldLagrangeMultiplierPressureBC1),
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idt(u),N()))(0,0)*id(M_fieldLagrangeMultiplierPressureBC1),
                               _geomap=this->geomap() );
            }
            else if ( nDim==3 )
            {
                auto alpha = 1./sqrt(1-Nz()*Nz());
                form2( _test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM1 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-timeSteppingScaling*trans(cross(id(u),N()))(0,2)*idt(M_fieldLagrangeMultiplierPressureBC1)*alpha,
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix+startBlockIndexPressureLM1,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr=-trans(cross(idt(u),N()))(0,2)*id(M_fieldLagrangeMultiplierPressureBC1)*alpha,
                               _geomap=this->geomap() );

                CHECK( this->hasStartSubBlockSpaceIndex("pressurelm2") ) << " start dof index for pressurelm2 is not present\n";
                size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");

                form2( _test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=J,_pattern=size_type(Pattern::COUPLED),
                       _rowstart=rowStartInMatrix,
                       _colstart=colStartInMatrix+startBlockIndexPressureLM2 ) +=
                    integrate( _range=markedfaces( this->mesh(),this->markerPressureBC() ),
                               _expr= -timeSteppingScaling*trans(cross(id(u),N()))(0,0)*alpha*idt(M_fieldLagrangeMultiplierPressureBC2)*Ny()
                               +timeSteppingScaling*trans(cross(id(u),N()))(0,1)*alpha*idt(M_fieldLagrangeMultiplierPressureBC2)*Nx(),
                               _geomap=this->geomap() );

                form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,_matrix=J,_pattern=size_type(Pattern::COUPLED),
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
    // windkessel implicit
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        CHECK( this->hasStartSubBlockSpaceIndex("windkessel") ) << " start dof index for windkessel is not present\n";
        size_type startBlockIndexWindkessel = this->startSubBlockSpaceIndex("windkessel");

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

                form2( _test=M_fluidOutletWindkesselSpace,_trial=XhV,_matrix=J,
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

                form2( _test=M_fluidOutletWindkesselSpace,_trial=XhV,_matrix=J,
                       _rowstart=blockStartWindkesselRow/*rowStartInMatrixWindkessel*/,
                       _colstart=colStartInMatrix )+=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr=-Rp*(trans(idt(u))*N())*id(presProximal),
                               _geomap=this->geomap() );
                //--------------------//
                // coupling with fluid model
                form2( _test=XhV, _trial=M_fluidOutletWindkesselSpace, _matrix=J,
                       _rowstart=rowStartInMatrix,
                       _colstart=blockStartWindkesselCol/*colStartInMatrixWindkessel*/ ) +=
                    integrate( _range=markedfaces(mesh,markerOutlet),
                               _expr= timeSteppingScaling*idt(presProximal)*trans(N())*id(v),
                               _geomap=this->geomap() );
            }
        }
    }
    //--------------------------------------------------------------------------------------------------//
    // slip bc
    if (BuildCstPart && !this->markerSlipBC().empty() )
    {
        // slip condition :
        auto P = Id-N()*trans(N());
        double gammaN = doption(_name="bc-slip-gammaN",_prefix=this->prefix());
        double gammaTau = doption(_name="bc-slip-gammaTau",_prefix=this->prefix());
        auto betaExtrapolate = M_bdfVelocity->poly();
        //auto beta = Beta.element<0>();
        auto beta = vf::project( _space=XhV,
                                 _range=boundaryfaces(mesh),
                                 _expr=idv(rho)*idv(betaExtrapolate) );
        auto Cn = val(gammaN*max(abs(trans(idv(beta))*N()),idv(mu)/vf::h()));
        auto Ctau = val(gammaTau*idv(mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) ));

        bilinearFormVV +=
            integrate( _range= markedfaces(mesh,this->markerSlipBC()),
                       _expr= Cn*(trans(idt(u))*N())*(trans(id(v))*N())+
                       Ctau*trans(idt(u))*id(v),
                       //+ trans(idt(p)*Id*N())*id(v)
                       //- trans(id(v))*N()* trans(2*idv(mu)*deft*N())*N()
                       _geomap=this->geomap()
                       );
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
            auto const& momentOfInertia = bpbc.body().momentOfInertia();

            bool hasActiveDofTranslationalVelocity = bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0;
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();
            bool hasActiveDofAngularVelocity = nLocalDofAngularVelocity > 0;
            if ( BuildCstPart )
            {
                if ( hasActiveDofTranslationalVelocity )
                {
                    auto const& basisToContainerGpTranslationalVelocityRow = J->mapRow().dofIdToContainerId( rowStartInMatrix+startBlockIndexTranslationalVelocity );
                    auto const& basisToContainerGpTranslationalVelocityCol = J->mapCol().dofIdToContainerId( colStartInMatrix+startBlockIndexTranslationalVelocity );
                    for (int d=0;d<nDim;++d)
                    {
                        J->add( basisToContainerGpTranslationalVelocityRow[d], basisToContainerGpTranslationalVelocityCol[d],
                                bpbc.bdfTranslationalVelocity()->polyDerivCoefficient(0)*massBody );
                    }
                }
                if ( hasActiveDofAngularVelocity )
                {
                    auto const& basisToContainerGpAngularVelocityRow = J->mapRow().dofIdToContainerId( rowStartInMatrix+startBlockIndexAngularVelocity );
                    auto const& basisToContainerGpAngularVelocityCol = J->mapCol().dofIdToContainerId( colStartInMatrix+startBlockIndexAngularVelocity );
                    for (int i=0;i<nLocalDofAngularVelocity;++i)
                    {
                        for (int j=0;j<nLocalDofAngularVelocity;++j)
                        {
                            J->add( basisToContainerGpAngularVelocityRow[i], basisToContainerGpAngularVelocityCol[j],
                                    bpbc.bdfAngularVelocity()->polyDerivCoefficient(0)*momentOfInertia(i,j) );
                        }
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------------------------------//
    double timeElapsed=t1.elapsed();
    this->log("FluidMechanics","updateJacobianWeakBC","finish in "+(boost::format("%1% s") % timeElapsed).str() );

} // updateJacobian



FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !this->hasStrongDirichletBC() ) return;
    this->log("FluidMechanics","updateJacobianDofElimination", "start" );

    this->timerTool("Solve").start();

    this->updateDofEliminationIds( "velocity", data );

    if ( this->hasMarkerPressureBC() )
    {
        this->updateDofEliminationIds( "pressurelm1", this->dofEliminationIds( "pressurebc-lm" ), data );
        if ( nDim == 3 )
            this->updateDofEliminationIds( "pressurelm2", this->dofEliminationIds( "pressurebc-lm" ), data );
    }

    for ( auto const& [bpname,bpbc] : M_bodySetBC )
    {
        if ( bpbc.hasTranslationalVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".translational-velocity";
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.translational-velocity" ), data );
        }
        if ( bpbc.hasAngularVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".angular-velocity";
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.angular-velocity" ), data );
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("FluidMechanics","updateJacobianDofElimination","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}



} // namespace FeelModels

} // namespace Feel
