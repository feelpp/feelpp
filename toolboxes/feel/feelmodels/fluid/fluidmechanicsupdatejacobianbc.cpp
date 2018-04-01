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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianWeakBC( DataUpdateJacobian & data, element_fluid_external_storage_type const& U ) const
{
    using namespace Feel::vf;

    this->log("FluidMechanics","updateJacobianWeakBC", "start" );
    boost::mpi::timer t1;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool BuildCstPart = data.buildCstPart();

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // density
    auto const& rho = this->densityViscosityModel()->fieldRho();
    // dynamic viscosity
    auto const& mu = this->densityViscosityModel()->fieldMu();

    //--------------------------------------------------------------------------------------------------//
    // Dirichlet bc by using Nitsche formulation
    if ( this->hasMarkerDirichletBCnitsche() && BuildCstPart )
    {
        // strain tensor
        auto const deft = sym(gradt(u));
        // stress tensor
        auto const Sigmat = -idt(p)*Id + 2*idv(mu)*deft;

        bilinearForm_PatternCoupled +=
            integrate( _range=markedfaces(mesh, this->markerDirichletBCnitsche()),
                       _expr= -trans(Sigmat*N())*id(v)
                       /**/   + this->dirichletBCnitscheGamma()*trans(idt(u))*id(v)/hFace(),
                       _geomap=this->geomap() );
    }
    //--------------------------------------------------------------------------------------------------//
    // Dirichlet bc by using Lagrange-multiplier
    if ( this->hasMarkerDirichletBClm() )
    {
        if ( BuildCstPart )
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
    // pressure bc
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
    // slip bc
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
    // fsi bc
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
    double timeElapsed=t1.elapsed();
    this->log("FluidMechanics","updateJacobianWeakBC","finish in "+(boost::format("%1% s") % timeElapsed).str() );

} // updateJacobian



FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianStrongDirichletBC( sparse_matrix_ptrtype& J,vector_ptrtype& RBis ) const
{
    if ( !this->hasDirichletBC() ) return;
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateJacobianStrongDirichletBC", "start",
                                               this->worldComm(),this->verboseAllProc());
    using namespace Feel::vf;

    boost::timer btimeStrongCL;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldVelocity();
    //auto RBis = this->backend()->newVector( J->mapRowPtr() );

    for( auto const& d : this->M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        auto exprUsed = vf::zero<nDim,1>();// 0*vf::one();
        if ( !listMarkerFaces.empty() )
            bilinearForm +=
                on( _range=markedfaces( mesh, listMarkerFaces ),
                    _element=u, _rhs=RBis, _expr=exprUsed );
        if ( !listMarkerEdges.empty() )
            bilinearForm +=
                on( _range=markededges(mesh, listMarkerEdges),
                    _element=u,_rhs=RBis,_expr=exprUsed );
        if ( !listMarkerPoints.empty() )
            bilinearForm +=
                on( _range=markedpoints(mesh, listMarkerPoints),
                    _element=u,_rhs=RBis,_expr=exprUsed );
    }
    for ( auto const& bcDirComp : this->M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto bilinearFormComp = form2( _test=this->functionSpaceVelocity(),_trial=this->functionSpaceVelocity(),_matrix=J,
                                           _rowstart=this->rowStartInMatrix(),
                                           _colstart=this->colStartInMatrix() );
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),comp ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            auto exprUsed = cst(0.);
            if ( !listMarkerFaces.empty() )
                bilinearFormComp +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=RBis,_expr=exprUsed );
            if ( !listMarkerEdges.empty() )
                bilinearFormComp +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=RBis,_expr=exprUsed );
            if ( !listMarkerPoints.empty() )
                bilinearFormComp +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=RBis,_expr=exprUsed );
        }
    }

    std::ostringstream ostr3;ostr3<<btimeStrongCL.elapsed()<<"s";
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateJacobianStrongDirichletBC", "finish in "+ostr3.str(),
                                               this->worldComm(),this->verboseAllProc());
}



} // namespace FeelModels

} // namespace Feel
