/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateOseenWeakBC( sparse_matrix_ptrtype& A , vector_ptrtype& F, bool _BuildCstPart ) const
{
#if defined(FEELMODELS_FLUID_BUILD_LINEAR_CODE)
    using namespace Feel::vf;

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("FluidMechanics","updateOseenWeakBC", "start"+sc );

    boost::timer thetimer;

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool BuildNonCstPart_robinFSI = BuildNonCstPart;
    if (this->useFSISemiImplicitScheme()) BuildNonCstPart_robinFSI=BuildCstPart;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    auto U = Xh->element("u");
    auto V = Xh->element("v");
    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

    //Deformations tensor (trial)
    auto deft = sym(gradt(u));
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // dynamic viscosity
    auto const& mu = this->densityViscosityModel()->fieldMu();
    auto const& rho = this->densityViscosityModel()->fieldRho();
    // Strain tensor (trial)
    auto Sigmat = -idt(p)*Id + 2*idv(mu)*deft;

    //--------------------------------------------------------------------------------------------------//

    if (BuildNonCstPart && !this->markerSlipBC().empty() )
    {
        auto P = Id-N()*trans(N());
        double gammaN = doption(_name="bc-slip-gammaN",_prefix=this->prefix());
        double gammaTau = doption(_name="bc-slip-gammaTau",_prefix=this->prefix());
        auto Beta = M_bdf_fluid->poly();
        auto beta = vf::project( _space=Beta.template element<0>().functionSpace(),
                                 _range=boundaryfaces(Beta.template element<0>().mesh()),
                                 _expr=idv(rho)*idv(Beta.template element<0>()) );
        //auto beta = Beta.element<0>();
        auto Cn = gammaN*max(abs(trans(idv(beta))*N()),idv(mu)/vf::h());
        auto Ctau = gammaTau*idv(mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) );

        bilinearForm_PatternCoupled +=
            integrate( _range= markedfaces(mesh,this->markerSlipBC()),
                       _expr= val(Cn)*(trans(idt(u))*N())*(trans(id(v))*N())+
                       val(Ctau)*trans(idt(u))*id(v),
                       //+ trans(idt(p)*Id*N())*id(v)
                       //- trans(id(v))*N()* trans(2*idv(*M_P0Mu)*deft*N())*N()
                       _geomap=this->geomap()
                       );

    }

    //--------------------------------------------------------------------------------------------------//

    if (this->hasMarkerDirichletBClm())
    {
        CHECK( this->startDofIndexFieldsInMatrix().find("dirichletlm") != this->startDofIndexFieldsInMatrix().end() )
            << " start dof index for dirichletlm is not present\n";
        size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;
        auto lambdaBC = this->XhDirichletLM()->element();
        if (BuildCstPart)
        {
            form2( _test=Xh,_trial=this->XhDirichletLM(),_matrix=A,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix,
                   _colstart=colStartInMatrix+startDofIndexDirichletLM ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(lambdaBC),id(u) ) );

            form2( _test=this->XhDirichletLM(),_trial=Xh,_matrix=A,_pattern=size_type(Pattern::COUPLED),
                   _rowstart=rowStartInMatrix+startDofIndexDirichletLM,
                   _colstart=colStartInMatrix ) +=
                integrate( _range=elements(this->meshDirichletLM()),
                           _expr= inner( idt(u),id(lambdaBC) ) );
        }
        if ( BuildNonCstPart)
        {
            this->updateBCDirichletLagMultLinearPDE( F );

#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet" )
            {
                std::list<std::string> movingBCmarkers = detail::intersectionList( this->markersNameMovingBoundary(),
                                                                                   this->markerDirichletBClm() );
                form1( _test=this->XhDirichletLM(),_vector=F,
                       _rowstart=rowStartInVector+startDofIndexDirichletLM ) +=
                    integrate( _range=markedfaces(mesh,movingBCmarkers ), //markedelements(this->meshDirichletLM(),movingBCmarkers),
                               _expr= inner( idv(this->meshVelocity2()),id(lambdaBC) ) );
            }
#endif
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // weak formulation of the boundaries conditions
    if ( this->hasMarkerDirichletBCnitsche() )
    {
        if ( BuildCstPart)
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerDirichletBCnitsche() ),
                           _expr= -trans(Sigmat*N())*id(v)
                           /**/   + this->dirichletBCnitscheGamma()*trans(idt(u))*id(v)/hFace(),
                           _geomap=this->geomap() );
        }
        if ( BuildNonCstPart)
        {
            this->updateBCDirichletNitscheLinearPDE( F );

#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet" )
            {
                std::list<std::string> movingBCmarkers = detail::intersectionList( this->markersNameMovingBoundary(),
                                                                                   this->markerDirichletBCnitsche() );
                form1( _test=Xh, _vector=F,
                       _rowstart=rowStartInVector) +=
                    integrate( _range=markedfaces(mesh,movingBCmarkers ),
                               _expr= this->dirichletBCnitscheGamma()*trans(idv(this->meshVelocity2()))*id(v)/hFace(),
                               _geomap=this->geomap() );
            }
#endif
        }

    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" )
    {
        if ( this->fluidOutletWindkesselCoupling() == "explicit" )
        {
            if (BuildNonCstPart)
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

                    std::string markerOutlet= this->fluidOutletMarkerName(k);//"Sortie";
                    auto outletQ = integrate(_range=markedfaces(mesh,markerOutlet),
                                             _expr=trans(idv(Beta.template element<0>()))*N() ).evaluate()(0,0);

                    double pressureDistalOld  = 0;
                    for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                        pressureDistalOld += Deltat*this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld()[k][i];

                    M_fluidOutletWindkesselPressureDistal[k] = alphaBF*pressureDistalOld + gammaBF*outletQ;
                    M_fluidOutletWindkesselPressureProximal[k] = kappaBF*outletQ + alphaBF*pressureDistalOld;

                    form1( _test=Xh, _vector=F,
                           _rowstart=rowStartInVector ) +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= -M_fluidOutletWindkesselPressureProximal[k]*trans(N())*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }
        else if ( this->fluidOutletWindkesselCoupling() == "implicit" )
        {
            CHECK( this->startDofIndexFieldsInMatrix().find("windkessel") != this->startDofIndexFieldsInMatrix().end() )
                << " start dof index for windkessel is not present\n";
            size_type startDofIndexWindkessel = this->startDofIndexFieldsInMatrix().find("windkessel")->second;

            if (BuildNonCstPart)
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

                    const size_type rowStartInMatrixWindkessel = rowStartInMatrix + startDofIndexWindkessel + 2*k;
                    const size_type colStartInMatrixWindkessel = colStartInMatrix + startDofIndexWindkessel + 2*k;
                    const size_type rowStartInVectorWindkessel = rowStartInVector + startDofIndexWindkessel + 2*k;
                    //--------------------//
                    // 1ere ligne
                    if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                    {
                        A->add( rowStartInMatrixWindkessel,colStartInMatrixWindkessel,
                                Cd*this->timeStepBDF()->polyDerivCoefficient(0)+1./Rd );
                    }

                    form2( _test=M_fluidOutletWindkesselSpace,_trial=Xh,_matrix=A,
                           _rowstart=rowStartInMatrixWindkessel,
                           _colstart=colStartInMatrix ) +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr=-(trans(idt(u))*N())*id(presDistal) );

                    if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                    {
                        double pressureDistalOld  = 0;
                        for ( uint8_type i = 0; i < this->timeStepBDF()->timeOrder(); ++i )
                            pressureDistalOld += this->timeStepBDF()->polyDerivCoefficient( i+1 )*this->fluidOutletWindkesselPressureDistalOld()[k][i];
                        // add in vector
                        F->add( rowStartInVectorWindkessel, Cd*pressureDistalOld);
                    }
                    //--------------------//
                    // 2eme ligne
                    if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
                    {
                        A->add( rowStartInMatrixWindkessel+1, colStartInMatrixWindkessel+1,  1.);

                        A->add( rowStartInMatrixWindkessel+1, colStartInMatrixWindkessel  , -1.);
                    }

                    form2( _test=M_fluidOutletWindkesselSpace,_trial=Xh,_matrix=A,
                           _rowstart=rowStartInMatrixWindkessel,
                           _colstart=colStartInMatrix )+=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr=-Rp*(trans(idt(u))*N())*id(presProximal) );
                    //--------------------//
                    // coupling with fluid model
                    form2( _test=Xh, _trial=M_fluidOutletWindkesselSpace, _matrix=A,
                           _rowstart=rowStartInMatrix,
                           _colstart=colStartInMatrixWindkessel ) +=
                        integrate( _range=markedfaces(mesh,markerOutlet),
                                   _expr= idt(presProximal)*trans(N())*id(v),
                                   _geomap=this->geomap() );



                }
            }

        }
    }

    //--------------------------------------------------------------------------------------------------//

    if ( this->isMoveDomain() && this->couplingFSIcondition() == "robin" )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        double gammaRobinFSI = this->gammaNitschFSI();//2500;//10;
        //double muFluid = this->mu();//0.03;
        if ( BuildNonCstPart_robinFSI )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= gammaRobinFSI*idv(mu)*inner(idt(u),id(u))/hFace(),
                           _geomap=this->geomap() );
        }

        if ( BuildNonCstPart )
        {
            form1( _test=Xh, _vector=F,
                   _rowstart=rowStartInVector ) +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= gammaRobinFSI*idv(mu)*inner(idv(this->meshVelocity2()),id(u))/hFace(),
                           _geomap=this->geomap() );

            //Deformations tensor (trial)
            auto const& uEval = this->fieldVelocity();
            auto const& pEval = this->fieldPressure();
            auto defv = sym(gradv(uEval));
            // Strain tensor (trial)
            auto Sigmav = -idv(pEval)*Id + 2*idv(mu)*defv;
            form1( _test=Xh, _vector=F,
                   _rowstart=rowStartInVector ) +=
                integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                           _expr= inner( Sigmav*N(),id(u)),
                           _geomap=this->geomap() );
        }
#endif
    }


    if ( UsePeriodicity && BuildNonCstPart )
    {
        std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
        double pressureJump = doption(_name="periodicity.pressure-jump",_prefix=this->prefix());
        form1( _test=Xh, _vector=F,
               _rowstart=rowStartInVector ) +=
            integrate( _range=markedfaces( this->mesh(),this->mesh()->markerName(marker1) ),
                       _expr=inner(pressureJump*N(),id(v) ) );
    }


    //--------------------------------------------------------------------------------------------------//

    std::ostringstream ostr;ostr<<thetimer.elapsed()<<"s";
    this->log("FluidMechanics","updateOseenWeakBC", "finish in "+ostr.str() );


#endif // defined(FEELMODELS_FLUID_BUILD_LINEAR_CODE)
} // updateOseenWeakBC

} // end namespace FeelModels
} // end namespace Feel
