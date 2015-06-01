/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels2/solid/solidmecbase.hpp>

#include <feel/feelmodels2/modelvf/solidmecstvenantkirchhoff.hpp>
#include <feel/feelmodels2/modelvf/solidmecfirstpiolakirchhoff.hpp>
#include <feel/feelmodels2/modelvf/solidmecincompressibility.hpp>

namespace Feel
{
namespace FeelModels
{


SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidual( const vector_ptrtype& X, vector_ptrtype& R,
                                                                         bool _buildCstPart, bool UseJacobianLinearTerms,
                                                                         bool _doClose, bool _doBCStrongDirichlet ) const
{
#if defined(FEELMODELS_SOLID_BUILD_RESIDUAL_CODE)
    using namespace Feel::vf;

    std::string sc=(_buildCstPart)?" (build cst part)":" (build non cst part)";
    if (this->verbose()) Feel::FeelModels::Log("--------------------------------------------------\n",
                                         this->prefix()+".SolidMechanics","updateResidual", "start"+sc,
                                         this->worldComm(),this->verboseAllProc());
    boost::mpi::timer thetimer,thetimerBis;

    bool BuildNonCstPart = !_buildCstPart;
    bool BuildCstPart = _buildCstPart;
    bool BuildCstPart_BoundaryParoiMobile = BuildCstPart;
    if (this->useFSISemiImplicitScheme() )
    {
        BuildCstPart_BoundaryParoiMobile=BuildNonCstPart;
    }

    //--------------------------------------------------------------------------------------------------//

    mesh_ptrtype mesh = M_Xh->mesh();

    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormDisplacement = form1( _test=M_Xh, _vector=R,_rowstart=rowStartInVector );

    //--------------------------------------------------------------------------------------------------//

    auto u = M_Xh->element();
    auto v = u;
    // copy vector values in fluid element
    for ( size_type k=0;k<M_Xh->nLocalDofWithGhost();++k )
        u(k) = X->operator()(rowStartInVector+k);
    //auto buzz1 = M_newmark_displ_struct->previousUnknown();

    //--------------------------------------------------------------------------------------------------//

    //#if defined(SOLIDMECHANICS_VOLUME_FORCE)
    //auto f = SOLIDMECHANICS_VOLUME_FORCE(this->shared_from_this());
    //#endif
    //auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());

    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    //std::cout << "alpha_f " << alpha_f << "alpha_m " <<alpha_m <<"\n";
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);


    //--------------------------------------------------------------------------------------------------//

    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    auto const& rho = this->mechanicalProperties()->fieldRho();
    //Identity Matrix
    auto Id = eye<nDim,nDim>();
    auto Fv = Id + alpha_f*gradv(u);
    auto Ev = alpha_f*sym(gradv(u)) + alpha_f*alpha_f*0.5*trans(gradv(u))*gradv(u);
    auto Sv = idv(coeffLame1)*trace(Ev)*Id + 2*idv(coeffLame2)*Ev;
    //case elastic
    auto Ev_elastic = alpha_f*sym(gradv(u));
    auto Sv_elastic = idv(coeffLame1)*trace(Ev_elastic)*Id + 2*idv(coeffLame2)*Ev_elastic;

    //--------------------------------------------------------------------------------------------------//
#if 0
#if (SOLIDMECHANICS_DIM==2)
    auto Fv11 = Fv(0,0);
    auto Fv12 = Fv(0,1);
    auto Fv21 = Fv(1,0);
    auto Fv22 = Fv(1,1);
    auto detFv = Fv11*Fv22-Fv12*Fv21;
    auto InvFv = (cst(1.)/detFv)*mat<2,2>( Fv22,-Fv12,-Fv21,Fv11);
    auto traceCv = pow(Fv11,2) + pow(Fv21,2) + pow(Fv12,2) + pow(Fv22,2);
#endif
#if (SOLIDMECHANICS_DIM==3)
    auto Fv11 = Fv(0,0);
    auto Fv12 = Fv(0,1);
    auto Fv13 = Fv(0,2);
    auto Fv21 = Fv(1,0);
    auto Fv22 = Fv(1,1);
    auto Fv23 = Fv(1,2);
    auto Fv31 = Fv(2,0);
    auto Fv32 = Fv(2,1);
    auto Fv33 = Fv(2,2);
    auto detFv = Fv11*(Fv22*Fv33-Fv23*Fv32) - Fv21*(Fv12*Fv33-Fv13*Fv32) + Fv31*(Fv12*Fv23 - Fv13*Fv22);
    auto InvFv = (cst(1.)/detFv)*mat<3,3>( Fv22*Fv33-Fv23*Fv32 , Fv13*Fv32-Fv12*Fv33 , Fv12*Fv23-Fv13*Fv22,
                                           Fv23*Fv31-Fv21*Fv33 , Fv11*Fv33-Fv13*Fv31 , Fv13*Fv21-Fv11*Fv23,
                                           Fv21*Fv32-Fv22*Fv31 , Fv12*Fv31-Fv11*Fv32 , Fv11*Fv22-Fv12*Fv21
                                           );
    auto traceCv = pow(Fv11,2) + pow(Fv21,2) + pow(Fv31,2) + pow(Fv12,2) + pow(Fv22,2) + pow(Fv32,2) + pow(Fv13,2) + pow(Fv23,2) + pow(Fv33,2);
#endif

    //--------------------------------------------------------------------------------------------------//

    double C10=this->mechanicalProperties()->cstYoungModulus()/6.;//approximation
    // 2(C01+C10) = M_coefflame2 <=> C01 = M_coefflame2/2 -C10
    double C01=this->mechanicalProperties()->cstCoeffLame2()/2.-C10;
    //auto FSv_mooneyrivlin = Fv*2(C10 + C01*traceCv) - 2*C01*Fv*trans(Fv)*Fv;
    auto Sv_mooneyrivlin = 2*(C10 + C01*traceCv)*Id - 2*C01*trans(Fv)*Fv;
    auto FSv_mooneyrivlin = Fv*Sv_mooneyrivlin;
#endif
    //--------------------------------------------------------------------------------------------------//
    // stress tensor terms
    thetimerBis.restart();
    if (M_pdeType=="Hyper-Elasticity")
    {
        if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
        {
            if (!BuildCstPart)
            {
                auto const FSv_neohookean = Feel::vf::FeelModels::solidMecFirstPiolaKirchhoffTensor<3*(nOrderDisplacement-1)>(u,*this->mechanicalProperties());
                linearFormDisplacement +=
                    integrate( _range=elements(mesh),
                               //_expr= trace(val(Fv*Sv)*trans(grad(v))),
                               //_expr=trace(Feel::vf::FeelModels::stressStVenantKirchhoff(u,coeffLame1,coeffLame2)*trans(grad(v))),
                               //_expr=Feel::vf::FeelModels::stressStVenantKirchhoffResidual(u,coeffLame1,coeffLame2),// le dernier 
                               _expr= inner(FSv_neohookean,grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else if (this->mechanicalProperties()->materialLaw() == "NeoHookean")
        {
            if (!BuildCstPart)
            {
                auto const FSv_neohookean = Feel::vf::FeelModels::solidMecFirstPiolaKirchhoffTensor<2*nOrderDisplacement>(u,*this->mechanicalProperties());
                linearFormDisplacement +=
                    integrate( _range=elements(mesh),
                               //_expr= trace(val(FSv_neohookean)*trans(grad(v))),
                               _expr= inner(FSv_neohookean,grad(v)),
                               _quad=_Q<2*nOrderDisplacement+1>(),
                               _geomap=this->geomap() );
            }
        }
#if 0
        else if (this->mechanicalProperties()->materialLaw() == "MooneyRivlin")
        {
            if (!BuildCstPart)
            {
                linearFormDisplacement +=
                    integrate( _range=elements(mesh),
                               _expr= trace(val(FSv_mooneyrivlin)*trans(grad(v))),
                               _geomap=this->geomap() );
            }
        }
#endif
    } // if (M_pdeType=="Hyper-Elasticity")
    else if (M_pdeType=="Elasticity-Large-Deformation")
    {
        if (!BuildCstPart)
            linearFormDisplacement +=
                integrate( _range=elements(mesh),
                           _expr= trace(Sv*trans(grad(v))),
                           _geomap=this->geomap() );
    }
    else if (M_pdeType=="Elasticity")
    {
        if (!BuildCstPart && !UseJacobianLinearTerms)
            linearFormDisplacement +=
                integrate( _range=elements(mesh),
                           _expr= trace( Sv_elastic*trans(grad(v))),
                           _geomap=this->geomap() );
    }

    double timeElapsedBis=thetimerBis.elapsed();
    this->log("SolidMechanics","updateResidual",
              "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//
    // incompressibility terms
    if ( this->useDisplacementPressureFormulation() && !BuildCstPart)
    {
        auto p = M_XhPressure->element();//*M_fieldPressure;
        // copy vector values in pressure element
        size_type startDofIndexPressure = this->startDofIndexFieldsInMatrix().find("pressure")->second;
        for ( size_type k=0;k<M_XhPressure->nLocalDofWithGhost();++k )
            p(k) = X->operator()(rowStartInVector+startDofIndexPressure+k);
        // assemble
        this->updateResidualIncompressibilityTerms(u,p,R);
    }
    //--------------------------------------------------------------------------------------------------//
#if 0
    // viscoelastic terms
    // VERY OLD! : must be fix and updated
    this->updateResidualViscoElasticityTerms(u,R);
#endif
    //--------------------------------------------------------------------------------------------------//
    // source term
    if (BuildCstPart)
    {
        this->updateSourceTermResidual( R );
    }
    //--------------------------------------------------------------------------------------------------//
    // discretisation acceleration term
    if (!this->isStationary())
    {
        if (!BuildCstPart && !UseJacobianLinearTerms)
        {
            linearFormDisplacement +=
                integrate( _range=elements(mesh),
                           _expr= M_newmark_displ_struct->polySecondDerivCoefficient()*idv(rho)*inner(idv(u),id(v)),
                           _geomap=this->geomap() );
        }
        if (BuildCstPart)
        {
            auto polySecondDerivDisp = M_newmark_displ_struct->polySecondDeriv();
            linearFormDisplacement +=
                integrate( _range=elements(mesh),
                           _expr= -idv(rho)*inner(idv(polySecondDerivDisp),id(v)),
                           _geomap=this->geomap() );
        }
    }
    //--------------------------------------------------------------------------------------------------//
    // neumann bc
    if (BuildCstPart)
    {
        this->updateBCNeumannResidual( R );
    }
    // follower pressure bc
    if (!BuildCstPart)
    {
        this->updateBCFollowerPressureResidual( u,R );
    }

    //--------------------------------------------------------------------------------------------------//
    // fsi bc
    if (this->getMarkerNameFSI().size()>0)
    {
        // neumann boundary condition with normal stress (fsi boundary condition)
        if (BuildCstPart_BoundaryParoiMobile)
        {
            linearFormDisplacement +=
                integrate( _range=markedfaces(mesh,this->getMarkerNameFSI()),
                           _expr= alpha_f*trans(idv(*M_normalStressFromFluid))*id(v),
                           _geomap=this->geomap() );
        }

        if ( this->couplingFSIcondition() == "robin" )
        {
            double gammaRobinFSI = this->gammaNitschFSI();//2500;//10;
            double muFluid = this->muFluidFSI();//0.03;
            if (!BuildCstPart && !UseJacobianLinearTerms)
            {
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->getMarkerNameFSI()),
                               _expr= gammaRobinFSI*muFluid*M_newmark_displ_struct->polyFirstDerivCoefficient()*inner(idv(u),id(v))/hFace(),
                               _geomap=this->geomap() );
            }

            if (BuildCstPart )
            {
                auto const& polyFirstDerivDisp = M_newmark_displ_struct->polyFirstDeriv();
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->getMarkerNameFSI()),
                               _expr= -gammaRobinFSI*muFluid*inner(idv( polyFirstDerivDisp/*M_newmark_displ_struct->polyFirstDeriv()*/ ),id(v))/hFace(),
                               _geomap=this->geomap() );
            }
            if (BuildCstPart_BoundaryParoiMobile)
            {
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->getMarkerNameFSI()),
                               _expr= -gammaRobinFSI*muFluid*inner(idv(this->velocityInterfaceFromFluid()),id(v))/hFace(),
                               _geomap=this->geomap() );
            }
        } // robin-robin fsi

    }

    //--------------------------------------------------------------------------------------------------//
    // robin boundary condition (used in wavePressure3d as external tissue for arterial wall)
    if ( M_markerNameBCRobin.size() > 0 && !BuildCstPart && !UseJacobianLinearTerms)
    {
        double alpha_robin = 1e4;
        linearFormDisplacement +=
            /**/ integrate( _range=markedfaces(mesh,M_markerNameBCRobin),
                            _expr= alpha_robin*trans(idv(u))*id(v),
                            _geomap=this->geomap() );

#if 0
        ForEachBC( bcDef,cl::robin_vec,
                   form1( _test=M_Xh, _vector=R) +=
                   /**/ integrate( _range=markedfaces(mesh,PhysicalName),
                                   _expr= alpha_robin*trans(idv(u))*id(v),
                                   _geomap=this->geomap() ) );
#endif
    }
    // TODO up second membre
    /*if (BuildCstPart)
     {
     //this->updateBCRobinResidual( R );

     double p0_robin = 115000;
     ForEachBC( bcDef,cl::robin_vec,
     form1( _test=M_Xh, _vector=R) +=
     integrate( _range=markedfaces(mesh,PhysicalName),
     _expr= p0_robin*trans(vf::N())*id(v),
     _geomap=this->geomap() ) );

     }
     */
    //--------------------------------------------------------------------------------------------------//
    // dirichlet condition by elimination
    if (this->hasMarkerDirichletBCelimination() && !BuildCstPart)
    {
        this->updateBCDirichletStrongResidual( R );
    }

    double timeElapsed=thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateResidual",
                                         "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str()+
                                         "\n--------------------------------------------------",
                                         this->worldComm(),this->verboseAllProc());
#endif //FEELMODELS_SOLID_BUILD_RESIDUAL_CODE
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidualIncompressibilityTerms( element_displacement_type const& u, element_pressure_type const& p, vector_ptrtype& R) const
{
#if defined(FEELMODELS_SOLID_BUILD_RESIDUAL_CODE)
    using namespace Feel::vf;

    boost::mpi::timer thetimer;
    this->log("SolidMechanics","updateResidualIncompressibilityTerms", "start" );

    auto mesh = M_Xh->mesh();
    auto v = u;
    auto q = p;
    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();

    /*double alpha_f=M_genAlpha_alpha_f;
     double alpha_m=M_genAlpha_alpha_m;
     double gamma=0.5+alpha_m-alpha_f;
     double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);*/

    size_type rowStartInVector = this->rowStartInVector();
    size_type startDofIndexPressure = this->startDofIndexFieldsInMatrix().find("pressure")->second;
    auto linearFormDisplacement = form1( _test=M_Xh, _vector=R,
                                         _rowstart=rowStartInVector );
    auto linearFormPressure = form1( _test=M_XhPressure, _vector=R,
                                     _rowstart=rowStartInVector+startDofIndexPressure );

    if (M_pdeType=="Hyper-Elasticity")
    {
        linearFormDisplacement +=
            integrate( _range=elements(mesh),
                       //_expr= trace(val(-idv(p)*trans(InvFv))*trans(grad(v))),
                       //_expr=inner( -idv(p)*Feel::vf::FeelModels::solidMecIncompressibilityPressure(u,p,*this->mechanicalProperties()),grad(v) ),
                       _expr=inner( /*-idv(p)**/Feel::vf::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()),grad(v) ),
                       _geomap=this->geomap() );
    }
    else if (M_pdeType=="Elasticity-Large-Deformation" || M_pdeType=="Elasticity")
    {
        //Identity Matrix
        auto Id = eye<nDim,nDim>();
        linearFormDisplacement +=
            integrate( _range= elements(mesh),
                       _expr= inner( -idv(p)*Id ,grad(v)),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    linearFormPressure +=
        integrate( _range=elements(mesh),
                   //_expr=val(Fav22+Fav11+Fav11*Fav22-Fav21*Fav12)*id(q),
                   //_expr= detFvM1*id(q),
                   _expr= Feel::vf::FeelModels::solidMecPressureFormulationConstraint(u,*this->mechanicalProperties())*id(q),
                   _geomap=this->geomap() );


    if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
    {
        linearFormPressure +=
            integrate(_range=elements(mesh),
                      _expr= -(cst(1.)/idv(coeffLame1))*idv(p)*id(q),
                      _geomap=this->geomap() );
    }
    else
    {
        auto kappa = idv(this->mechanicalProperties()->fieldBulkModulus());
        linearFormPressure +=
            integrate( _range=elements(mesh),
                       _expr= -(cst(1.)/kappa)*idv(p)*id(q),
                       _geomap=this->geomap() );
    }


    double timeElapsed=thetimer.elapsed();
    this->log("SolidMechanics","updateResidualIncompressibilityTerms",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );

#endif //FEELMODELS_SOLID_BUILD_RESIDUAL_CODE
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidualViscoElasticityTerms( /*const*/ element_displacement_type& U, vector_ptrtype& R) const
{
#if 0
#if defined(FEELMODELS_SOLID_BUILD_RESIDUAL_CODE)
    using namespace Feel::vf;

    if (this->verbose()) std::cout << "[SolidMechanics] : updateResidualViscoElasticityTerms start\n";

    auto mesh = M_Xh->mesh();

    auto u = U.element<0>();
    auto v = U.element<0>();

    auto Ev2 = 0.5*(gradv(u)+trans(gradv(u)) );// + 0.5*trans(gradv(u))*gradv(u);


    double gammav=0.01;

    auto Buzz1bis = M_bdf_displ_struct->polyDeriv();
    auto buzz1bis = Buzz1bis.element<0>();


    auto Evbis = 0.5*(gradv(buzz1bis)+trans(gradv(buzz1bis)) );

    form1( _test=M_Xh, _vector=R ) +=
        integrate( _range=elements(mesh),
                   _expr= gammav*M_bdf_displ_struct->polyDerivCoefficient(0)*trace( Ev2*trans(grad(v))),
                   _geomap=this->geomap() );

    form1( _test=M_Xh, _vector=R ) +=
        integrate( _range=elements(mesh),
                   _expr= -gammav*trace( Evbis*trans(grad(v))),
                   _geomap=this->geomap() );

    if (this->verbose()) std::cout << "[SolidMechanics] : updateResidualViscoElasticityTerms finish\n";
#endif //FEELMODELS_SOLID_BUILD_RESIDUAL_CODE
#endif
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//


} // FeelModels
} // Feel



