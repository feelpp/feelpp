
#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{


#define jumpgradt( u ) ::Feel::vf::leftfacet(dnt(u)) + ::Feel::vf::rightfacet(dnt(u))
#define jumpgrad( u ) ::Feel::vf::leftface(dn(u)) + ::Feel::vf::rightface(dn(u))


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilisation( DataUpdateLinear & data ) const
{
    if ( M_stabilizationGLS && M_stabilizationGLSDoAssembly )
    {
        auto rho = idv(this->materialProperties()->fieldRho());
        //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
        auto mu = idv(this->materialProperties()->fieldMu());

        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            this->updateLinearPDEStabilisationGLS( data, rho, mu, rangeData.first );
    }

    //using namespace Feel::vf;
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();
    sparse_matrix_ptrtype& A_extended = data.matrixExtended();
    bool _BuildExtendedPart = data.buildExtendedPart();

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("FluidMechanics","updateLinearPDEStabilisation", "start"+sc );
    boost::mpi::timer thetimer;

    //----------------------------------------------------------------------------------------------------//

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;
    bool BuildExtendedPart = _BuildExtendedPart;

    bool BuildTermStabCIP = BuildNonCstPart;
    if (this->isMoveDomain() /*this->useFSISemiImplicitScheme()*/)
    {
        BuildTermStabCIP = BuildCstPart;
    }

    //----------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    auto U = Xh->element("u");
    auto V = Xh->element("v");
    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    /*auto bilinearForm_PatternDefault = form2( _test=Xh,_trial=Xh,_matrix=A,
     _pattern=size_type(Pattern::DEFAULT),
     _rowstart=rowStartInMatrix,
     _colstart=colStartInMatrix );*/
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    // a detail very important (A is cst or not)
    sparse_matrix_ptrtype A_extendedReal = (this->buildMatrixPrecond())?A_extended:A;
    auto bilinearForm_PatternExtended = form2( _test=Xh,_trial=Xh,_matrix=A_extendedReal/*A_extended*/,
                                               _pattern=size_type(Pattern::EXTENDED),
                                               _rowstart=rowStartInMatrix,
                                               _colstart=colStartInMatrix );

    //----------------------------------------------------------------------------------------------------//
    // stabilisation-div-div
    if ( this->doStabDivDiv() && BuildCstPart  )
    {
        double beta = doption(_name="stabilisation-div-div-beta",_prefix=this->prefix());
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(Xh->mesh()),
                       _expr=/*vf::h()**/beta*divt(u)*div(v),
                       _geomap=this->geomap() );
    }
    //----------------------------------------------------------------------------------------------------//
    // stabilisation-cip-convection
    if ( this->doCIPStabConvection() && BuildTermStabCIP  )
    {
        double gamma = this->stabCIPConvectionGamma();
        double order_scaling = math::pow( double(nOrderVelocity), 3.5 );
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();

        if (M_isMoveDomain)
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            auto cip_stab_coeff_ale_expr = (gamma/order_scaling)*abs( inner(idv(u_extrapoled)-idv(this->meshVelocity()),N() ) )*pow(hFace(),2.0);
            bilinearForm_PatternExtended +=
                integrate( _range=marked3faces(Xh->mesh(),1),
                           //_expr=val(cip_stab_coeff_ale_expr)*inner( jumpt(gradt(u)),jump(grad(v)) ),
                           _expr=val(cip_stab_coeff_ale_expr)*inner( jumpgradt(u), jumpgrad(v) ),
                           _geomap=this->geomap() );
#endif
        }
        else
        {
            auto cip_stab_coeff_expr = (gamma/order_scaling)*abs( trans(idv(u_extrapoled))*N() )*pow(hFace(),2.0);
            bilinearForm_PatternExtended +=
                integrate( _range=marked3faces(Xh->mesh(),1),
                           //_expr=val(cip_stab_coeff_expr)*inner( jumpt(gradt(u)),jump(grad(v)) ),
                           _expr=val(cip_stab_coeff_expr)*inner( jumpgradt(u), jumpgrad(v) ),
                           _geomap=this->geomap() );
        }
    }
    //----------------------------------------------------------------------------------------------------//
    // stabilisation-cip-divergence
    if ( this->doCIPStabDivergence() && BuildTermStabCIP  )
    {
        double gamma = this->stabCIPDivergenceGamma();
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto Re = beta_abs*hFace()/idv(mu);
        auto cip_div_coeff_expr = gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)*beta_abs;
        bilinearForm_PatternExtended +=
            integrate( _range=marked3faces(Xh->mesh(),1),
                       _expr=val(cip_div_coeff_expr)*inner( jumpt(divt(u)),jump(div(v)) ),
                       _geomap=this->geomap() );
    }
    //----------------------------------------------------------------------------------------------------//
    // stabilisation-cip-pressure
    if ( this->doCIPStabPressure() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPPressureGamma();
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
#if 1
        //double order_scaling = math::pow( double(nOrderPressure), 3.5 );
        //auto beta_abs = vf::sqrt( trans(idv(u_extrapoled))*idv(u_extrapoled));
        //auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, idv(mu)) )/cst(order_scaling);
        auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, idv(mu)) );
#else
        auto Re = beta_abs*hFace()/idv(mu);
        auto cip_stab_coeff_expr =  gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)/max(beta_abs,cst(1e-6)); //last max in order to not divide by 0
#endif
        bilinearForm_PatternExtended +=
            integrate( _range=marked3faces(Xh->mesh(),1),
                       _expr=val(cip_stab_coeff_expr)*( jumpt(gradt(p))*jump(grad(q)) ),
                       _geomap=this->geomap() );
    }

    //----------------------------------------------------------------------------------------------------//


#if 0
    if (this->application()->vm()[prefixvm(this->prefix(),"stabilisation-divergence")].as<bool>() ) //stab div
    {
        double alphaPSPG = 0.1;
        auto coeffPSPG = alphaPSPG*h()*h()/(2*idv(M_P0Nu));
#if 0
        form2( Xh_Loc_p, Xh_Loc_p, A_Lp_Lp ) += integrate( elements(meshLoc), coeffPSPG*grad(qLoc)*trans(gradt(qLoc)) );
        form2( Xh_Loc_p, Xh_Loc_ux, A_Lp_Lux ) += integrate( elements(meshLoc), - coeffPSPG*dx(qLoc)*trace(hesst(uxLoc)) );
        form2( Xh_Loc_p, Xh_Loc_uy, A_Lp_Luy ) += integrate( elements(meshLoc), - coeffPSPG*dy(qLoc)*trace(hesst(uyLoc)) );

        form2( Xh_Loc_ux, Xh_Loc_p, A_Lux_Lp ) += integrate( elements(meshLoc), -coeffPSPG*trace(hess(vxLoc))*dxt(qLoc));
        form2( Xh_Loc_uy, Xh_Loc_p, A_Lux_Lp ) += integrate( elements(meshLoc), -coeffPSPG*trace(hess(vyLoc))*dyt(qLoc));

        form2( Xh_Loc_ux, Xh_Loc_ux, A_Lux_Lux ) += integrate( elements(meshLoc), coeffPSPG*trace(hess(vxLoc))*trace(hesst(uxLoc)) );
        form2( Xh_Loc_uy, Xh_Loc_uy, A_Luy_Luy ) += integrate( elements(meshLoc), coeffPSPG*trace(hess(vyLoc))*trace(hesst(uyLoc)) );
#else

        bilinearForm_PatternCoupled +=
            integrate( elements(mesh), coeffPSPG*div(v)*divt(u) );

    }
#endif
#endif


    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateLinearPDEStabilisation",(boost::format("finish in %1% s") % timeElapsed).str() );

} // updateLinearPDEStabilisation

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualStabilisation( DataUpdateResidual & data, element_fluid_external_storage_type const& U ) const
{
    if ( M_stabilizationGLS && M_stabilizationGLSDoAssembly )
    {
        auto rho = idv(this->materialProperties()->fieldRho());
        //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
        auto mu = idv(this->materialProperties()->fieldMu());

        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            this->updateResidualStabilisationGLS( data, U, rho, mu, rangeData.first );
    }

    using namespace Feel::vf;

    this->log("FluidMechanics","updateResidualStabilisation", "start" );
    boost::mpi::timer thetimer;

    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool BuildTermStabCIP = !BuildCstPart;
    if ( this->isMoveDomain() )
        BuildTermStabCIP = !BuildCstPart && !UseJacobianLinearTerms;

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    //auto U = Xh->element("u"); U = *XVec;
    //auto V = Xh->element("v");
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();

    //--------------------------------------------------------------------------------------------------//

    size_type rowStartInVector = this->rowStartInVector();
    /*auto linearForm_PatternDefault = form1( _test=Xh, _vector=R,
     _pattern=size_type(Pattern::DEFAULT),
     _rowstart=rowStartInVector );*/
    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=rowStartInVector );
    auto linearForm_PatternExtended = form1( _test=Xh, _vector=R,
                                             _pattern=size_type(Pattern::EXTENDED),
                                             _rowstart=rowStartInVector );

    if (!BuildCstPart && !UseJacobianLinearTerms )
    {
        if ( this->doStabDivDiv()  )
        {
            double beta = doption(_name="stabilisation-div-div-beta",_prefix=this->prefix());
            linearForm_PatternCoupled +=
                integrate( _range=elements(Xh->mesh()),
                           _expr=/*vf::h()**/beta*divv(u)*div(v),
                           _geomap=this->geomap() );
        }
    }

    // stabilisation-cip-convection
    if ( this->doCIPStabConvection() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPConvectionGamma();
        auto order_scaling = std::pow( double(nOrderVelocity), 3.5 );
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();

        if (M_isMoveDomain)
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            auto cip_stab_coeff_ale_expr = gamma*abs( trans(idv(u_extrapoled)-idv(this->meshVelocity()) )*N() )*pow(hFace(),2.0)/cst(order_scaling);
            linearForm_PatternExtended +=
                integrate( _range=marked3faces(mesh,1),
                           //_expr= val( cip_stab_coeff_ale_expr*trans(jumpv(gradv(u))) )*jump(grad(v)), // this line not work!
                           _expr= val( cip_stab_coeff_ale_expr )*inner( jumpv(gradv(u)),jump(grad(v)) ),
                           _geomap=this->geomap() );
#endif
        }
        else
        {
            auto cip_stab_coeff_expr = gamma*abs( trans(idv(u_extrapoled) )*N() )*pow(hFace(),2.0)/cst(order_scaling);
            linearForm_PatternExtended +=
                integrate( _range=marked3faces(mesh,1),
                           //_expr= val( cip_stab_coeff_expr*trans(jumpv(gradv(u))) )*jump(grad(v)), // this line not work!
                           _expr= val( cip_stab_coeff_expr)*inner( jumpv(gradv(u)),jump(grad(v)) ),
                           _geomap=this->geomap() );
        }
    }
    // stabilisation-cip-divergence
    if ( this->doCIPStabDivergence() && BuildTermStabCIP  )
    {
        double gamma = this->stabCIPDivergenceGamma();
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto Re = beta_abs*hFace()/idv(mu);
        auto cip_div_coeff_expr = gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)*beta_abs;
        linearForm_PatternExtended +=
            integrate( _range=marked3faces(Xh->mesh(),1),
                       _expr=val(cip_div_coeff_expr)*inner( jumpv(divv(u)),jump(div(v)) ),
                       _geomap=this->geomap() );
    }
    // stabilisation-cip-pressure
    if ( this->doCIPStabPressure() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPPressureGamma();
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, idv(mu)) );
        linearForm_PatternExtended +=
            integrate( _range=marked3faces(Xh->mesh(),1),
                       _expr=val(cip_stab_coeff_expr)*( jumpv(gradv(p))*jump(grad(q)) ),
                       _geomap=this->geomap() );
    }


    //------------------------------------------------------------------------------------//
#if 0 // STAB_SUPG_PSPG

    auto BBuzz = M_bdf_fluid->polyDeriv();
    auto Bbuzz = BBuzz.element<0>();

    auto supg_c1 = cst(30.);
    auto pspg_c = cst(100.);

    auto u_norm = max(cst(1e-16),sqrt(trans(idv(u))*idv(u)));
    auto Reh = u_norm*vf::h()*idv(M_P0Rho)/(2*idv(mu));
    auto chiSupg = supg_c1/(idv(M_P0Rho)*sqrt( cst(1.) + (cst(3.)/Reh)*(cst(3.)/Reh) ) );
    //auto chiSupg =max(cst(0.),min(Reh/cst(6.),cst(1.) ) );
    //auto chiSupg =chi(Reh>=cst(0.) && Reh <=cst(3.))*Reh/cst(3.) + chi(Reh>cst(3.));
    auto chiPspg = cst(1.)/(idv(M_P0Rho)*sqrt( cst(1.) + (cst(3.)*pspg_c/Reh)*(cst(3.)*pspg_c/Reh) ) );
    //auto chiPspg =max(cst(0.),min(Reh/cst(6.),cst(1.) ) );
    //auto chiPspg =chi(Reh>=cst(0.) && Reh <=cst(3.))*Reh/cst(3.) + chi(Reh>cst(3.));
    auto tau_supg = chi(u_norm>cst(1.e-16) )*chiSupg*vf::h()/(2*u_norm);
    auto tau_pspg = chi(u_norm>cst(1.e-16) )*chiPspg*vf::h()/(2*u_norm);

    auto Sigmav_11 = trans(Sigmav*oneX())*oneX();
    auto Sigmav_12 = trans(Sigmav*oneY())*oneX();
    auto Sigmav_21 = trans(Sigmav*oneX())*oneY();
    auto Sigmav_22 = trans(Sigmav*oneY())*oneY();
    auto SigmaProj1 = vf::project(M_Xh->functionSpace<0>(),elements(M_mesh), vec(Sigmav_11,Sigmav_12) );
    auto SigmaProj2 = vf::project(M_Xh->functionSpace<0>(),elements(M_mesh), vec(Sigmav_21,Sigmav_22) );
    auto divSigmav = trans(vec( divv(SigmaProj1),divv(SigmaProj2)) );
    auto time_supg = idv(M_P0Rho)*(trans(idv(u))*M_bdf_fluid->polyDerivCoefficient(0) - trans(idv(Bbuzz)));
    auto convec_supg = idv(M_P0Rho)*trans( gradv(u)*idv(u) );
    auto force_supg = -divSigmav;

    if (!BuildCstPart)
    {
        form1( Xh, R ) +=
            integrate (elements(mesh),
                       val(time_supg + convec_supg + force_supg)*(val(tau_supg)*grad(v)*idv(u)+val(tau_pspg)*trans(grad(q))) );
    }
#endif
    //------------------------------------------------------------------------------------//
#if 0 //STAB_FIC
    element_fluid_2_type c = U.element<2>();
    element_fluid_2_type d = V.element<2>();
    element_fluid_3_type pp = U.element<3>();
    element_fluid_3_type qq = V.element<3>();
    //auto ToTo=cst(1.)/(8*M_mu/(3*vf::h()*vf::h())+ 2*trans(idv(M_bdf_fluid->poly().element<0>()))*idv(M_bdf_fluid->poly().element<0>())/vf::h());
    auto ToTo=cst(1.)/(8*M_mu/(3*vf::h()*vf::h())+ 0./vf::h());
    auto hij=vec(vf::h()/2.,vf::h()/2.);
    form1( Xh, R ) += integrate (elements(mesh),M_rho*(trans(idv(u))*gradv(u)-trans(idv(c)))*grad(v)*hij);
    form1( Xh, R ) += integrate (elements(mesh),ToTo*grad(p)*(trans(gradv(p))-idv(pp)) );
    form1( Xh, R ) += integrate (elements(mesh),+M_rho*trans(id(c))*(idv(c)-gradv(u)*idv(u)) );
    form1( Xh, R ) += integrate (elements(mesh),ToTo*trans(id(pp))*(idv(pp)-trans(gradv(p))) );
#endif
    //------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateResidualStabilisation",
              (boost::format("finish in %1% s") % timeElapsed).str() );

} // updateResidualStabilisation

//--------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianStabilisation( DataUpdateJacobian & data, element_fluid_external_storage_type const& U ) const
{
    using namespace Feel::vf;

    if ( M_stabilizationGLS && M_stabilizationGLSDoAssembly )
    {
        auto rho = idv(this->materialProperties()->fieldRho());
        //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
        auto mu = idv(this->materialProperties()->fieldMu());
        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            this->updateJacobianStabilisationGLS( data, U, rho, mu, rangeData.first );
    }

    this->log("FluidMechanics","updateJacobianStabilisation", "start" );
    boost::mpi::timer thetimer;

    //--------------------------------------------------------------------------------------------------//
    sparse_matrix_ptrtype& J = data.jacobian();
    bool BuildCstPart = data.buildCstPart();

    bool BuildNonCstPart = !BuildCstPart;
    bool BuildTermStabCIP = BuildNonCstPart;
    if (this->isMoveDomain() /*this->useFSISemiImplicitScheme()*/)
    {
        BuildTermStabCIP = BuildCstPart;
    }

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    //auto U = Xh->element("u"); U = *XVec;
    //auto V = Xh->element("v");
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto bilinearForm_PatternExtended = form2( _test=Xh,_trial=Xh,_matrix=J,
                                               _pattern=size_type(Pattern::EXTENDED),
                                               _rowstart=rowStartInMatrix,
                                               _colstart=colStartInMatrix );

    //--------------------------------------------------------------------------------------------------//
    if ( BuildCstPart && this->doStabDivDiv() )
    {
        double beta = doption(_name="stabilisation-div-div-beta",_prefix=this->prefix());
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(Xh->mesh()),
                       _expr=/*vf::h()**/beta*divt(u)*div(v),
                       _geomap=this->geomap() );
    }

    // stabilisation-cip-convection
    if ( this->doCIPStabConvection() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPConvectionGamma();
        auto order_scaling = std::pow( double(nOrderVelocity), 3.5 );
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
#if 0
        auto betaField = Xh->template functionSpace<0>()->element();
        if ( this->isMoveDomain() )
        {
            // warning here !! extended element maybe!!
            betaField.on(_range=elements(mesh),_expr=idv(u_extrapoled)-idv(this->meshVelocity()) );
        }
        else
            betaField.add(1.,u_extrapoled );
#endif

        if ( this->isMoveDomain() )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            auto cip_stab_coeff_ale_expr = (gamma/order_scaling)*abs( trans(idv(u_extrapoled)-idv(this->meshVelocity()))*N() )*pow(hFace(),2.0);
            bilinearForm_PatternExtended +=
                integrate( _range=marked3faces(Xh->mesh(),1),
                           //_expr= val(cip_stab_coeff_ale_expr)*inner( jumpt(gradt(u)),jump(grad(v)) ),
                           _expr=val(cip_stab_coeff_ale_expr)*inner( jumpgradt(u), jumpgrad(v) ),
                           _geomap=this->geomap() );
#endif
        }
        else
        {
#if 1
            auto cip_stab_coeff_expr = (gamma/order_scaling)*abs( trans(idv(u_extrapoled))*N() )*pow(hFace(),2.0);
#else
            auto beta_abs = vf::norm2(idv(u_extrapoled));
            auto Re = beta_abs*hFace()/idv(mu);
            auto cip_stab_coeff_expr = gamma*min(cst(1.),Re)*abs( trans(idv(u_extrapoled))*N() )*pow(hFace(),2.0)/cst(order_scaling);
#endif
            bilinearForm_PatternExtended +=
                integrate( _range=marked3faces(Xh->mesh(),1),
                           //_expr= val(cip_stab_coeff_expr)*inner( jumpt(gradt(u)),jump(grad(v)) ),
                           _expr=val(cip_stab_coeff_expr)*inner( jumpgradt(u), jumpgrad(v) ),
                           _geomap=this->geomap() );
        }
    }

    // stabilisation-cip-divergence
    if ( this->doCIPStabDivergence() && BuildTermStabCIP  )
    {
        double gamma = this->stabCIPDivergenceGamma();
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto Re = beta_abs*hFace()/idv(mu);
        auto cip_div_coeff_expr = gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)*beta_abs;
        bilinearForm_PatternExtended +=
            integrate( _range=marked3faces(Xh->mesh(),1),
                       _expr=val(cip_div_coeff_expr)*inner( jumpt(divt(u)),jump(div(v)) ),
                       _geomap=this->geomap() );
    }
    // stabilisation-cip-pressure
    if ( this->doCIPStabPressure() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPPressureGamma();
        auto U_extrapoled = this->timeStepBDF()->poly();
        auto u_extrapoled = U_extrapoled.template element<0>();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, idv(mu)) );

        bilinearForm_PatternExtended +=
            integrate( _range=marked3faces(Xh->mesh(),1),
                       _expr=val(cip_stab_coeff_expr)*( jumpt(gradt(p))*jump(grad(q)) ),
                       _geomap=this->geomap() );
    }
    //----------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateJacobianStabilisation",
              (boost::format("finish in %1% s") % timeElapsed).str() );

} // updateJacobianStabilisation

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//


} // end namespace FeelModels
} // end namespace Feel
