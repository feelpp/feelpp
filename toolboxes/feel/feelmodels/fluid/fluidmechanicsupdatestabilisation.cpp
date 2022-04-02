
#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>
// #include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

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
#if 0 // VINCENT
        auto rho = idv(this->materialProperties()->fieldRho());
        //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
        auto mu = idv(this->materialProperties()->fieldMu());

        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            this->updateLinearPDEStabilisationGLS( data, rho, mu, rangeData.first );
#endif // VINCENT
    }

    //using namespace Feel::vf;
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("FluidMechanics","updateLinearPDEStabilisation", "start"+sc );
    boost::mpi::timer thetimer;

    //----------------------------------------------------------------------------------------------------//

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool BuildTermStabCIP = BuildNonCstPart;
    if (this->hasMeshMotion() /*this->useFSISemiImplicitScheme()*/)
    {
        BuildTermStabCIP = BuildCstPart;
    }

    //----------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto const& u = this->fieldVelocity();
    auto const& v = u;
    auto const& p = this->fieldPressure();
    auto const& q = p;

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    /*auto bilinearForm_PatternDefault = form2( _test=Xh,_trial=Xh,_matrix=A,
     _pattern=size_type(Pattern::DEFAULT),
     _rowstart=rowStartInMatrix,
     _colstart=colStartInMatrix );*/
    auto bilinearFormVV_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                                _pattern=size_type(Pattern::COUPLED),
                                                _rowstart=rowStartInMatrix,
                                                _colstart=colStartInMatrix );
    // a detail very important (A is cst or not)
    auto bilinearFormVV_PatternExtended = form2( _test=XhV,_trial=XhV,_matrix=A,
                                               _pattern=size_type(Pattern::EXTENDED),
                                               _rowstart=rowStartInMatrix,
                                               _colstart=colStartInMatrix );

    //----------------------------------------------------------------------------------------------------//
    // // stabilisation-div-div
    // if ( this->doStabDivDiv() && BuildCstPart  )
    // {
    //     double beta = doption(_name="stabilisation-div-div-beta",_prefix=this->prefix());
    //     bilinearFormVV_PatternCoupled +=
    //         integrate( _range=M_rangeMeshElements,
    //                    _expr=/*vf::h()**/beta*divt(u)*div(v),
    //                    _geomap=this->geomap() );
    // }
    //----------------------------------------------------------------------------------------------------//
    // stabilisation-cip-convection
    if ( this->doCIPStabConvection() && BuildTermStabCIP  )
    {
        double gamma = this->stabCIPConvectionGamma();
        double order_scaling = math::pow( double(nOrderVelocity), 3.5 );
        auto const& u_extrapoled = this->timeStepBDF()->poly();

        if ( this->hasMeshMotion() )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            auto cip_stab_coeff_ale_expr = (gamma/order_scaling)*abs( inner(idv(u_extrapoled)-idv(M_fieldMeshVelocityUsedWithStabCIP/*this->meshVelocity()*/),N() ) )*pow(hFace(),2.0);
            bilinearFormVV_PatternExtended +=
                integrate( _range=marked3faces(XhV->mesh(),1),
                           //_expr=val(cip_stab_coeff_ale_expr)*inner( jumpt(gradt(u)),jump(grad(v)) ),
                           _expr=val(cip_stab_coeff_ale_expr)*inner( jumpgradt(u), jumpgrad(v) ),
                           _geomap=this->geomap() );
#endif
        }
        else
        {
            auto cip_stab_coeff_expr = (gamma/order_scaling)*abs( trans(idv(u_extrapoled))*N() )*pow(hFace(),2.0);
            bilinearFormVV_PatternExtended +=
                integrate( _range=marked3faces(XhV->mesh(),1),
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
        auto const& u_extrapoled = this->timeStepBDF()->poly();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto muExpr = this->dynamicViscosityExpr( u_extrapoled/*, se*/ );
        auto Re = beta_abs*hFace()/muExpr;
        auto cip_div_coeff_expr = gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)*beta_abs;
        bilinearFormVV_PatternExtended +=
            integrate( _range=marked3faces(XhV->mesh(),1),
                       _expr=val(cip_div_coeff_expr)*inner( jumpt(divt(u)),jump(div(v)) ),
                       _geomap=this->geomap() );
    }
    //----------------------------------------------------------------------------------------------------//
    // stabilisation-cip-pressure
    if ( this->doCIPStabPressure() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPPressureGamma();
        auto const& u_extrapoled = this->timeStepBDF()->poly();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto muExpr = this->dynamicViscosityExpr( u_extrapoled/*, se*/ );
#if 1
        //double order_scaling = math::pow( double(nOrderPressure), 3.5 );
        //auto beta_abs = vf::sqrt( trans(idv(u_extrapoled))*idv(u_extrapoled));
        //auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, idv(mu)) )/cst(order_scaling);
        auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, muExpr) );
#else
        auto Re = beta_abs*hFace()/muExpr;
        auto cip_stab_coeff_expr =  gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)/max(beta_abs,cst(1e-6)); //last max in order to not divide by 0
#endif
        auto bilinearFormPP_PatternExtended = form2( _test=XhP,_trial=XhP,_matrix=A,
                                                     _pattern=size_type(Pattern::EXTENDED),
                                                     _rowstart=rowStartInMatrix+1,
                                                     _colstart=colStartInMatrix+1 );
        bilinearFormPP_PatternExtended +=
            integrate( _range=marked3faces(XhV->mesh(),1),
                       _expr=val(cip_stab_coeff_expr)*( jumpt(gradt(p))*jump(grad(q)) ),
                       _geomap=this->geomap() );
    }

    //----------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateLinearPDEStabilisation",(boost::format("finish in %1% s") % timeElapsed).str() );

} // updateLinearPDEStabilisation

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualStabilisation( DataUpdateResidual & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const
{
    if ( M_stabilizationGLS && M_stabilizationGLSDoAssembly )
    {
#if 0 // VINCENT
        auto rho = idv(this->materialProperties()->fieldRho());
        //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
        auto mu = idv(this->materialProperties()->fieldMu());

        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            this->updateResidualStabilisationGLS( data, u, p, rho, mu, rangeData.first );
#endif // VINCENT
    }

    this->log("FluidMechanics","updateResidualStabilisation", "start" );
    boost::mpi::timer thetimer;

    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool BuildTermStabCIP = !BuildCstPart;
    if ( this->hasMeshMotion() )
        BuildTermStabCIP = !BuildCstPart && !UseJacobianLinearTerms;

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto const& v = u;
    auto const& q = p;

    //--------------------------------------------------------------------------------------------------//

    size_type rowStartInVector = this->rowStartInVector();
    /*auto linearForm_PatternDefault = form1( _test=Xh, _vector=R,
     _pattern=size_type(Pattern::DEFAULT),
     _rowstart=rowStartInVector );*/
    auto linearFormV_PatternCoupled = form1( _test=XhV, _vector=R,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=rowStartInVector );
    auto linearFormV_PatternExtended = form1( _test=XhV, _vector=R,
                                              _pattern=size_type(Pattern::EXTENDED),
                                              _rowstart=rowStartInVector );

    // if (!BuildCstPart && !UseJacobianLinearTerms )
    // {
    //     if ( this->doStabDivDiv()  )
    //     {
    //         double beta = doption(_name="stabilisation-div-div-beta",_prefix=this->prefix());
    //         linearFormV_PatternCoupled +=
    //             integrate( _range=M_rangeMeshElements,
    //                        _expr=/*vf::h()**/beta*divv(u)*div(v),
    //                        _geomap=this->geomap() );
    //     }
    // }

    // stabilisation-cip-convection
    if ( this->doCIPStabConvection() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPConvectionGamma();
        auto order_scaling = std::pow( double(nOrderVelocity), 3.5 );
        auto const& u_extrapoled = this->timeStepBDF()->poly();

        if ( this->hasMeshMotion() )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            auto cip_stab_coeff_ale_expr = gamma*abs( trans(idv(u_extrapoled)-idv(M_fieldMeshVelocityUsedWithStabCIP/*this->meshVelocity()*/) )*N() )*pow(hFace(),2.0)/cst(order_scaling);
            linearFormV_PatternExtended +=
                integrate( _range=marked3faces(mesh,1),
                           //_expr= val( cip_stab_coeff_ale_expr*trans(jumpv(gradv(u))) )*jump(grad(v)), // this line not work!
                           _expr= val( cip_stab_coeff_ale_expr )*inner( jumpv(gradv(u)),jump(grad(v)) ),
                           _geomap=this->geomap() );
#endif
        }
        else
        {
            auto cip_stab_coeff_expr = gamma*abs( trans(idv(u_extrapoled) )*N() )*pow(hFace(),2.0)/cst(order_scaling);
            linearFormV_PatternExtended +=
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
        auto const& u_extrapoled = this->timeStepBDF()->poly();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto muExpr = this->dynamicViscosityExpr( u/*, se*/ );
        auto Re = beta_abs*hFace()/muExpr;
        auto cip_div_coeff_expr = gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)*beta_abs;
        linearFormV_PatternExtended +=
            integrate( _range=marked3faces(XhV->mesh(),1),
                       _expr=val(cip_div_coeff_expr)*inner( jumpv(divv(u)),jump(div(v)) ),
                       _geomap=this->geomap() );
    }
    // stabilisation-cip-pressure
    if ( this->doCIPStabPressure() && BuildTermStabCIP )
    {
        auto linearFormP_PatternExtended = form1( _test=XhP, _vector=R,
                                                  _pattern=size_type(Pattern::EXTENDED),
                                                  _rowstart=rowStartInVector+1 );
        double gamma = this->stabCIPPressureGamma();
        auto const& u_extrapoled = this->timeStepBDF()->poly();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto muExpr = this->dynamicViscosityExpr( u/*, se*/ );
        auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, muExpr) );
        linearFormP_PatternExtended +=
            integrate( _range=marked3faces(XhV->mesh(),1),
                       _expr=val(cip_stab_coeff_expr)*( jumpv(gradv(p))*jump(grad(q)) ),
                       _geomap=this->geomap() );
    }

    //------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    this->log("FluidMechanics","updateResidualStabilisation",
              (boost::format("finish in %1% s") % timeElapsed).str() );
} // updateResidualStabilisation

//--------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianStabilisation( DataUpdateJacobian & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const
{
    if ( M_stabilizationGLS && M_stabilizationGLSDoAssembly )
    {
 #if 0 // VINCENT
        auto rho = idv(this->materialProperties()->fieldRho());
        //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
        auto mu = idv(this->materialProperties()->fieldMu());
        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            this->updateJacobianStabilisationGLS( data, u, p, rho, mu, rangeData.first );
#endif // VINCENT
    }

    this->log("FluidMechanics","updateJacobianStabilisation", "start" );
    boost::mpi::timer thetimer;

    //--------------------------------------------------------------------------------------------------//
    sparse_matrix_ptrtype& J = data.jacobian();
    bool BuildCstPart = data.buildCstPart();

    bool BuildNonCstPart = !BuildCstPart;
    bool BuildTermStabCIP = BuildNonCstPart;
    if ( this->hasMeshMotion() /*this->useFSISemiImplicitScheme()*/)
    {
        BuildTermStabCIP = BuildCstPart;
    }

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();
    auto const& v = u;
    auto const& q = p;

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    auto bilinearFormVV_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                                _pattern=size_type(Pattern::COUPLED),
                                                _rowstart=rowStartInMatrix,
                                                _colstart=colStartInMatrix );
    auto bilinearFormVV_PatternExtended = form2( _test=XhV,_trial=XhV,_matrix=J,
                                                 _pattern=size_type(Pattern::EXTENDED),
                                                 _rowstart=rowStartInMatrix,
                                                 _colstart=colStartInMatrix );

    //--------------------------------------------------------------------------------------------------//
    // if ( BuildCstPart && this->doStabDivDiv() )
    // {
    //     double beta = doption(_name="stabilisation-div-div-beta",_prefix=this->prefix());
    //     bilinearFormVV_PatternCoupled +=
    //         integrate( _range=elements(XhV->mesh()),
    //                    _expr=/*vf::h()**/beta*divt(u)*div(v),
    //                    _geomap=this->geomap() );
    // }

    // stabilisation-cip-convection
    if ( this->doCIPStabConvection() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPConvectionGamma();
        auto order_scaling = std::pow( double(nOrderVelocity), 3.5 );
        auto const& u_extrapoled = this->timeStepBDF()->poly();
#if 0
        auto betaField = Xh->template functionSpace<0>()->element();
        if ( this->hasMeshMotion() )
        {
            // warning here !! extended element maybe!!
            betaField.on(_range=elements(mesh),_expr=idv(u_extrapoled)-idv(M_fieldMeshVelocityUsedWithStabCIP/*this->meshVelocity()*/) );
        }
        else
            betaField.add(1.,u_extrapoled );
#endif

        if ( this->hasMeshMotion() )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            auto cip_stab_coeff_ale_expr = (gamma/order_scaling)*abs( trans(idv(u_extrapoled)-idv(M_fieldMeshVelocityUsedWithStabCIP/*this->meshVelocity()*/))*N() )*pow(hFace(),2.0);
            bilinearFormVV_PatternExtended +=
                integrate( _range=marked3faces(XhV->mesh(),1),
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
            bilinearFormVV_PatternExtended +=
                integrate( _range=marked3faces(XhV->mesh(),1),
                           //_expr= val(cip_stab_coeff_expr)*inner( jumpt(gradt(u)),jump(grad(v)) ),
                           _expr=val(cip_stab_coeff_expr)*inner( jumpgradt(u), jumpgrad(v) ),
                           _geomap=this->geomap() );
        }
    }

    // stabilisation-cip-divergence
    if ( this->doCIPStabDivergence() && BuildTermStabCIP  )
    {
        double gamma = this->stabCIPDivergenceGamma();
        auto const& u_extrapoled = this->timeStepBDF()->poly();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto muExpr = this->dynamicViscosityExpr( u/*, se*/ );
        auto Re = beta_abs*hFace()/muExpr;
        auto cip_div_coeff_expr = gamma*min(cst(1.),Re)*vf::pow(hFace(),2.0)*beta_abs;
        bilinearFormVV_PatternExtended +=
            integrate( _range=marked3faces(XhV->mesh(),1),
                       _expr=val(cip_div_coeff_expr)*inner( jumpt(divt(u)),jump(div(v)) ),
                       _geomap=this->geomap() );
    }
    // stabilisation-cip-pressure
    if ( this->doCIPStabPressure() && BuildTermStabCIP )
    {
        double gamma = this->stabCIPPressureGamma();
        auto const& u_extrapoled = this->timeStepBDF()->poly();
        auto beta_abs = vf::norm2(idv(u_extrapoled));
        auto muExpr = this->dynamicViscosityExpr( u/*, se*/ );
        auto cip_stab_coeff_expr =  gamma*( vf::pow(hFace(),3.0) / vf::max( hFace()*beta_abs, muExpr) );

        auto bilinearFormPP_PatternExtended = form2( _test=XhP,_trial=XhP,_matrix=J,
                                                     _pattern=size_type(Pattern::EXTENDED),
                                                     _rowstart=rowStartInMatrix+1,
                                                     _colstart=colStartInMatrix+1 );
        bilinearFormPP_PatternExtended +=
            integrate( _range=marked3faces(XhV->mesh(),1),
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
