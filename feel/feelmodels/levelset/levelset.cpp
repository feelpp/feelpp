#include "levelset.hpp"

namespace Feel {
namespace FeelModels {

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
LEVELSETBASE_CLASS_TEMPLATE_TYPE::LevelSet(
        mesh_ptrtype mesh, 
        std::string const& prefix, 
        PeriodicityType periodicityLS )
   :
    M_prefix(prefix),
    M_periodicity(periodicityLS),
    M_doUpdateMarkers(true)
{
    M_iterSinceReinit = 0;
    M_mass = 0.;

    // +++++++++++++  informations from options +++++++++
    reinitevery = ioption(prefixvm(prefix,"reinitevery"));
    enable_reinit = boption(prefixvm(prefix,"enable-reinit"));
    M_useMarker2AsMarkerDoneFmm = boption(prefixvm(prefix,"fm-use-markerdelta"));
    hj_max_iter = ioption(prefixvm(prefix,"hj-max-iter"));
    hj_dtau = doption(prefixvm(prefix,"hj-dtau"));
    hj_tol = doption(prefixvm(prefix,"hj-tol"));
    impose_inflow = ioption(prefixvm(prefix,"impose-inflow"));
    stabStrategy = ioption(prefixvm(prefix,"stabilization-strategy"));
    useRegularPhi = boption(_name=prefixvm(prefix,"use-regularized-phi"));
    hdNodalProj = boption(_name=prefixvm(prefix,"h-d-nodal-proj"));
    M_thicknessInterface=doption(prefixvm(prefix,"thickness-interface"));

    // -------------- Advection -------------
    const std::string _stabmeth = soption(_name=prefixvm(prefix,"advec-stab-method"));
    std::map<std::string, AdvectionStabMethod> m_stabmeth;
    m_stabmeth["NO"]=NO;
    m_stabmeth["CIP"]=CIP;
    m_stabmeth["GALS"]=GALS;
    m_stabmeth["SUPG"]=SUPG;
    m_stabmeth["SGS"]=SGS;
    CHECK(m_stabmeth.count(_stabmeth)) << _stabmeth <<" is not in the list of possible stabilization methods\n";
    M_advecstabmethod = m_stabmeth.find( _stabmeth )->second;


    // ----------- time discretization scheme ----------
    const std::string _tds = soption(prefixvm(prefix,"time-discr-scheme"));
    std::map<std::string, LevelSetTimeDiscretization> m_discr_method;
    m_discr_method["BDF2"]=BDF2;
#if defined(LEVELSET_CONSERVATIVE_ADVECTION)
    m_discr_method["CN_CONSERVATIVE"]=CN_CONSERVATIVE;
#endif
    m_discr_method["EU"]=EU;
    CHECK(m_discr_method.count(_tds)) << _tds << " is not in the list of possible time discretization\n";
    M_discrMethod = m_discr_method.find( _tds )->second;

    // --------------- mesh adaptation -----------------
#if defined (MESH_ADAPTATION)
    auto backend_mesh_adapt = backend_type::build(Environment::vm(), "mesh-adapt-backend");
    mesh_adapt.reset( new mesh_adaptation_type ( backend_mesh_adapt ));
#endif


    // ------------ reinitialization method -----------
    const std::string _rm = soption(prefixvm(prefix,"reinit-method"));
    CHECK( _rm == "hj" || _rm == "fm" ) << _rm << " is not a possible reinitialization method\n";
    if (_rm == "hj" )
        M_reinitmethod = HJ;
    else
        M_reinitmethod = FM;

    // in the case of fast marching method, strategy to set the first elements
    strategyBeforeFm = (strategyBeforeFm_type) ioption(prefixvm(prefix,"fm-init-first-elts-strategy"));

    initWithMesh(mesh);

    __iter=0;

    //M_phi = M_spaceLevelSet->elementPtr();
    //M_phio = M_spaceLevelSet->elementPtr();
    //M_phinl=M_spaceLevelSet->elementPtr();
    M_heaviside = M_spaceLevelSet->elementPtr();
    M_dirac = M_spaceLevelSet->elementPtr();

#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
    if (M_discrMethod==CN_CONSERVATIVE)
        {
            phic = M_spaceLSCorr->elementPtr();
        }
#endif

    reinitializerUpdated = false;

    this->levelsetInfos(true);

} //constructor

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::initWithMesh(mesh_ptrtype mesh)
{
    // +++++++++ initialize every quantities which need the mesh +++++++++++
    M_mesh = mesh;

    //  -------  spaces -------
    M_spaceP0 = spaceP0_type::New(_mesh=M_mesh, _periodicity=periodicity(NoPeriodicity()) );
    M_spaceLSVec = spaceLSVec_type::New(_mesh=M_mesh, _periodicity= periodicity(M_periodicity) );
    if (M_advecstabmethod == CIP)
        M_spaceLevelSet = space_levelset_type::New(_mesh=M_mesh, _periodicity= periodicity(M_periodicity), _extended_doftable=std::vector<bool>(1,true) );
    else
        M_spaceLevelSet = space_levelset_type::New(_mesh=M_mesh, _periodicity= periodicity(M_periodicity) );

    //  -------  backends ------
    M_backend_advrea = backend(_name="ls-advec");
    M_backend_smooth = backend(_name="ls-smooth");

    // ---------- advection -----------
    LOG(INFO)<<"start creating advection for level set"<<std::endl;
    M_advection = advection_type::New( M_spaceLevelSet, M_backend_advrea, M_advecstabmethod, option(prefixvm(M_prefix,"coeff-stab")).template as<double>() );
    LOG(INFO)<<"end creating advection for level set"<<std::endl;

    if ( (M_reinitmethod == HJ) || (strategyBeforeFm==HJ_EQ) )
    {
        LOG(INFO)<<"start creating advection for hamilton jacobi"<<std::endl;
        auto backend_reinit_hj = backend(_name="ls-hj-reinit");
        M_advection_hj = advection_type::New( M_spaceLevelSet, backend_reinit_hj, SUPG, option( prefixvm(M_prefix,"hj-coeff-stab")).template as<double>() );
        LOG(INFO)<<"end creating advection for hamilton jacobi"<<std::endl;
    }

    if (M_thicknessInterface > 1.)
        std::cout<<"\n\n------------ WARNING : thickness_interface > 1 ----------\n\n";

    //------------- initialize matrices for conservative advection-------------
#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
    if (M_discrMethod==CN_CONSERVATIVE)
        {
#if (LEVELSET_CONSERVATIVE_ADVECTION == 1)
            M_spaceLSCorr = spaceLSCorr_type::New(M_mesh);
#else
            M_spaceLSCorr = M_spaceLevelSet;
#endif

            // to solve non linear problem for correction
            backend_nl = backend(_name="ls-nl-advec");

            J = backend_nl->newMatrix(M_spaceLSCorr, M_spaceLSCorr);
            form2(M_spaceLSCorr, M_spaceLSCorr, J);

            R = backend_nl->newVector(M_spaceLSCorr);
            form1(M_spaceLSCorr, R);

            // to compute conservative heavyside function
            backend_h = backend(_name="ls-nl-advec-h");
            D_h = backend_h->newMatrix(M_spaceLevelSet, M_spaceLevelSet );
            form2( M_spaceLevelSet, M_spaceLevelSet, D_h);

            F_h = backend_h->newVector(M_spaceLevelSet);
            form1(M_spaceLevelSet, F_h);

            // might be an option (one more !)
            k_correction = (LEVELSET_CONSERVATIVE_ADVECTION==1) ? 0. : M_thicknessInterface;
        }
#endif

    // --------------- projectors ---------
    LOG(INFO)<<"start creating projectors for level set"<<std::endl;
    M_l2p = projector(M_spaceLevelSet, M_spaceLevelSet, backend(_name="ls-l2p") );
    M_l2pVec = projector(M_spaceLSVec, M_spaceLSVec, backend(_name="ls-l2pVec") );
    LOG(INFO)<<"end creating projectors for level set"<<std::endl;

} // init


LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void LEVELSETBASE_CLASS_TEMPLATE_TYPE::imposePhi( elementLS_ptrtype phi )
{
    using namespace Feel::vf;

    // impose phi and save the old one in phio
    M_phio= M_phi;
    M_phi = phi;
    updateHeavyside();
    updateDirac();
    updateMass();
}//imposePhi


LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void LEVELSETBASE_CLASS_TEMPLATE_TYPE::initialize(elementLS_ptrtype phi, bool doFirstReinit, ReinitMethod method, int max_iter, double dtau, double tol)
    {
        using namespace Feel;
        using namespace Feel::vf;
        using namespace Feel::levelset;

        //if phi is not a distance function, force to do a reinitialization at first
        *M_phi=*phi;

        if (doFirstReinit==true)
            {

              if ( (method == HJ) && (max_iter == -1) )
                {
                  // if default value, use the one given in levelset options for the first reinit (this choice is done automatically for the other reinit)
                      max_iter = hj_max_iter;
                      dtau=hj_dtau;
                      tol = hj_tol;
                }

              this->reinitialize(max_iter, dtau, tol, false, method);

            }

        *M_phio=*M_phi;

        updateHeavyside();
        updateDirac();
        updateMass();

    } //initialize


LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::elementLS_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::circleShape(double r0, double x0, double y0, double z0)
{
    auto shape = M_spaceLevelSet->elementPtr();

    *shape = vf::project(M_spaceLevelSet, elements(M_mesh),
                         sqrt((Px()-x0)*(Px()-x0)
                              + (Py()-y0)*(Py()-y0)
                              + (Pz()-z0)*(Pz()-z0)) - r0 );
    return shape;
} //circleShape

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::elementLS_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::ellipseShape(double a_ell, double b_ell, double x0, double y0, double z0)
{
    auto shape = M_spaceLevelSet->elementPtr();
    *shape = vf::project(M_spaceLevelSet, elements(M_mesh),
                         sqrt((Px()-x0)*(Px()-x0) / (a_ell * a_ell)
                              + (Py()-y0)*(Py()-y0) / (b_ell * b_ell)
                              + (Pz()-z0)*(Pz()-z0) / (b_ell * b_ell) ) - 1 );
    return shape;
}//ellipseShape

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::self_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::New( mesh_ptrtype mesh, std::string const& prefix, double TimeStep, PeriodicityType periodocity )
{
    self_ptrtype new_ls( new self_type( mesh, prefix, TimeStep, periodocity) );
    return new_ls;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::setStrategyBeforeFm( int strat )
{
    if (reinitializerUpdated)
        LOG(INFO)<<" !!!  WARNING !!! : setStrategyBeforeFm set after the fast marching has been actually initialized ! \n";
    strategyBeforeFm = (strategyBeforeFm_type) strat;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
Feel::levelset::strategyBeforeFm_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::getStrategyBeforeFm()
{
    return strategyBeforeFm;
}


LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::setUseMarker2AsMarkerDoneFmm(bool value)
{
    M_useMarker2AsMarkerDoneFmm = value;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::setThicknessInterface( double value )
{
    M_thicknessInterface = value;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update levelset-dependent functions
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateDirac()
{
    // derivative of Heavyside function
    auto eps = this->thicknessInterface();

    if (M_useRegularPhi)
    {
        auto psi = idv(M_phi) / sqrt( gradv(M_phi) * trans(gradv(M_phi)) );
        auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/(2*eps) *( 1 + cos(pi*psi/eps) )
            +
            vf::chi(psi>eps)*vf::constant(0.0);

        if (boption(_name=prefixvm(M_prefix,"h-d-nodal-proj")))
            *M_dirac = vf::project( M_spaceLevelSet, elements(M_mesh), D_expr );
        else
            *M_dirac = M_l2p->project(D_expr);
    }
    else
    {
        auto psi = idv(M_phi) ;
        auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/(2*eps) *( 1 + cos(pi*psi/eps) )
            +
            vf::chi(psi>eps)*vf::constant(0.0);

        if (boption(_name=prefixvm(M_prefix,"h-d-nodal-proj")))
            *M_dirac = vf::project( M_spaceLevelSet, elements(M_mesh), D_expr );
        else
            *M_dirac = M_l2p->project(D_expr);
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateHeaviside()
{ 
    auto eps = this->thicknessInterface();

    if (useRegularPhi)
    {
        auto psi = idv(M_phi) / vf::sqrt( gradv(M_phi) * trans(gradv(M_phi)) );
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/2*(1 + psi/eps + 1/pi*vf::sin( pi*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (boption(_name=prefixvm(M_prefix,"h-d-nodal-proj")))
            *M_heaviside = vf::project(M_spaceLevelSet, elements(M_mesh), H_expr);
        else
            *M_heaviside = M_l2p->project(H_expr);
    }
    else
    {
        auto psi = idv(M_phi);
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/2*(1 + psi/eps + 1/pi*vf::sin( pi*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (boption(_name=prefixvm(M_prefix,"h-d-nodal-proj")))
            *M_heaviside = vf::project(M_spaceLevelSet, elements(M_mesh), H_expr);
        else
            *M_heaviside = M_l2p->project(H_expr);
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMass()
{
    M_mass = integrate(
            _range=elements(M_mesh),
            _expr=(1-idv(this->heaviside())) ).evaluate()(0,0);
            //_expr=vf::chi( idv(M_phi)<0.0) ).evaluate()(0,0); // gives very noisy results
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Markers accessors
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerInterface()
{
    if( !M_markerInterface )
        M_markerInterface.reset( new element_markers_type(M_spaceMarkers, "MarkerInterface") );

    if( M_doUpdateMarkers )
       this->updateMarkerInterface(); 

    return M_markerInterface;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerDirac()
{
    if( !M_markerDirac )
        M_markerDirac.reset( new element_markers_type(M_spaceMarkers, "MarkerDirac") );

    if( M_doUpdateMarkers )
       this->updateMarkerDirac(); 

    return M_markerDirac;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerHeaviside()
{
    if( !M_markerHeaviside )
        M_markerHeaviside.reset( new element_markers_type(M_spaceMarkers, "MarkerHeaviside") );

    if( M_doUpdateMarkers )
       this->updateMarkerHeaviside(); 

    return M_markerHeaviside;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerCrossedElements()
{
    if( !M_markerCrossedElements )
        M_markerCrossedElements.reset( new element_markers_type(M_spaceMarkers, "MarkerCrossedElements") );

    if( M_doUpdateMarkers )
       this->updateMarkerCrossedElements(); 

    return M_markerCrossedElements;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Advection
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
template<typename ExprT>
bool 
LEVELSETBASE_CLASS_TEMPLATE_TYPE::advect(vf::Expr<ExprT> const& velocity)
{
    //output : true if reinitialized
    bool didReinit=false;

    /*
       update stabilization strategy
       1 - update stab every iteration
       2 - update stab when reinitialize
       3 - update stab every n iterations (not implemented yet)
       */

    if( M_iterSinceReinit >= M_advection->timeSchemeOrder()-1 )
    {
        M_advection->updateAdvectionVelocity(velocity);
        M_advection->solve();

        M_iterSinceReinit++;
    }
    else
    {
        M_advection->setTimeOrder( M_iterSinceReinit + 1 );
        M_advection->updateAdvectionVelocity(velocity);
        M_advection->solve();

        M_iterSinceReinit++;
    }

    /*//        if ( enable_reinit && (ForceReinit || doReinit() ))
    if (M_discrMethod != CN_CONSERVATIVE)
    {
        bool timeToReinit;
        if (reinitevery > 0)
            timeToReinit = (M_iterSinceReinit == 0) ? false : (M_iterSinceReinit%reinitevery)==0 ;
        else
        {
            double dtd = distToDist();
            std::cout<<"dtd = "<<dtd<<std::endl;
            double reinitif = option(prefixvm(M_prefix,"reinit-if-dist-smaller")).template as<double>();
            timeToReinit = dtd > reinitif ;
        }


        if ( enable_reinit && (ForceReinit || timeToReinit ) && (updateTime) )
        {

            if (!updateTime)
                M_phinl.swap(this->phi());
            this->reinitialize(hj_max_iter, hj_dtau, hj_tol, true);
            didReinit=true;
            if (!updateTime)
                M_phinl.swap(this->phi());
        }
    }*/
    updateDirac();
    updateHeavyside();
    updateMass();
    M_doUpdateMarkers = true;
}

//----------------------------------------------------------------------------//
// Reinitialization
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::reinitialize()
{ 
    if( !M_reinitializerIsUpdatedForUse )
        this->createReinitializer();

    if ( M_reinitMethod == FM )
    {
        bool useMarker2AsMarkerDoneFmm = M_useMarker2AsMarkerDoneFmm;

        elementP0_ptrtype markerdelta;
        if ( useMarker2AsMarkerDoneFmm )
        {
            markerdelta = markerDelta();
            M_mesh->updateMarker2( *markerdelta );
        }

        mesh_ptrtype mesh_reinit;
        auto phi_reinit = M_spaceReinitP1->element();
        auto phi = this->phi();

        switch (strategyBeforeFm)
        {
            case ILP :
            {
                // save the smoothed gradient magnitude of phi
                auto modgradphi = M_smooth->project( vf::min(vf::max(vf::sqrt(inner(gradv(phi), gradv(phi))), 0.92), 2.) );

                if (Order > 1)
                {
                    auto modgradphiP1 = M_spaceReinitP1->element();
                    M_opInterpolationLStoP1->apply( *phi, phi_reinit );
                    M_opInterpolationLStoP1->apply( modgradphi, modgradphiP1 );

                    mesh_reinit = M_opLagrangeP1->mesh();

                    phi_reinit = vf::project(M_spaceReinitP1, elements(mesh_reinit),
                                               idv(phi_reinit) / idv(modgradphiP1) );

                    if (useMarker2AsMarkerDoneFmm)
                    {
                        auto spaceP0OpLag = spaceP0_type::New( mesh_reinit );
                        auto mark2opLag = vf::project(spaceP0OpLag, elements(mesh_reinit), idv( markerdelta ) > 1e-6 );
                        mesh_reinit->updateMarker2( mark2opLag );
                    }
                }
                else
                {
                    phi_reinit = vf::project(M_spaceReinitP1, elements(M_mesh), idv(phi) / idv(modgradphi) );
                    mesh_reinit = M_mesh;
                }
            }
            break;

            case HJ_EQ :
            case NONE :
            {
                if (strategyBeforeFm == HJ_EQ)
                {
                    CHECK(false) << "TODO\n";
                    //*phi = *explicitHJ(max_iter, dtau, tol);
                }

                if (Order > 1)
                {
                    M_opInterpolationLStoP1->apply( *phi, phi_reinit );
                    mesh_reinit = M_opLagrangeP1->mesh();

                    if (useMarker2AsMarkerDoneFmm)
                    {
                        auto spaceP0OpLag = spaceP0_type::New( mesh_reinit );
                        auto mark2opLag = vf::project(spaceP0OpLag, elements(mesh_reinit), idv( markerdelta ) > 1e-6 );
                        mesh_reinit->updateMarker2( mark2opLag );
                    }
                }
                // if P1 periodic, project on non periodic space for reinit
                else if ( this->M_periodicity.isPeriodic() )
                    phi_reinit = vf::project(M_spaceReinitP1, elements(M_mesh), idv(phi) );
            }
            break;

            default:
            {
                CHECK(false)<<"no strategy chosen to initialize first elements before fast marching\n"
                            <<"please, consider setting the option fm-init-first-elts-strategy to 0, 1 or 2\n";
            }
            break;
        } // switch strategyBeforeFm


        // Fast Marching Method
        phi_reinit = M_reinitializerFMS->march(phi_reinit, useMarker2AsMarkerDoneFmm);

        if (Order > 1)
            M_opInterpolationP1toLS->apply( phi_reinit, *phi );
        else
            *phi = vf::project(this->functionSpace(), elements(M_mesh), idv(phi_reinit));

        LOG(INFO)<< "reinit with FMM done"<<std::endl;

    } // Fast Marching


    else if ( M_reinitMethod == HJ )
    {
        CHECK(false) << "TODO\n";
        //ch.restart();
        //*phi = *explicitHJ(max_iter, dtau, tol);
        //LOG(INFO)<<"reinit done in "<<ch.elapsed()<<" s\n";
    } //HJ explicit

    M_iterSinceReinit=0;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update markers
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMarkerInterface()
{
    /* returns a marker (P0_type) on the elements crossed by the levelset
       ie :
       - if the element has not phi on all the dof of the same sign -> the element is mark as 1
       - if phi on the dof of the element are of the same sign -> mark as 0
      */

    const int ndofv = space_levelset_type::fe_type::nDof;

    auto phi = this->phi();

    auto it_elt = M_mesh->beginElementWithProcessId(M_mesh->worldComm().localRank());
    auto en_elt = M_mesh->endElementWithProcessId(M_mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    for (; it_elt!=en_elt; it_elt++)
    {
        int nbplus = 0;
        int nbminus = 0;

        for (int j=0; j<ndofv ; j++)
        {
            if (phi->localToGlobal(it_elt->id(), j, 0) >= 0.)
                nbplus++;
            else
                nbminus++;
        }

        //if elt crossed by interface
        if ( (nbminus != ndofv) && (nbplus!=ndofv) )
            M_markerInterface->assign(it_elt->id(), 0, 0, 1);
        else
            M_markerInterface->assign(it_elt->id(), 0, 0, 0);
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMarkerDirac()
{
    const int ndofv = space_levelset_type::fe_type::nDof;

    auto it_elt = M_mesh->beginElementWithProcessId(M_mesh->worldComm().localRank());
    auto en_elt = M_mesh->endElementWithProcessId(M_mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    double dirac_cut = M_dirac->max() / 10.;

    for (; it_elt!=en_elt; it_elt++)
    {
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs( M_dirac->localToGlobal(it_elt->id(), j, 0) ) > dirac_cut )
            {
                mark_elt = true;
                break; //don't need to do the others dof
            }
        }
        if( mark_elt )
            M_markerDirac->assign(it_elt->id(), 0, 0, 1);
        else
            M_markerDirac->assign(it_elt->id(), 0, 0, 0);
    }
}//markerDelta


LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMarkerHeaviside(bool invert, bool cut_at_half)
{
    /* returns P0 element having :
    if invert == true : 1 on elements inside Heaviside function (where H is smaller than epsilon on at least 1 dof)
    if invert == false : 1 on elements where Heaviside function greater than epsilon
    if cut_in_out == true : marker = 1 where H>0.5 and 0 where H<0.5
    */

    const int ndofv = space_levelset_type::fe_type::nDof;

    auto it_elt = M_mesh->beginElementWithProcessId(M_mesh->worldComm().localRank());
    auto en_elt = M_mesh->endElementWithProcessId(M_mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    double cut;
    if (cut_at_half) cut = 0.5;
    else             cut = invert ? 1e-3 : 0.999;

    for (; it_elt!=en_elt; it_elt++)
    {
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs( M_heaviside->localToGlobal(it_elt->id(), j, 0) ) > cut )
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            M_markerHeaviside->assign(it_elt->id(), 0, 0, (invert)?0:1);
        else
            M_markerHeaviside->assign(it_elt->id(), 0, 0, (invert)?1:0);
    }
} 

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMarkerCrossedElements()
{
    // return a "marker" on the elements traversed by the interface between phio and phi
    // ie : mark as 1 is one of the dof of the element has  (phi * phio < 0)
    const int ndofv = space_levelset_type::fe_type::nDof;

    auto it_elt = M_mesh->beginElementWithProcessId(M_mesh->worldComm().localRank());
    auto en_elt = M_mesh->endElementWithProcessId(M_mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    auto phi = this->phi();
    auto phio = this->phio();

    auto prod = vf::project(M_spaceLevelSet, elements(M_mesh),
                            idv(phio) * idv(phi) );


    for (; it_elt!=en_elt; it_elt++)
    {
        bool mark_elt = false;
        for (int j=0; j<ndofv ; j++)
        {
            if (prod.localToGlobal(it_elt->id(), j, 0) <= 0.)
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            M_markerCrossedElements->assign(it_elt->id(), 0, 0, 1);
        else
            M_markerCrossedElements->assign(it_elt->id(), 0, 0, 0);
    }
}


} // FeelModels
} // Feel

//instantiation
#define FEELPP_INSTANTIATE_LEVELSET 1
#include "levelset_instance.hpp"
