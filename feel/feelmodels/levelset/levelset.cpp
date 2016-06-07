#include <feel/feelmodels/levelset/levelset.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel {
namespace FeelModels {

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
LEVELSET_CLASS_TEMPLATE_TYPE::LevelSet( 
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository ) 
:
    super_type( prefix, worldComm, subPrefix, rootRepository ),
    M_mass(0.),
    //M_periodicity(periodicityLS),
    M_doUpdateMarkers(true),
    M_reinitializerIsUpdatedForUse(false),
    M_iterSinceReinit(0)
{
    this->loadParametersFromOptionsVm();

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".LevelSetConstructor.data";
    this->addTimerTool("Constructor", nameFileConstructor);

    /*// --------------- mesh adaptation -----------------
#if defined (MESH_ADAPTATION)
    auto backend_mesh_adapt = backend_type::build(Environment::vm(), "mesh-adapt-backend");
    mesh_adapt.reset( new mesh_adaptation_type ( backend_mesh_adapt ));
#endif*/

    __iter=0;

/*#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
    if (M_discrMethod==CN_CONSERVATIVE)
    {
        phic = M_spaceLSCorr->elementPtr();
    }
#endif*/

    this->levelsetInfos(true);
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::self_ptrtype
LEVELSET_CLASS_TEMPLATE_TYPE::New(
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
{
    self_ptrtype new_ls( new self_type(prefix, worldComm, subPrefix, rootRepository) );
    return new_ls;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::build()
{
    this->createAdvection();
    this->createFunctionSpaces();
    this->createReinitialization();
    this->createExporters();
    this->createOthers();
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::build( mesh_ptrtype const& mesh )
{
    this->createAdvection( mesh );
    this->createFunctionSpaces();
    this->createReinitialization();
    this->createExporters();
    this->createOthers();
}

//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::init()
{
    M_advection->init();
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createAdvection()
{
    M_advection.reset( 
            new advection_type(
                this->prefix(), 
                this->worldComm(), 
                "advection", 
                this->rootRepositoryWithoutNumProc()
                ) );

    M_advection->build();
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createAdvection( mesh_ptrtype const& mesh )
{
    M_advection.reset( 
            new advection_type(
                this->prefix(), 
                this->worldComm(), 
                this->subPrefix(), 
                this->rootRepository()
                ) );

    M_advection->build( mesh );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    M_spaceLevelSetVec = space_levelset_vectorial_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    M_spaceMarkers = space_markers_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );

    // Init P1 reinitialization if required
    if( M_reinitMethod == LevelSetReinitMethod::FM )
    {
        if( Order > 1 )
        {
            M_opLagrangeP1 = lagrangeP1( this->functionSpace() );
            M_spaceReinitP1 = space_levelset_reinitP1_type::New( 
                    _mesh = M_opLagrangeP1->mesh(), 
                    _periodicity = periodicity(NoPeriodicity()) 
                    );
        }
        else
        {
            M_spaceReinitP1 = space_levelset_reinitP1_type::New( 
                    _mesh = this->mesh(), 
                    _periodicity = periodicity(NoPeriodicity()) 
                    );
        }
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createReinitialization()
{
    switch( M_reinitMethod )
    { 
        case LevelSetReinitMethod::FM :
        {
            //M_reinitializer.reset( 
                    //new ReinitializerFMS<space_levelset_reinitP1_type, periodicity_type>( M_spaceReinitP1, M_periodicity ) 
                    //);
            M_reinitializerFMS.reset( 
                    new ReinitializerFMS<space_levelset_reinitP1_type, periodicity_type>( M_spaceReinitP1, M_periodicity ) 
                    );
            if( Order > 1)
            {
                M_opInterpolationLStoP1 = opInterpolation(
                        _domainSpace = this->functionSpace(),
                        _imageSpace = this->functionSpaceReinitP1(),
                        _type = InterpolationNonConforme(false)
                        );
                M_opInterpolationP1toLS = opInterpolation(
                        _domainSpace = this->functionSpaceReinitP1(),
                        _imageSpace = this->functionSpace(),
                        _type = InterpolationNonConforme(false)
                        );
            }
        }
        break;
        case LevelSetReinitMethod::HJ :
        {
            // TODO
        }
        break;
    }

    M_reinitializerIsUpdatedForUse = true;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createExporters()
{
    this->log("LevelSet","createExporters", "start");
    this->timerTool("Constructor").start();

    std::string geoExportType="static";//change_coords_only, change, static
    std::string exporterPath = this->rootRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports"));
    M_exporter = exporter( _mesh=this->mesh(),
                           _name="Export",
                           _geo=geoExportType,
                           _path=exporterPath );

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("LevelSet","createExporters",(boost::format("finish in %1% s")%tElpased).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createOthers()
{
    M_heaviside.reset( new element_levelset_type(this->functionSpace(), "Heaviside") );
    M_dirac.reset( new element_levelset_type(this->functionSpace(), "Dirac") );
}

/*LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initWithMesh(mesh_ptrtype mesh)
{
    // +++++++++ initialize every quantities which need the mesh +++++++++++
    this->mesh() = mesh;

    //  -------  spaces -------
    M_spaceP0 = space_markers_type::New(_mesh=this->mesh(), _periodicity=periodicity(NoPeriodicity()) );
    M_spaceLSVec = spaceLSVec_type::New(_mesh=this->mesh(), _periodicity= periodicity(M_periodicity) );
    if (M_advecstabmethod == CIP)
        this->functionSpace() = space_levelset_type::New(_mesh=this->mesh(), _periodicity= periodicity(M_periodicity), _extended_doftable=std::vector<bool>(1,true) );
    else
        this->functionSpace() = space_levelset_type::New(_mesh=this->mesh(), _periodicity= periodicity(M_periodicity) );

    //  -------  backends ------
    M_backend_advrea = backend(_name="ls-advec");
    M_backend_smooth = backend(_name="ls-smooth");

    // ---------- advection -----------
    LOG(INFO)<<"start creating advection for level set"<<std::endl;
    M_advection = advection_type::New( this->functionSpace(), M_backend_advrea, M_advecstabmethod, option(prefixvm(this->prefix(),"coeff-stab")).template as<double>() );
    LOG(INFO)<<"end creating advection for level set"<<std::endl;

    if ( (M_reinitMethod == HJ) || (M_strategyBeforeFM==HJ_EQ) )
    {
        LOG(INFO)<<"start creating advection for hamilton jacobi"<<std::endl;
        auto backend_reinit_hj = backend(_name="ls-hj-reinit");
        M_advection_hj = advection_type::New( this->functionSpace(), backend_reinit_hj, SUPG, option( prefixvm(this->prefix(),"hj-coeff-stab")).template as<double>() );
        LOG(INFO)<<"end creating advection for hamilton jacobi"<<std::endl;
    }

    if (M_thicknessInterface > 1.)
        std::cout<<"\n\n------------ WARNING : thickness_interface > 1 ----------\n\n";

    //------------- initialize matrices for conservative advection-------------
#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
    if (M_discrMethod==CN_CONSERVATIVE)
        {
#if (LEVELSET_CONSERVATIVE_ADVECTION == 1)
            M_spaceLSCorr = spaceLSCorr_type::New(this->mesh());
#else
            M_spaceLSCorr = this->functionSpace();
#endif

            // to solve non linear problem for correction
            backend_nl = backend(_name="ls-nl-advec");

            J = backend_nl->newMatrix(M_spaceLSCorr, M_spaceLSCorr);
            form2(M_spaceLSCorr, M_spaceLSCorr, J);

            R = backend_nl->newVector(M_spaceLSCorr);
            form1(M_spaceLSCorr, R);

            // to compute conservative heavyside function
            backend_h = backend(_name="ls-nl-advec-h");
            D_h = backend_h->newMatrix(this->functionSpace(), this->functionSpace() );
            form2( this->functionSpace(), this->functionSpace(), D_h);

            F_h = backend_h->newVector(this->functionSpace());
            form1(this->functionSpace(), F_h);

            // might be an option (one more !)
            k_correction = (LEVELSET_CONSERVATIVE_ADVECTION==1) ? 0. : M_thicknessInterface;
        }
#endif

    // --------------- projectors ---------
    LOG(INFO)<<"start creating projectors for level set"<<std::endl;
    M_l2p = projector(this->functionSpace(), this->functionSpace(), backend(_name="ls-l2p") );
    M_l2pVec = projector(M_spaceLSVec, M_spaceLSVec, backend(_name="ls-l2pVec") );
    LOG(INFO)<<"end creating projectors for level set"<<std::endl;

} // init*/


/*LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void 
LEVELSET_CLASS_TEMPLATE_TYPE::imposePhi( elementLS_ptrtype phi )
{
    using namespace Feel::vf;

    // impose phi and save the old one in phio
    M_phio= this->phi();
    M_phi = phi;
    updateHeaviside();
    updateDirac();
    updateMass();
}//imposePhi*/


/*LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void LEVELSET_CLASS_TEMPLATE_TYPE::initialize(elementLS_ptrtype phi, bool doFirstReinit, ReinitMethod method, int max_iter, double dtau, double tol)
{
    using namespace Feel;
    using namespace Feel::vf;

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

    *M_phio=*this->phi();

    updateHeaviside();
    updateDirac();
    updateMass();

} //initialize*/

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::setStrategyBeforeFm( int strat )
{
    if (M_reinitializerIsUpdatedForUse)
        LOG(INFO)<<" !!!  WARNING !!! : setStrategyBeforeFm set after the fast marching has been actually initialized ! \n";
    M_strategyBeforeFM = (strategy_before_FM_type) strat;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    //M_enableReinit = boption(prefixvm(this->prefix(),"enable-reinit"));
    //M_reinitEvery = ioption(prefixvm(this->prefix(),"reinit-every"));
    M_useMarker2AsMarkerDoneFmm = boption(prefixvm(this->prefix(),"fm-use-markerdirac"));
    //hj_max_iter = ioption(prefixvm(this->prefix(),"hj-max-iter"));
    //hj_dtau = doption(prefixvm(this->prefix(),"hj-dtau"));
    //hj_tol = doption(prefixvm(this->prefix(),"hj-tol"));
    //impose_inflow = ioption(prefixvm(this->prefix(),"impose-inflow"));
    //stabStrategy = ioption(prefixvm(this->prefix(),"stabilization-strategy"));
    M_thicknessInterface=doption(prefixvm(this->prefix(),"thickness-interface"));
    M_useRegularPhi = boption(_name=prefixvm(this->prefix(),"use-regularized-phi"));
    M_useHeavisideDiracNodalProj = boption(_name=prefixvm(this->prefix(),"h-d-nodal-proj"));

    std::string reinitmethod = soption( _name="reinit-method", _prefix=this->prefix() );
    if( reinitmethod == "fm" )
        M_reinitMethod = LevelSetReinitMethod::FM;
    else if( reinitmethod == "hj" )
        M_reinitMethod = LevelSetReinitMethod::HJ;
    else
        CHECK( false ) << reinitmethod << " is not a valid reinitialization method\n";

    M_strategyBeforeFM = (strategy_before_FM_type) ioption(prefixvm(this->prefix(),"fm-init-first-elts-strategy"));

    M_doExportAdvection = boption(_name="export-advection", _prefix=this->prefix());
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update levelset-dependent functions
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateDirac()
{
    // derivative of Heaviside function
    auto eps = this->thicknessInterface();

    if (M_useRegularPhi)
    {
        auto psi = idv(this->phi()) / sqrt( gradv(this->phi()) * trans(gradv(this->phi())) );
        auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/(2*eps) *( 1 + cos(pi*psi/eps) )
            +
            vf::chi(psi>eps)*vf::constant(0.0);

        if ( M_useHeavisideDiracNodalProj )
            *M_dirac = vf::project( this->functionSpace(), elements(this->mesh()), D_expr );
        else
            *M_dirac = M_l2p->project(D_expr);
    }
    else
    {
        auto psi = idv(this->phi()) ;
        auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/(2*eps) *( 1 + cos(pi*psi/eps) )
            +
            vf::chi(psi>eps)*vf::constant(0.0);

        if (M_useHeavisideDiracNodalProj)
            *M_dirac = vf::project( this->functionSpace(), elements(this->mesh()), D_expr );
        else
            *M_dirac = M_l2p->project(D_expr);
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateHeaviside()
{ 
    auto eps = this->thicknessInterface();

    if (M_useRegularPhi)
    {
        auto psi = idv(this->phi()) / vf::sqrt( gradv(this->phi()) * trans(gradv(this->phi())) );
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/2*(1 + psi/eps + 1/pi*vf::sin( pi*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (M_useHeavisideDiracNodalProj)
            *M_heaviside = vf::project(this->functionSpace(), elements(this->mesh()), H_expr);
        else
            *M_heaviside = M_l2p->project(H_expr);
    }
    else
    {
        auto psi = idv(this->phi());
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/2*(1 + psi/eps + 1/pi*vf::sin( pi*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (M_useHeavisideDiracNodalProj)
            *M_heaviside = vf::project(this->functionSpace(), elements(this->mesh()), H_expr);
        else
            *M_heaviside = M_l2p->project(H_expr);
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMass()
{
    M_mass = integrate(
            _range=elements(this->mesh()),
            _expr=(1-idv(this->heaviside())) ).evaluate()(0,0);
            //_expr=vf::chi( idv(this->phi())<0.0) ).evaluate()(0,0); // gives very noisy results
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Markers accessors
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerInterface()
{
    if( !M_markerInterface )
        M_markerInterface.reset( new element_markers_type(M_spaceMarkers, "MarkerInterface") );

    if( M_doUpdateMarkers )
       this->updateMarkerInterface(); 

    return M_markerInterface;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerDirac()
{
    if( !M_markerDirac )
        M_markerDirac.reset( new element_markers_type(M_spaceMarkers, "MarkerDirac") );

    if( M_doUpdateMarkers )
       this->updateMarkerDirac(); 

    return M_markerDirac;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerHeaviside(bool invert, bool cut_at_half)
{
    if( !M_markerHeaviside )
        M_markerHeaviside.reset( new element_markers_type(M_spaceMarkers, "MarkerHeaviside") );

    //if( M_doUpdateMarkers )
       this->updateMarkerHeaviside( invert, cut_at_half );

    return M_markerHeaviside;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerCrossedElements()
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
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
template<typename ExprT>
void 
LEVELSET_CLASS_TEMPLATE_TYPE::advect(vf::Expr<ExprT> const& velocity)
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
            double reinitif = option(prefixvm(this->prefix(),"reinit-if-dist-smaller")).template as<double>();
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
    updateHeaviside();
    updateMass();
    M_doUpdateMarkers = true;
}

//----------------------------------------------------------------------------//
// Reinitialization
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::reinitialize()
{ 
    if( !M_reinitializerIsUpdatedForUse )
        this->createReinitialization();

    if ( M_reinitMethod == LevelSetReinitMethod::FM )
    {
        bool useMarker2AsMarkerDoneFmm = M_useMarker2AsMarkerDoneFmm;

        if ( useMarker2AsMarkerDoneFmm )
        {
            this->mesh()->updateMarker2( *this->markerDirac() );
        }

        mesh_ptrtype mesh_reinit;
        auto phi_reinit = M_spaceReinitP1->element();
        auto phi = this->phi();

        switch (M_strategyBeforeFM)
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
                        auto spaceP0OpLag = space_markers_type::New( mesh_reinit );
                        auto mark2opLag = vf::project(spaceP0OpLag, elements(mesh_reinit), idv( this->markerDirac() ) > 1e-6 );
                        mesh_reinit->updateMarker2( mark2opLag );
                    }
                }
                else
                {
                    phi_reinit = vf::project(M_spaceReinitP1, elements(this->mesh()), idv(phi) / idv(modgradphi) );
                    mesh_reinit = this->mesh();
                }
            }
            break;

            case HJ_EQ :
            case NONE :
            {
                if (M_strategyBeforeFM == HJ_EQ)
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
                        auto spaceP0OpLag = space_markers_type::New( mesh_reinit );
                        auto mark2opLag = vf::project(spaceP0OpLag, elements(mesh_reinit), idv( this->markerDirac() ) > 1e-6 );
                        mesh_reinit->updateMarker2( mark2opLag );
                    }
                }
                // if P1 periodic, project on non periodic space for reinit
                else if ( this->M_periodicity.isPeriodic() )
                    phi_reinit = vf::project(M_spaceReinitP1, elements(this->mesh()), idv(phi) );
            }
            break;

            default:
            {
                CHECK(false)<<"no strategy chosen to initialize first elements before fast marching\n"
                            <<"please, consider setting the option fm-init-first-elts-strategy to 0, 1 or 2\n";
            }
            break;
        } // switch M_strategyBeforeFM


        // Fast Marching Method
        phi_reinit = M_reinitializerFMS->march(phi_reinit, useMarker2AsMarkerDoneFmm);

        if (Order > 1)
            M_opInterpolationP1toLS->apply( phi_reinit, *phi );
        else
            *phi = vf::project(this->functionSpace(), elements(this->mesh()), idv(phi_reinit));

        LOG(INFO)<< "reinit with FMM done"<<std::endl;

    } // Fast Marching


    else if ( M_reinitMethod == LevelSetReinitMethod::HJ )
    {
        CHECK(false) << "TODO\n";
        //ch.restart();
        //*phi = *explicitHJ(max_iter, dtau, tol);
        //LOG(INFO)<<"reinit done in "<<ch.elapsed()<<" s\n";
    } //HJ explicit

    M_iterSinceReinit=0;
}

//----------------------------------------------------------------------------//
// Initial value
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::setInitialValue(element_levelset_ptrtype const& phiv, bool doReinitialize)
{
    *this->phi() = *phiv;

    if (doReinitialize)
        this->reinitialize();

    updateHeaviside();
    updateDirac();
    updateMass();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
std::string
LEVELSET_CLASS_TEMPLATE_TYPE::levelsetInfos( bool show )
{
    /* print some info when creating levelset instance */
    std::ostringstream infos;
    infos<< "\n||==============================================||"
         << "\n||----------Info :   LEVELSET    ---------------||"
         << "\n||==============================================||"
         << "\n   Prefix : " << this->prefix()
         << "\n   Dim : " << nDim
         << "\n   Order : " << Order
         << "\n   Periodicity : " << M_periodicity.isPeriodic()
         << "\n   Nb proc : " << this->worldComm().globalSize()
         << "\n\n  Level set parameters :"
         << "\n     -- thickness interface : " << this->thicknessInterface()
         << "\n     -- use regular phi (phi / |grad(phi)|) " << this->M_useRegularPhi
#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
         << "\n levelset built with conservative advection mode : "<<LEVELSET_CONSERVATIVE_ADVECTION
#else
         << "\n levelset built without conservative advection mode"
#endif
         << "\n     -- stabilization of advection : "<< soption( _name="advec-stab-method", _prefix=this->prefix() )
         //<< "\n     -- impose dirichlet at inflow : " << impose_inflow
         << "\n\n  Reinitialization parameters :"
         //<< "\n     -- enable reinitialization : "<<enable_reinit
         ;
    //if (enable_reinit)
    //{
        const std::string reinitmethod = soption( _name="reinit-method", _prefix=this->prefix() );
        infos << "\n     -- reinitialization method : " << reinitmethod;
        //if (reinitmethod == "hj")
        //{
            //infos << "\n      * hj maximum iteration per reinit : " << hj_max_iter
                  //<< "\n      * hj pseudo time step dtau : " << hj_dtau
                  //<< "\n      * hj stabilization : SUPG"
                  //<< "\n      * hj coeff stab : " << option( prefixvm(M_prefix,"hj-coeff-stab")).template as<double>()
                  //<< "\n      * hj tolerence on dist to dist error : "<<hj_tol;
        //}
        //else
        //{
            //infos << "\n      * fm smoothing coefficient for ILP : " << Environment::vm()[prefixvm(M_prefix,"fm-smooth-coeff")].template as< double >();
        //}
    //}
    //infos << "\n\n  Level set spaces :"
          //<< "\n     -- scalar LS space ndof : "<< this->functionSpace()->nDof()
          //<< "\n     -- vectorial LS ndof : "<< this->functionSpaceVectorial()->nDof()
          //<< "\n     -- scalar P0 space ndof : "<< this->functionSpaceMarkers()->nDof()
          //<<"\n||==============================================||\n\n";

    if (show)
        LOG(INFO)<<infos.str();

    return infos.str();
}

//----------------------------------------------------------------------------//
// Export results
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("LevelSet","exportResults", "start");
    this->timerTool("PostProcessing").start();
 
    if ( !M_exporter->doExport() ) return;

    M_exporter->step( time )->add( prefixvm(this->prefix(),"phi"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"phi")),
                                   this->phi() );
    M_exporter->step( time )->add( prefixvm(this->prefix(),"Dirac"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Dirac")),
                                   this->dirac() );
    M_exporter->save();

    if( M_doExportAdvection )
        M_advection->exportResults( time );

    this->timerTool("PostProcessing").stop("exportResults");
    this->log("LevelSet","exportResults", "finish");
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update markers
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerInterface()
{
    /* returns a marker (P0_type) on the elements crossed by the levelset
       ie :
       - if the element has not phi on all the dof of the same sign -> the element is mark as 1
       - if phi on the dof of the element are of the same sign -> mark as 0
      */

    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();
    auto phi = this->phi();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
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

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerDirac()
{
    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
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
}//markerDirac


LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerHeaviside(bool invert, bool cut_at_half)
{
    /* returns P0 element having :
    if invert == true : 1 on elements inside Heaviside function (where H is smaller than epsilon on at least 1 dof)
    if invert == false : 1 on elements where Heaviside function greater than epsilon
    if cut_in_out == true : marker = 1 where H>0.5 and 0 where H<0.5
    */

    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
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

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerCrossedElements()
{
    // return a "marker" on the elements traversed by the interface between phio and phi
    // ie : mark as 1 is one of the dof of the element has  (phi * phio < 0)
    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    auto phi = this->phi();
    auto phio = this->phio();

    auto prod = vf::project(this->functionSpace(), elements(mesh),
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
