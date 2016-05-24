#include "levelset.hpp"

template<int Order, int Dim, typename PeriodicityType>
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::LevelSet(mesh_ptrtype mesh, std::string const& prefix, double TimeStep, PeriodicityType periodicityLS)
   :
    M_prefix(prefix),
    M_periodicity(periodicityLS),
    pi(M_PI),
    Dt(TimeStep)
{
    itersincereinit.reset( new int(0) );
    M_mass.reset( new double(0) );

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
    M_epsilon=doption(prefixvm(prefix,"thickness-interface"));

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

    M_phi = M_spaceLS->elementPtr();
    M_phio = M_spaceLS->elementPtr();
    M_phinl=M_spaceLS->elementPtr();
    M_H = M_spaceLS->elementPtr();
    M_delta = M_spaceLS->elementPtr();

#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
    if (M_discrMethod==CN_CONSERVATIVE)
        {
            phic = M_spaceLSCorr->elementPtr();
        }
#endif

    reinitializerUpdated = false;

    this->levelsetInfos(true);

} //constructor

template<int Order, int Dim, typename PeriodicityType>
void
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::initWithMesh(mesh_ptrtype mesh)
{
    // +++++++++ initialize every quantities which need the mesh +++++++++++
    M_mesh = mesh;

    //  -------  spaces -------
    M_spaceP0 = spaceP0_type::New(_mesh=M_mesh, _periodicity=periodicity(NoPeriodicity()) );
    M_spaceLSVec = spaceLSVec_type::New(_mesh=M_mesh, _periodicity= periodicity(M_periodicity) );
    if (M_advecstabmethod == CIP)
        M_spaceLS = spaceLS_type::New(_mesh=M_mesh, _periodicity= periodicity(M_periodicity), _extended_doftable=std::vector<bool>(1,true) );
    else
        M_spaceLS = spaceLS_type::New(_mesh=M_mesh, _periodicity= periodicity(M_periodicity) );

    //  -------  backends ------
    M_backend_advrea = backend(_name="ls-advec");
    M_backend_smooth = backend(_name="ls-smooth");

    // ---------- advection -----------
    LOG(INFO)<<"start creating advection for level set"<<std::endl;
    M_advection = advection_type::New( M_spaceLS, M_backend_advrea, M_advecstabmethod, option(prefixvm(M_prefix,"coeff-stab")).template as<double>() );
    LOG(INFO)<<"end creating advection for level set"<<std::endl;

    if ( (M_reinitmethod == HJ) || (strategyBeforeFm==HJ_EQ) )
    {
        LOG(INFO)<<"start creating advection for hamilton jacobi"<<std::endl;
        auto backend_reinit_hj = backend(_name="ls-hj-reinit");
        M_advection_hj = advection_type::New( M_spaceLS, backend_reinit_hj, SUPG, option( prefixvm(M_prefix,"hj-coeff-stab")).template as<double>() );
        LOG(INFO)<<"end creating advection for hamilton jacobi"<<std::endl;
    }

    if (M_epsilon > 1.)
        std::cout<<"\n\n------------ WARNING : thickness_interface > 1 ----------\n\n";

    //------------- initialize matrices for conservative advection-------------
#if defined (LEVELSET_CONSERVATIVE_ADVECTION)
    if (M_discrMethod==CN_CONSERVATIVE)
        {
#if (LEVELSET_CONSERVATIVE_ADVECTION == 1)
            M_spaceLSCorr = spaceLSCorr_type::New(M_mesh);
#else
            M_spaceLSCorr = M_spaceLS;
#endif

            // to solve non linear problem for correction
            backend_nl = backend(_name="ls-nl-advec");

            J = backend_nl->newMatrix(M_spaceLSCorr, M_spaceLSCorr);
            form2(M_spaceLSCorr, M_spaceLSCorr, J);

            R = backend_nl->newVector(M_spaceLSCorr);
            form1(M_spaceLSCorr, R);

            // to compute conservative heavyside function
            backend_h = backend(_name="ls-nl-advec-h");
            D_h = backend_h->newMatrix(M_spaceLS, M_spaceLS );
            form2( M_spaceLS, M_spaceLS, D_h);

            F_h = backend_h->newVector(M_spaceLS);
            form1(M_spaceLS, F_h);

            // might be an option (one more !)
            k_correction = (LEVELSET_CONSERVATIVE_ADVECTION==1) ? 0. : M_epsilon;
        }
#endif

    // --------------- projectors ---------
    LOG(INFO)<<"start creating projectors for level set"<<std::endl;
    M_l2p = projector(M_spaceLS, M_spaceLS, backend(_name="ls-l2p") );
    M_l2pVec = projector(M_spaceLSVec, M_spaceLSVec, backend(_name="ls-l2pVec") );
    LOG(INFO)<<"end creating projectors for level set"<<std::endl;

} // init


template<int Order, int Dim, typename PeriodicityType>
void Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::imposePhi( elementLS_ptrtype phi )
{
    using namespace Feel::vf;

    // impose phi and save the old one in phio
    M_phio= M_phi;
    M_phi = phi;
    updateHeavyside();
    updateDirac();
    updateMass();
}//imposePhi


template<int Order, int Dim, typename PeriodicityType>
void Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::initialize(elementLS_ptrtype phi, bool doFirstReinit, ReinitMethod method, int max_iter, double dtau, double tol)
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


template<int Order, int Dim, typename PeriodicityType>
typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::elementLS_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::circleShape(double r0, double x0, double y0, double z0)
{
    auto shape = M_spaceLS->elementPtr();

    *shape = vf::project(M_spaceLS, elements(M_mesh),
                         sqrt((Px()-x0)*(Px()-x0)
                              + (Py()-y0)*(Py()-y0)
                              + (Pz()-z0)*(Pz()-z0)) - r0 );
    return shape;
} //circleShape

template<int Order, int Dim, typename PeriodicityType>
typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::elementLS_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::ellipseShape(double a_ell, double b_ell, double x0, double y0, double z0)
{
    auto shape = M_spaceLS->elementPtr();
    *shape = vf::project(M_spaceLS, elements(M_mesh),
                         sqrt((Px()-x0)*(Px()-x0) / (a_ell * a_ell)
                              + (Py()-y0)*(Py()-y0) / (b_ell * b_ell)
                              + (Pz()-z0)*(Pz()-z0) / (b_ell * b_ell) ) - 1 );
    return shape;
}//ellipseShape


/* accessors defined here. Mac doesn't like it to be in the .hpp */

template<int Order, int Dim, typename PeriodicityType>
inline const typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::elementLS_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::phi()
{
    return M_phi;
}

template<int Order, int Dim, typename PeriodicityType>
inline const typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::elementLS_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::phinl()
{
    return M_phinl;
}

template<int Order, int Dim, typename PeriodicityType>
inline const typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::elementLS_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::H()
{
    return M_H;
}

template<int Order, int Dim, typename PeriodicityType>
inline const typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::elementLS_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::D()
{
    return M_delta;
}

template<int Order, int Dim, typename PeriodicityType>
inline double
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::mass()
{
    return *M_mass;
}

template<int Order, int Dim, typename PeriodicityType>
inline int
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::iterSinceReinit()
{
    return *itersincereinit;
}

/*delta and heavyside thickness*/
template<int Order, int Dim, typename PeriodicityType>
inline double
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::thicknessInterface()
{
    return M_epsilon;
}




template<int Order, int Dim, typename PeriodicityType>
typename Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::self_ptrtype
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::New( mesh_ptrtype mesh, std::string const& prefix, double TimeStep, PeriodicityType periodocity )
{
    self_ptrtype new_ls( new self_type( mesh, prefix, TimeStep, periodocity) );
    return new_ls;
}

template<int Order, int Dim, typename PeriodicityType>
void
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::setDt( double time )
{ Dt = time; }


template<int Order, int Dim, typename PeriodicityType>
void
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::setStrategyBeforeFm( int strat )
{
    if (reinitializerUpdated)
        LOG(INFO)<<" !!!  WARNING !!! : setStrategyBeforeFm set after the fast marching has been actually initialized ! \n";
    strategyBeforeFm = (strategyBeforeFm_type) strat;
}

template<int Order, int Dim, typename PeriodicityType>
Feel::levelset::strategyBeforeFm_type
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::getStrategyBeforeFm()
{
    return strategyBeforeFm;
}



template<int Order, int Dim, typename PeriodicityType>
void
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::setUseMarker2AsMarkerDoneFmm(bool value)
{
    M_useMarker2AsMarkerDoneFmm = value;
}

template<int Order, int Dim, typename PeriodicityType>
void
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::setThicknessInterface( double value )
{
    M_epsilon = value;
}


//instantiation
#define FEELPP_INSTANTIATE_LEVELSET 1
#include "levelset_instance.hpp"
