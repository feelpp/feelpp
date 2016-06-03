#ifndef _LEVELSET_HPP
#define _LEVELSET_HPP 1

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/projector.hpp>

#include <feel/feelmodels/advection/advection.hpp>

#include <feel/feells/reinitializer.hpp>
#include <feel/feells/reinit_fms.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>

#if defined (MESH_ADAPTATION_LS)
 #include <levelsetmesh/meshadaptation.hpp>
// #warning MESH_ADAPTATION_LS is defined in levelset. Need to be defined identically in the application
#endif


/*
about time discretization and saving :
- BDF2 is implemented but the save() option saves only one iteration back in time which is incompatible with BDF2. If needed, do not implement the saving of earlier levelset but use the BDF2() class of Feel++ which does it already (and even at higher order).
- In practice, higher order time discretization is not compatible with reinitialization, CN sheme should be prefered. (and this explains why I didn't add Feel::BDF2 in levelset)
*/

/*
  about reinitialization :
Fast Marching Method is faster and more precise than solving Hamilton Jacobi equation :

- for periodic boundary conditions, FMM takes care of periodic boundary condition but needs a non periodic space. Thus, one has to give to the reinitialization a non periodic mesh and separately give the appropriate periodicity to apply on it.
(in reinitialization, the periodicity is given has a template parameter and has to be instantiated properly)


- at order > 1, since reinitFM accepts only P1 spaces, one pass by OperatorLagrange to reinit a P1 levelset and reset it to a Pn space. So this method is exact at quadrature points only.
-------> the reinitialization space is always P1 NON PERIODIC !

*/

namespace Feel {
namespace FeelModels {

// time discretization of the advection equation
enum LevelSetTimeDiscretization {BDF2, /*CN,*/ EU, CN_CONSERVATIVE};

/* Levelset reinitialization strategy
 * FM -> Fast-Marching
 * HJ -> Hamilton-Jacobi
 */
enum class LevelSetReinitMethod {FM, HJ};

template<typename ConvexType, int Order=1, typename PeriodicityType = NoPeriodicity>
class LevelSet : 
    public ModelBase
{
public:
    typedef ModelBase super_type;
    typedef LevelSet<ConvexType, Order, PeriodicityType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef double value_type;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Periodicity
    typedef PeriodicityType periodicity_type;

    //--------------------------------------------------------------------//
    // Advection
    typedef Lagrange<Order, Scalar> basis_levelset_type;
    //typedef FunctionSpace<mesh_type, bases<basis_levelset_type>, value_type, Periodicity<periodicity_type>> space_levelset_type;

    typedef Advection<convex_type, basis_levelset_type, periodicity_type> advection_type;

    //--------------------------------------------------------------------//
    // Space levelset
    typedef typename advection_type::space_advection_type space_levelset_type;
    typedef boost::shared_ptr<space_levelset_type> space_levelset_ptrtype;
    typedef typename space_levelset_type::element_type element_levelset_type;
    typedef boost::shared_ptr<element_levelset_type> element_levelset_ptrtype;

    //--------------------------------------------------------------------//
    // Space vectorial levelset
    typedef Lagrange<Order, Vectorial> basis_levelset_vectorial_type;
    typedef FunctionSpace<mesh_type, bases<basis_levelset_vectorial_type>, value_type, Periodicity<periodicity_type> > space_levelset_vectorial_type;
    typedef boost::shared_ptr<space_levelset_vectorial_type> space_levelset_vectorial_ptrtype;
    typedef typename space_levelset_vectorial_type::element_type element_levelset_vectorial_type;
    typedef boost::shared_ptr< element_levelset_vectorial_type > element_levelset_vectorial_ptrtype;

    // ---------------- Correction levelset space -------
#if (LEVELSET_CONSERVATIVE_ADVECTION == 1)
    // correction only one ddl
    typedef bases<Lagrange<0, Scalar, Continuous> > basisLSCorr_type;
    typedef FunctionSpace<mesh_type, basisLSCorr_type > spaceLSCorr_type;
    typedef boost::shared_ptr<spaceLSCorr_type> spaceLSCorr_ptrtype;
    typedef typename spaceLSCorr_type::element_type elementLSCorr_type;
    typedef boost::shared_ptr< elementLSCorr_type > elementLSCorr_ptrtype;
#elif (LEVELSET_CONSERVATIVE_ADVECTION == 2)
    // correction in the same space than phi
    typedef space_levelset_type spaceLSCorr_type;
    typedef space_levelset_ptrtype spaceLSCorr_ptrtype;
    typedef element_levelset_type elementLSCorr_type;
    typedef element_levelset_ptrtype elementLSCorr_ptrtype;
#endif

    //--------------------------------------------------------------------//
    // Space levelset reinitP1
    // used for FM reinitialization, always P1 non-periodic !
    typedef Lagrange<1, Scalar> basis_levelset_reinitP1_type;
    typedef FunctionSpace<mesh_type, bases<basis_levelset_reinitP1_type>, value_type, Periodicity<NoPeriodicity>, mortars<NoMortar> > space_levelset_reinitP1_type;
    typedef boost::shared_ptr<space_levelset_reinitP1_type> space_levelset_reinitP1_ptrtype;
    typedef typename space_levelset_reinitP1_type::element_type element_levelset_reinitP1_type;
    typedef boost::shared_ptr< element_levelset_reinitP1_type > element_levelset_reinitP1_ptrtype;

    //--------------------------------------------------------------------//
    // Space markers P0
    typedef Lagrange<0, Scalar, Discontinuous> basis_markers_type;
    typedef FunctionSpace<mesh_type, bases<basis_markers_type>, value_type, Periodicity<NoPeriodicity> > space_markers_type;
    typedef boost::shared_ptr<space_markers_type> space_markers_ptrtype;
    typedef typename space_markers_type::element_type element_markers_type;
    typedef boost::shared_ptr< element_markers_type > element_markers_ptrtype;

#if defined (MESH_ADAPTATION_LS)
    typedef MeshAdaptation<Dim, Order, 1, periodicity_type > mesh_adaptation_type;
    typedef boost::shared_ptr< mesh_adaptation_type > mesh_adaptation_ptrtype;
#endif

    //--------------------------------------------------------------------//
    // Interpolation operators
    typedef boost::tuple<
        boost::mpl::size_t<MESH_ELEMENTS>,
        typename MeshTraits<mesh_type>::element_const_iterator,
        typename MeshTraits<mesh_type>::element_const_iterator> range_visu_ho_type;

    typedef OperatorInterpolation<
        space_levelset_type, // from space
        space_levelset_reinitP1_type, // to space
        range_visu_ho_type> op_interpolation_LS_to_P1_type;

    typedef OperatorInterpolation<
        space_levelset_reinitP1_type, // from space
        space_levelset_type, // to space
        range_visu_ho_type> op_interpolation_P1_to_LS_type;

    typedef boost::shared_ptr<op_interpolation_LS_to_P1_type> op_interpolation_LS_to_P1_ptrtype;
    typedef boost::shared_ptr<op_interpolation_P1_to_LS_type> op_interpolation_P1_to_LS_ptrtype;

    typedef OperatorLagrangeP1<space_levelset_type> op_lagrangeP1_type;
    typedef boost::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    //--------------------------------------------------------------------//
    // Reinitialization
    typedef Reinitializer<space_levelset_type> reinitializer_type;
    typedef boost::shared_ptr<reinitializer_type> reinitializer_ptrtype;

    enum strategy_before_FM_type {NONE=0, ILP=1, HJ_EQ=2, IL_HJ_EQ=3};

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    LevelSet(
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );


    LevelSet( self_type const& L ) = default;

    static self_ptrtype New( 
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    //--------------------------------------------------------------------//
    // Initialization
    void init();
    void initFromMesh( mesh_ptrtype const& mesh );

    virtual void loadParametersFromOptionsVm();

    void createFunctionSpaces();
    void createAdvection();
    void createReinitialization();
    void createOthers();

    //--------------------------------------------------------------------//
    space_levelset_ptrtype const& functionSpace() const { return M_advection->functionSpace(); }
    space_markers_ptrtype const& functionSpaceMarkers() const { return M_spaceMarkers; }
    space_levelset_vectorial_ptrtype const& functionsSpaceVectorial() const { return M_spaceLevelSetVec; }
    space_levelset_reinitP1_ptrtype const& functionSpaceReinitP1() const { return M_spaceReinitP1; }

    space_levelset_ptrtype const& functionSubspace() const { return M_subspaceLevelSet; }
    space_levelset_vectorial_ptrtype const& functionSubspaceVectorial() const { return M_subspaceLevelSetVec; }

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"LevelsetMesh.path"); }

    mesh_ptrtype const& mesh() const { return M_mesh; }
    mesh_ptrtype const& submesh() const { return M_submesh; }

    //--------------------------------------------------------------------//
    // Mesh adaptation
#if defined (MESH_ADAPTATION_LS)
    mesh_ptrtype adaptMesh(double time, element_levelset_ptrtype elt);
    mesh_adaptation_ptrtype mesh_adapt;
    template<typename ExprT> mesh_ptrtype adaptedMeshFromExp(ExprT expr);
#endif

    //--------------------------------------------------------------------//
    // Levelset
    element_levelset_ptrtype const& phi() const { return M_advection->fieldSolutionPtr(); }
    //element_levelset_ptrtype const& phinl() const { return M_phinl; }
    element_levelset_ptrtype const& heaviside() const { return M_heaviside; }
    element_levelset_ptrtype const& H() const { return this->heaviside(); }
    element_levelset_ptrtype const& dirac() const { return M_dirac; }
    element_levelset_ptrtype const& D() const { return this->dirac(); }

    double mass() const { return M_mass; }

    double thicknessInterface() const { return M_thicknessInterface; }

    int iterSinceReinit() const { return M_iterSinceReinit; }

    //--------------------------------------------------------------------//
    // Markers
    element_markers_ptrtype const& markerInterface();
    element_markers_ptrtype const& markerDirac();
    element_markers_ptrtype const& markerHeaviside(bool invert = false, bool cut_at_half = false);
    element_markers_ptrtype const& markerCrossedElements();

    //--------------------------------------------------------------------//
    // Advection
    template<typename ExprT>
    void advect(vf::Expr<ExprT> const& velocity);

    /* returns phi after advection by Velocity
      do not change the value of this->M_phi or any other variable (delta heavyside ...) */
    template < typename TVeloc >
    element_levelset_ptrtype phiAdvected(TVeloc const& Velocity, bool stabilization = true)
    { return this->advReactUpdate(Velocity, stabilization); }
    
    //--------------------------------------------------------------------//
    // Reinitialization
    void reinitialize();

    // template < typename TVeloc >
    // void updateE(TVeloc& Velocity);

    //void initialize(element_levelset_ptrtype, bool doFirstReinit = false, LevelSetReinitMethod method=FM, int max_iter=-1, double dtau=0.01, double tol=0.1);
    //element_levelset_ptrtype circleShape(double r0, double x0, double y0, double z0=0);
    //element_levelset_ptrtype ellipseShape(double a_ell, double b_ell, double x0, double y0, double z0=0);
    //void imposePhi( element_levelset_ptrtype );

    //element_levelset_ptrtype makeDistFieldFromParametrizedCurve(std::function<double(double)> xexpr, std::function<double(double)> yexpr,
                                                         //double tStart, double tEnd, double dt,
                                                         //bool useRandomPt=true, double randomness=doption("gmsh.hsize") / 5.,
                                                         //bool export_points=false, std::string export_name="");


    //element_levelset_ptrtype makeDistFieldFromParametrizedCurve(std::tuple< std::function<double(double)>, std::function<double(double)>, double, double >  paramCurve,
                                                         //double dt, bool useRandomPt=true, double randomness=doption("gmsh.hsize") / 5.,
                                                         //bool export_points=false, std::string export_name="");

    //template< typename Elt1, typename Elt2 >
    //std::vector<double> getStatReinit(Elt1 __phio, Elt2 __phi);

    /* update the submesh and subspaces*/
    //void updateSubMeshSubSpace(element_markers_ptrtype marker);
    //void updateSubMeshSubSpace();

    std::string levelsetInfos( bool show = false );


    /*// ----------- serialization, save, restart
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {ar & *itersincereinit & *M_phi & *M_phio;}

    void restart(double restart_time);
    void save(double time);
    static double getRestartTimeFromExporter( std::string name );*/

    // ------------ setters -------------
    void setStrategyBeforeFm( int strat = 1 );
    strategy_before_FM_type strategyBeforeFm() { return M_strategyBeforeFM; }
    void setUseMarker2AsMarkerDoneFmm( bool val=true ) { M_useMarker2AsMarkerDoneFmm = val; }
    void setThicknessInterface( double value ) { M_thicknessInterface = value; }

    // ------------------- projectors ----------------
    //M_l2p : for projection (update H, D ...), M_smooth : only for reinit !
    boost::shared_ptr< Projector<space_levelset_type, space_levelset_type> >  M_l2p;
    boost::shared_ptr< Projector<space_levelset_vectorial_type, space_levelset_vectorial_type> >  M_l2pVec;

protected:
    //--------------------------------------------------------------------//
    // Levelset data update functions
    void updateDirac();
    void updateHeaviside();
    void updateMass();

    void updateMarkerDirac();
    void updateMarkerHeaviside(bool invert, bool cut_at_half);
    void updateMarkerCrossedElements();
    void updateMarkerInterface();

private:
    void initWithMesh(mesh_ptrtype mesh);
    void initFastMarching(mesh_ptrtype const& mesh);

    //--------------------------------------------------------------------//
    // Levelset previous step phi
    element_levelset_ptrtype const& phio() const { return M_advection->timeStepBDF()->unknowns()[1]; }

    template < typename TVeloc >
    element_levelset_ptrtype advReactUpdate(TVeloc& Velocity, bool updateStabilization=true, bool updateBilinearForm=true);

#if defined(LEVELSET_CONSERVATIVE_ADVECTION)
    template < typename TVeloc >
    void conservativeH(TVeloc& Velocity, element_levelset_ptrtype Hc);

    void updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    void updateResidual(const vector_ptrtype& X, vector_ptrtype& R, element_levelset_ptrtype Hc);
    void computeCorrection(element_levelset_ptrtype Hc);
#endif

    element_levelset_ptrtype explicitHJ(int, double, double);
    element_levelset_ptrtype explicitILHJ(int, double, double);
    double distToDist();
    //        element_levelset_reinitP1_ptrtype makeShape(std::function<double(double)> xexpr, std::function<double(double)> yexpr, double tStart, double tEnd, double dt, space_markers_ptrtype mspaceP0, mesh_ptrtype msmesh, bool useRandomPt=false);



protected:
    //--------------------------------------------------------------------//
    // Levelset data
    //element_levelset_ptrtype M_phi;
    //element_levelset_ptrtype M_phio;
    //element_levelset_ptrtype M_phinl;

    element_levelset_ptrtype M_heaviside;
    element_levelset_ptrtype M_dirac;
    double M_mass;

#if defined(LEVELSET_CONSERVATIVE_ADVECTION)
    elementLSCorr_ptrtype phic;
#endif


private:
    //--------------------------------------------------------------------//
    // Mesh 
    mesh_ptrtype M_mesh;
    mesh_ptrtype M_submesh;

    //--------------------------------------------------------------------//
    periodicity_type M_periodicity;

    //--------------------------------------------------------------------//
    // Spaces
    //space_levelset_ptrtype M_spaceLevelSet;
    space_levelset_vectorial_ptrtype M_spaceLevelSetVec;
    space_markers_ptrtype M_spaceMarkers;

    space_levelset_ptrtype M_subspaceLevelSet;
    space_levelset_vectorial_ptrtype M_subspaceLevelSetVec;

    space_levelset_reinitP1_ptrtype M_spaceReinitP1;
    
    //--------------------------------------------------------------------//
    // Markers
    element_markers_ptrtype M_markerDirac;
    element_markers_ptrtype M_markerHeaviside;
    element_markers_ptrtype M_markerCrossedElements;
    element_markers_ptrtype M_markerInterface;
    bool M_doUpdateMarkers;

    //--------------------------------------------------------------------//
    // Reinitialization
    op_interpolation_LS_to_P1_ptrtype M_opInterpolationLStoP1;
    op_interpolation_P1_to_LS_ptrtype M_opInterpolationP1toLS;
    op_lagrangeP1_ptrtype M_opLagrangeP1;

    //reinitializer_ptrtype M_reinitializer;
    boost::shared_ptr< ReinitializerFMS<space_levelset_reinitP1_type, PeriodicityType> > M_reinitializerFMS;
    bool M_reinitializerIsUpdatedForUse;

    boost::shared_ptr<Projector<space_levelset_type, space_levelset_type>> M_smooth;

    int M_iterSinceReinit;

    //--------------------------------------------------------------------//
    // Backends
    backend_ptrtype M_backend_smooth;

    //--------------------------------------------------------------------//
    // Advection
    boost::shared_ptr<advection_type> M_advection;
    //boost::shared_ptr<advection_type> M_advection_hj;

#if defined(LEVELSET_CONSERVATIVE_ADVECTION)
    spaceLSCorr_ptrtype M_spaceLSCorr;

    backend_ptrtype backend_nl;
    backend_ptrtype backend_h;

    sparse_matrix_ptrtype J;
    vector_ptrtype R;

    sparse_matrix_ptrtype D_h;
    vector_ptrtype F_h;
#endif

    //--------------------------------------------------------------------//
    // Parameters
    double M_thicknessInterface;
    bool M_useRegularPhi;

    int impose_inflow;
    bool hdNodalProj;
    double k_correction;
    //--------------------------------------------------------------------//
    // Reinitialization
    //bool M_enableReinit;
    //int M_reinitEvery;
    LevelSetReinitMethod M_reinitMethod;
    strategy_before_FM_type M_strategyBeforeFM;
    bool M_useMarker2AsMarkerDoneFmm;
    //int M_hjMaxIter;
    //double M_hjDtau;
    //double M_hjTol;

    //LevelSetTimeDiscretization M_discrMethod;


    // -------------- variables -----------
    boost::timer ch;
    std::ofstream statReinitFile;
    int __iter;

}; //class LevelSet



/*template<int Order, int Dim, typename PeriodicityType>
template< typename Elt1, typename Elt2 >
std::vector<double>
Feel::levelset::LevelSet<Order, Dim, PeriodicityType>::getStatReinit(Elt1 __phio,  Elt2 __phi)
{

    using namespace Feel;
    using namespace Feel::vf;

    // double& MassError, double& SignChangeError, double& distDist;
    std::vector< double > stat(3);

    //Mass_Error=mass(phi)/mass(phio)-1;

    double mass_phi = integrate(elements(M_mesh),
                                idv(__phi) < 0 ).evaluate()(0,0);
    double mass_phio= integrate(elements(M_mesh),
                                idv(__phio) < 0 ).evaluate()(0,0);

    stat[0] = mass_phi / mass_phio - 1.0 ; // MassError

    stat[1] = integrate(elements(M_mesh),
                        idv(__phi)*idv(__phio) < 0. ).evaluate()(0,0); // SignChangeError

    stat[2] = integrate(elements(M_mesh),
                vf::abs(sqrt(gradv(__phi)*trans(gradv(__phi)))-1.0)).evaluate()(0,0); // distDist

    stat[2] /= integrate(elements(M_mesh), vf::cst(1.)).evaluate()(0,0);

    return stat;
}//GetStatReinit*/


} //namespace FeelModels
} //namespace Feel

// #if CONSERVATIVE_ADVECTION
// #include "nonlinearcorrection.cpp"
// #endif


#endif
