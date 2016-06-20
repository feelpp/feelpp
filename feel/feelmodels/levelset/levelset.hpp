/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2016-05-20

 Copyright (C) 2016 Université Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 \file levelset.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2016-05-20
 */
#ifndef _LEVELSET_HPP
#define _LEVELSET_HPP 1

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/projector.hpp>

#include <feel/feelmodels/advection/advection.hpp>

#include <feel/feelmodels/levelset/reinitializer.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>

#if defined (MESH_ADAPTATION_LS)
 #include <levelsetmesh/meshadaptation.hpp>
// #warning MESH_ADAPTATION_LS is defined in levelset. Need to be defined identically in the application
#endif


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
class LevelSet 
    : public AdvectionBase< ConvexType, Lagrange<Order, Scalar>, PeriodicityType >
    , public boost::enable_shared_from_this< LevelSet<ConvexType, Order, PeriodicityType> >
{
public:
    typedef AdvectionBase< ConvexType, Lagrange<Order, Scalar>, PeriodicityType > super_type;
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
    // Space levelset
    typedef typename super_type::space_advection_type space_levelset_type;
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
    // Exporter
    typedef Exporter<mesh_type, nOrderGeo> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;

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

    void build();
    void build( mesh_ptrtype const& mesh );

    //--------------------------------------------------------------------//
    // Initialization
    void init();

    virtual void loadParametersFromOptionsVm();
    virtual void loadConfigICFile();

    void createFunctionSpaces();
    void createReinitialization();
    void createOthers();

    //--------------------------------------------------------------------//
    space_markers_ptrtype const& functionSpaceMarkers() const { return M_spaceMarkers; }
    space_levelset_vectorial_ptrtype const& functionSpaceVectorial() const { return M_spaceLevelSetVec; }

    space_levelset_ptrtype const& functionSubspace() const { return M_subspaceLevelSet; }
    space_levelset_vectorial_ptrtype const& functionSubspaceVectorial() const { return M_subspaceLevelSetVec; }

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"LevelsetMesh.path"); }

    //mesh_ptrtype const& mesh() const { return M_advection->mesh(); }
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
    element_levelset_ptrtype const& phi() const { return this->fieldSolutionPtr(); }
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
    bool hasSourceTerm() const { return false; }
    void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const {}
    void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const {}

    template<typename ExprT>
    void advect(vf::Expr<ExprT> const& velocity);
    void solve();

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

    //--------------------------------------------------------------------//
    // Initial value
    void setInitialValue(element_levelset_type const& phiv, bool doReinitialize);
    void setInitialValue(element_levelset_type const& phiv)
    {
        this->setInitialValue(phiv, M_reinitInitialValue);
    }
    void setInitialValue(element_levelset_ptrtype const& phiv, bool doReinitialize) 
    { 
        this->setInitialValue(*phiv, doReinitialize);
    }
    void setInitialValue(element_levelset_ptrtype const& phiv)
    {
        this->setInitialValue(*phiv, M_reinitInitialValue);
    }
    template<typename ExprT>
    void setInitialValue(vf::Expr<ExprT> const& expr, bool doReinitialize)
    {
        this->phi()->on( 
                _range=elements(this->mesh()),
                _expr=expr
                );
        if(doReinitialize)
            this->reinitialize();

        updateHeaviside();
        updateDirac();
        updateMass();
    }
    template<typename ExprT>
    void setInitialValue(vf::Expr<ExprT> const& expr)
    {
        this->setInitialValue(expr, M_reinitInitialValue );
    }
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

    //--------------------------------------------------------------------//
    // Export results
    using super_type::exportResults;
    void exportResults( double time );

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
    void setUseMarkerDiracAsMarkerDoneFM( bool val = true ) { M_useMarkerDiracAsMarkerDoneFM  = val; }
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
    // Levelset
    element_levelset_ptrtype & phi() { return this->fieldSolutionPtr(); }
    element_levelset_ptrtype const& phio() const { return this->timeStepBDF()->unknowns()[1]; }

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

    //--------------------------------------------------------------------//
    // Levelset initial value
    map_scalar_field<2> M_icDirichlet;

#if defined(LEVELSET_CONSERVATIVE_ADVECTION)
    elementLSCorr_ptrtype phic;
#endif


private:
    //--------------------------------------------------------------------//
    // Mesh 
    //mesh_ptrtype M_mesh;
    mesh_ptrtype M_submesh;

    //--------------------------------------------------------------------//
    // Periodicity
    periodicity_type M_periodicity;

    //--------------------------------------------------------------------//
    // Spaces
    //space_levelset_ptrtype M_spaceLevelSet;
    space_levelset_vectorial_ptrtype M_spaceLevelSetVec;
    space_markers_ptrtype M_spaceMarkers;

    space_levelset_ptrtype M_subspaceLevelSet;
    space_levelset_vectorial_ptrtype M_subspaceLevelSetVec;

    //--------------------------------------------------------------------//
    // Markers
    element_markers_ptrtype M_markerDirac;
    element_markers_ptrtype M_markerHeaviside;
    element_markers_ptrtype M_markerCrossedElements;
    element_markers_ptrtype M_markerInterface;
    bool M_doUpdateMarkers;

    //--------------------------------------------------------------------//
    // Reinitialization
    reinitializer_ptrtype M_reinitializer;
    bool M_reinitializerIsUpdatedForUse;

    boost::shared_ptr<Projector<space_levelset_type, space_levelset_type>> M_smooth;

    int M_iterSinceReinit;

    //--------------------------------------------------------------------//
    // Backends
    backend_ptrtype M_backend_smooth;

    //--------------------------------------------------------------------//
    // Advection
    //advection_ptrtype M_advection;
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
    // Export
    //exporter_ptrtype M_exporter;
    //bool M_doExportAdvection;
    //--------------------------------------------------------------------//
    // Parameters
    double M_thicknessInterface;
    bool M_useRegularPhi;
    bool M_useHeavisideDiracNodalProj;

    //int impose_inflow;
    double k_correction;
    //--------------------------------------------------------------------//
    // Reinitialization
    //bool M_enableReinit;
    //int M_reinitEvery;
    LevelSetReinitMethod M_reinitMethod;
    strategy_before_FM_type M_strategyBeforeFM;
    bool M_useMarkerDiracAsMarkerDoneFM;
    //int M_hjMaxIter;
    //double M_hjDtau;
    //double M_hjTol;
    bool M_reinitInitialValue;

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
