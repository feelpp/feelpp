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

#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelmodels/levelset/reinitializer.hpp>
#include <feel/feelmodels/levelset/reinitializer_fm.hpp>
#include <feel/feelmodels/levelset/reinitializer_hj.hpp>
#include <feel/feelfilters/straightenmesh.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>

#include <boost/parameter/preprocessor.hpp>

#include <feel/feelmodels/levelset/parameter_map.hpp>

#if defined (MESH_ADAPTATION_LS)
 #include <levelsetmesh/meshadaptation.hpp>
// #warning MESH_ADAPTATION_LS is defined in levelset. Need to be defined identically in the application
#endif


namespace Feel {

    // LevelSet::build parameters
    BOOST_PARAMETER_NAME(space_vectorial)
    BOOST_PARAMETER_NAME(space_markers)
    BOOST_PARAMETER_NAME(reinitializer)
    BOOST_PARAMETER_NAME(projectorL2)
    BOOST_PARAMETER_NAME(projectorL2_vectorial)
    BOOST_PARAMETER_NAME(smoother)
    BOOST_PARAMETER_NAME(smoother_vectorial)

namespace FeelModels {

// time discretization of the advection equation
enum LevelSetTimeDiscretization {BDF2, /*CN,*/ EU, CN_CONSERVATIVE};

/* Levelset reinitialization strategy
 * FM -> Fast-Marching
 * HJ -> Hamilton-Jacobi
 */
enum class LevelSetReinitMethod {FM, HJ};

enum class LevelSetMeasuresExported
{
    Volume, Perimeter
};
enum class LevelSetFieldsExported
{
    GradPhi, ModGradPhi
};

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
    // Projectors
    typedef Projector<space_levelset_type, space_levelset_type> projector_levelset_type;
    typedef boost::shared_ptr<projector_levelset_type> projector_levelset_ptrtype;
    
    typedef Projector<space_levelset_vectorial_type, space_levelset_vectorial_type> projector_levelset_vectorial_type;
    typedef boost::shared_ptr<projector_levelset_vectorial_type> projector_levelset_vectorial_ptrtype;

    //--------------------------------------------------------------------//
    // Reinitialization
    typedef Reinitializer<space_levelset_type> reinitializer_type;
    typedef boost::shared_ptr<reinitializer_type> reinitializer_ptrtype;
    typedef ReinitializerFM<space_levelset_type> reinitializerFM_type;
    typedef boost::shared_ptr<reinitializerFM_type> reinitializerFM_ptrtype;
    typedef ReinitializerHJ<space_levelset_type> reinitializerHJ_type;
    typedef boost::shared_ptr<reinitializerHJ_type> reinitializerHJ_ptrtype;

    enum strategy_before_FM_type {NONE=0, ILP=1, HJ_EQ=2, IL_HJ_EQ=3};

    //--------------------------------------------------------------------//
    // Initial value
    enum class ShapeType {
        SPHERE, ELLIPSE
    };
    static std::map<std::string, ShapeType> ShapeTypeMap;

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //--------------------------------------------------------------------//
    // ModGradPhi advection
    typedef Advection<ConvexType, Lagrange<Order, Scalar>, PeriodicityType> modgradphi_advection_type;
    typedef boost::shared_ptr<modgradphi_advection_type> modgradphi_advection_ptrtype;

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

    BOOST_PARAMETER_MEMBER_FUNCTION( 
            (void), build, tag,
            ( required
              ( space, * )
            )
            ( optional
              ( space_vectorial, (space_levelset_vectorial_ptrtype), space_levelset_vectorial_type::New(_mesh=space->mesh(), _worldscomm=this->worldsComm()) )
              ( space_markers, (space_markers_ptrtype), space_markers_type::New(_mesh=space->mesh(), _worldscomm=this->worldsComm()) )
              ( reinitializer, *( boost::is_convertible<mpl::_, reinitializer_ptrtype> ), reinitializer_ptrtype() )
              ( projectorL2, (projector_levelset_ptrtype), Feel::projector(space, space, backend(_name=prefixvm(this->prefix(),"projector-l2"))) )
              ( projectorL2_vectorial, (projector_levelset_vectorial_ptrtype), Feel::projector(space_vectorial, space_vectorial, backend(_name=prefixvm(this->prefix(),"projector-l2-vec"))) )
              ( smoother, (projector_levelset_ptrtype), Feel::projector(space, space, backend(_name=prefixvm(this->prefix(),"smoother")), DIFF, space->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=this->prefix())/Order, 30) )
              ( smoother_vectorial, (projector_levelset_vectorial_ptrtype), Feel::projector(space_vectorial, space_vectorial, backend(_name=prefixvm(this->prefix(),"smoother-vec")), DIFF, space->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=this->prefix())/Order, 30) )
            )
            )
    {
        super_type::build( space );
        // createFunctionSpaces
        M_spaceLevelSetVec = space_vectorial;
        M_spaceMarkers = space_markers;
        // createInterfaceQuantities
        this->createInterfaceQuantities();
        // createReinitialization
        if( reinitializer )
        {
            M_reinitializer = reinitializer;
        }
        this->createReinitialization();
        // createOthers
        M_projectorL2 = projectorL2;
        M_projectorL2Vec = projectorL2_vectorial;
        M_smoother = smoother;
        M_smootherVectorial = smoother_vectorial;

    }

    //--------------------------------------------------------------------//
    // Initialization
    void init();
    void initLevelsetValue();
    void initPostProcess();

    virtual void loadParametersFromOptionsVm();
    virtual void loadConfigICFile();
    virtual void loadConfigBCFile();
    virtual void loadConfigPostProcess();

    void createFunctionSpaces();
    void createInterfaceQuantities();
    void createReinitialization();
    void createOthers();

    boost::shared_ptr<std::ostringstream> getInfo() const;

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
    element_levelset_ptrtype & phi() { return this->fieldSolutionPtr(); }
    element_levelset_ptrtype const& phi() const { return this->fieldSolutionPtr(); }
    //element_levelset_ptrtype const& phinl() const { return M_phinl; }
    element_levelset_vectorial_ptrtype const& gradPhi() const;
    element_levelset_ptrtype const& modGradPhi() const;
    element_levelset_ptrtype const& stretch() const;
    element_levelset_ptrtype const& heaviside() const;
    element_levelset_ptrtype const& H() const { return this->heaviside(); }
    element_levelset_ptrtype const& dirac() const;
    element_levelset_ptrtype const& D() const { return this->dirac(); }

    element_levelset_vectorial_ptrtype const& normal() const;
    element_levelset_vectorial_ptrtype const& N() const { return this->normal(); }
    element_levelset_ptrtype const& curvature() const;
    element_levelset_ptrtype const& K() const { return this->curvature(); }

    double thicknessInterface() const { return M_thicknessInterface; }
    void setThicknessInterface( double value ) { M_thicknessInterface = value; }

    int iterSinceReinit() const { return M_iterSinceReinit; }

    projector_levelset_ptrtype const& projectorL2() const { return M_projectorL2; }
    projector_levelset_vectorial_ptrtype const& projectorL2Vectorial() const { return M_projectorL2Vec; }
    projector_levelset_ptrtype const& smoother();
    projector_levelset_vectorial_ptrtype const& smootherVectorial();

    void updateInterfaceQuantities();

    //--------------------------------------------------------------------//
    // Markers
    element_markers_ptrtype const& markerInterface();
    element_markers_ptrtype const& markerDirac();
    element_markers_ptrtype const& markerHeaviside(bool invert = false, bool cut_at_half = false);
    element_markers_ptrtype const& markerCrossedElements();

    //--------------------------------------------------------------------//
    // Utility distances
    element_levelset_ptrtype distToBoundary();
    element_levelset_ptrtype distToMarkedFaces( boost::any const& marker );
    element_levelset_ptrtype distToMarkedFaces( std::initializer_list<boost::any> marker );

    //--------------------------------------------------------------------//
    // Advection
    bool hasSourceTerm() const { return false; }
    void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const {}
    void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;

    template<typename ExprT>
    void advect(vf::Expr<ExprT> const& velocity);
    void solve();

    void updateTimeStep();

    //--------------------------------------------------------------------//
    // Reinitialization
    void reinitialize( bool useSmoothReinit = false );

    void setStrategyBeforeFm( int strat = 1 );
    strategy_before_FM_type strategyBeforeFm() { return M_strategyBeforeFM; }
    void setUseMarkerDiracAsMarkerDoneFM( bool val = true ) { M_useMarkerDiracAsMarkerDoneFM  = val; }

    reinitializer_ptrtype const& reinitializer() const { return M_reinitializer; }
    reinitializerFM_ptrtype const& reinitializerFM( bool buildOnTheFly = true );
    reinitializerHJ_ptrtype const& reinitializerHJ( bool buildOnTheFly = true );

    bool hasReinitialized() const { return M_hasReinitialized; }

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
        auto phi_init = this->functionSpace()->element();
        phi_init.on( 
                _range=elements(this->mesh()),
                _expr=expr
                );
        this->setInitialValue( phi_init, doReinitialize );
    }
    template<typename ExprT>
    void setInitialValue(vf::Expr<ExprT> const& expr)
    {
        this->setInitialValue(expr, M_reinitInitialValue );
    }

    //--------------------------------------------------------------------//
    // Export results
    bool hasPostProcessMeasureExported( LevelSetMeasuresExported const& measure) const;
    bool hasPostProcessFieldExported( LevelSetFieldsExported const& field) const;

    //--------------------------------------------------------------------//
    // Physical quantities
    double volume() const;
    double perimeter() const;

    //--------------------------------------------------------------------//
    // Utility functions
    static reinitializer_ptrtype buildReinitializer( 
            LevelSetReinitMethod method, 
            space_levelset_ptrtype const& space,
            std::string const& prefix = "" );

protected:
    //--------------------------------------------------------------------//
    // Levelset data update functions
    void updateGradPhi();
    void updateModGradPhi();
    void updateDirac();
    void updateHeaviside();

    void updateNormal();
    void updateCurvature();

    void updateMarkerDirac();
    void updateMarkerHeaviside(bool invert, bool cut_at_half);
    void updateMarkerCrossedElements();
    void updateMarkerInterface();

private:
    void initWithMesh(mesh_ptrtype mesh);
    void initFastMarching(mesh_ptrtype const& mesh);

    //--------------------------------------------------------------------//
    void addShape( 
            std::pair<ShapeType, parameter_map> const& shape, 
            element_levelset_type & phi );

    //--------------------------------------------------------------------//
    element_levelset_ptrtype const& phio() const { return this->timeStepBDF()->unknowns()[1]; }

    //--------------------------------------------------------------------//
    // Interface rectangular function
    element_levelset_type interfaceRectangularFunction() const { return this->interfaceRectangularFunction(this->phi()); }
    element_levelset_type interfaceRectangularFunction( element_levelset_ptrtype const& p ) const;
    //--------------------------------------------------------------------//
    // Export
    void exportResultsImpl( double time );
    void exportMeasuresImpl( double time );
    // Save
    void saveCurrent() const;


protected:
    //--------------------------------------------------------------------//
    // Interface quantities update flags
    mutable bool M_doUpdateDirac;
    mutable bool M_doUpdateHeaviside;
    mutable bool M_doUpdateNormal;
    mutable bool M_doUpdateCurvature;
    mutable bool M_doUpdateGradPhi;
    mutable bool M_doUpdateModGradPhi;

    //--------------------------------------------------------------------//
    // Levelset initial value
    map_scalar_field<2> M_icDirichlet;
    std::vector<std::pair<ShapeType, parameter_map>> M_icShapes;

    //--------------------------------------------------------------------//
    // Boundary conditions
    std::list<std::string> M_bcMarkersInflow;

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
    // Projectors
    projector_levelset_ptrtype M_projectorL2;
    projector_levelset_vectorial_ptrtype M_projectorL2Vec;
    projector_levelset_ptrtype M_smoother;
    projector_levelset_vectorial_ptrtype M_smootherVectorial;
    //--------------------------------------------------------------------//
    // Levelset data
    mutable element_levelset_vectorial_ptrtype M_levelsetGradPhi;
    mutable element_levelset_ptrtype M_levelsetModGradPhi;
    mutable element_levelset_ptrtype M_heaviside;
    mutable element_levelset_ptrtype M_dirac;
    //--------------------------------------------------------------------//
    // Normal, curvature
    mutable element_levelset_vectorial_ptrtype M_levelsetNormal;
    mutable element_levelset_ptrtype M_levelsetCurvature;
    bool M_doSmoothCurvature;

    //--------------------------------------------------------------------//
    // Reinitialization
    reinitializer_ptrtype M_reinitializer;
    reinitializerFM_ptrtype M_reinitializerFM;
    reinitializerHJ_ptrtype M_reinitializerHJ;
    bool M_reinitializerIsUpdatedForUse;

    boost::shared_ptr<Projector<space_levelset_type, space_levelset_type>> M_smootherFM;

    bool M_hasReinitialized;
    bool M_hasReinitializedSmooth;
    int M_iterSinceReinit;
    // Vector that stores the iterSinceReinit of each time-step
    std::vector<int> M_vecIterSinceReinit;
    //bool M_useSmoothReinitialization;

    //--------------------------------------------------------------------//
    // Backends
    backend_ptrtype M_backend_smooth;

    //--------------------------------------------------------------------//
    // Advection
    int M_timeOrder;

    //--------------------------------------------------------------------//
    // ModGradPhi advection
    bool M_useGradientAugmented;
    bool M_reinitGradientAugmented;
    bool M_reinitStretchAugmented;
    modgradphi_advection_ptrtype M_modGradPhiAdvection;

    //--------------------------------------------------------------------//
    // Stretch advection
    bool M_useStretchAugmented;
    modgradphi_advection_ptrtype M_stretchAdvection;
    mutable element_levelset_ptrtype M_levelsetStretch;

    //--------------------------------------------------------------------//
    // Export
    std::set<LevelSetMeasuresExported> M_postProcessMeasuresExported;
    std::set<LevelSetFieldsExported> M_postProcessFieldsExported;
    //--------------------------------------------------------------------//
    // Parameters
    double M_thicknessInterface;
    bool M_useAdaptiveThicknessInterface;
    bool M_useRegularPhi;
    bool M_useHeavisideDiracNodalProj;

    double k_correction;
    //--------------------------------------------------------------------//
    // Reinitialization
    LevelSetReinitMethod M_reinitMethod;
    strategy_before_FM_type M_strategyBeforeFM;
    bool M_useMarkerDiracAsMarkerDoneFM;
    //int M_hjMaxIter;
    //double M_hjDtau;
    //double M_hjTol;
    bool M_reinitInitialValue;

    //LevelSetTimeDiscretization M_discrMethod;

}; //class LevelSet

template<typename ConvexType, int Order, typename PeriodicityType>
template<typename ExprT>
void 
LevelSet<ConvexType, Order, PeriodicityType>::advect(vf::Expr<ExprT> const& velocity)
{
    this->updateAdvectionVelocity(velocity);
    this->solve();
}

} // namespace FeelModels
} // namespace Feel

#endif
