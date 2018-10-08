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

#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>
#include <feel/feelmodels/levelset/levelsettoolmanager.hpp>
#include <feel/feelmodels/levelset/reinitializer.hpp>
#include <feel/feelmodels/levelset/reinitializer_fm.hpp>
#include <feel/feelmodels/levelset/reinitializer_hj.hpp>

#include <feel/feelfilters/straightenmesh.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>

#include <boost/parameter/preprocessor.hpp>

#include <feel/feelmodels/levelset/parameter_map.hpp>

#include <boost/bimap.hpp>

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
enum class LevelSetDistanceMethod { NONE, FASTMARCHING, HAMILTONJACOBI, RENORMALISATION };

enum class LevelSetMeasuresExported
{
    Volume, Perimeter, Position_COM, Velocity_COM
};
enum class LevelSetFieldsExported
{
    GradPhi, ModGradPhi, 
    AdvectionVelocity,
    BackwardCharacteristics, CauchyGreenInvariant1, CauchyGreenInvariant2
};

template<
    typename ConvexType, typename BasisType, typename PeriodicityType = NoPeriodicity, 
    typename FunctionSpaceAdvectionVelocityType = FunctionSpace< Mesh<ConvexType>, bases<typename detail::ChangeBasisPolySet<Vectorial,BasisType>::type>, Periodicity<PeriodicityType> >,
    typename BasisPnType = BasisType
    >
class LevelSet : public ModelNumerical,
                 public std::enable_shared_from_this< LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType> >
{
    typedef ModelNumerical super_type;
public:
    typedef LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    static const uint16_type Order = BasisType::nOrder;
    typedef double value_type;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type mymesh_type;

    //--------------------------------------------------------------------//
    // Periodicity
    typedef PeriodicityType periodicity_type;
    //--------------------------------------------------------------------//
    // Function space manager
    typedef LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_space_manager_type;
    typedef std::shared_ptr<levelset_space_manager_type> levelset_space_manager_ptrtype;
    //--------------------------------------------------------------------//
    // Function spaces and elements
    // levelset
    typedef typename levelset_space_manager_type::basis_scalar_type basis_levelset_type;
    typedef typename levelset_space_manager_type::space_scalar_type space_levelset_type;
    typedef typename levelset_space_manager_type::space_scalar_ptrtype space_levelset_ptrtype;
    typedef typename space_levelset_type::element_type element_levelset_type;
    typedef std::shared_ptr<element_levelset_type> element_levelset_ptrtype;
    // levelset PN
    typedef typename levelset_space_manager_type::basis_scalar_PN_type basis_levelset_PN_type;
    typedef typename levelset_space_manager_type::space_scalar_PN_type space_levelset_PN_type;
    typedef typename levelset_space_manager_type::space_scalar_PN_ptrtype space_levelset_PN_ptrtype;
    typedef typename space_levelset_PN_type::element_type element_levelset_PN_type;
    typedef std::shared_ptr<element_levelset_PN_type> element_levelset_PN_ptrtype;
    // vectorial
    typedef typename levelset_space_manager_type::basis_vectorial_type basis_vectorial_type;
    typedef typename levelset_space_manager_type::space_vectorial_type space_vectorial_type;
    typedef typename levelset_space_manager_type::space_vectorial_ptrtype space_vectorial_ptrtype;
    typedef typename space_vectorial_type::element_type element_vectorial_type;
    typedef std::shared_ptr< element_vectorial_type > element_vectorial_ptrtype;
    // markers P0
    typedef typename levelset_space_manager_type::basis_markers_type basis_markers_type;
    typedef typename levelset_space_manager_type::space_markers_type space_markers_type;
    typedef typename levelset_space_manager_type::space_markers_ptrtype space_markers_ptrtype;
    typedef typename space_markers_type::element_type element_markers_type;
    typedef std::shared_ptr<element_markers_type> element_markers_ptrtype;
    // tensor2symm
    typedef typename levelset_space_manager_type::basis_tensor2symm_type basis_tensor2symm_type;
    typedef typename levelset_space_manager_type::space_tensor2symm_type space_tensor2symm_type;
    typedef typename levelset_space_manager_type::space_tensor2symm_ptrtype space_tensor2symm_ptrtype;
    typedef typename space_tensor2symm_type::element_type element_tensor2symm_type;
    typedef std::shared_ptr< element_tensor2symm_type > element_tensor2symm_ptrtype;

#if defined (MESH_ADAPTATION_LS)
    typedef MeshAdaptation<Dim, Order, 1, periodicity_type > mesh_adaptation_type;
    typedef std::shared_ptr< mesh_adaptation_type > mesh_adaptation_ptrtype;
#endif
    //--------------------------------------------------------------------//
    // Advection toolbox
    typedef AdvDiffReac< space_levelset_type, FunctionSpaceAdvectionVelocityType > advection_toolbox_type;
    typedef std::shared_ptr<advection_toolbox_type> advection_toolbox_ptrtype;
    static_assert( advection_toolbox_type::is_scalar, "LevelSet function basis must be scalar" );

    typedef typename advection_toolbox_type::space_advection_velocity_type space_advection_velocity_type;
    typedef typename advection_toolbox_type::space_advection_velocity_ptrtype space_advection_velocity_ptrtype;
    typedef typename advection_toolbox_type::element_advection_velocity_type element_advection_velocity_type;
    typedef typename advection_toolbox_type::element_advection_velocity_ptrtype element_advection_velocity_ptrtype;

    //--------------------------------------------------------------------//
    // Range types
    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    typedef typename MeshTraits<mesh_type>::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef faces_reference_wrapper_t<mesh_type> range_faces_type;

    //--------------------------------------------------------------------//
    // Stretch and shear types
    typedef element_levelset_type element_stretch_type;
    typedef element_levelset_ptrtype element_stretch_ptrtype;

    //--------------------------------------------------------------------//
    // Tool manager
    typedef LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_tool_manager_type;
    typedef std::shared_ptr<levelset_tool_manager_type> levelset_tool_manager_ptrtype;
    //--------------------------------------------------------------------//
    // Projectors
    typedef Projector<space_levelset_type, space_levelset_type> projector_levelset_type;
    typedef std::shared_ptr<projector_levelset_type> projector_levelset_ptrtype;
    
    typedef Projector<space_vectorial_type, space_vectorial_type> projector_levelset_vectorial_type;
    typedef std::shared_ptr<projector_levelset_vectorial_type> projector_levelset_vectorial_ptrtype;

    typedef Projector<space_tensor2symm_type, space_tensor2symm_type> projector_tensor2symm_type;
    typedef std::shared_ptr<projector_tensor2symm_type> projector_tensor2symm_ptrtype;

    //--------------------------------------------------------------------//
    // Reinitialization
    typedef Reinitializer<space_levelset_type> reinitializer_type;
    typedef std::shared_ptr<reinitializer_type> reinitializer_ptrtype;
    typedef ReinitializerFM<space_levelset_type> reinitializerFM_type;
    typedef std::shared_ptr<reinitializerFM_type> reinitializerFM_ptrtype;
    typedef ReinitializerHJ<space_levelset_type> reinitializerHJ_type;
    typedef std::shared_ptr<reinitializerHJ_type> reinitializerHJ_ptrtype;

    enum class FastMarchingInitializationMethod { 
        NONE=0, ILP, SMOOTHED_ILP, HJ_EQ, IL_HJ_EQ
    };

    typedef boost::bimap<std::string, FastMarchingInitializationMethod> fastmarchinginitializationmethodidmap_type;

    //--------------------------------------------------------------------//
    // Initial value
    enum class ShapeType {
        SPHERE, ELLIPSE
    };
    static std::map<std::string, ShapeType> ShapeTypeMap;

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //--------------------------------------------------------------------//
    // Derivation methods
    enum class DerivationMethod { 
        NODAL_PROJECTION, L2_PROJECTION, SMOOTH_PROJECTION, PN_NODAL_PROJECTION
    };
    typedef boost::bimap<std::string, DerivationMethod> derivationmethod_maptype;
    static const derivationmethod_maptype DerivationMethodMap;

    //--------------------------------------------------------------------//
    // ModGradPhi advection
    typedef basis_levelset_type basis_modgradphi_advection_type;
    typedef AdvDiffReac< space_levelset_type, space_advection_velocity_type > modgradphi_advection_type;
    typedef std::shared_ptr<modgradphi_advection_type> modgradphi_advection_ptrtype;
    // Stretch advection
    typedef basis_levelset_type basis_stretch_advection_type;
    typedef AdvDiffReac< space_levelset_type, space_advection_velocity_type > stretch_advection_type;
    typedef std::shared_ptr<stretch_advection_type> stretch_advection_ptrtype;

    //--------------------------------------------------------------------//
    // Backward characteristics advection
    typedef basis_vectorial_type basis_backwardcharacteristics_advection_type;
    typedef AdvDiffReac< space_vectorial_type, space_advection_velocity_type > backwardcharacteristics_advection_type;
    typedef std::shared_ptr<backwardcharacteristics_advection_type> backwardcharacteristics_advection_ptrtype;
    typedef typename backwardcharacteristics_advection_type::element_advection_type element_backwardcharacteristics_type;
    typedef std::shared_ptr<element_backwardcharacteristics_type> element_backwardcharacteristics_ptrtype;
    // Cauchy-Green tensor invariants types
    typedef element_levelset_type element_cauchygreen_invariant_type;
    typedef std::shared_ptr<element_cauchygreen_invariant_type> element_cauchygreen_invariant_ptrtype;

    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mymesh_type, nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;
    typedef exporter_ptrtype exporter_manager_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    LevelSet(
            std::string const& prefix,
            worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );


    LevelSet( self_type const& L ) = default;

    static self_ptrtype New( 
            std::string const& prefix,
            worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );

    //--------------------------------------------------------------------//
    // Initialization
    void init();
    void initLevelsetValue();
    void initPostProcess();

    std::shared_ptr<std::ostringstream> getInfo() const;

    //--------------------------------------------------------------------//
    // Advection data
    typename advection_toolbox_type::bdf_ptrtype timeStepBDF() { return M_advectionToolbox->timeStepBDF(); }
    typename advection_toolbox_type::bdf_ptrtype /*const&*/ timeStepBDF() const { return M_advectionToolbox->timeStepBDF(); }
    std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    template<typename ExprT>
    void
    updateAdvectionVelocity(vf::Expr<ExprT> const& v_expr) { M_advectionToolbox->updateAdvectionVelocity( v_expr ); }
    void updateAdvectionVelocity( element_advection_velocity_ptrtype const& velocity ) { return M_advectionToolbox->updateAdvectionVelocity( velocity ); }
    void updateAdvectionVelocity( element_advection_velocity_type const& velocity ) { return M_advectionToolbox->updateAdvectionVelocity( velocity ); }
    //--------------------------------------------------------------------//
    // Spaces
    levelset_space_manager_ptrtype const& functionSpaceManager() const { return M_spaceManager; }
    void setFunctionSpaceManager( levelset_space_manager_ptrtype const& manager ) { M_spaceManager = manager; }

    space_levelset_ptrtype const& functionSpace() const { return M_spaceLevelset; }
    space_markers_ptrtype const& functionSpaceMarkers() const { return M_spaceMarkers; }
    space_vectorial_ptrtype const& functionSpaceVectorial() const { return M_spaceVectorial; }
    space_tensor2symm_ptrtype const& functionSpaceTensor2Symm() const { return M_spaceTensor2Symm; }

    space_advection_velocity_ptrtype const& functionSpaceAdvectionVelocity() const { return M_spaceAdvectionVelocity; }
    void setFunctionSpaceAdvectionVelocity( space_advection_velocity_ptrtype const& space ) { M_spaceAdvectionVelocity = space; }

    mesh_ptrtype const& mesh() const { return this->functionSpace()->mesh(); }

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"LevelsetMesh.path"); }

    mesh_ptrtype const& submeshDirac() const;
    mesh_ptrtype const& submeshOuter( double cut = 0.999 ) const;
    mesh_ptrtype const& submeshInner( double cut = 1e-3 ) const;

    //--------------------------------------------------------------------//
    // Mesh adaptation
#if defined (MESH_ADAPTATION_LS)
    mesh_ptrtype adaptMesh(double time, element_levelset_ptrtype elt);
    mesh_adaptation_ptrtype mesh_adapt;
    template<typename ExprT> mesh_ptrtype adaptedMeshFromExp(ExprT expr);
#endif

    //--------------------------------------------------------------------//
    // Levelset
    element_levelset_ptrtype & phi() { return M_advectionToolbox->fieldSolutionPtr(); }
    element_levelset_ptrtype const& phi() const { return M_advectionToolbox->fieldSolutionPtr(); }
    //element_levelset_ptrtype const& phinl() const { return M_phinl; }
    element_levelset_PN_ptrtype const& phiPN() const;
    element_vectorial_ptrtype const& gradPhi() const;
    element_levelset_ptrtype const& modGradPhi() const;
    element_levelset_ptrtype const& distance() const;

    element_stretch_ptrtype const& stretch() const;
    element_backwardcharacteristics_ptrtype const& backwardCharacteristics() const;

    element_levelset_ptrtype const& heaviside() const;
    element_levelset_ptrtype const& H() const { return this->heaviside(); }
    element_levelset_ptrtype const& dirac() const;
    element_levelset_ptrtype const& D() const { return this->dirac(); }

    element_vectorial_ptrtype const& normal() const;
    element_vectorial_ptrtype const& N() const { return this->normal(); }
    element_levelset_ptrtype const& curvature() const;
    element_levelset_ptrtype const& K() const { return this->curvature(); }
    element_vectorial_ptrtype const& distanceNormal() const;
    element_levelset_ptrtype const& distanceCurvature() const;

    void updateInterfaceQuantities();

    double thicknessInterface() const { return M_thicknessInterface; }
    void setThicknessInterface( double value ) { M_thicknessInterface = value; }

    int iterSinceReinit() const { return M_iterSinceReinit; }

    //--------------------------------------------------------------------//
    // Tools
    levelset_tool_manager_ptrtype const& toolManager() const { return M_toolManager; }
    void setToolManager( levelset_tool_manager_ptrtype const& manager ) { M_toolManager = manager; }

    projector_levelset_ptrtype const& projectorL2() const { return M_projectorL2Scalar; }
    projector_levelset_vectorial_ptrtype const& projectorL2Vectorial() const { return M_projectorL2Vectorial; }
    projector_tensor2symm_ptrtype const& projectorL2Tensor2Symm() const { return M_projectorL2Tensor2Symm; }

    projector_levelset_ptrtype const& smoother() const { return M_projectorSMScalar; }
    projector_levelset_vectorial_ptrtype const& smootherVectorial() const { return M_projectorSMVectorial; }
    projector_levelset_ptrtype const& smootherInterface() const;
    projector_levelset_vectorial_ptrtype const& smootherInterfaceVectorial() const;

    //--------------------------------------------------------------------//
    // Cauchy-Green tensor related quantities
    auto leftCauchyGreenTensorExpr() const;
    element_tensor2symm_ptrtype const& leftCauchyGreenTensor() const;
    auto cauchyGreenInvariant1Expr() const;
    element_cauchygreen_invariant_ptrtype const& cauchyGreenInvariant1() const;
    auto cauchyGreenInvariant2Expr() const;
    element_cauchygreen_invariant_ptrtype const& cauchyGreenInvariant2() const;

    //--------------------------------------------------------------------//
    // Markers
    element_markers_ptrtype const& markerInterface() const;
    element_markers_ptrtype const& markerDirac() const;
    element_markers_ptrtype const& markerOuter( double cut = 0.999 ) const;
    element_markers_ptrtype const& markerInner( double cut = 1e-3 ) const;
    element_markers_ptrtype const& markerHeaviside( double cut = 0.999 ) const { return this->markerOuter(cut); }
    element_markers_ptrtype const& markerCrossedElements() const;

    //--------------------------------------------------------------------//
    // Ranges
    range_elements_type const& rangeMeshElements() const { return M_advectionToolbox->rangeMeshElements(); }

    range_elements_type interfaceElements() const;
    range_elements_type outerElementsRange( double cut );
    range_elements_type outerElementsRange() { return outerElementsRange( -this->thicknessInterface() ); }
    range_elements_type const& rangeDiracElements() const;

    range_faces_type rangeInterfaceFaces() const;

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
    void advect( element_advection_velocity_ptrtype const& velocity );
    void advect( element_advection_velocity_type const& velocity );
    void solve();

    void updateTimeStep();

    //--------------------------------------------------------------------//
    // Reinitialization
    void reinitialize();
    element_levelset_type redistantiate( element_levelset_type const& phi, LevelSetDistanceMethod method ) const;

    void setFastMarchingInitializationMethod( FastMarchingInitializationMethod m );
    FastMarchingInitializationMethod fastMarchingInitializationMethod() { return M_fastMarchingInitializationMethod; }
    void setUseMarkerDiracAsMarkerDoneFM( bool val = true ) { M_useMarkerDiracAsMarkerDoneFM  = val; }

    reinitializer_ptrtype const& reinitializer() const { return M_reinitializer; }
    reinitializerFM_ptrtype const& reinitializerFM( bool buildOnTheFly = true );
    reinitializerHJ_ptrtype const& reinitializerHJ( bool buildOnTheFly = true );

    bool hasReinitialized() const { return M_hasReinitialized; }

    //--------------------------------------------------------------------//
    // Extension velocity
    template<typename ExprT>
    element_advection_velocity_type extensionVelocity( vf::Expr<ExprT> const& u ) const;

    //--------------------------------------------------------------------//
    // Initial value
    void setInitialValue(element_levelset_ptrtype const& phiv, bool doReinitialize);
    void setInitialValue(element_levelset_ptrtype const& phiv)
    {
        this->setInitialValue(phiv, M_reinitInitialValue);
    }
    template<typename ExprT>
    void setInitialValue(vf::Expr<ExprT> const& expr, bool doReinitialize)
    {
        auto phi_init = this->functionSpace()->elementPtr();
        phi_init->on( 
                _range=this->rangeMeshElements(),
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
    void exportResults( bool save = true ) { this->exportResults( this->currentTime(), save ); }
    void exportResults( double time, bool save = true );
    bool hasPostProcessMeasureExported( LevelSetMeasuresExported const& measure) const;
    bool hasPostProcessFieldExported( LevelSetFieldsExported const& field) const;
    exporter_ptrtype & getExporter() { return M_exporter; }
    exporter_ptrtype const& getExporter() const { return M_exporter; }
    //--------------------------------------------------------------------//
    // Physical quantities
    double volume() const;
    double perimeter() const;
    auto positionCOM() const;
    auto velocityCOM() const;

protected:
    //--------------------------------------------------------------------//
    void buildImpl();
    //--------------------------------------------------------------------//
    // Levelset data update functions
    void updateGradPhi();
    void updateModGradPhi();
    void updateDirac();
    void updateHeaviside();

    void updateNormal();
    void updateCurvature();

    void updatePhiPN();

    void updateDistance();
    void updateDistanceNormal();
    void updateDistanceCurvature();

    void updateMarkerDirac();
    void markerHeavisideImpl( element_markers_ptrtype const& marker, bool invert, double cut );
    void updateMarkerCrossedElements();
    void updateMarkerInterface();

    void updateLeftCauchyGreenTensor();

private:
    void loadParametersFromOptionsVm();
    void loadConfigICFile();
    void loadConfigBCFile();
    void loadConfigPostProcess();

    void createFunctionSpaces();
    void createInterfaceQuantities();
    void createReinitialization();
    void createTools();
    void createExporters();

    //--------------------------------------------------------------------//
    void addShape( 
            std::pair<ShapeType, parameter_map> const& shape, 
            element_levelset_type & phi );

    //--------------------------------------------------------------------//
    element_levelset_ptrtype const& phio() const { return this->timeStepBDF()->unknowns()[1]; }

    //--------------------------------------------------------------------//
    // Interface rectangular function
    element_levelset_type interfaceRectangularFunction() const { return this->interfaceRectangularFunction(this->phi()); }
    element_levelset_type interfaceRectangularFunction( element_levelset_ptrtype const& p ) const { return this->interfaceRectangularFunction(*p); }
    element_levelset_type interfaceRectangularFunction( element_levelset_type const& p ) const;
    //--------------------------------------------------------------------//
    // Export
    void exportResultsImpl( double time, bool save );
    void exportMeasuresImpl( double time, bool save );
    // Save
    void saveCurrent() const;


protected:
    //--------------------------------------------------------------------//
    // Interface quantities update flags
    mutable bool M_doUpdateDirac;
    mutable bool M_doUpdateHeaviside;
    mutable bool M_doUpdateInterfaceElements;
    mutable bool M_doUpdateRangeDiracElements;
    mutable bool M_doUpdateInterfaceFaces;
    mutable bool M_doUpdateSmootherInterface;
    mutable bool M_doUpdateSmootherInterfaceVectorial;
    mutable bool M_doUpdateNormal;
    mutable bool M_doUpdateCurvature;
    mutable bool M_doUpdateGradPhi;
    mutable bool M_doUpdateModGradPhi;
    mutable bool M_doUpdatePhiPN;
    mutable bool M_doUpdateDistance;
    mutable bool M_doUpdateDistanceNormal;
    mutable bool M_doUpdateDistanceCurvature;

    //--------------------------------------------------------------------//
    // Levelset initial value
    map_scalar_field<2> M_icDirichlet;
    std::vector<std::pair<ShapeType, parameter_map>> M_icShapes;

    //--------------------------------------------------------------------//
    // Boundary conditions
    std::list<std::string> M_bcMarkersInflow;

private:
    //--------------------------------------------------------------------//
    // Meshes 
    mutable mesh_ptrtype M_submeshDirac;
    mutable bool M_doUpdateSubmeshDirac;
    mutable mesh_ptrtype M_submeshOuter;
    mutable bool M_doUpdateSubmeshOuter;
    mutable mesh_ptrtype M_submeshInner;
    mutable bool M_doUpdateSubmeshInner;

    //--------------------------------------------------------------------//
    // Periodicity
    periodicity_type M_periodicity;
    // Advection toolbox
    advection_toolbox_ptrtype M_advectionToolbox;
    bool M_doExportAdvection;
    //--------------------------------------------------------------------//
    // Spaces
    levelset_space_manager_ptrtype M_spaceManager;
    bool M_useSpaceIsoPN;
    space_levelset_ptrtype M_spaceLevelset;
    space_vectorial_ptrtype M_spaceVectorial;
    space_markers_ptrtype M_spaceMarkers;
    space_tensor2symm_ptrtype M_spaceTensor2Symm;
    space_advection_velocity_ptrtype M_spaceAdvectionVelocity;

    //--------------------------------------------------------------------//
    // Markers
    mutable element_markers_ptrtype M_markerDirac;
    mutable element_markers_ptrtype M_markerOuter;
    mutable double M_markerOuterCut;
    mutable element_markers_ptrtype M_markerInner;
    mutable double M_markerInnerCut;
    mutable element_markers_ptrtype M_markerCrossedElements;
    mutable element_markers_ptrtype M_markerInterface;
    bool M_doUpdateMarkers;

    //--------------------------------------------------------------------//
    // Ranges
    mutable range_elements_type M_rangeDiracElements;
    //--------------------------------------------------------------------//
    // Tools (projectors)
    levelset_tool_manager_ptrtype M_toolManager;

    projector_levelset_ptrtype M_projectorL2Scalar;
    projector_levelset_vectorial_ptrtype M_projectorL2Vectorial;
    projector_tensor2symm_ptrtype M_projectorL2Tensor2Symm;

    projector_levelset_ptrtype M_projectorSMScalar;
    projector_levelset_vectorial_ptrtype M_projectorSMVectorial;
    mutable projector_levelset_ptrtype M_smootherInterface;
    mutable projector_levelset_vectorial_ptrtype M_smootherInterfaceVectorial;
    //--------------------------------------------------------------------//
    // Levelset data
    mutable element_levelset_PN_ptrtype M_levelsetPhiPN;
    mutable element_vectorial_ptrtype M_levelsetGradPhi;
    mutable element_levelset_ptrtype M_levelsetModGradPhi;
    mutable element_levelset_ptrtype M_heaviside;
    mutable element_levelset_ptrtype M_dirac;
    mutable element_levelset_ptrtype M_distance;

    mutable range_elements_type M_interfaceElements;
    mutable range_faces_type M_interfaceFaces;
    //--------------------------------------------------------------------//
    // Normal, curvature
    mutable element_vectorial_ptrtype M_levelsetNormal;
    mutable element_levelset_ptrtype M_levelsetCurvature;
    mutable element_vectorial_ptrtype M_distanceNormal;
    mutable element_levelset_ptrtype M_distanceCurvature;
    //--------------------------------------------------------------------//
    // Advection

    //--------------------------------------------------------------------//
    // Derivation methods
    DerivationMethod M_gradPhiMethod;
    DerivationMethod M_curvatureMethod;

    //--------------------------------------------------------------------//
    // ModGradPhi advection
    bool M_useGradientAugmented;
    bool M_reinitGradientAugmented;
    bool M_reinitStretchAugmented;
    modgradphi_advection_ptrtype M_modGradPhiAdvection;

    //--------------------------------------------------------------------//
    // Stretch advection
    bool M_useStretchAugmented;
    stretch_advection_ptrtype M_stretchAdvection;
    mutable element_stretch_ptrtype M_levelsetStretch;

    //--------------------------------------------------------------------//
    // Backward characteristics advection
    bool M_useCauchyAugmented;
    backwardcharacteristics_advection_ptrtype M_backwardCharacteristicsAdvection;
    vector_field_expression<nDim> M_initialBackwardCharacteristics;
    bool M_hasInitialBackwardCharacteristics;
    // Cauchy-Green tensor
    element_tensor2symm_ptrtype M_leftCauchyGreenTensor;
    mutable bool M_doUpdateCauchyGreenTensor;
    // Cauchy-Green tensor invariants
    mutable element_cauchygreen_invariant_ptrtype M_cauchyGreenInvariant1;
    mutable bool M_doUpdateCauchyGreenInvariant1;
    mutable element_cauchygreen_invariant_ptrtype M_cauchyGreenInvariant2;
    mutable bool M_doUpdateCauchyGreenInvariant2;

    //--------------------------------------------------------------------//
    // Redistantiation
    static const std::map<std::string, LevelSetDistanceMethod> LevelSetDistanceMethodIdMap;

    LevelSetDistanceMethod M_distanceMethod;

    //--------------------------------------------------------------------//
    // Reinitialization
    reinitializer_ptrtype M_reinitializer;
    reinitializerFM_ptrtype M_reinitializerFM;
    reinitializerHJ_ptrtype M_reinitializerHJ;
    bool M_reinitializerIsUpdatedForUse;

    LevelSetDistanceMethod M_reinitMethod;
    FastMarchingInitializationMethod M_fastMarchingInitializationMethod;
    static const fastmarchinginitializationmethodidmap_type FastMarchingInitializationMethodIdMap;
    bool M_useMarkerDiracAsMarkerDoneFM;

    bool M_reinitInitialValue;

    bool M_hasReinitialized;
    int M_iterSinceReinit;
    // Vector that stores the iterSinceReinit of each time-step
    std::vector<int> M_vecIterSinceReinit;

    //--------------------------------------------------------------------//
    // Extension velocity
    bool M_useExtensionVelocity;
    double M_extensionVelocityNitscheGamma;
    mutable sparse_matrix_ptrtype M_extensionVelocityLHSMatrix;
    mutable vector_ptrtype M_extensionVelocityRHSVector;
    
    //--------------------------------------------------------------------//
    // Export
    exporter_ptrtype M_exporter;
    std::set<LevelSetMeasuresExported> M_postProcessMeasuresExported;
    std::set<LevelSetFieldsExported> M_postProcessFieldsExported;
    //--------------------------------------------------------------------//
    // Parameters
    double M_thicknessInterface;
    double M_thicknessInterfaceRectangularFunction;
    bool M_useAdaptiveThicknessInterface;
    bool M_useRegularPhi;
    bool M_useHeavisideDiracNodalProj;

    bool M_fixVolume;
    double M_initialVolume;

    //LevelSetTimeDiscretization M_discrMethod;

}; //class LevelSet

//----------------------------------------------------------------------------//
// Advection
template<typename ConvexType, typename BasisType, typename PeriodicityType, typename FunctionSpaceAdvectionVelocityType, typename BasisPnType>
template<typename ExprT>
void 
LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType>::advect(vf::Expr<ExprT> const& velocity)
{
    this->updateAdvectionVelocity(velocity);
    this->solve();
}
//----------------------------------------------------------------------------//
// Extension velocity
template<typename ConvexType, typename BasisType, typename PeriodicityType, typename FunctionSpaceAdvectionVelocityType, typename BasisPnType>
template<typename ExprT>
typename LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType>::element_advection_velocity_type
LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType>::extensionVelocity( vf::Expr<ExprT> const& u) const
{
    this->log("LevelSet", "extensionVelocity", "start");
    this->timerTool("Solve").start();

    auto const& spaceVelocity = this->functionSpaceAdvectionVelocity();
    auto const& spaceVelocityComp = spaceVelocity->compSpace();
    auto const& backendExtensionVelocity = backend(
                _name=prefixvm(this->prefix(),"extension-velocity"),
                _worldcomm=spaceVelocityComp->worldCommPtr()
            );

    if( !M_extensionVelocityLHSMatrix )
    {
        M_extensionVelocityLHSMatrix = backendExtensionVelocity->newMatrix(
                _trial=spaceVelocityComp,
                _test=spaceVelocityComp,
                _pattern=size_type(Pattern::COUPLED)
                );
    }
    if( !M_extensionVelocityRHSVector )
    {
        M_extensionVelocityRHSVector = backendExtensionVelocity->newVector(
                _test=spaceVelocityComp
                );
    }

    auto bilinearForm = form2( _trial=spaceVelocityComp, _test=spaceVelocityComp, _matrix=M_extensionVelocityLHSMatrix );
    auto linearForm = form1( _test=spaceVelocityComp, _vector=M_extensionVelocityRHSVector );
    auto Fext = spaceVelocityComp->element();

    auto const& phi = this->phi();
    //auto gradPhiExpr = trans(gradv(phi));
    auto gradPhiExpr = idv(this->gradPhi());;
    auto NExpr = gradPhiExpr / sqrt(trans(gradPhiExpr)*gradPhiExpr);

    bilinearForm = integrate(
            _range=this->rangeMeshElements(),
            _expr=(gradt(Fext)*gradPhiExpr)*id(Fext)
            );
    // BC at interface with penalty method
    bilinearForm += integrate(
            _range=this->rangeDiracElements(),
            _expr=M_extensionVelocityNitscheGamma/h() * idt(Fext)*id(Fext) * idv(this->dirac())
            );
    linearForm = integrate(
            _range=this->rangeDiracElements(),
            _expr=M_extensionVelocityNitscheGamma/h() * (trans(u)*NExpr)*id(Fext) * idv(this->dirac())
            );
    double timeElapsedAssembly = this->timerTool("Solve").stop();
    this->log("LevelSet", "extensionVelocity", "assembly finish in "+(boost::format("%1% s") %timeElapsedAssembly).str() );
    this->timerTool("Solve").start();

    //bilinearForm.solve( _rhs=linearForm, _solution=Fext, _name=prefixvm(this->prefix(),"extension-velocity") );
    backendExtensionVelocity->solve(
            _matrix=M_extensionVelocityLHSMatrix,
            _solution=Fext,
            _rhs=M_extensionVelocityRHSVector
            );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("LevelSet", "extensionVelocity", "finish in "+(boost::format("%1% s") %timeElapsed).str() );

    return vf::project(
            _space=spaceVelocity,
            _range=this->rangeMeshElements(),
            _expr=(idv(Fext)*NExpr - u)*idv(this->heaviside()) + u
            );
}

} // namespace FeelModels
} // namespace Feel

#endif
