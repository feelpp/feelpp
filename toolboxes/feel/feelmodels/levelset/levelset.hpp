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

#include <feel/feelmodels/levelset/levelsetbase.hpp>
#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>
#include <feel/feelmodels/levelset/levelsettoolmanager.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_fm.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_hj.hpp>
#include <feel/feelmodels/levelset/levelsetparticleinjector.hpp>

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

template<
    typename ConvexType, typename BasisType, typename PeriodicityType = NoPeriodicity, 
    typename FunctionSpaceAdvectionVelocityType = FunctionSpace< Mesh<ConvexType>, bases<typename detail::ChangeBasisPolySet<Vectorial,BasisType>::type>/*, Periodicity<PeriodicityType>*/ >,
    typename BasisPnType = BasisType
    >
class LevelSet : public LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType>,
                 public std::enable_shared_from_this< LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType> >
{
    typedef LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType> super_type;
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
    // scalar
    typedef typename levelset_space_manager_type::basis_scalar_type basis_scalar_type;
    typedef typename levelset_space_manager_type::space_scalar_type space_scalar_type;
    typedef typename levelset_space_manager_type::space_scalar_ptrtype space_scalar_ptrtype;
    typedef typename space_scalar_type::element_type element_scalar_type;
    typedef std::shared_ptr<element_scalar_type> element_scalar_ptrtype;
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
    // Heaviside and Dirac expressions
    typedef Expr< LevelsetDeltaExpr<element_levelset_type> > levelset_delta_expr_type;
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
    // Redistanciation
    typedef LevelSetRedistanciation<space_levelset_type> redistanciation_type;
    typedef std::shared_ptr<redistanciation_type> redistanciation_ptrtype;
    typedef LevelSetRedistanciationFM<space_levelset_type> redistanciationFM_type;
    typedef std::shared_ptr<redistanciationFM_type> redistanciationFM_ptrtype;
    typedef LevelSetRedistanciationHJ<space_levelset_type> redistanciationHJ_type;
    typedef std::shared_ptr<redistanciationHJ_type> redistanciationHJ_ptrtype;

    enum class FastMarchingInitializationMethod { 
        NONE=0, ILP_NODAL, ILP_L2, ILP_SMOOTH, HJ_EQ, IL_HJ_EQ
    };

    typedef boost::bimap<std::string, FastMarchingInitializationMethod> fastmarchinginitializationmethodidmap_type;

    //--------------------------------------------------------------------//
    // Particle injector
    typedef LevelSetParticleInjector<self_type> levelsetparticleinjector_type;
    typedef std::shared_ptr<levelsetparticleinjector_type> levelsetparticleinjector_ptrtype;

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

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
    using exporter_type = typename super_type::exporter_type;
    using exporter_ptrtype = typename super_type::exporter_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    LevelSet( std::string const& prefix,
              std::string const& keyword = "levelset",
              worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() );
    LevelSet( std::string const& prefix,
              worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() )
        : LevelSet( prefix, prefix, _worldComm, subPrefix, modelRep )
    {}

    LevelSet( self_type const& L ) = default;

    static self_ptrtype New( std::string const& prefix,
                             std::string const& keyword = "levelset",
                             worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }
    std::shared_ptr<self_type const> shared_from_this() const { return std::dynamic_pointer_cast<self_type const>( super_type::shared_from_this() ); }

    //--------------------------------------------------------------------//
    // Initialization
    void init();
    void initInitialValues();
    void initPostProcess() override;

    std::shared_ptr<std::ostringstream> getInfo() const override;

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
    space_advection_velocity_ptrtype const& functionSpaceAdvectionVelocity() const { return M_spaceAdvectionVelocity; }
    void setFunctionSpaceAdvectionVelocity( space_advection_velocity_ptrtype const& space ) { M_spaceAdvectionVelocity = space; }
    space_tensor2symm_ptrtype const& functionSpaceTensor2Symm() const { return M_spaceTensor2Symm; }

    //std::string fileNameMeshPath() const override { return prefixvm(this->prefix(),"LevelsetMesh.path"); }

    //--------------------------------------------------------------------//
    // Levelset
    element_levelset_ptrtype & phi() { return this->phiPtr(); }
    element_levelset_ptrtype const& phi() const { return this->phiPtr(); }
    element_stretch_ptrtype const& stretch() const;
    element_backwardcharacteristics_ptrtype const& backwardCharacteristics() const;

    void updateInterfaceQuantities() override;

    int iterSinceRedistanciation() const { return M_iterSinceRedistanciation; }

    //--------------------------------------------------------------------//
    // Tools
    projector_tensor2symm_ptrtype const& projectorL2Tensor2Symm() const { return M_projectorL2Tensor2Symm; }
    //--------------------------------------------------------------------//
    // Redistanciation
    void redistanciate() override;

    bool useOrder1AfterRedist() const { return M_useOrder1AfterRedist; }
    void setUseOrder1AfterRedist( bool b ) { M_useOrder1AfterRedist = b; }

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
    element_markers_ptrtype const& markerCrossedElements() const;

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
    // Extension velocity
    template<typename ExprT>
    element_advection_velocity_type extensionVelocity( vf::Expr<ExprT> const& u ) const;

    //--------------------------------------------------------------------//
    // Symbols expressions
    using super_type::symbolsExpr;
    // Fields
#if 0
    auto allFields( std::string const& prefix = "" ) const
    {
        return hana::concat( super_type::allFields( prefix ), hana::make_tuple(
                std::make_pair( prefixvm( prefix, "advection-velocity" ), M_advectionToolbox->fieldAdvectionVelocityPtr() ),
                this->optionalScalarFields( prefix ),
                this->optionalVectorialFields( prefix )
                ) );
    }
#else
    auto modelFields( std::string const& prefix = "" ) const
        {
            return super_type::modelFields( prefix );
        }

#endif
    auto optionalScalarFields( std::string const& prefix = "" ) const
    {
        std::map<std::string, element_scalar_ptrtype> fields;
        if( M_useStretchAugmented )
        {
            fields[prefixvm( prefix, "stretch" )] = this->stretch();
        }
        if( M_useCauchyAugmented )
        {
            fields[prefixvm( prefix, "cauchygreeninvariant1" )] = this->cauchyGreenInvariant1();
            fields[prefixvm( prefix, "cauchygreeninvariant2" )] = this->cauchyGreenInvariant2();
        }
        return fields;
    }
    auto optionalVectorialFields( std::string const& prefix = "" ) const
    {
        std::map<std::string, element_vectorial_ptrtype> fields;
        if( M_useCauchyAugmented )
        {
            fields[prefixvm( prefix, "backwardcharacteristics" )] = M_backwardCharacteristicsAdvection->fieldSolutionPtr();
        }
        return fields;
    }
    // Measures quantities
    auto allMeasuresQuantities( std::string const& prefix = "" ) const
    {
        // TODO : we need to explicitly convert Eigen::Matrix to std::vector as
        // the begin and end iterators used in ModelNumerical::588 are only implemented from
        // Eigen v3.4 on (not released yet).
        auto eigenToVec = []( Eigen::Matrix<value_type, nDim, 1> const& m )
        {
            std::vector<value_type> v( m.size() );
            Eigen::Matrix<value_type, nDim, 1>::Map( &v[0], v.size() ) = m;
            return v;
        };
        return hana::concat( super_type::allMeasuresQuantities( prefix ), hana::make_tuple(
                    ModelMeasuresQuantity( prefix, "velocity-com", std::bind( eigenToVec, std::bind( &self_type::velocityCOM, this ) ) )
                    ) );
    }
    //--------------------------------------------------------------------//
    // Export
    std::set<std::string> postProcessSaveAllFieldsAvailable() const override;
    std::set<std::string> postProcessExportsAllFieldsAvailable() const override;

    using super_type::exportResults;
    void exportResults( double time ) override;
    template<typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr ) { this->exportResults( time, symbolsExpr, this->modelFields(), this->allMeasuresQuantities() ); }
    //--------------------------------------------------------------------//
    // Physical quantities
    Eigen::Matrix<value_type, nDim, 1> velocityCOM() const { 
        return integrate( 
            _range=this->rangeMeshElements(), 
            _expr=idv(M_advectionToolbox->fieldAdvectionVelocity()) * (1.-idv(this->H()))
            ).evaluate() / this->volume();
    }

protected:
    //--------------------------------------------------------------------//
    void buildImpl();
    //--------------------------------------------------------------------//
    // Levelset data update functions
    void updateMarkerCrossedElements();

    void updateLeftCauchyGreenTensor();

private:
    void loadParametersFromOptionsVm();
    void initBoundaryConditions();

    void createFunctionSpaces();
    void createInterfaceQuantities();
    void createTools();
    void createExporters();

    //--------------------------------------------------------------------//
    element_levelset_ptrtype const& phiPreviousTimeStepPtr() const { return this->timeStepBDF()->unknowns()[1]; }

    //--------------------------------------------------------------------//
    // Save
    void saveCurrent() const;


protected:
    //--------------------------------------------------------------------//
    // Boundary conditions
    std::list<std::string> M_bcMarkersInflow;

private:
    //--------------------------------------------------------------------//
    // Advection toolbox
    advection_toolbox_ptrtype M_advectionToolbox;
    bool M_doExportAdvection;
    //--------------------------------------------------------------------//
    // Spaces
    space_advection_velocity_ptrtype M_spaceAdvectionVelocity;
    space_tensor2symm_ptrtype M_spaceTensor2Symm;

    //--------------------------------------------------------------------//
    // Markers
    mutable element_markers_ptrtype M_markerCrossedElements;

    //--------------------------------------------------------------------//
    // Ranges
    //--------------------------------------------------------------------//
    // Tools (projectors)
    projector_tensor2symm_ptrtype M_projectorL2Tensor2Symm;
    //--------------------------------------------------------------------//
    // Levelset data
    //--------------------------------------------------------------------//
    // Normal, curvature
    //--------------------------------------------------------------------//
    // Advection

    //--------------------------------------------------------------------//
    // Derivation methods

    //--------------------------------------------------------------------//
    // Curvature diffusion

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
    int M_iterSinceRedistanciation;
    // Vector that stores the iterSinceRedistanciation of each time-step
    std::vector<int> M_vecIterSinceRedistanciation;

    bool M_useOrder1AfterRedist;

    //--------------------------------------------------------------------//
    // Extension velocity
    bool M_useExtensionVelocity;
    double M_extensionVelocityNitscheGamma;
    mutable sparse_matrix_ptrtype M_extensionVelocityLHSMatrix;
    mutable vector_ptrtype M_extensionVelocityRHSVector;

    //--------------------------------------------------------------------//
    // Particle injectors
    std::vector<levelsetparticleinjector_ptrtype> M_levelsetParticleInjectors;
    
    //--------------------------------------------------------------------//
    // Export
    //--------------------------------------------------------------------//
    // Parameters
    bool M_fixVolume;
    bool M_fixArea;

    //LevelSetTimeDiscretization M_discrMethod;

}; //class LevelSet

#ifndef LEVELSET_CLASS_TEMPLATE_DECLARATIONS
#define LEVELSET_CLASS_TEMPLATE_DECLARATIONS \
    template< typename ConvexType, typename BasisType, typename PeriodicityType, typename FunctionSpaceAdvectionVelocityType, typename BasisPnType > \
        /**/
#endif
#ifndef LEVELSET_CLASS_TEMPLATE_TYPE
#define LEVELSET_CLASS_TEMPLATE_TYPE \
    LevelSet<ConvexType, BasisType, PeriodicityType, FunctionSpaceAdvectionVelocityType, BasisPnType> \
        /**/
#endif
//----------------------------------------------------------------------------//
// Advection
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
template<typename ExprT>
void 
LEVELSET_CLASS_TEMPLATE_TYPE::advect(vf::Expr<ExprT> const& velocity)
{
    this->updateAdvectionVelocity(velocity);
    this->solve();
}
//----------------------------------------------------------------------------//
// Extension velocity
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
template<typename ExprT>
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_advection_velocity_type
LEVELSET_CLASS_TEMPLATE_TYPE::extensionVelocity( vf::Expr<ExprT> const& u) const
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

    //auto const& phi = this->phi();
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
