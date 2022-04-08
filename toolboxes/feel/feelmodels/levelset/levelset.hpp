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
#ifndef FEELPP_TOOLBOXES_LEVELSET_HPP
#define FEELPP_TOOLBOXES_LEVELSET_HPP 1

#include <fmt/core.h>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/projector.hpp>

#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

#include <feel/feelmodels/levelset/levelsetbase.hpp>
#include <feel/feelmodels/levelset/levelsetboundaryconditions.hpp>
#include <feel/feelmodels/levelset/levelsetparticleinjector.hpp>

#include <feel/feelfilters/straightenmesh.hpp>

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
    typename BasisPnType = BasisType
    >
class LevelSet: 
    public LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType>,
    public ModelPhysics<ConvexType::nDim>
{
    typedef LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType> super_type;
public:
    typedef LevelSet<ConvexType, BasisType, PeriodicityType, BasisPnType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    static constexpr uint16_type Order = BasisType::nOrder;
    typedef double value_type;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static constexpr uint16_type nDim = convex_type::nDim;
    static constexpr uint16_type nOrderGeo = convex_type::nOrder;
    static constexpr uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // Materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;
    // Periodicity
    typedef PeriodicityType periodicity_type;
    // Function space manager
    typedef LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_space_manager_type;
    typedef std::shared_ptr<levelset_space_manager_type> levelset_space_manager_ptrtype;
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

    static_assert( space_levelset_type::is_scalar, "LevelSetBase function basis must be scalar" );

    // Cached fields
    using cached_scalar_field_type = CachedModelField<element_scalar_type, InplaceUpdatePolicy>;
    using cached_vectorial_field_type = CachedModelField<element_vectorial_type, InplaceUpdatePolicy>;

    // Heaviside and Dirac expressions
    typedef Expr< LevelsetDeltaExpr<element_levelset_type> > levelset_delta_expr_type;

    // Advection toolbox
    typedef CoefficientFormPDE<convex_type, basis_scalar_type> cfpde_toolbox_type;
    typedef std::shared_ptr<cfpde_toolbox_type> cfpde_toolbox_ptrtype;
    using cfpde_infos_type = typename ModelPhysicCoefficientFormPDE<nDim>::infos_type;
    using CFPDECoefficient = typename cfpde_toolbox_type::Coefficient;

    // Range types
    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    typedef typename MeshTraits<mesh_type>::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef faces_reference_wrapper_t<mesh_type> range_faces_type;

    // Stretch and shear types
    typedef element_levelset_type element_stretch_type;
    typedef element_levelset_ptrtype element_stretch_ptrtype;

    // Tool manager
    typedef LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_tool_manager_type;
    typedef std::shared_ptr<levelset_tool_manager_type> levelset_tool_manager_ptrtype;
    // Projectors
    typedef Projector<space_levelset_type, space_levelset_type> projector_levelset_type;
    typedef std::shared_ptr<projector_levelset_type> projector_levelset_ptrtype;
    
    typedef Projector<space_vectorial_type, space_vectorial_type> projector_levelset_vectorial_type;
    typedef std::shared_ptr<projector_levelset_vectorial_type> projector_levelset_vectorial_ptrtype;

    typedef Projector<space_tensor2symm_type, space_tensor2symm_type> projector_tensor2symm_type;
    typedef std::shared_ptr<projector_tensor2symm_type> projector_tensor2symm_ptrtype;

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

    // Particle injector
    typedef LevelSetParticleInjector<self_type> levelsetparticleinjector_type;
    typedef std::shared_ptr<levelsetparticleinjector_type> levelsetparticleinjector_ptrtype;

    // Backend
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    // Algebraic data
    using ModelAlgebraic::DataUpdateLinear;
    using ModelAlgebraic::DataNewtonInitialGuess;
    using ModelAlgebraic::DataUpdateJacobian;
    using ModelAlgebraic::DataUpdateResidual;

    // ModGradPhi advection
    typedef CoefficientFormPDE<convex_type, basis_scalar_type> modgradphi_advection_type;
    typedef std::shared_ptr<modgradphi_advection_type> modgradphi_advection_ptrtype;
    // Stretch advection
    typedef CoefficientFormPDE<convex_type, basis_scalar_type> stretch_advection_type;
    typedef std::shared_ptr<stretch_advection_type> stretch_advection_ptrtype;

    // Backward characteristics advection
    typedef basis_vectorial_type basis_backwardcharacteristics_advection_type;
    typedef CoefficientFormPDE<convex_type, basis_vectorial_type> backwardcharacteristics_advection_type;
    typedef std::shared_ptr<backwardcharacteristics_advection_type> backwardcharacteristics_advection_ptrtype;
    typedef typename backwardcharacteristics_advection_type::element_unknown_type element_backwardcharacteristics_type;
    typedef std::shared_ptr<element_backwardcharacteristics_type> element_backwardcharacteristics_ptrtype;
    // Cauchy-Green tensor invariants types
    typedef element_scalar_type element_cauchygreen_invariant_type;
    typedef std::shared_ptr<element_cauchygreen_invariant_type> element_cauchygreen_invariant_ptrtype;

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

    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory() { M_advectionToolbox->initAlgebraicFactory(); this->setAlgebraicFactory( M_advectionToolbox->algebraicFactory() ); }
    void initPostProcess() override;

    // Infos
    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    //--------------------------------------------------------------------//
    // Initial condition
    void updateInitialConditions();

    //--------------------------------------------------------------------//
    // Physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    void updateParameterValues() override;
    void setParameterValues( std::map<std::string,double> const& paramValues ) override;

    //--------------------------------------------------------------------//
    // Advection data
    typename cfpde_toolbox_type::bdf_unknown_ptrtype /*const&*/ timeStepBDF() const { return M_advectionToolbox->timeStepBdfUnknown(); }
    std::shared_ptr<TSBase> timeStepBase() { return M_advectionToolbox->timeStepBase(); }
    std::shared_ptr<TSBase> timeStepBase() const { return M_advectionToolbox->timeStepBase(); }
    template<typename ExprT>
    void setAdvectionVelocityExpr( vf::Expr<ExprT> const& v_expr )
    { 
        ModelExpression velocityExpr;
        if constexpr ( nDim == 2 )
            velocityExpr.setExprVectorial2( v_expr );
        else if constexpr ( nDim == 3 )
            velocityExpr.setExprVectorial3( v_expr );

        for( std::string const& matName : M_advectionToolbox->materialsProperties()->physicToMaterials( M_advectionToolbox->physicDefault() ) )
        {
            M_advectionToolbox->materialsProperties()->materialProperties( matName ).add( M_advectionToolbox->convectionCoefficientName(), velocityExpr );
        }
    }
    //void updateAdvectionVelocity( element_advection_velocity_ptrtype const& velocity ) { [>return M_advectionToolbox->updateAdvectionVelocity( velocity );<] }
    //void updateAdvectionVelocity( element_advection_velocity_type const& velocity ) { [>return M_advectionToolbox->updateAdvectionVelocity( velocity );<] }
    //--------------------------------------------------------------------//
    // Spaces
    space_vectorial_ptrtype const& functionSpaceAdvectionVelocity() const { return M_spaceAdvectionVelocity; }
    void setFunctionSpaceAdvectionVelocity( space_vectorial_ptrtype const& space ) { M_spaceAdvectionVelocity = space; }
    space_tensor2symm_ptrtype const& functionSpaceTensor2Symm() const { return M_spaceTensor2Symm; }

    //std::string fileNameMeshPath() const override { return prefixvm(this->prefix(),"LevelsetMesh.path"); }

    //--------------------------------------------------------------------//
    // Levelset
    element_levelset_ptrtype phi() const { return this->phiPtr(); }
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
    // Assembly and solve
    int nBlockMatrixGraph() const { return 1; }

    bool hasSourceTerm() const { return false; }
    void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const {}
    void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;

    void solve();

    // Vector solution: forward to CFPDE
    ModelAlgebraic::block_vector_ptrtype algebraicBlockVectorSolution( std::string const& name = "" ) const { return M_advectionToolbox->algebraicBlockVectorSolution( name ); }

    // Assembly: forward CFPDE assembly
    void updateLinearPDE( DataUpdateLinear & data ) const override { M_advectionToolbox->updateLinearPDE( data ); }
    template <typename ModelContextType>
    void updateLinearPDE( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const { M_advectionToolbox->updateLinearPDE( data, mctx ); }
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override { M_advectionToolbox->updateLinearPDEDofElimination( data ); }
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const { M_advectionToolbox->updateLinearPDEDofElimination( data, mctx ); }
    template <typename ModelContextType,typename RangeType>
    void updateLinearPDEStabilizationGLS( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const { M_advectionToolbox->updateLinearPDEStabilizationGLS( data, mctx, matName, range ); }

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override { M_advectionToolbox->updateNewtonInitialGuess( data ); }
    template <typename ModelContextType>
    void updateNewtonInitialGuess( ModelAlgebraic::DataNewtonInitialGuess & data, ModelContextType const& mctx ) const { M_advectionToolbox->updateNewtonInitialGuess( data, mctx ); }
    void updateJacobian( DataUpdateJacobian & data ) const override { M_advectionToolbox->updateJacobian( data ); }
    template <typename ModelContextType>
    void updateJacobian( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx ) const { M_advectionToolbox->updateJacobian( data, mctx ); }
    template <typename ModelContextType,typename RangeType>
    void updateJacobianStabilizationGLS( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const { M_advectionToolbox->updateJacobianStabilizationGLS( data, mctx, matName, range ); }
    void updateJacobianDofElimination( ModelAlgebraic::DataUpdateJacobian & data ) const override { M_advectionToolbox->updateJacobianDofElimination( data ); }
    void updateResidual( DataUpdateResidual & data ) const override { M_advectionToolbox->updateResidual( data ); }
    template <typename ModelContextType>
    void updateResidual( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx ) const { M_advectionToolbox->updateResidual( data, mctx ); }
    template <typename ModelContextType,typename RangeType>
    void updateResidualStabilizationGLS( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const { M_advectionToolbox->updateResidualStabilizationGLS( data, mctx, matName, range ); }
    void updateResidualDofElimination( ModelAlgebraic::DataUpdateResidual & data ) const override { M_advectionToolbox->updateResidualDofElimination( data ); }

    //--------------------------------------------------------------------//
    // Time stepping
    void startTimeStep();
    void updateTimeStep();

    //--------------------------------------------------------------------//
    // Extension velocity
    template<typename ExprT>
    element_vectorial_type extensionVelocity( vf::Expr<ExprT> const& u ) const;

    //--------------------------------------------------------------------//
    // Model fields
    auto modelFields( std::string const& prefix = "" ) const
    {
        return Feel::FeelModels::modelFields( 
                super_type::modelFields( prefix ),
                M_advectionToolbox->modelFields( prefix )
                );
    }
    template <typename LevelsetFieldType>
    auto modelFields( LevelsetFieldType const& phi, std::string const& prefix = "" ) const
    {
        return Feel::FeelModels::modelFields( 
                super_type::modelFields( phi, prefix ),
                M_advectionToolbox->modelFields( phi, prefix )
                );
    }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
    {
        return Feel::FeelModels::modelFields( 
                super_type::modelFields( sol, rowStartInVector, prefix ),
                M_advectionToolbox->modelFields( sol, rowStartInVector, prefix )
                );
    }
    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
    {
        return Feel::FeelModels::selectorModelFields( 
                super_type::trialSelectorModelFields( startBlockSpaceIndex ),
                M_advectionToolbox->trialSelectorModelFields( startBlockSpaceIndex )
                );
    }

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
            fields[prefixvm( prefix, "backwardcharacteristics" )] = M_backwardCharacteristicsAdvection->fieldUnknownPtr();
        }
        return fields;
    }
    // Symbols expressions
    using super_type::symbolsExpr;
    // Measures quantities
    auto allMeasuresQuantities( std::string const& prefix = "" ) const
    {
#if 0
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
#endif
        return Feel::FeelModels::modelMeasuresQuantities( super_type::allMeasuresQuantities( prefix ),
                                                          modelMeasuresQuantity(prefix, "velocity-com", std::bind( &self_type::velocityCOM, this ) )
                                                          );
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
            //_expr=idv(M_advectionToolbox->fieldAdvectionVelocity()) * (1.-idv(this->H()))
            _expr=vec( cst(0), cst(0) ) * (1.-idv(this->H()))
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

    void initFunctionSpaces();
    void initInterfaceQuantities();
    void initTools();
    void initExporters();

    //--------------------------------------------------------------------//
    element_levelset_ptrtype const& phiPreviousTimeStepPtr() const { return this->timeStepBDF()->unknowns()[1]; }

    //--------------------------------------------------------------------//
    // Save
    void saveCurrent() const;


protected:
    //--------------------------------------------------------------------//
    // Boundary conditions
    using boundary_conditions_type = LevelSetBoundaryConditions;
    std::shared_ptr<boundary_conditions_type> M_boundaryConditions;

private:
    //--------------------------------------------------------------------//
    // Advection toolbox
    cfpde_toolbox_ptrtype M_advectionToolbox;
    bool M_doExportAdvection;
    //--------------------------------------------------------------------//
    // Spaces
    space_vectorial_ptrtype M_spaceAdvectionVelocity;
    space_tensor2symm_ptrtype M_spaceTensor2Symm;
    //--------------------------------------------------------------------//
    // Materials properties
    materialsproperties_ptrtype M_materialsProperties;

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
    template< typename ConvexType, typename BasisType, typename PeriodicityType, typename BasisPnType > \
        /**/
#endif
#ifndef LEVELSET_CLASS_TEMPLATE_TYPE
#define LEVELSET_CLASS_TEMPLATE_TYPE \
    LevelSet<ConvexType, BasisType, PeriodicityType, BasisPnType> \
        /**/
#endif
//----------------------------------------------------------------------------//
// Extension velocity
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
template<typename ExprT>
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_vectorial_type
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
