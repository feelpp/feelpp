/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-07-05

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file fsi.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-05
 */

#ifndef FEELPP_MODELS_FSI_H
#define FEELPP_MODELS_FSI_H 1

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/solid/solidmechanics.hpp>

#include <feel/feelmodels/fsi/aitkenrelaxationfsi.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
//#include <feel/feelts/tsbase.hpp>

namespace Feel
{
namespace FeelModels
{

template< class FluidType, class SolidType >
class FSI : public ModelNumerical,
            public ModelPhysics<FluidType::convex_type::nRealDim>,
            public std::enable_shared_from_this< FSI<FluidType,SolidType> >
{
    using super_physics_type = ModelPhysics<FluidType::convex_type::nRealDim>;
public :
    typedef ModelNumerical super_type;
    typedef FSI<FluidType,SolidType> self_type;

    typedef FluidType fluid_type;
    typedef SolidType solid_type;
    typedef std::shared_ptr<fluid_type> fluid_ptrtype;
    typedef std::shared_ptr<solid_type> solid_ptrtype;

    typedef typename fluid_type::mesh_type mesh_fluid_type;
    typedef typename fluid_type::trace_mesh_type trace_mesh_fluid_type;
    typedef typename solid_type::mesh_type mesh_solid_type;
    typedef typename solid_type::trace_mesh_type trace_mesh_solid_type;
    typedef typename solid_type::solid_1dreduced_type::mesh_type mesh_solid_1dreduced_type;

    // materials properties
    typedef MaterialsProperties<fluid_type::mesh_type::nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;


    // mesh velocity on FSI boundary
    typedef FunctionSpace<trace_mesh_fluid_type, bases<typename fluid_type::basis_fluid_u_type> > space_fluid_meshvelocityonboundary_type;
    typedef std::shared_ptr<space_fluid_meshvelocityonboundary_type> space_fluid_meshvelocityonboundary_ptrtype;
    typedef typename space_fluid_meshvelocityonboundary_type::element_type element_fluid_meshvelocityonboundary_type;
    typedef std::shared_ptr<element_fluid_meshvelocityonboundary_type> element_fluid_meshvelocityonboundary_ptrtype;
    //___________________________________________________________________________________//
    // normal stress from fluid into solid model
    static const uint16_type nOrder_solid_normalStressFromFluid = fluid_type::space_normalstress_type::basis_type::nOrder;
    typedef bases<Lagrange< nOrder_solid_normalStressFromFluid, Vectorial,Discontinuous,PointSetFekete>> basis_solid_normalstressfromfluid_type;
    //typedef FunctionSpace<mesh_solid_type,basis_solid_normalstressfromfluid_type> space_solid_normalstressfromfluid_type;
    typedef FunctionSpace<trace_mesh_solid_type,basis_solid_normalstressfromfluid_type> space_solid_normalstressfromfluid_type;
    //typedef typename solid_type::space_normal_stress_type space_solid_normalstressfromfluid_type;
    typedef std::shared_ptr<space_solid_normalstressfromfluid_type> space_solid_normalstressfromfluid_ptrtype;
    typedef typename space_solid_normalstressfromfluid_type::element_type element_solid_normalstressfromfluid_type;
    typedef std::shared_ptr<element_solid_normalstressfromfluid_type> element_solid_normalstressfromfluid_ptrtype;
    //___________________________________________________________________________________//
    // normal stress from fluid into solid 1dreduced model
    static const uint16_type nOrder_solid1dReduced_normalStressFromFluid = fluid_type::space_normalstress_type::basis_type::nOrder;
    typedef bases<Lagrange< nOrder_solid1dReduced_normalStressFromFluid, Vectorial,Discontinuous,PointSetFekete>> basis_solid1dreduced_normalstressfromfluid_vect_type;
    typedef FunctionSpace<mesh_solid_1dreduced_type, basis_solid1dreduced_normalstressfromfluid_vect_type> space_solid1dreduced_normalstressfromfluid_vect_type;
    typedef std::shared_ptr<space_solid1dreduced_normalstressfromfluid_vect_type> space_solid1dreduced_normalstressfromfluid_vect_ptrtype;
    typedef typename space_solid1dreduced_normalstressfromfluid_vect_type::element_type element_solid1dreduced_normalstressfromfluid_vect_type;
    typedef std::shared_ptr<element_solid1dreduced_normalstressfromfluid_vect_type> element_solid1dreduced_normalstressfromfluid_vect_ptrtype;

    typedef typename space_solid1dreduced_normalstressfromfluid_vect_type::component_functionspace_type space_solid1dreduced_normalstressfromfluid_scal_type;
    typedef typename space_solid1dreduced_normalstressfromfluid_vect_type::component_functionspace_ptrtype space_solid1dreduced_normalstressfromfluid_scal_ptrtype;
    typedef typename space_solid1dreduced_normalstressfromfluid_scal_type::element_type element_solid1dreduced_normalstressfromfluid_scal_type;
    typedef std::shared_ptr<element_solid1dreduced_normalstressfromfluid_scal_type> element_solid1dreduced_normalstressfromfluid_scal_ptrtype;
    //___________________________________________________________________________________//

    typedef faces_reference_wrapper_t<mesh_fluid_type> range_fluid_face_type;
    typedef elements_reference_wrapper_t<trace_mesh_fluid_type> range_fluid_trace_elt_type;
    typedef faces_reference_wrapper_t<mesh_solid_type> range_solid_face_type;
    typedef elements_reference_wrapper_t<trace_mesh_solid_type> range_solid_trace_elt_type;
    typedef elements_reference_wrapper_t<mesh_solid_1dreduced_type> range_solid_elt_1dreduced_type;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //typedef InterpolationFSI<fluid_type,solid_type> interpolationFSI_type;
    //typedef std::shared_ptr<interpolationFSI_type> interpolationFSI_ptrtype;
        // space and element displacement with interaction 2d/2d or 3d/3d
    typedef typename fluid_type::mesh_ale_type::ale_map_functionspace_type space_fluid_disp_type;
    typedef typename fluid_type::mesh_ale_type::ale_map_element_type element_fluid_disp_type;
    typedef typename fluid_type::mesh_ale_type::ale_map_element_ptrtype element_fluid_disp_ptrtype;

    typedef typename solid_type::space_displacement_type space_struct_disp_type;
    typedef typename solid_type::element_displacement_type element_struct_disp_type;

    //operator interpolation for this displacement
    typedef OperatorInterpolation<space_struct_disp_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation2dTo2dnonconf_disp_type;
    typedef std::shared_ptr<op_interpolation2dTo2dnonconf_disp_type> op_interpolation2dTo2dnonconf_disp_ptrtype;
    typedef OperatorInterpolation<space_struct_disp_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation2dTo2dconf_disp_type;
    typedef std::shared_ptr<op_interpolation2dTo2dconf_disp_type> op_interpolation2dTo2dconf_disp_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element velocity with interaction 2d/2d or 3d/3d
    //typedef typename fluid_type::element_meshvelocityonboundary_type::functionspace_type space_fluid_velocity_type;
    //typedef typename fluid_type::element_meshvelocityonboundary_type element_fluid_velocity_type;

    typedef typename solid_type::space_displacement_type space_struct_velocity_type;
    typedef typename solid_type::element_vectorial_type element_struct_velocity_type;

    //operator interpolation for this velocity
    typedef OperatorInterpolation<space_struct_velocity_type, space_fluid_meshvelocityonboundary_type/*space_fluid_velocity_type*/,
                                  range_fluid_trace_elt_type/*range_fluid_face_type*/,InterpolationNonConforme> op_interpolation2dTo2dnonconf_velocity_type;
    typedef std::shared_ptr<op_interpolation2dTo2dnonconf_velocity_type> op_interpolation2dTo2dnonconf_velocity_ptrtype;
    typedef OperatorInterpolation<space_struct_velocity_type, space_fluid_meshvelocityonboundary_type/*space_fluid_velocity_type*/,
                                  range_fluid_trace_elt_type/*range_fluid_face_type*/,InterpolationConforme> op_interpolation2dTo2dconf_velocity_type;
    typedef std::shared_ptr<op_interpolation2dTo2dconf_velocity_type> op_interpolation2dTo2dconf_velocity_ptrtype;
    
    typedef OperatorInterpolation<space_struct_velocity_type, typename fluid_type::space_velocity_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation2dTo2dnonconf_velocityBis_type;
    typedef std::shared_ptr<op_interpolation2dTo2dnonconf_velocityBis_type> op_interpolation2dTo2dnonconf_velocityBis_ptrtype;
    typedef OperatorInterpolation<space_struct_velocity_type, typename fluid_type::space_velocity_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation2dTo2dconf_velocityBis_type;
    typedef std::shared_ptr<op_interpolation2dTo2dconf_velocityBis_type> op_interpolation2dTo2dconf_velocityBis_ptrtype;
    
    //-----------------------------------------------------------------------------------//
    // space and element displacement with interaction 2d/1d (disp is scalar)
    typedef typename solid_type::solid_1dreduced_type::space_displacement_type/*space_vect_1dreduced_type*/ space_struct_vect_disp_1dreduced_type;
    typedef typename solid_type::solid_1dreduced_type::element_displacement_type/*element_vect_1dreduced_type*/ element_struct_vect_disp_1dreduced_type;

    //operator interpolation for this displacement
    typedef OperatorInterpolation<space_struct_vect_disp_1dreduced_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation1dToNdnonconf_disp_type;
    typedef std::shared_ptr<op_interpolation1dToNdnonconf_disp_type> op_interpolation1dToNdnonconf_disp_ptrtype;
    typedef OperatorInterpolation<space_struct_vect_disp_1dreduced_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation1dToNdconf_disp_type;
    typedef std::shared_ptr<op_interpolation1dToNdconf_disp_type> op_interpolation1dToNdconf_disp_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element velocity with interaction 2d/1d (velocity is scalar)
    typedef typename solid_type::solid_1dreduced_type::space_displacement_type/*space_vect_1dreduced_type*/ space_struct_vect_velocity_1dreduced_type;
    typedef typename solid_type::solid_1dreduced_type::element_displacement_component_type/*element_vect_1dreduced_type*/ element_struct_vect_velocity_1dreduced_type;

    //operator interpolation for this velocity
    typedef OperatorInterpolation<space_struct_vect_velocity_1dreduced_type, space_fluid_meshvelocityonboundary_type/*space_fluid_velocity_type*/,
                                  range_fluid_trace_elt_type/*range_fluid_face_type*/, InterpolationNonConforme> op_interpolation1dToNdnonconf_velocity_type;
    typedef std::shared_ptr<op_interpolation1dToNdnonconf_velocity_type> op_interpolation1dToNdnonconf_velocity_ptrtype;
    typedef OperatorInterpolation<space_struct_vect_velocity_1dreduced_type, space_fluid_meshvelocityonboundary_type/*space_fluid_velocity_type*/,
                                  range_fluid_trace_elt_type/*range_fluid_face_type*/, InterpolationConforme> op_interpolation1dToNdconf_velocity_type;
    typedef std::shared_ptr<op_interpolation1dToNdconf_velocity_type> op_interpolation1dToNdconf_velocity_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element stress with interaction 2d/2d or 3d/3d
    typedef typename fluid_type::space_normalstress_type space_fluid_stress_type;
    typedef typename fluid_type::element_normalstress_type element_fluid_stress_type;

    typedef typename solid_type::space_normal_stress_type space_struct_stress_type;
    //typedef typename solid_type::element_normal_stress_type element_struct_stress_type;

    //operator interpolation for this stress
    typedef OperatorInterpolation<space_fluid_stress_type, space_solid_normalstressfromfluid_type,
                                  range_solid_trace_elt_type/*range_solid_face_type*/,InterpolationNonConforme> op_interpolation2dTo2dnonconf_stress_type;
    typedef std::shared_ptr<op_interpolation2dTo2dnonconf_stress_type> op_interpolation2dTo2dnonconf_stress_ptrtype;
    typedef OperatorInterpolation<space_fluid_stress_type, space_solid_normalstressfromfluid_type,
                                  range_solid_trace_elt_type/*range_solid_face_type*/,InterpolationConforme> op_interpolation2dTo2dconf_stress_type;
    typedef std::shared_ptr<op_interpolation2dTo2dconf_stress_type> op_interpolation2dTo2dconf_stress_ptrtype;

    // NEW ( robin-neumann)
    typedef OperatorInterpolation<space_struct_stress_type, space_fluid_stress_type,
                                  range_fluid_face_type,InterpolationConforme> op_s2f_interpolation2dTo2dconf_stress_type;
    typedef std::shared_ptr<op_s2f_interpolation2dTo2dconf_stress_type> op_s2f_interpolation2dTo2dconf_stress_ptrtype;

    //typedef typename fluid_type::space_fluid_velocity_type space_fluid_velocity_type;
    typedef OperatorInterpolation<typename fluid_type::space_velocity_type/*space_fluid_velocity_type*/,space_struct_velocity_type,
                                  range_solid_face_type,InterpolationConforme> op_f2s_interpolation2dTo2dconf_velocity_type;
    typedef std::shared_ptr<op_f2s_interpolation2dTo2dconf_velocity_type> op_f2s_interpolation2dTo2dconf_velocity_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element stress with interaction 2d/1d
    //typedef typename solid_type::space_stress_vect_1dreduced_type space_struct_stress_vect_1dreduced_type;
    //typedef typename solid_type::element_stress_vect_1dreduced_type element_struct_stress_vect_1dreduced_type;
    //operator interpolation for this stress
    typedef OperatorInterpolation<space_fluid_stress_type, space_solid1dreduced_normalstressfromfluid_vect_type,
                                  range_solid_elt_1dreduced_type,InterpolationNonConforme> op_interpolation1dToNdnonconf_stress_type;
    typedef std::shared_ptr<op_interpolation1dToNdnonconf_stress_type> op_interpolation1dToNdnonconf_stress_ptrtype;
    typedef OperatorInterpolation<space_fluid_stress_type, space_solid1dreduced_normalstressfromfluid_vect_type,
                                  range_solid_elt_1dreduced_type,InterpolationConforme> op_interpolation1dToNdconf_stress_type;
    typedef std::shared_ptr<op_interpolation1dToNdconf_stress_type> op_interpolation1dToNdconf_stress_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//

    typedef AitkenRelaxationFSI<solid_type> aitkenrelaxationFSI_type;
    typedef std::shared_ptr<aitkenrelaxationFSI_type> aitkenrelaxationFSI_ptrtype;
    typedef FixPointConvergenceFSI<solid_type> fixpointconvergenceFSI_type;
    typedef std::shared_ptr<fixpointconvergenceFSI_type> fixpointconvergenceFSI_ptrtype;

    //---------------------------------------------------------------------------------------------------------//

    FSI( std::string const& prefix,
         std::string const& keyword = "fsi",
         worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
         ModelBaseRepository const& modelRep = ModelBaseRepository() );
    FSI( self_type const & M ) = default;

    //static std::string expandStringFromSpec( std::string const& expr );

    //---------------------------------------------------------------------------------------------------------//

    double meshSize() const { return M_meshSize; }

    fluid_ptrtype const& fluidModel() const { return M_fluidModel; }
    solid_ptrtype const& solidModel() const { return M_solidModel; }
    void setFluidModel( fluid_ptrtype const& fm ) { M_fluidModel=fm; }
    void setSolidModel( solid_ptrtype const& sm ) { M_solidModel=sm; }

    std::string fsiCouplingType() const { return M_fsiCouplingType; }
    std::string fsiCouplingBoundaryCondition() const { return M_fsiCouplingBoundaryCondition; }
    bool useFSISemiImplicitScheme() const { return ( this->fsiCouplingType() == "Semi-Implicit" ); }
    bool interfaceFSIisConforme() const { return M_interfaceFSIisConforme; }
    double fixPointTolerance() const { return M_fixPointTolerance; }
    double fixPointInitialTheta() const { return M_fixPointInitialTheta; }
    double fixPointMinTheta() const { return M_fixPointMinTheta; }
    int fixPointMaxIt() const { return M_fixPointMaxIt; }
    int fixPointMinItConvergence() const { return M_fixPointMinItConvergence; }

    //interpolationFSI_ptrtype interpolationTool() { return M_interpolationFSI; }
    //interpolationFSI_ptrtype const& interpolationTool() const { return M_interpolationFSI; }
    aitkenrelaxationFSI_ptrtype aitkenRelaxTool() { return M_aitkenFSI; }
    aitkenrelaxationFSI_ptrtype const& aitkenRelaxTool() const { return M_aitkenFSI; }

    //---------------------------------------------------------------------------------------------------------//

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    //---------------------------------------------------------------------------------------------------------//

    void createMesh();
    void init();
    void solve();
private :
    void updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models ) override;

    void initCouplingRobinNeumannGeneralized();

    void initInterpolation();
    void initDispInterpolation();
    void initDisp1dToNdInterpolation();
    void initStressInterpolation();
    void initStress1dToNdInterpolation();
    void initVelocityInterpolation();
    void initVelocity1dToNdInterpolation();
    void initVelocityInterpolationF2S();
    void initStressInterpolationS2F();

    void transfertDisplacement();
    void transfertDisplacementAndApplyMeshMoving();
    void transfertStress();
    void transfertVelocity( bool useExtrap=false);
    void transfertRobinNeumannGeneralizedS2F( int iterationFSI );
    void transfertRobinNeumannGeneralizedS2F_BdfNewmark( int iterationFSI );
    void transfertRobinNeumannGeneralizedS2F_BdfBdf( int iterationFSI );


    void transfertStressS2F();
    void transfertVelocityF2S( int iterationFSI, bool _useExtrapolation );

    void transfertGradVelocityF2S();


    double couplingRNG_coeffForm2() const { return M_couplingRNG_coeffForm2; }
    typename fluid_type::space_velocity_type::element_ptrtype/*element_meshvelocityonboundary_ptrtype*/ const& couplingRNG_evalForm1() const { return M_couplingRNG_evalForm1; }

    auto
    couplingRNG_operatorExpr( mpl::int_<2> /**/ ) const
        {
            return mat<2,2>(idv(M_coulingRNG_operatorDiagonalOnFluid)(0),cst(0.),cst(0.),idv(M_coulingRNG_operatorDiagonalOnFluid)(1));
        }
    auto
    couplingRNG_operatorExpr( mpl::int_<3> /**/ ) const
        {
            return mat<3,3>(idv(M_coulingRNG_operatorDiagonalOnFluid)(0),cst(0.),cst(0.),
                            cst(0.),idv(M_coulingRNG_operatorDiagonalOnFluid)(1),cst(0.),
                            cst(0.),cst(0.),idv(M_coulingRNG_operatorDiagonalOnFluid)(2) );
        }
public :
    //---------------------------------------------------------------------------------------------------------//

    void updateTime(double time);

    std::shared_ptr<TSBase> timeStepBase() const { return this->fluidTimeStepBase(); }
    std::shared_ptr<TSBase> fluidTimeStepBase() const { return this->fluidModel()->timeStepBase(); }
    std::shared_ptr<TSBase> solidTimeStepBase() const { return this->solidModel()->timeStepBase(); }
    void startTimeStep();
    void updateTimeStep();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time )
    {
        this->fluidModel()->exportResults(time);
        this->solidModel()->exportResults(time);
    }


    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //---------------------------------------------------------------------------------------------------------//
    void updateLinearPDE_Fluid( DataUpdateLinear & data ) const;
    void updateJacobian_Fluid( DataUpdateJacobian & data ) const;
    void updateResidual_Fluid( DataUpdateResidual & data ) const;
    void updateLinearPDE_Solid( DataUpdateLinear & data ) const;
    void updateJacobian_Solid( DataUpdateJacobian & data ) const;
    void updateResidual_Solid( DataUpdateResidual & data ) const;
    void updateLinearPDEDofElimination_Fluid( DataUpdateLinear & data ) const;
    void updateNewtonInitialGuess_Fluid( DataNewtonInitialGuess & data ) const;
    void updateJacobianDofElimination_Fluid( DataUpdateJacobian & data ) const;
    void updateResidualDofElimination_Fluid( DataUpdateResidual & data ) const;

    void updateLinearPDE_Solid1dReduced( DataUpdateLinear & data ) const;

private :
    void updateBackendOptimisation( int iterationFSI, double lastErrorRelative );
    void solveImpl1();
    void solveImpl2();
    void solveImpl3();

    //---------------------------------------------------------------------------------------------------------//
    typedef typename fluid_type::operatorpcdbase_type operatorpcdbase_fluid_type;
    void initInHousePreconditionerPCD_fluid( operatorpcdbase_fluid_type & opPCD ) const;
    void updateInHousePreconditionerPCD_fluid( operatorpcdbase_fluid_type & opPCD, DataUpdateBase & data ) const;
private :

    fluid_ptrtype M_fluidModel;
    solid_ptrtype M_solidModel;

    materialsproperties_ptrtype M_materialsProperties;

    double M_meshSize;
    fs::path M_mshfilepathFluidPart1,M_mshfilepathSolidPart1;
    fs::path M_mshfilepathFluidPartN,M_mshfilepathSolidPartN;
    std::set<std::string> M_markersNameFluid,M_markersNameSolid;
    std::string M_tagFileNameMeshGenerated;

    range_fluid_face_type M_rangeFSI_fluid;
    range_solid_face_type M_rangeFSI_solid;
    std::map<std::string,range_fluid_face_type> M_rangeMeshFacesByMaterial_fluid;

    std::string M_fsiCouplingType; // implicit,semi-implicit
    std::string M_fsiCouplingBoundaryCondition; // dirichlet-neumann, robin-robin, ...
    bool M_interfaceFSIisConforme;
    double M_fixPointTolerance, M_fixPointInitialTheta, M_fixPointMinTheta;
    int M_fixPointMaxIt, M_fixPointMinItConvergence;

    //interpolationFSI_ptrtype M_interpolationFSI;
    aitkenrelaxationFSI_ptrtype M_aitkenFSI;
    fixpointconvergenceFSI_ptrtype M_fixPointConvergenceFSI;

    int M_previousTimeOrder,M_currentTimeOrder;
    bool M_reusePrecOptFluid,M_reusePrecRebuildAtFirstFSIStepOptFluid,M_reuseJacOptFluid,M_reuseJacRebuildAtFirstNewtonStepOptFluid,M_reuseJacRebuildAtFirstFSIStepOptFluid;
    bool M_reusePrecOptSolid,M_reusePrecRebuildAtFirstFSIStepOptSolid,M_reuseJacOptSolid,M_reuseJacRebuildAtFirstNewtonStepOptSolid,M_reuseJacRebuildAtFirstFSIStepOptSolid;
    int M_reusePrecActivatedAfterNbFsiIterationFluid,M_reusePrecActivatedAfterNbFsiIterationSolid;
    double M_reusePrecActivatedToleranceFluid,M_reusePrecActivatedToleranceSolid;

    double M_couplingNitscheFamily_gamma, M_couplingNitscheFamily_gamma0, M_couplingNitscheFamily_alpha;

    double M_couplingRNG_coeffForm2;
    typename fluid_type::space_velocity_type::element_ptrtype/*element_meshvelocityonboundary_ptrtype*/ M_couplingRNG_evalForm1;
    typename fluid_type::space_velocity_type::element_ptrtype/*element_meshvelocityonboundary_ptrtype*/ M_coulingRNG_operatorDiagonalOnFluid;
    sparse_matrix_ptrtype M_coulingRNG_matrixTimeDerivative, M_coulingRNG_matrixStress;
    vector_ptrtype M_coulingRNG_vectorTimeDerivative,  M_coulingRNG_vectorStress;
    std::string M_coulingRNG_strategyTimeStepCompatibility;

    op_interpolation2dTo2dnonconf_disp_ptrtype M_opDisp2dTo2dnonconf;
    op_interpolation2dTo2dconf_disp_ptrtype M_opDisp2dTo2dconf;

    op_interpolation2dTo2dconf_stress_ptrtype M_opStress2dTo2dconf;
    op_interpolation2dTo2dnonconf_stress_ptrtype M_opStress2dTo2dnonconf;

    op_interpolation2dTo2dconf_velocity_ptrtype M_opVelocity2dTo2dconf;
    op_interpolation2dTo2dnonconf_velocity_ptrtype M_opVelocity2dTo2dnonconf;

    op_interpolation2dTo2dconf_velocityBis_ptrtype M_opVelocityBis2dTo2dconf;
    op_interpolation2dTo2dnonconf_velocityBis_ptrtype M_opVelocityBis2dTo2dnonconf;

    // 1d reduced model
    op_interpolation1dToNdnonconf_disp_ptrtype M_opDisp1dToNdnonconf;
    op_interpolation1dToNdconf_disp_ptrtype M_opDisp1dToNdconf;

    op_interpolation1dToNdnonconf_stress_ptrtype M_opStress1dToNdnonconf;
    op_interpolation1dToNdconf_stress_ptrtype M_opStress1dToNdconf;

    op_interpolation1dToNdnonconf_velocity_ptrtype M_opVelocity1dToNdnonconf;
    op_interpolation1dToNdconf_velocity_ptrtype M_opVelocity1dToNdconf;

    op_s2f_interpolation2dTo2dconf_stress_ptrtype M_opStress2dTo2dconfS2F;

    op_f2s_interpolation2dTo2dconf_velocity_ptrtype M_opVelocity2dTo2dconfF2S;


    // normal stress of fluid at fsi interface and computed on reference domain
    typename fluid_type::space_normalstress_ptrtype M_spaceNormalStress_fluid;
    typename fluid_type::element_normalstress_ptrtype M_fieldNormalStressRefMesh_fluid;

    // normal stress from fluid into solid model
    space_solid_normalstressfromfluid_ptrtype M_spaceNormalStressFromFluid_solid;
    element_solid_normalstressfromfluid_ptrtype M_fieldNormalStressFromFluid_solid;
    typename solid_type::element_vectorial_ptrtype M_fieldVelocityInterfaceFromFluid_solid;
    element_solid_normalstressfromfluid_ptrtype fieldNormalStressFromFluidPtr_solid() const { return M_fieldNormalStressFromFluid_solid; }
    typename solid_type::element_vectorial_ptrtype fieldVelocityInterfaceFromFluidPtr_solid() const { return M_fieldVelocityInterfaceFromFluid_solid; }
    element_solid_normalstressfromfluid_type const& fieldNormalStressFromFluid_solid() const { return *M_fieldNormalStressFromFluid_solid; }
    typename solid_type::element_vectorial_type const& fieldVelocityInterfaceFromFluid_solid() const { return *M_fieldVelocityInterfaceFromFluid_solid; }

    // normal stress from fluid into solid 1dreduced model
    space_solid1dreduced_normalstressfromfluid_vect_ptrtype M_spaceNormalStressFromFluid_solid1dReduced;
    element_solid1dreduced_normalstressfromfluid_scal_ptrtype M_fieldNormalStressFromFluidScalar_solid1dReduced;
    element_solid1dreduced_normalstressfromfluid_vect_ptrtype  M_fieldNormalStressFromFluidVectorial_solid1dReduced;

    // gradient of velocity from fluid into solid model (use normal spress function space and store the matric field by a vector of vectorial field)
    std::vector<typename fluid_type::element_normalstress_ptrtype> M_fieldsGradVelocity_fluid;
    std::vector<typename space_solid_normalstressfromfluid_type::/*component_functionspace_type::*/element_ptrtype> M_fieldsGradVelocity_solid;

    auto
    gradVelocityExpr_fluid2solid(  hana::int_<2> /**/ ) const
        {
            return mat<2,2>( idv(M_fieldsGradVelocity_solid[0])(0),idv(M_fieldsGradVelocity_solid[0])(1),
                             idv(M_fieldsGradVelocity_solid[1])(0),idv(M_fieldsGradVelocity_solid[1])(1) );
        }
    auto
    gradVelocityExpr_fluid2solid(  hana::int_<3> /**/ ) const
        {
            return mat<3,3>( idv(M_fieldsGradVelocity_solid[0])(0),idv(M_fieldsGradVelocity_solid[0])(1),idv(M_fieldsGradVelocity_solid[0])(2),
                             idv(M_fieldsGradVelocity_solid[1])(0),idv(M_fieldsGradVelocity_solid[1])(1),idv(M_fieldsGradVelocity_solid[1])(2),
                             idv(M_fieldsGradVelocity_solid[2])(0),idv(M_fieldsGradVelocity_solid[2])(1),idv(M_fieldsGradVelocity_solid[2])(2) );
        }



    element_fluid_meshvelocityonboundary_type & meshVelocity2() { return *M_meshVelocityInterface; }
    element_fluid_meshvelocityonboundary_ptrtype meshVelocity2Ptr() { return M_meshVelocityInterface; }
    element_fluid_meshvelocityonboundary_type const & meshVelocity2() const { return *M_meshVelocityInterface; }
    element_fluid_meshvelocityonboundary_ptrtype const & meshVelocity2Ptr() const { return M_meshVelocityInterface; }
    space_fluid_meshvelocityonboundary_ptrtype M_XhMeshVelocityInterface;
    element_fluid_meshvelocityonboundary_ptrtype M_meshVelocityInterface;
    //std::set<size_type> M_dofsVelocityInterfaceOnMovingBoundary;

    std::set<size_type> M_dofsMultiProcessVelocitySpaceOnFSI_fluid;

    element_fluid_disp_ptrtype M_meshDisplacementOnInterface_fluid;
};

} // namespace FeelModels
} // namespace Feel



#endif // FEELPP_MODELS_FSI_H
