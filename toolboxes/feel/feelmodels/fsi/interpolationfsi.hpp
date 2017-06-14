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
 \file interpolationfsi.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-05
 */

#ifndef FEELPP_MODELS_INTERPOLATION_FSI_H
#define FEELPP_MODELS_INTERPOLATION_FSI_H 1


#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

#include <feel/feelmodels/modelcore/log.hpp>


namespace Feel
{
namespace FeelModels
{


template< class FluidType, class SolidType >
class InterpolationFSI
{
public :
    typedef InterpolationFSI<FluidType,SolidType> self_type;
    typedef FluidType fluid_type;
    typedef SolidType solid_type;

    typedef boost::shared_ptr<fluid_type> fluid_ptrtype;
    typedef boost::shared_ptr<solid_type> solid_ptrtype;

    typedef typename fluid_type::mesh_type mesh_fluid_type;
    typedef typename solid_type::mesh_type mesh_solid_type;
    typedef typename solid_type::mesh_1dreduced_type mesh_solid_1dreduced_type;

    typedef faces_reference_wrapper_t<mesh_fluid_type> range_fluid_face_type;
    typedef faces_reference_wrapper_t<mesh_solid_type> range_solid_face_type;
    typedef elements_reference_wrapper_t<mesh_solid_1dreduced_type> range_solid_elt_1dreduced_type;


    //-----------------------------------------------------------------------------------//
    // space and element displacement with interaction 2d/2d or 3d/3d
    typedef typename fluid_type::mesh_ale_type::ale_map_functionspace_type space_fluid_disp_type;
    typedef typename fluid_type::mesh_ale_type::ale_map_element_type element_fluid_disp_type;

    typedef typename solid_type::space_displacement_type space_struct_disp_type;
    typedef typename solid_type::element_displacement_type element_struct_disp_type;

    //operator interpolation for this displacement
    typedef OperatorInterpolation<space_struct_disp_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation2dTo2dnonconf_disp_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dnonconf_disp_type> op_interpolation2dTo2dnonconf_disp_ptrtype;
    typedef OperatorInterpolation<space_struct_disp_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation2dTo2dconf_disp_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dconf_disp_type> op_interpolation2dTo2dconf_disp_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element velocity with interaction 2d/2d or 3d/3d
    typedef typename fluid_type::element_meshvelocityonboundary_type::functionspace_type space_fluid_velocity_type;
    typedef typename fluid_type::element_meshvelocityonboundary_type element_fluid_velocity_type;

    typedef typename solid_type::space_displacement_type space_struct_velocity_type;
    typedef typename solid_type::element_vectorial_type element_struct_velocity_type;

    //operator interpolation for this velocity
    typedef OperatorInterpolation<space_struct_velocity_type, space_fluid_velocity_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation2dTo2dnonconf_velocity_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dnonconf_velocity_type> op_interpolation2dTo2dnonconf_velocity_ptrtype;
    typedef OperatorInterpolation<space_struct_velocity_type, space_fluid_velocity_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation2dTo2dconf_velocity_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dconf_velocity_type> op_interpolation2dTo2dconf_velocity_ptrtype;
    
    typedef OperatorInterpolation<space_struct_velocity_type, typename fluid_type::space_fluid_velocity_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation2dTo2dnonconf_velocityBis_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dnonconf_velocityBis_type> op_interpolation2dTo2dnonconf_velocityBis_ptrtype;
    typedef OperatorInterpolation<space_struct_velocity_type, typename fluid_type::space_fluid_velocity_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation2dTo2dconf_velocityBis_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dconf_velocityBis_type> op_interpolation2dTo2dconf_velocityBis_ptrtype;
    
    //-----------------------------------------------------------------------------------//
    // space and element displacement with interaction 2d/1d (disp is scalar)
    typedef typename solid_type::space_vect_1dreduced_type space_struct_vect_disp_1dreduced_type;
    typedef typename solid_type::element_vect_1dreduced_type element_struct_vect_disp_1dreduced_type;

    //operator interpolation for this displacement
    typedef OperatorInterpolation<space_struct_vect_disp_1dreduced_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationNonConforme> op_interpolation1dToNdnonconf_disp_type;
    typedef boost::shared_ptr<op_interpolation1dToNdnonconf_disp_type> op_interpolation1dToNdnonconf_disp_ptrtype;
    typedef OperatorInterpolation<space_struct_vect_disp_1dreduced_type, space_fluid_disp_type,
                                  range_fluid_face_type,InterpolationConforme> op_interpolation1dToNdconf_disp_type;
    typedef boost::shared_ptr<op_interpolation1dToNdconf_disp_type> op_interpolation1dToNdconf_disp_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element velocity with interaction 2d/1d (velocity is scalar)
    typedef typename solid_type::space_vect_1dreduced_type space_struct_vect_velocity_1dreduced_type;
    typedef typename solid_type::element_vect_1dreduced_type element_struct_vect_velocity_1dreduced_type;

    //operator interpolation for this velocity
    typedef OperatorInterpolation<space_struct_vect_velocity_1dreduced_type, space_fluid_velocity_type,
                                  range_fluid_face_type, InterpolationNonConforme> op_interpolation1dToNdnonconf_velocity_type;
    typedef boost::shared_ptr<op_interpolation1dToNdnonconf_velocity_type> op_interpolation1dToNdnonconf_velocity_ptrtype;
    typedef OperatorInterpolation<space_struct_vect_velocity_1dreduced_type, space_fluid_velocity_type,
                                  range_fluid_face_type, InterpolationConforme> op_interpolation1dToNdconf_velocity_type;
    typedef boost::shared_ptr<op_interpolation1dToNdconf_velocity_type> op_interpolation1dToNdconf_velocity_ptrtype;

    //-----------------------------------------------------------------------------------//
    // space and element stress with interaction 2d/2d or 3d/3d
    typedef typename fluid_type::space_stress_type space_fluid_stress_type;
    typedef typename fluid_type::element_stress_type element_fluid_stress_type;

    typedef typename solid_type::space_normal_stress_type space_struct_stress_type;
    typedef typename solid_type::element_normal_stress_type element_struct_stress_type;

    //operator interpolation for this stress
    typedef OperatorInterpolation<space_fluid_stress_type, space_struct_stress_type,
                                  range_solid_face_type,InterpolationNonConforme> op_interpolation2dTo2dnonconf_stress_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dnonconf_stress_type> op_interpolation2dTo2dnonconf_stress_ptrtype;
    typedef OperatorInterpolation<space_fluid_stress_type, space_struct_stress_type,
                                  range_solid_face_type,InterpolationConforme> op_interpolation2dTo2dconf_stress_type;
    typedef boost::shared_ptr<op_interpolation2dTo2dconf_stress_type> op_interpolation2dTo2dconf_stress_ptrtype;

    // NEW ( robin-neumann)
    typedef OperatorInterpolation<space_struct_stress_type, space_fluid_stress_type,
                                  range_fluid_face_type,InterpolationConforme> op_s2f_interpolation2dTo2dconf_stress_type;
    typedef boost::shared_ptr<op_s2f_interpolation2dTo2dconf_stress_type> op_s2f_interpolation2dTo2dconf_stress_ptrtype;

    //typedef typename fluid_type::space_fluid_velocity_type space_fluid_velocity_type;
    typedef OperatorInterpolation<typename fluid_type::space_fluid_velocity_type/*space_fluid_velocity_type*/,space_struct_velocity_type,
                                  range_solid_face_type,InterpolationConforme> op_f2s_interpolation2dTo2dconf_velocity_type;
    typedef boost::shared_ptr<op_f2s_interpolation2dTo2dconf_velocity_type> op_f2s_interpolation2dTo2dconf_velocity_ptrtype;


    //-----------------------------------------------------------------------------------//
    // space and element stress with interaction 2d/1d
    typedef typename solid_type::space_stress_vect_1dreduced_type space_struct_stress_vect_1dreduced_type;
    typedef typename solid_type::element_stress_vect_1dreduced_type element_struct_stress_vect_1dreduced_type;

    //operator interpolation for this stress
    typedef OperatorInterpolation<space_fluid_stress_type, space_struct_stress_vect_1dreduced_type,
                                  range_solid_elt_1dreduced_type,InterpolationNonConforme> op_interpolation1dToNdnonconf_stress_type;
    typedef boost::shared_ptr<op_interpolation1dToNdnonconf_stress_type> op_interpolation1dToNdnonconf_stress_ptrtype;
    typedef OperatorInterpolation<space_fluid_stress_type, space_struct_stress_vect_1dreduced_type,
                                  range_solid_elt_1dreduced_type,InterpolationConforme> op_interpolation1dToNdconf_stress_type;
    typedef boost::shared_ptr<op_interpolation1dToNdconf_stress_type> op_interpolation1dToNdconf_stress_ptrtype;

    //-----------------------------------------------------------------------------------//

    InterpolationFSI(fluid_ptrtype fluid, solid_ptrtype solid, bool buildOperators=true);
    InterpolationFSI( self_type const & M ) = default;

    void buildOperatorDispStructToFluid();
    void buildOperatorNormalStressFluidToStruct();

    //-----------------------------------------------------------------------------------//

    fluid_ptrtype const& fluid() const { return M_fluid; }
    fluid_ptrtype fluid() { return M_fluid; }
    solid_ptrtype const& solid() const { return M_solid; }
    solid_ptrtype solid(){ return M_solid; }

    bool verbose() const { return M_verbose; }
    bool verboseAllProc() const { return M_verboseAllProc; }
    WorldComm worldComm() const { return M_worldComm; }

    std::string fsiCouplingType() const { return M_fsiCouplingType; }
    std::string fsiCouplingBoundaryCondition() const { return M_fsiCouplingBoundaryCondition; }

    vector_ptrtype const& robinNeumannInterfaceOperator() const { return M_robinNeumannInterfaceOperator; }
    void setRobinNeumannInterfaceOperator( vector_ptrtype const& vec ) { M_robinNeumannInterfaceOperator = vec; }

    //-----------------------------------------------------------------------------------//

    void initDispInterpolation();
    void initDisp1dToNdInterpolation();

    void initStressInterpolation();
    void initStress1dToNdInterpolation();

    void initVelocityInterpolation();
    void initVelocity1dToNdInterpolation();

    void initVelocityInterpolationF2S();
    void initStressInterpolationS2F();

    //-----------------------------------------------------------------------------------//

    void transfertDisplacement();
    void transfertStress();
    void transfertVelocity( bool useExtrap=false);
    void transfertRobinNeumannGeneralizedS2F( int iterationFSI, double manualScaling = 1 );

    void transfertStressS2F();
    void transfertVelocityF2S( int iterationFSI, bool _useExtrapolation );

    void transfertRobinNeumannInterfaceOperatorS2F();

    //void updateFieldVelocitySolidPreviousPrevious( typename solid_type::element_displacement_type const& vel ) {  *M_fieldVelocitySolidPreviousPrevious = vel; }
    //void updateFieldVelocitySolid1dReducedPreviousPrevious( typename solid_type::element_1dreduced_type const& vel ) {  *M_fieldVelocitySolid1dReducedPreviousPrevious = vel; }
private :

    fluid_ptrtype M_fluid;
    solid_ptrtype M_solid;

    WorldComm M_worldComm;

    bool M_interfaceFSIisConforme;
    std::string M_fsiCouplingType;
    std::string M_fsiCouplingBoundaryCondition;

    op_interpolation2dTo2dnonconf_disp_ptrtype M_opDisp2dTo2dnonconf;
    op_interpolation2dTo2dconf_disp_ptrtype M_opDisp2dTo2dconf;
    bool M_opDispIsInit;

    op_interpolation2dTo2dconf_stress_ptrtype M_opStress2dTo2dconf;
    op_interpolation2dTo2dnonconf_stress_ptrtype M_opStress2dTo2dnonconf;
    bool M_opStressIsInit;

    op_interpolation2dTo2dconf_velocity_ptrtype M_opVelocity2dTo2dconf;
    op_interpolation2dTo2dnonconf_velocity_ptrtype M_opVelocity2dTo2dnonconf;
    bool M_opVelocityIsInit;
    
    op_interpolation2dTo2dconf_velocityBis_ptrtype M_opVelocityBis2dTo2dconf;
    op_interpolation2dTo2dnonconf_velocityBis_ptrtype M_opVelocityBis2dTo2dnonconf;
    
    // 1d reduced model
    op_interpolation1dToNdnonconf_disp_ptrtype M_opDisp1dToNdnonconf;
    op_interpolation1dToNdconf_disp_ptrtype M_opDisp1dToNdconf;
    bool M_opDisp1dToNdIsInit;

    op_interpolation1dToNdnonconf_stress_ptrtype M_opStress1dToNdnonconf;
    op_interpolation1dToNdconf_stress_ptrtype M_opStress1dToNdconf;
    bool M_opStress1dToNdIsInit;

    //op_interpolation1dTo2d_velocity_scalar_ptrtype M_opVelocityScal1dTo2d;
    //bool M_opVelocityScalIsInit;

    op_interpolation1dToNdnonconf_velocity_ptrtype M_opVelocity1dToNdnonconf;
    op_interpolation1dToNdconf_velocity_ptrtype M_opVelocity1dToNdconf;

    op_s2f_interpolation2dTo2dconf_stress_ptrtype M_opStress2dTo2dconfS2F;

    op_f2s_interpolation2dTo2dconf_velocity_ptrtype M_opVelocity2dTo2dconfF2S;

    bool M_verbose,M_verboseAllProc;

    vector_ptrtype M_robinNeumannInterfaceOperator;
    //typename solid_type::element_displacement_ptrtype M_fieldVelocitySolidPreviousPrevious;
    //typename solid_type::element_1dreduced_ptrtype M_fieldVelocitySolid1dReducedPreviousPrevious;

};


} // namespace FeelModels
} // namespace Feel



#endif // FEELPP_MODELS_INTERPOLATION_FSI_H
