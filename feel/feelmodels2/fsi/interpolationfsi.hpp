/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
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
 \author Vincent Chabannes <vincent.chabannes@imag.fr>
 \date 2011-07-05
 */

#ifndef FEELPP_MODELS_INTERPOLATION_FSI_H
#define FEELPP_MODELS_INTERPOLATION_FSI_H 1


#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/projector.hpp>

#include <feel/feelmodels2/modelcore/log.hpp>


namespace Feel
{
namespace FeelModels
{


template< class FluidMech, class SolidMech >
class InterpolationFSI
{
public :

    typedef FluidMech fluid_type;
    typedef SolidMech solid_type;

    typedef boost::shared_ptr<fluid_type> fluid_ptrtype;
    typedef boost::shared_ptr<solid_type> solid_ptrtype;

    typedef typename fluid_type::mesh_type mesh_fluid_type;
    typedef typename solid_type::mesh_type mesh_solid_type;
    typedef typename solid_type::mesh_1dreduced_type mesh_solid_1dreduced_type;

    typedef boost::tuple<boost::mpl::size_t<MESH_FACES>,
                         typename MeshTraits<mesh_fluid_type>::marker_face_const_iterator,
                         typename MeshTraits<mesh_fluid_type>::marker_face_const_iterator> range_fluid_face_type;

    typedef boost::tuple<boost::mpl::size_t<MESH_FACES>,
                         typename MeshTraits<mesh_solid_type>::marker_face_const_iterator,
                         typename MeshTraits<mesh_solid_type>::marker_face_const_iterator> range_solid_face_type;

    typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<mesh_solid_1dreduced_type>::element_const_iterator,
                         typename MeshTraits<mesh_solid_1dreduced_type>::element_const_iterator> range_solid_elt_1dreduced_type;


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

    typedef typename solid_type::space_stress_type space_struct_stress_type;
    typedef typename solid_type::element_stress_type element_struct_stress_type;

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

    InterpolationFSI(fluid_ptrtype fluid, solid_ptrtype solid, bool buildOperators=true)
        :
        M_fluid(fluid),
        M_solid(solid),
        M_worldComm(fluid->worldComm()),
        M_interfaceFSIisConforme( option(_name="fsi.conforming-interface").template as<bool>()),
        M_fsiCouplingType( option(_name="fsi.coupling-type").template as<std::string>()),
        M_fsiCouplingBoundaryCondition( option(_name="fsi.coupling-bc").template as<std::string>() ),
        M_opDispIsInit(false),
        M_opStressIsInit(false),
        M_opVelocityIsInit(false),
        M_verbose(fluid->verbose() || solid->verbose()),
        M_verboseAllProc(fluid->verboseAllProc() || solid->verboseAllProc())
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","constructor", "start",
                                                       this->worldComm(),this->verboseAllProc());

            bool doBuild = buildOperators && !this->fluid()->markersNameMovingBoundary().empty();
            if (M_solid->isStandardModel())
                doBuild = doBuild && !this->solid()->getMarkerNameFSI().empty();
            if ( doBuild )
            {
                boost::mpi::timer btime;double thet;
                if (M_solid->isStandardModel())
                {
                    //---------------------------------------------//
                    this->initDispInterpolation();
                    std::ostringstream ostr1;ostr1<<btime.elapsed()<<"s";
                    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initDispInterpolation","done in "+ostr1.str(),
                                                        this->worldComm(),this->verboseAllProc());
                    btime.restart();
                    //---------------------------------------------//
                    this->initStressInterpolation();
                    std::ostringstream ostr2;ostr2<<btime.elapsed()<<"s";
                    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initStressInterpolation","done in "+ostr2.str(),
                                                        this->worldComm(),this->verboseAllProc());
                    btime.restart();
                    //---------------------------------------------//
                    if (this->fsiCouplingType()=="Semi-Implicit" || this->fsiCouplingBoundaryCondition()=="robin-neumann")
                    {
                        this->initVelocityInterpolation();
                        std::ostringstream ostr3;ostr3<<btime.elapsed()<<"s";
                        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initVelocityInterpolation","done in "+ostr3.str(),
                                                            this->worldComm(),this->verboseAllProc());
                    }
                    btime.restart();
                    //---------------------------------------------//
#if 0
                    if (this->fsiCouplingBoundaryCondition()=="robin-neumann")
                    {
                        this->initStressInterpolationS2F();
                        std::ostringstream ostr4;ostr4<<btime.elapsed()<<"s";
                        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initStressInterpolationS2F","done in "+ostr4.str(),
                                                            this->worldComm(),this->verboseAllProc());

                    }
#endif
                    //---------------------------------------------//
                    if (this->fsiCouplingBoundaryCondition()=="robin-robin")
                    {
                        this->initVelocityInterpolationF2S();
                    }
                }
                else if ( M_solid->is1dReducedModel() )
                {
                    //---------------------------------------------//
                    this->initDisp1dToNdInterpolation();
                    thet = btime.elapsed();
                    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initDisp1dToNdInterpolation",
                                                        (boost::format("done in %1%s")% thet).str(),
                                                        this->worldComm(),this->verboseAllProc());
                    btime.restart();
                    //---------------------------------------------//
                    this->initStress1dToNdInterpolation();
                    thet = btime.elapsed();
                    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initStress1dToNdInterpolation",
                                                        (boost::format("done in %1%s")% thet).str(),
                                                        this->worldComm(),this->verboseAllProc());
                    btime.restart();
                    //---------------------------------------------//
                    if (this->fsiCouplingType()=="Semi-Implicit" || this->fsiCouplingBoundaryCondition()=="robin-neumann")
                    {
                        this->initVelocity1dToNdInterpolation();
                        thet = btime.elapsed();
                        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initVelocity1dToNdInterpolation",
                                                            (boost::format("done in %1%s")% thet).str(),
                                                            this->worldComm(),this->verboseAllProc());
                    }
                    //---------------------------------------------//
                }
            }
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","constructor", "finish",
                                                this->worldComm(),this->verboseAllProc());

        }

    //-----------------------------------------------------------------------------------//

    void buildOperatorDispStructToFluid()
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","buildDispStructToFluid", "start",
                                                this->worldComm(),this->verboseAllProc());
            boost::timer btime;
            if (M_solid->isStandardModel())
            {
                this->initDispInterpolation();
            }
            else if ( M_solid->is1dReducedModel() )
            {
                this->initDisp1dToNdInterpolation();
            }
            else std::cout << "[InterpolationFSI] : BUG " << std::endl;
            std::ostringstream ostr1;ostr1<<btime.elapsed()<<"s";
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","buildDispStructToFluid", "finish in "+ostr1.str(),
                                                this->worldComm(),this->verboseAllProc());
        }

    //-----------------------------------------------------------------------------------//

    void buildOperatorNormalStressFluidToStruct()
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","buildOperatorNormalStressFluidToStruct", "start",
                                                this->worldComm(),this->verboseAllProc());
            boost::timer btime;
            if ( M_solid->isStandardModel() )
            {
                this->initStressInterpolation();
            }
            else if ( M_solid->is1dReducedModel() )
            {
                this->initStress1dToNdInterpolation();
            }
            else std::cout << "[InterpolationFSI] : BUG " << std::endl;
            std::ostringstream ostr1;ostr1<<btime.elapsed()<<"s";
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","buildOperatorNormalStressFluidToStruct", "finish in "+ostr1.str(),
                                                this->worldComm(),this->verboseAllProc());
        }

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//



    void initDispInterpolation()
        {
            if (M_interfaceFSIisConforme)
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank() )
                    std::cout << "initDispInterpolation() CONFORME"  << std::endl;
                M_opDisp2dTo2dconf = opInterpolation(_domainSpace=this->solid()->functionSpaceDisplacement(),
                                                     _imageSpace=this->fluid()->getMeshALE()->displacement()->functionSpace(),
                                                     _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                     _type=InterpolationConforme(),
                                                     _backend=this->fluid()->backend() );
            }
            else
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initDispInterpolation() NONCONFORME" << std::endl;
                M_opDisp2dTo2dnonconf = opInterpolation(_domainSpace=this->solid()->functionSpaceDisplacement(),
                                                        _imageSpace=this->fluid()->getMeshALE()->displacement()->functionSpace(),
                                                        _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                        _type=InterpolationNonConforme(),
                                                        _backend=this->fluid()->backend() );
            }
        }


    //-----------------------------------------------------------------------------------//

    void initDisp1dToNdInterpolation()
    {
        if ( this->fluid()->markersNameMovingBoundary().empty() ) return;
        for ( auto const& s : this->fluid()->markersNameMovingBoundary() )
            std::cout << " ssss = " << s << "\n";

        if (M_interfaceFSIisConforme)
        {
            if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                std::cout << "initDisp1dToNdInterpolation() CONFORME"  << std::endl;
            M_opDisp1dToNdconf = opInterpolation(_domainSpace=this->solid()->fieldDisplacementVect1dReduced().functionSpace(),
                                                 _imageSpace=this->fluid()->getMeshALE()->displacement()->functionSpace(),//M_fluid->meshDisplacementOnInterface().functionSpace(),
                                                 _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                 _type=InterpolationConforme(),
                                                 _backend=M_fluid->backend() );
        }
        else
        {
            if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                std::cout << "initDisp1dToNdInterpolation() NONCONFORME" << std::endl;
            M_opDisp1dToNdnonconf = opInterpolation(_domainSpace=this->solid()->fieldDisplacementVect1dReduced().functionSpace(),
                                                    _imageSpace=this->fluid()->getMeshALE()->displacement()->functionSpace(),//M_fluid->meshDisplacementOnInterface().functionSpace(),
                                                    _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                    _type=InterpolationNonConforme(),
                                                    _backend=M_fluid->backend() );
        }
    }

    //-----------------------------------------------------------------------------------//

    void initStressInterpolation()
        {
            std::vector<int> saveActivities_stress = M_fluid->getNormalStress()->functionSpace()->worldComm().activityOnWorld();
            if (M_fluid->worldComm().globalSize()>1 && !M_fluid->functionSpace()->hasEntriesForAllSpaces() )
                M_fluid->getNormalStress()->functionSpace()->worldComm().applyActivityOnlyOn(0/*VelocityWorld*/);

            if (M_interfaceFSIisConforme)
            {
                if (this->verbose() && M_fluid->worldComm().globalRank()==M_fluid->worldComm().masterRank())
                    std::cout << "initStressInterpolation() CONFORME" << std::endl;
                M_opStress2dTo2dconf = opInterpolation(_domainSpace=this->fluid()->getNormalStress()->functionSpace(),
                                                       _imageSpace=this->solid()->normalStressFromFluid()->functionSpace(),
                                                       _range=markedfaces(this->solid()->mesh(),this->solid()->getMarkerNameFSI()),
                                                       _type=InterpolationConforme(),
                                                       _backend=M_fluid->backend() );
            }
            else
            {
                if (this->verbose() && M_fluid->worldComm().globalRank()==M_fluid->worldComm().masterRank())
                    std::cout << "initStressInterpolation() NONCONFORME" << std::endl;
                M_opStress2dTo2dnonconf = opInterpolation(_domainSpace=this->fluid()->getNormalStress()->functionSpace(),
                                                          _imageSpace=this->solid()->normalStressFromFluid()->functionSpace(),
                                                          _range=markedfaces(this->solid()->mesh(),this->solid()->getMarkerNameFSI()),
                                                          _type=InterpolationNonConforme(true,true,true,15),
                                                          _backend=M_fluid->backend() );
            }
            // revert initial activities
            if (M_fluid->worldComm().globalSize()>1  && !M_fluid->functionSpace()->hasEntriesForAllSpaces() )
                M_fluid->getNormalStress()->functionSpace()->worldComm().setIsActive(saveActivities_stress);
        }

    //-----------------------------------------------------------------------------------//

    void initStress1dToNdInterpolation()
    {
        if (M_interfaceFSIisConforme)
        {
            if (this->verbose() && M_fluid->worldComm().globalRank()==M_fluid->worldComm().masterRank())
                std::cout << "initStress1dToNdInterpolation() CONFORME" << std::endl;
            M_opStress1dToNdconf = opInterpolation(_domainSpace=M_fluid->getNormalStress()->functionSpace(),
                                                   _imageSpace=M_solid->getStressVect1dReduced().functionSpace(),
                                                   _range=elements(M_solid->mesh1dReduced()),
                                                   _type=InterpolationConforme(),
                                                   _backend=M_fluid->backend() );

        }
        else
        {
            if (this->verbose() && M_fluid->worldComm().globalRank()==M_fluid->worldComm().masterRank())
                std::cout << "initStress1dToNdInterpolation() NONCONFORME" << std::endl;
            M_opStress1dToNdnonconf = opInterpolation(_domainSpace=M_fluid->getNormalStress()->functionSpace(),
                                                      _imageSpace=M_solid->getStressVect1dReduced().functionSpace(),
                                                      _range=elements(M_solid->mesh1dReduced()),
                                                      _type=InterpolationNonConforme(),
                                                      _backend=M_fluid->backend() );
        }
    }

    //-----------------------------------------------------------------------------------//

    void initVelocityInterpolation()
        {
            if (M_interfaceFSIisConforme)
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initVelocityInterpolation() CONFORME" << std::endl;
                M_opVelocity2dTo2dconf = opInterpolation(_domainSpace=this->solid()->fieldVelocity().functionSpace(),
                                                         _imageSpace=this->fluid()->meshVelocity2().functionSpace(),
                                                         _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                         _type=InterpolationConforme(),
                                                         _backend=this->fluid()->backend() );

            }
            else
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initVelocityInterpolation() NONCONFORME" << std::endl;
                M_opVelocity2dTo2dnonconf = opInterpolation(_domainSpace=this->solid()->fieldVelocity().functionSpace(),
                                                            _imageSpace=this->fluid()->meshVelocity2().functionSpace(),
                                                            _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                            _type=InterpolationNonConforme(),
                                                            _backend=this->fluid()->backend() );
            }
        }

    //-----------------------------------------------------------------------------------//

    void initVelocity1dToNdInterpolation()
        {
            if (M_interfaceFSIisConforme)
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initVelocity1dToNdInterpolation() CONFORME" << std::endl;
                M_opVelocity1dToNdconf = opInterpolation(_domainSpace=this->solid()->fieldVelocityVect1dReduced().functionSpace(),
                                                         _imageSpace=this->fluid()->meshVelocity2().functionSpace(),
                                                         _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                         _type=InterpolationConforme(),
                                                         _backend=this->fluid()->backend() );
            }
            else
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initVelocity1dToNdInterpolation() NONCONFORME" << std::endl;
                M_opVelocity1dToNdnonconf = opInterpolation(_domainSpace=this->solid()->fieldVelocityVect1dReduced().functionSpace(),
                                                            _imageSpace=this->fluid()->meshVelocity2().functionSpace(),
                                                            _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                            _type=InterpolationNonConforme(),
                                                            _backend=this->fluid()->backend() );
            }
        }

    //-----------------------------------------------------------------------------------//

    void initVelocityInterpolationF2S()
        {

            if (M_interfaceFSIisConforme)
            {
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initVelocityInterpolationF2S() CONFORME" << std::endl;
                M_opVelocity2dTo2dconfF2S = opInterpolation(_domainSpace=this->fluid()->functionSpaceVelocity(),
                                                            _imageSpace=this->solid()->velocityInterfaceFromFluid()->functionSpace(),
                                                            _range=markedfaces(this->solid()->mesh(),this->solid()->getMarkerNameFSI()),
                                                            _type=InterpolationConforme(),
                                                            _backend=this->fluid()->backend() );

            }
            else
            {
#if 0
                if (this->verbose() && this->fluid()->worldComm().isMasterRank())
                    std::cout << "initVelocityInterpolation() NONCONFORME" << std::endl;
                M_opVelocity2dTo2dnonconfF2S = opInterpolation(_domainSpace=this->fluid()->meshVelocity2().functionSpace(),
                                                               _imageSpace=this->solid()->fieldVelocity()->functionSpace(),
                                                               _range=markedfaces(this->solid()->mesh(),this->solid()->getMarkerNameFSI()),
                                                               _type=InterpolationNonConforme(),
                                                               _backend=this->fluid()->backend() );
#endif
            }
        }

    //-----------------------------------------------------------------------------------//

    void initStressInterpolationS2F()
        {
            M_opStress2dTo2dconfS2F = opInterpolation(_domainSpace=this->solid()->normalStressFromStruct()->functionSpace(),
                                                      _imageSpace=this->fluid()->normalStressFromStruct()->functionSpace(),
                                                      _range=markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()),
                                                      _type=InterpolationConforme(),
                                                      _backend=this->fluid()->backend() );
        }

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//

    void transfertDisplacement()
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertDisplacement", "start",
                                                this->worldComm(),this->verboseAllProc());

            if (M_solid->isStandardModel())
            {
                if (M_interfaceFSIisConforme)
                {
                    if ( M_opDisp2dTo2dconf )
                        M_opDisp2dTo2dconf->apply(M_solid->fieldDisplacement(),
                                                  *(M_fluid->meshDisplacementOnInterface()) );
                }
                else
                {
                    if ( M_opDisp2dTo2dnonconf )
                        M_opDisp2dTo2dnonconf->apply(M_solid->fieldDisplacement(),
                                                     *(M_fluid->meshDisplacementOnInterface()) );
                }
            }
            else if ( M_solid->is1dReducedModel() )
            {
                M_solid->updateInterfaceDispFrom1dDisp();
                if (M_interfaceFSIisConforme)
                {
                    if ( M_opDisp1dToNdconf )
                        M_opDisp1dToNdconf->apply(M_solid->fieldDisplacementVect1dReduced(),
                                                  *(M_fluid->meshDisplacementOnInterface() ) );
                }
                else
                {
                    if ( M_opDisp1dToNdnonconf )
                        M_opDisp1dToNdnonconf->apply(M_solid->fieldDisplacementVect1dReduced(),
                                                     *(M_fluid->meshDisplacementOnInterface() ) );
                }
            }
            else std::cout << "[InterpolationFSI] : BUG " << std::endl;

            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertDisplacement", "finish",
                                                       this->worldComm(),this->verboseAllProc());
        }

    //-----------------------------------------------------------------------------------//

    void transfertStress()
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStress", "start",
                                                       this->worldComm(),this->verboseAllProc());

            M_fluid->updateNormalStress();

            std::vector<int> saveActivities_stress = M_fluid->getNormalStress()->map().worldComm().activityOnWorld();
            if (M_fluid->worldComm().globalSize()>1  && !M_fluid->functionSpace()->hasEntriesForAllSpaces()) M_fluid->getNormalStress()->map().worldComm().applyActivityOnlyOn(0/*VelocityWorld*/);

            if (M_solid->isStandardModel())
            {
                if (M_interfaceFSIisConforme)
                {
                    M_opStress2dTo2dconf->apply( *(M_fluid->getNormalStress()), *(M_solid->normalStressFromFluid()) );
                }
                else
                {
#if 1
                    M_opStress2dTo2dnonconf->apply( *(M_fluid->getNormalStress()), *(M_solid->normalStressFromFluid()) );
#else
                    //auto FluidPhysicalName = M_fluid->getMarkerNameFSI().front();
                    auto mysubmesh = createSubmesh(this->fluid()->mesh(),markedfaces(this->fluid()->mesh(),this->fluid()->markersNameMovingBoundary()/*FluidPhysicalName*/));
                    typedef Mesh<Simplex<1,1,2> > mymesh_type;
                    typedef bases<Lagrange<1, Vectorial,Continuous,PointSetFekete> > mybasis_stress_type;
                    typedef FunctionSpace<mymesh_type, mybasis_stress_type> myspace_stress_type;
                    auto myXh = myspace_stress_type::New(mysubmesh);

                    typedef bases<Lagrange<1, Vectorial,Discontinuous,PointSetFekete> > mybasis_disc_stress_type;
                    typedef FunctionSpace<mymesh_type, mybasis_disc_stress_type> myspace_disc_stress_type;
                    auto myXhdisc = myspace_disc_stress_type::New(mysubmesh);
                    auto myudisc = myXhdisc->element();

                    auto myopStressConv = opInterpolation(_domainSpace=M_fluid->getNormalStress()->functionSpace(),
                                                          _imageSpace=myXhdisc,
                                                          _range=elements(mysubmesh),
                                                          _type=InterpolationConforme(),
                                                          _backend=M_fluid->backend() );

                    myopStressConv->apply(*M_fluid->getNormalStress(),myudisc);


                    auto opProj = opProjection(_domainSpace=myXh,
                                               _imageSpace=myXh,
                                               //_backend=backend_type::build( this->vm(), "", worldcomm ),
                                               _type=Feel::L2);//L2;
                    auto proj_beta = opProj->operator()(vf::trans(vf::idv(myudisc/*M_fluid->getNormalStress()*/)));

                    auto SolidPhysicalName = M_solid->getMarkerNameFSI().front();
                    auto myopStress2dTo2dnonconf = opInterpolation(_domainSpace=myXh,//M_fluid->getNormalStress()->functionSpace(),
                                                                   _imageSpace=M_solid->getStress()->functionSpace(),
                                                                   _range=markedfaces(M_solid->mesh(),SolidPhysicalName),
                                                                   _type=InterpolationNonConforme(),
                                                                   _backend=M_fluid->backend() );

                    myopStress2dTo2dnonconf->apply(proj_beta,*(M_solid->getStress()));
#endif
                }




            }
            else if ( M_solid->is1dReducedModel() )
            {
                if (M_interfaceFSIisConforme)
                {
                    if ( M_opStress1dToNdconf )
                        M_opStress1dToNdconf->apply(*(M_fluid->getNormalStress()), M_solid->getStressVect1dReduced() );
                }
                else
                {
                    if ( M_opStress1dToNdnonconf )
                        M_opStress1dToNdnonconf->apply(*(M_fluid->getNormalStress()), M_solid->getStressVect1dReduced() );
                }
                M_solid->updateInterfaceScalStressDispFromVectStress();
            }

            // revert initial activities
            if (M_fluid->worldComm().globalSize()>1  && !M_fluid->functionSpace()->hasEntriesForAllSpaces())
                M_fluid->getNormalStress()->map().worldComm().setIsActive(saveActivities_stress);

            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStress", "finish",
                                                this->worldComm(),this->verboseAllProc());

        }

    //-----------------------------------------------------------------------------------//

    void transfertVelocity()
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertVelocity", "start",
                                                this->worldComm(),this->verboseAllProc());

            if (M_solid->isStandardModel())
            {
                if (M_interfaceFSIisConforme)
                {
                    M_opVelocity2dTo2dconf->apply( M_solid->fieldVelocity(), M_fluid->meshVelocity2() );
                }
                else
                {
                    M_opVelocity2dTo2dnonconf->apply( M_solid->fieldVelocity(), M_fluid->meshVelocity2() );
                }
            }
            else if ( M_solid->is1dReducedModel() )
            {
                M_solid->updateInterfaceVelocityFrom1dVelocity();
                if (M_interfaceFSIisConforme)
                {
                    if ( M_opVelocity1dToNdconf )
                        M_opVelocity1dToNdconf->apply(M_solid->fieldVelocityVect1dReduced(),
                                                      M_fluid->meshVelocity2() );
                }
                else
                {
                    if ( M_opVelocity1dToNdnonconf )
                        M_opVelocity1dToNdnonconf->apply(M_solid->fieldVelocityVect1dReduced(),
                                                         M_fluid->meshVelocity2() );
                }
            }

            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertVelocity", "finish",
                                                this->worldComm(),this->verboseAllProc());
        }

    //-----------------------------------------------------------------------------------//

    void transfertStressS2F()
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStressS2F", "start",
                                                this->worldComm(),this->verboseAllProc());

            this->solid()->updateNormalStressFromStruct();

            M_opStress2dTo2dconfS2F->apply(*this->solid()->normalStressFromStruct(),
                                           *this->fluid()->normalStressFromStruct());

            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStressS2F", "finish",
                                                this->worldComm(),this->verboseAllProc());

        }

    void transfertVelocityF2S(bool useExtrapolation)
        {
            if ( useExtrapolation )
            {
                if (false)
                {
                    // bdf extrapolation
                    auto solExtrap = this->fluid()->timeStepBDF()->polyDeriv();
                    auto velExtrap = solExtrap.template element<0>();
                    if (M_interfaceFSIisConforme)
                    {
                        M_opVelocity2dTo2dconfF2S->apply( velExtrap,*this->solid()->velocityInterfaceFromFluid() );
                    }
                    else
                    {
                        CHECK(false) << "TODO\n";
                    }
                }
                else
                {
                    // extrap in Explicit strategies for incompressible fluid-structure interaction problems: Nitche ....
#if 0
                    auto velExtrap = this->fluid()->functionSpace()->template functionSpace<0>()->element();
                    velExtrap.add(  2.0, this->fluid()->getSolution()->template element<0>() );
                    velExtrap.add( -1.0, this->fluid()->timeStepBDF()->unknown(0).template element<0>() );
#else
                    auto velExtrap = this->fluid()->functionSpace()->template functionSpace<0>()->element();
                    velExtrap.add(  2.0, this->fluid()->timeStepBDF()->unknown(0).template element<0>() );
                    velExtrap.add( -1.0, this->fluid()->timeStepBDF()->unknown(1).template element<0>() );
#endif
                    if (M_interfaceFSIisConforme)
                    {
                        M_opVelocity2dTo2dconfF2S->apply( velExtrap,*this->solid()->velocityInterfaceFromFluid() );
                    }
                    else
                    {
                        CHECK(false) << "TODO\n";
                    }
                }
            }
            else // no extrap : take current solution
            {
                if (M_interfaceFSIisConforme)
                {
                    M_opVelocity2dTo2dconfF2S->apply( this->fluid()->timeStepBDF()->unknown(0).template element<0>(),//this->fluid()->getSolution()->template element<0>(),
                                                      *this->solid()->velocityInterfaceFromFluid() );
                }
                else
                {
                    CHECK(false) << "TODO\n";
                }
            }
        }

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
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

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//


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

};


} // namespace FeelModels
} // namespace Feel



#endif // FEELPP_MODELS_INTERPOLATION_FSI_H
