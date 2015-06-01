/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-07-17

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
   \file solidmecbase.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-17
 */

#ifndef FEELPP_SOLIDMECHANICSBASE_HPP
#define FEELPP_SOLIDMECHANICSBASE_HPP 1


#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelts/bdf.hpp>
#include <feel/feelts/newmark.hpp>

#include <feel/feeldiscr/operatorinterpolation.hpp>

#include <feel/feelmodels2/modelcore/modelnumerical.hpp>
#include <feel/feelmodels2/modelcore/markermanagement.hpp>
#include <feel/feelmodels2/solid/mechanicalpropertiesdescription.hpp>
#include <feel/feelmodels2/modelcore/options.hpp>
#include <feel/feelmodels2/modelalg/modelalgebraicfactory.hpp>

namespace Feel
{
namespace FeelModels
{
template< typename ConvexType, typename BasisDisplacementType,bool UseCstMechProp >
class SolidMechanicsBase : public ModelNumerical,
                           public MarkerManagementDirichletBC,
                           public MarkerManagementNeumannBC
{
public:
    typedef ModelNumerical super_type;

    typedef SolidMechanicsBase<ConvexType,BasisDisplacementType,UseCstMechProp> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // Standart Model
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //___________________________________________________________________________________//
    // basis
    static const uint16_type nOrder = BasisDisplacementType::nOrder;//OrderDisp;
    static const uint16_type nOrderDisplacement = nOrder;
    static const uint16_type nOrderPressure = (nOrder>1)? nOrder-1:1;
    typedef BasisDisplacementType basis_u_type;
    //typedef Lagrange<nOrderDisplacement, Vectorial,Continuous,PointSetFekete> basis_u_type;
    typedef Lagrange<nOrderPressure, Scalar,Continuous,PointSetFekete> basis_l_type;
    //typedef Lagrange<nOrder, Vectorial,Discontinuous,PointSetFekete> basis_stress_type;
    typedef Lagrange<nOrder+1, Vectorial,Discontinuous,PointSetFekete> basis_stress_type;
    typedef Lagrange<0, Vectorial,Continuous> basis_constraint_vec_type;
    //___________________________________________________________________________________//
    // displacement space
    typedef FunctionSpace<mesh_type, bases<basis_u_type>/*,double,NoPeriodicity*/> space_displacement_type;
    typedef boost::shared_ptr<space_displacement_type> space_displacement_ptrtype;
    typedef typename space_displacement_type::element_type element_displacement_type;
    typedef boost::shared_ptr<element_displacement_type> element_displacement_ptrtype;
    typedef typename space_displacement_type::element_type element_vectorial_type;
    typedef boost::shared_ptr<element_vectorial_type> element_vectorial_ptrtype;
    //___________________________________________________________________________________//
    // pressure space
    typedef FunctionSpace<mesh_type, bases<basis_l_type> > space_pressure_type;
    typedef boost::shared_ptr<space_pressure_type> space_pressure_ptrtype;
    typedef typename space_pressure_type::element_type element_pressure_type;
    typedef boost::shared_ptr<element_pressure_type> element_pressure_ptrtype;
    //___________________________________________________________________________________//
    // vectorial constraint space
    typedef FunctionSpace<mesh_type, bases<basis_constraint_vec_type> > space_constraint_vec_type;
    typedef boost::shared_ptr<space_constraint_vec_type> space_constraint_vec_ptrtype;
    typedef typename space_constraint_vec_type::element_type element_constraint_vec_type;
    typedef boost::shared_ptr<element_constraint_vec_type> element_constraint_vec_ptrtype;
    //___________________________________________________________________________________//
    // scalar stress space
    typedef Lagrange<nOrder, Scalar,Continuous,PointSetFekete> basis_stress_scal_type;
    typedef FunctionSpace<mesh_type, bases<basis_stress_scal_type>/*, double, NoPeriodicity*/> space_stress_scal_type;
    typedef boost::shared_ptr<space_stress_scal_type> space_stress_scal_ptrtype;
    typedef typename space_stress_scal_type::element_type element_stress_scal_type;
    typedef boost::shared_ptr<element_stress_scal_type> element_stress_scal_ptrtype;
    // vectorial stress space
    typedef FunctionSpace<mesh_type, bases<basis_stress_type>/*, double, NoPeriodicity*/> space_stress_type;
    typedef boost::shared_ptr<space_stress_type> space_stress_ptrtype;
    typedef typename space_stress_type::element_type element_stress_type;
    typedef boost::shared_ptr<element_stress_type> element_stress_ptrtype;
    //___________________________________________________________________________________//
    // methodsnum tool
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;
    typedef typename model_algebraic_factory_type::graph_type graph_type;
    typedef typename model_algebraic_factory_type::graph_ptrtype graph_ptrtype;
    //___________________________________________________________________________________//
    // newmark or savets(bdf) class
    typedef Newmark<space_displacement_type> newmark_displacement_type;
    typedef boost::shared_ptr<newmark_displacement_type> newmark_displacement_ptrtype;
    typedef Bdf<space_pressure_type>  savets_pressure_type;
    typedef boost::shared_ptr<savets_pressure_type> savets_pressure_ptrtype;
    //___________________________________________________________________________________//
    // functionspace for rho,coefflame1,coefflame2
#if 0
#if SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME
    typedef bases<Lagrange<0, Scalar,Continuous> > basis_scalar_P0_type;
#else
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
#endif
#endif
    typedef bases<Lagrange<0, Scalar,Continuous> > basis_scalar_P0_continuous_type;
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_discontinuous_type;
    typedef typename mpl::if_< mpl::bool_<UseCstMechProp>,
                               basis_scalar_P0_continuous_type,
                               basis_scalar_P0_discontinuous_type >::type basis_scalar_P0_type;

    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef boost::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef typename space_scalar_P0_type::element_type element_scalar_P0_type;
    typedef boost::shared_ptr<element_scalar_P0_type> element_scalar_P0_ptrtype;
    // mechanical properties desc
    typedef MechanicalPropertiesDescription<space_scalar_P0_type> mechanicalproperties_type;
    typedef boost::shared_ptr<mechanicalproperties_type> mechanicalproperties_ptrtype;
    //___________________________________________________________________________________//
    // exporter
    typedef Exporter<mesh_type,nOrderGeo> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;
    //typedef Exporter<mesh_type,nOrderGeo> gmsh_export_type;
    //typedef boost::shared_ptr<gmsh_export_type> gmsh_export_ptrtype;
#if defined(FEELPP_HAS_VTK)
    //fais comme ca car bug dans opeartorlagrangeP1 pour les champs vectorielles
    typedef FunctionSpace<mesh_type,bases<Lagrange<nOrder,Scalar,Continuous,PointSetFekete> > > space_create_ho_type;

    typedef Mesh<Simplex<nDim,1,nDim> > mesh_visu_ho_type;

    typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > > space_vectorial_visu_ho_type;
    typedef boost::shared_ptr<space_vectorial_visu_ho_type> space_vectorial_visu_ho_ptrtype;
    typedef typename space_vectorial_visu_ho_type::element_type element_vectorial_visu_ho_type;
    typedef boost::shared_ptr<element_vectorial_visu_ho_type> element_vectorial_visu_ho_ptrtype;


    typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<mesh_visu_ho_type>::element_const_iterator,
                         typename MeshTraits<mesh_visu_ho_type>::element_const_iterator> range_visu_ho_type;
    typedef boost::tuple<boost::mpl::size_t<MESH_FACES>,
                         typename MeshTraits<mesh_visu_ho_type>::location_face_const_iterator,
                         typename MeshTraits<mesh_visu_ho_type>::location_face_const_iterator> range_visu_boundaryfaces_ho_type;

    typedef OperatorInterpolation<space_displacement_type,
                                  space_vectorial_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_disp_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_disp_type> op_interpolation_visu_ho_disp_ptrtype;

    typedef OperatorInterpolation<space_stress_type,
                                  space_vectorial_visu_ho_type,
                                  range_visu_boundaryfaces_ho_type/*range_visu_ho_type*/> op_interpolation_visu_ho_normalstress_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_normalstress_type> op_interpolation_visu_ho_normalstress_ptrtype;

    //typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > > space_scalar_visu_ho_type;
    typedef typename space_vectorial_visu_ho_type::component_functionspace_type space_scalar_visu_ho_type;
    typedef boost::shared_ptr<space_scalar_visu_ho_type> space_scalar_visu_ho_ptrtype;
    typedef typename space_scalar_visu_ho_type::element_type element_scalar_visu_ho_type;
    typedef boost::shared_ptr<element_scalar_visu_ho_type> element_scalar_visu_ho_ptrtype;

    typedef OperatorInterpolation<space_pressure_type,
                                  space_scalar_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_pressure_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_pressure_type> op_interpolation_visu_ho_pressure_ptrtype;

    typedef Exporter<mesh_visu_ho_type> export_ho_type;
    typedef boost::shared_ptr<export_ho_type> export_ho_ptrtype;
#endif
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //1d reduced_part
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // mesh
    typedef Simplex<1,1/*nOrderGeo*/,/*2*/nDim> convex_1d_reduced_type;
    typedef Mesh<convex_1d_reduced_type> mesh_1d_reduced_type;
    typedef boost::shared_ptr<mesh_1d_reduced_type> mesh_1d_reduced_ptrtype;
    //___________________________________________________________________________________//
    // basis
    //typedef bases<Lagrange<nOrder, Scalar,Continuous,PointSetFekete> > basis_1d_reduced_type;
    typedef bases<Lagrange<nOrder, Vectorial,Continuous,PointSetFekete> > basis_vect_1d_reduced_type;
    typedef bases<Lagrange<nOrder+1, Vectorial,Discontinuous,PointSetFekete> > basis_stress_vect_1d_reduced_type;
    //___________________________________________________________________________________//
    // function space
    typedef FunctionSpace<mesh_1d_reduced_type, basis_vect_1d_reduced_type> space_vect_1d_reduced_type;
    typedef boost::shared_ptr<space_vect_1d_reduced_type> space_vect_1d_reduced_ptrtype;
    typedef typename space_vect_1d_reduced_type::element_type element_vect_1d_reduced_type;
    typedef boost::shared_ptr<element_vect_1d_reduced_type> element_vect_1d_reduced_ptrtype;

    typedef typename space_vect_1d_reduced_type::component_functionspace_type space_1d_reduced_type;
    typedef typename space_vect_1d_reduced_type::component_functionspace_ptrtype space_1d_reduced_ptrtype;
    typedef typename space_1d_reduced_type::element_type element_1d_reduced_type;
    typedef boost::shared_ptr<element_1d_reduced_type> element_1d_reduced_ptrtype;

    typedef FunctionSpace<mesh_1d_reduced_type, basis_stress_vect_1d_reduced_type> space_stress_vect_1d_reduced_type;
    typedef boost::shared_ptr<space_stress_vect_1d_reduced_type> space_stress_vect_1d_reduced_ptrtype;
    typedef typename space_stress_vect_1d_reduced_type::element_type element_stress_vect_1d_reduced_type;
    typedef boost::shared_ptr<element_stress_vect_1d_reduced_type> element_stress_vect_1d_reduced_ptrtype;

    typedef typename space_stress_vect_1d_reduced_type::component_functionspace_type space_stress_scal_1d_reduced_type;
    typedef typename space_stress_vect_1d_reduced_type::component_functionspace_ptrtype space_stress_scal_1d_reduced_ptrtype;
    typedef typename space_stress_scal_1d_reduced_type::element_type element_stress_scal_1d_reduced_type;
    typedef boost::shared_ptr<element_stress_scal_1d_reduced_type> element_stress_scal_1d_reduced_ptrtype;
    //___________________________________________________________________________________//
    // bdf class
    typedef Newmark<space_1d_reduced_type>  newmark_1d_reduced_type;
    typedef boost::shared_ptr<newmark_1d_reduced_type> newmark_1d_reduced_ptrtype;
    //___________________________________________________________________________________//
    // exporter
    typedef Exporter<mesh_1d_reduced_type,mesh_1d_reduced_type::nOrder/*nOrderGeo*/> exporter_1d_reduced_type;
    typedef boost::shared_ptr<exporter_1d_reduced_type> exporter_1d_reduced_ptrtype;
    //___________________________________________________________________________________//

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // interpolation between 1d and 2d
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//

    typedef boost::tuple<boost::mpl::size_t<MESH_FACES>,
                         typename MeshTraits<mesh_type>::marker_face_const_iterator,
                         typename MeshTraits<mesh_type>::marker_face_const_iterator> range_face_type;


    typedef OperatorInterpolation<space_vect_1d_reduced_type, space_displacement_type, range_face_type> op_interpolation1dTo2d_disp_type;
    typedef boost::shared_ptr<op_interpolation1dTo2d_disp_type> op_interpolation1dTo2d_disp_ptrtype;

    typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<mesh_1d_reduced_type>::element_const_iterator,
                         typename MeshTraits<mesh_1d_reduced_type>::element_const_iterator> range_elt1d_reduced_type;

    typedef OperatorInterpolation<space_stress_scal_type, space_1d_reduced_type ,range_elt1d_reduced_type> op_interpolation2dTo1d_normalstress_type;
    typedef boost::shared_ptr<op_interpolation2dTo1d_normalstress_type> op_interpolation2dTo1d_normalstress_ptrtype;

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//


    SolidMechanicsBase( bool __isStationary,
                        std::string prefix,
                        WorldComm const& _worldComm=WorldComm(),
                        bool __buildMesh=true,
                        std::string subPrefix="",
                        std::string appliShortRepository=option(_name="exporter.directory").as<std::string>() );

    void build();
    void build(mesh_ptrtype mesh );
    void build( mesh_1d_reduced_ptrtype mesh );
    void buildStandardModel( mesh_ptrtype mesh = mesh_ptrtype() );
    void build1dReducedModel( mesh_1d_reduced_ptrtype mesh = mesh_1d_reduced_ptrtype() );

    //virtual void loadConfigBCFile() = 0;
    virtual void loadConfigMeshFile(std::string const& geofilename) = 0;
    virtual void loadConfigMeshFile1dReduced(std::string const& geofilename) = 0;
    //void loadConfigBCFile();

    void loadParameterFromOptionsVm();
    void createWorldsComm();
    //void createModel();

    void createMesh();
    void createMesh1dReduced();
    void createFunctionSpaces();
    void createFunctionSpaces1dReduced();
    void createTimeDiscretisation();
    void createTimeDiscretisation1dReduced();
    void createExporters();
    void createExporters1dReduced();
    void createOthers();

    void createAdditionalFunctionSpacesFSI();
    void createAdditionalFunctionSpacesFSIStandard();
    void createAdditionalFunctionSpacesFSI1dReduced();

    void restartExporters();
    void restartExporters1dReduced();

    void loadMesh(mesh_ptrtype __mesh );


    void changeRepository(std::string pathLoc="");

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"SolidMechanicsMesh.path"); }

    //-----------------------------------------------------------------------------------//

    void init( bool buildMethodNum, typename model_algebraic_factory_type::appli_ptrtype const& app );


    boost::shared_ptr<std::ostringstream> getInfo() const;
    void printInfo() const;
    void saveInfo() const;
    void printAndSaveInfo() const { this->printInfo();this->saveInfo(); }

    void pdeType(std::string __type);
    void pdeSolver(std::string __type);
    //void materialLaw(std::string materialLaw);
    //std::string materialLaw() const;

    bool useDisplacementPressureFormulation() const { return M_useDisplacementPressureFormulation; }
    void setUseDisplacementPressureFormulation( bool b ) {
        M_useDisplacementPressureFormulation = b;
        if (M_mechanicalProperties) M_mechanicalProperties->setUseDisplacementPressureFormulation(b);
    }

    // fsi coupling
    bool useFSISemiImplicitScheme() const { return M_useFSISemiImplicitScheme; }
    void useFSISemiImplicitScheme(bool b) { M_useFSISemiImplicitScheme=b; }
    std::string couplingFSIcondition() const { return M_couplingFSIcondition; }
    void couplingFSIcondition(std::string s) { M_couplingFSIcondition=s; }

    std::list<std::string> getMarkerNameFSI() const { return M_markerNameFSI; }
    double gammaNitschFSI() const { return M_gammaNitschFSI; }
    void gammaNitschFSI(double d) { M_gammaNitschFSI=d; }
    double muFluidFSI() const { return M_muFluidFSI; }
    void muFluidFSI(double d) { M_muFluidFSI=d; }


    bool is1dReducedModel() const { return M_is1dReducedModel; }
    bool isStandardModel() const { return M_isStandardModel; }

    //-----------------------------------------------------------------------------------//
    // generic access
    //-----------------------------------------------------------------------------------//

    block_pattern_type blockPattern() const;
    bool hasExtendedPattern() const { return false; }

    backend_ptrtype backend() { return M_backend; }
    backend_ptrtype const& backend() const { return M_backend; }


    //-----------------------------------------------------------------------------------//
    // all model
    //-----------------------------------------------------------------------------------//

    boost::shared_ptr<TSBase> timeStepBase()
    {
        if (this->isStandardModel())
            return this->timeStepNewmark();
        else// if (this->is1dReducedModel())
            return this->timeStepNewmark1dReduced();
    }
    boost::shared_ptr<TSBase> timeStepBase() const
    {
        if (this->isStandardModel())
            return this->timeStepNewmark();
        else// if (this->is1dReducedModel())
            return this->timeStepNewmark1dReduced();
    }

    void updateTimeStep();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
private :
    void exportResultsImpl( double time );
    void exportResultsImplHO( double time );
public :

    void predictorDispl();

    //-----------------------------------------------------------------------------------//
    // standart model
    //-----------------------------------------------------------------------------------//

    //mesh_ptrtype mesh() { return M_mesh; }
    mesh_ptrtype const& mesh() const { return M_mesh; }

    //space_displacement_ptrtype functionSpace() { return M_Xh; }
    space_displacement_ptrtype const& functionSpace() const { return M_Xh; }
    //space_displacement_ptrtype functionSpaceDisplacement() { return M_Xh; }
    space_displacement_ptrtype const& functionSpaceDisplacement() const { return M_Xh; }
    space_pressure_ptrtype const& functionSpacePressure() const { CHECK( M_XhPressure ) << "space pressure not define";return M_XhPressure; }

    FEELPP_DEPRECATED element_displacement_ptrtype getSolution() { return U_displ_struct; }
    FEELPP_DEPRECATED element_displacement_ptrtype const& getSolution() const { return U_displ_struct; }
    FEELPP_DEPRECATED element_displacement_type & getDisplacement() { return *U_displ_struct; }
    FEELPP_DEPRECATED element_displacement_type const& getDisplacement() const { return *U_displ_struct; }
    element_displacement_type & fieldDisplacement() { return *U_displ_struct; }
    element_displacement_type const& fieldDisplacement() const { return *U_displ_struct; }
    element_displacement_ptrtype & fieldDisplacementPtr() { return U_displ_struct; }
    element_displacement_ptrtype const& fieldDisplacementPtr() const { return U_displ_struct; }
    element_pressure_type & fieldPressure() { CHECK( M_fieldPressure ) << "field pressure not define"; return *M_fieldPressure; }
    element_pressure_type const& fieldPressure() const { CHECK( M_fieldPressure ) << "field pressure not define"; return *M_fieldPressure; }
    element_pressure_ptrtype & fieldPressurePtr() { CHECK( M_fieldPressure ) << "field pressure not define"; return M_fieldPressure; }
    element_pressure_ptrtype const& fieldPressurePtr() const { CHECK( M_fieldPressure ) << "field pressure not define"; return M_fieldPressure; }

    newmark_displacement_ptrtype & timeStepNewmark() { return M_newmark_displ_struct; }
    newmark_displacement_ptrtype const& timeStepNewmark() const { return M_newmark_displ_struct; }
    savets_pressure_ptrtype const& timeStepSavetsPressure() const { CHECK( M_savetsPressure ) << "savets pressure not define"; return M_savetsPressure; }

    element_displacement_ptrtype & getVelocityPtr() { return M_newmark_displ_struct->currentVelocityPtr(); }
    element_displacement_ptrtype const& getVelocityPtr() const { return M_newmark_displ_struct->currentVelocityPtr(); }
    element_displacement_type & getVelocity() { return M_newmark_displ_struct->currentVelocity(); }
    element_displacement_type const& getVelocity() const { return M_newmark_displ_struct->currentVelocity(); }

    element_stress_ptrtype normalStressFromFluid() { return M_normalStressFromFluid; }
    element_stress_ptrtype const& normalStressFromFluid() const { return M_normalStressFromFluid; }
    element_stress_ptrtype normalStressFromStruct() { return M_normalStressFromStruct; }
    element_stress_ptrtype const& normalStressFromStruct() const { return M_normalStressFromStruct; }
    element_vectorial_ptrtype velocityInterfaceFromFluid() { return M_velocityInterfaceFromFluid; }
    element_vectorial_ptrtype const& velocityInterfaceFromFluid() const { return M_velocityInterfaceFromFluid; }

    //----------------------------------//

    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    graph_ptrtype buildMatrixGraph() const;
    int nBlockMatrixGraph() const;
    std::map<std::string,size_type> const& startDofIndexFieldsInMatrix() const { return M_startDofIndexFieldsInMatrix; }
    BlocksBaseVector<double> blockVectorSolution() { return M_blockVectorSolution; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }


    //-----------------------------------------------------------------------------------//
    // 1d reduced model
    //-----------------------------------------------------------------------------------//

    mesh_1d_reduced_ptrtype mesh1dReduced() { return M_mesh_1d_reduced; }
    mesh_1d_reduced_ptrtype const& mesh1dReduced() const { return M_mesh_1d_reduced; }

    space_1d_reduced_ptrtype functionSpace1dReduced() { return M_Xh_1d_reduced; }
    space_1d_reduced_ptrtype const& functionSpace1dReduced() const { return M_Xh_1d_reduced; }

    element_1d_reduced_type & getDisplacementScal1dReduced() { return *M_disp_1d_reduced; }
    element_1d_reduced_type const & getDisplacementScal1dReduced() const { return *M_disp_1d_reduced; }
    element_vect_1d_reduced_type & getDisplacementVect1dReduced() { return *M_disp_vect_1d_reduced; }
    element_vect_1d_reduced_type const & getDisplacementVect1dReduced() const { return *M_disp_vect_1d_reduced; }

    element_stress_scal_1d_reduced_type & getStressScal1dReduced() { return *M_stress_1d_reduced; }
    element_stress_vect_1d_reduced_type & getStressVect1dReduced() { return *M_stress_vect_1d_reduced; }

    element_1d_reduced_type & getVelocityScal1dReduced() { return *M_velocity_1d_reduced; }
    element_vect_1d_reduced_type & getVelocityVect1dReduced() { return *M_velocity_vect_1d_reduced; }

    newmark_1d_reduced_ptrtype timeStepNewmark1dReduced() { return M_newmark_displ_1d_reduced; }
    newmark_1d_reduced_ptrtype const& timeStepNewmark1dReduced() const { return M_newmark_displ_1d_reduced; }

    //-----------------------------------------------------------------------------------//

    void TransfertDisp1dTo2d(std::string __rangename);
    void TransfertStress2dTo1d(std::string __rangename);
    void precomputeTransfertStress2dTo1d(std::string __rangename);
    void precomputeTransfertDisp1dTo2d(std::string __rangename);
#if 0
    void precomputeDisp1dTo2dWithInterpolation();
    void precomputeNormalStress2dTo1dWithInterpolation();
    void transfertDisp1dTo2dWithInterpolation();
    void transfertNormalStress2dTo1dWithInterpolation();
#endif


    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//

    // assemble methods
    virtual void updateCLDirichlet(vector_ptrtype& U) const = 0;

    void updateLinearPDE(const vector_ptrtype& X, sparse_matrix_ptrtype& A, vector_ptrtype& F,
                         bool _buildCstPart,
                         sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                         bool _doClose=true,
                         bool _doBCStrongDirichlet=true ) const;

    void updateLinearGeneralisedString(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F) const;
    void updateLinearGeneralisedStringGeneralisedAlpha(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                                       bool _buildCstPart,
                                                       bool _doClose=true,
                                                       bool _doBCStrongDirichlet=true) const;

    //void updateLinearElasticity(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F);
    void updateLinearElasticityGeneralisedAlpha(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                                bool _buildCstPart,
                                                bool _doClose=true,
                                                bool _doBCStrongDirichlet=true) const;
    void updateLinearElasticityAxiSymGeneralisedAlpha(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                                      bool _buildCstPart,
                                                      bool _doClose=true,
                                                      bool _doBCStrongDirichlet=true) const;

    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J, vector_ptrtype& R,
                         bool BuildCstPart, sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                         bool _doClose=true, bool _doBCStrongDirichlet=true ) const;
    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R,
                         bool BuildCstPart, bool UseJacobianLinearTerms,
                         bool _doClose=true, bool _doBCStrongDirichlet=true ) const;

    //void updateJacobianIncompressibilityTerms( const vector_ptrtype& X /*const*/ /*element_displacement_type& U*/, sparse_matrix_ptrtype& J) const;
    //void updateResidualIncompressibilityTerms( const vector_ptrtype& X /*const*/ /*element_displacement_type& U*/, vector_ptrtype& R) const;
    void updateJacobianIncompressibilityTerms( element_displacement_type const& u, element_pressure_type const& p, sparse_matrix_ptrtype& J) const;
    void updateResidualIncompressibilityTerms( element_displacement_type const& u, element_pressure_type const& p, vector_ptrtype& R) const;

    void updateJacobianViscoElasticityTerms( /*const*/ element_displacement_type& U, sparse_matrix_ptrtype& J) const;
    void updateResidualViscoElasticityTerms( /*const*/ element_displacement_type& U, vector_ptrtype& R) const;


    virtual void updateBCDirichletStrongResidual(vector_ptrtype& R) const = 0;
    virtual void updateBCNeumannResidual( vector_ptrtype& R ) const = 0;
    virtual void updateBCRobinResidual( vector_ptrtype& R ) const = 0;
    virtual void updateBCFollowerPressureResidual(element_displacement_type const& u, vector_ptrtype& R ) const = 0;

    virtual void updateBCDirichletStrongJacobian(sparse_matrix_ptrtype& J) const = 0;
    virtual void updateBCFollowerPressureJacobian(element_displacement_type const& u, sparse_matrix_ptrtype& J) const = 0;

    virtual void updateSourceTermResidual( vector_ptrtype& R ) const = 0;

    virtual void updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const = 0;
    virtual void updateSourceTermLinearPDE( vector_ptrtype& F ) const = 0;
    virtual void updateBCNeumannLinearPDE( vector_ptrtype& F ) const = 0;
    //virtual void updateBCRobinLinearPDE( vector_ptrtype& F ) const = 0;



#if 0
    void updatePreconditioner(const vector_ptrtype& X,
                              sparse_matrix_ptrtype& J,
                              sparse_matrix_ptrtype& A_extended,
                              sparse_matrix_ptrtype& Prec) const;
#endif
    //void solve() { this->solve(true); }
    virtual void solve( bool upVelAcc=true );


    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    // physical parameters
#if 0
    double youngModulus() const { return M_youngmodulus; } //E
    void setYoungModulus( double val ) { M_youngmodulus=val;this->updateLameCoeffFromYoungPoisson(); }
    double coeffPoisson() const { return M_coeffpoisson; } //sigma;
    void setCoeffPoisson( double val ) { M_coeffpoisson=val;this->updateLameCoeffFromYoungPoisson(); }
    double rho() const { return M_rho; }
    void setRho( double val ) { M_rho=val; this->updateRho( cst(M_rho) ); }

    element_scalar_P0_ptrtype const& rhoP0() const { return M_P0Rho;}
    //element_scalar_P0_ptrtype const& coefflame1P0() const { return M_P0Coefflame1; }
    //element_scalar_P0_ptrtype const& coefflame2P0() const { return M_P0Coefflame2; }
    element_scalar_P0_ptrtype const& coefflame1Ptr() const { return this->mechanicalProperties()->coefflame1Ptr(); }
    element_scalar_P0_ptrtype const& coefflame2Ptr() const { return this->mechanicalProperties()->coefflame2Ptr(); }
#endif

    mechanicalproperties_ptrtype const& mechanicalProperties() const { return M_mechanicalProperties; }

#if 0
    void updateLameCoeffFromYoungPoisson();

    template < typename ExprT >
    void updateRho(vf::Expr<ExprT> const& __expr)
    {
        *M_P0Rho= vf::project(M_XhScalarP0,elements( M_mesh), __expr );
    }
    template < typename ExprT >
    void updateCoefflame1(vf::Expr<ExprT> const& __expr)
    {
        //*M_P0Coefflame1 = vf::project(M_XhScalarP0,elements( M_mesh), __expr );
        this->mechanicalProperties()->updateCoefflame1( __expr);
    }
    template < typename ExprT >
    void updateCoefflame2(vf::Expr<ExprT> const& __expr)
    {
        //*M_P0Coefflame2 = vf::project(M_XhScalarP0,elements( M_mesh), __expr );
        this->mechanicalProperties()->updateCoefflame2( __expr);
    }
#endif

    //-----------------------------------------------------------------------------------//

    void updateNormalStressFromStruct();

    void updatePreStress() { *U_displ_struct_prestress=*U_displ_struct; }

    void updateVelocity();

    template <typename element_fluid_ptrtype>
    void updateStressTensor(element_fluid_ptrtype fluidSol);

    void updateStressTensorBis(element_stress_ptrtype stressN);

    template <typename element_fluid_ptrtype>
    void updateVelocityInterface(element_fluid_ptrtype fluidSol);

    //usefull for 1d reduced model
    void updateInterfaceDispFrom1dDisp();
    void updateInterfaceVelocityFrom1dVelocity();
    void updateInterfaceScalStressDispFromVectStress();

    //-----------------------------------------------------------------------------------//
    // post processing computation
    void saveDataInPoint(const std::list<boost::tuple<std::string,typename mesh_type::node_type> > & __listPt, bool extrapolate=false);

    double computeMaxDisp() const;
    double computeMaxDispOnBoundary(std::string __marker) const;

    double computeIncompressibility() const;



protected:

    backend_ptrtype M_backend;

    double M_meshSize;

    // weak boundary
    //coef de penalisation
    double M_penalbc;
    // ho visu
    bool M_isHOVisu;

    // model
    std::string M_pdeType;
    std::string M_pdeSolver;
    //ASUP std::string M_materialLaw;
    bool M_useDisplacementPressureFormulation;
    bool M_is1dReducedModel;
    bool M_isStandardModel;

    //model parameters
    //ASUP double M_youngmodulus;//EE;
    //ASUP double M_coeffpoisson;//sigma;
    //ASUP double M_coefflame1;//lambda;
    //ASUP double M_coefflame2;//mu;
    //ASUP double M_rho;

    space_scalar_P0_ptrtype M_XhScalarP0;
    //ASUP element_scalar_P0_ptrtype M_P0Rho;//densite
#if 0
    element_scalar_P0_ptrtype M_P0Coefflame1;
    element_scalar_P0_ptrtype M_P0Coefflame2;
#endif

    mechanicalproperties_ptrtype M_mechanicalProperties;

    //generalised-alpha parameters
    double M_genAlpha_rho;
    double M_genAlpha_alpha_m;
    double M_genAlpha_alpha_f;
    double M_genAlpha_gamma;
    double M_genAlpha_beta;

    std::list<std::string> M_markerNameBCRobin;
    // fsi
    bool M_useFSISemiImplicitScheme;
    std::string M_couplingFSIcondition;
    std::list<std::string> M_markerNameFSI;
    double M_gammaNitschFSI;
    double M_muFluidFSI;

    std::set<std::string> M_nameFilesData;

    //-------------------------------------------//
    // standart model
    //-------------------------------------------//

    // mesh
    mesh_ptrtype M_mesh;
    MeshMover<mesh_type> M_mesh_mover;
    // function space
    space_displacement_ptrtype M_Xh;
    element_displacement_ptrtype U_displ_struct;
    space_pressure_ptrtype M_XhPressure;
    element_pressure_ptrtype M_fieldPressure;
    space_constraint_vec_ptrtype M_XhConstraintVec;
    element_constraint_vec_ptrtype M_fieldConstraintVec;
    space_displacement_ptrtype M_XhVectorial;
    element_vectorial_ptrtype U_displ_struct_prestress;
    element_vectorial_ptrtype M_velocityInterfaceFromFluid;
    // normal stress space
    space_stress_ptrtype M_XhStress;
    space_stress_scal_ptrtype M_XhScalarStress;
    element_stress_ptrtype M_normalStressFromFluid;
    element_stress_ptrtype M_normalStressFromStruct;
    // time discretisation
    newmark_displacement_ptrtype M_newmark_displ_struct;
    savets_pressure_ptrtype M_savetsPressure;
    // tool solver ( assembly+solver )
    model_algebraic_factory_ptrtype M_methodNum;
    // start dof index fields in matrix (lm,windkessel,...)
    std::map<std::string,size_type> M_startDofIndexFieldsInMatrix;
    // block vector solution
    BlocksBaseVector<double> M_blockVectorSolution;
    // exporter
    exporter_ptrtype M_exporter;
    bool M_doExportVelocity;
    bool M_doExportAcceleration;
    bool M_doExportNormalStress;
    bool M_doExportVelocityInterfaceFromFluid;
#if defined(FEELPP_HAS_VTK)
    export_ho_ptrtype M_exporter_ho;
    space_vectorial_visu_ho_ptrtype M_XhVectorialVisuHO;
    element_vectorial_visu_ho_ptrtype M_displacementVisuHO;
    op_interpolation_visu_ho_disp_ptrtype M_opIdisplacement;
    op_interpolation_visu_ho_normalstress_ptrtype M_opInormalstress;
    space_scalar_visu_ho_ptrtype M_XhScalarVisuHO;
    element_scalar_visu_ho_ptrtype M_pressureVisuHO;
    op_interpolation_visu_ho_pressure_ptrtype M_opIpressure;
#endif

    //-------------------------------------------//
    // 1d_reduced model
    //-------------------------------------------//

    // mesh
    mesh_1d_reduced_ptrtype M_mesh_1d_reduced;
    // function space
    space_1d_reduced_ptrtype M_Xh_1d_reduced;
    space_stress_vect_1d_reduced_ptrtype M_XhStressVect_1d_reduced;
    //element disp,vel,acc
    element_1d_reduced_ptrtype M_disp_1d_reduced;
    element_1d_reduced_ptrtype M_velocity_1d_reduced;
    element_1d_reduced_ptrtype M_acceleration_1d_reduced;
    // normal stress space
    element_stress_scal_1d_reduced_ptrtype M_stress_1d_reduced;
    // vectorial 1d_reduced space
    space_vect_1d_reduced_ptrtype M_Xh_vect_1d_reduced;
    element_vect_1d_reduced_ptrtype M_disp_vect_1d_reduced;
    element_vect_1d_reduced_ptrtype M_velocity_vect_1d_reduced;
    element_stress_vect_1d_reduced_ptrtype M_stress_vect_1d_reduced;
    // time discretisation
    newmark_1d_reduced_ptrtype M_newmark_displ_1d_reduced;
    // method num
    model_algebraic_factory_ptrtype M_methodNum_1d_reduced;
    // exporter
    exporter_1d_reduced_ptrtype M_exporter_1d_reduced;



#if 0
    // dof define in 1d -> pts define on the surface stress 2d (transfert stress)
    boost::shared_ptr<std::vector<typename mesh_type::node_type> > M_map_stress_1d_reduced;
    // dof define in 2d -> pts define on the 1d mesh (transfert diplacement)
    boost::shared_ptr<std::vector<typename mesh_type::node_type> > M_map_disp_1d_reduced;
    // inutile!!!
    op_interpolation1dTo2d_disp_ptrtype M_op1dTo2d_disp;
    op_interpolation2dTo1d_normalstress_ptrtype M_op2dTo1d_normalStress;
#endif


}; // SolidMechanics


#if 0
    template <typename element_stressN_ptrtype>
    void
    SOLIDMECHANICSBASE_CLASS_NAME::updateStressTensor( element_stressN_ptrtype __stressN)
    {

        if (this->verbose()) std::cout << "[SolidMechanics] : updateStressTensor start\n";

        boost::timer btime; btime.restart();

        auto opI=opInterpolation( _domainSpace=__stressN->functionSpace(),
                                  _imageSpace=M_normalStressFromFluid->functionSpace(),//M_Xh,
                                  _range=markedfaces(M_mesh,M_markerNameFSI),
                                  _backend=M_backend );
        opI->apply(*__stressN,*M_normalStressFromFluid);


#if 0
        if (M_is1dReduced)
            {
                ForEachBC( bcDef,cl::paroi_mobile,
                           TransfertStress2dTo1d(PhysicalName); )
            }
#endif
        if (this->verbose()) std::cout << "[SolidMechanics] : updateStressTensor finish in "<<btime.elapsed()<<"\n";

    }
#endif

} // FeelModels

} // Feel
#endif /* FEELPP_SOLIDMECHANICSBASE_HPP */
