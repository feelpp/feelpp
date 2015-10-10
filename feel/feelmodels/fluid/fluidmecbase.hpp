/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
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
 \file fluidmecbase.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#ifndef FEELPP_FLUIDMECHANICSBASE_HPP
#define FEELPP_FLUIDMECHANICSBASE_HPP 1


#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/projectors.hpp>


#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

#include <feel/feelmodels/fluid/densityviscositymodel.hpp>

#if defined( FEELPP_MODELS_HAS_MESHALE )
#include <feel/feelmodels/modelmesh/meshale.hpp>
#endif




namespace Feel
{
namespace FeelModels
{
template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity=false>
class FluidMechanicsBase : public ModelNumerical,
                           public MarkerManagementDirichletBC,
                           public MarkerManagementNeumannBC,
                           public MarkerManagementALEMeshBC,
                           public MarkerManagementSlipBC,
                           public MarkerManagementPressureBC
{
public:
    typedef ModelNumerical super_type;

    typedef FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity > self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // trace mesh
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // basis fluid
    static const bool useMixedBasis = true;
    static const uint16_type nOrderVelocity = BasisVelocityType::nOrder;
    static const uint16_type nOrderPressure = BasisPressureType::nOrder;
    typedef BasisVelocityType basis_fluid_u_type;
    typedef BasisPressureType basis_fluid_p_type;
    typedef Lagrange<0, Scalar,Continuous> basis_l_type;
    //___________________________________________________________________________________//
    // mixed basis
    typedef bases<basis_fluid_u_type,basis_fluid_p_type> basis_fluid_type;
    //___________________________________________________________________________________//
    // function space
    typedef FunctionSpace<mesh_type, basis_fluid_type, Periodicity< Periodic<>, NoPeriodicity > > space_fluid_periodic_type;
    typedef FunctionSpace<mesh_type, basis_fluid_type> space_fluid_nonperiodic_type;
    typedef typename mpl::if_< mpl::bool_<UsePeriodicity>,
                               space_fluid_periodic_type,
                               space_fluid_nonperiodic_type >::type space_fluid_type;
    typedef boost::shared_ptr<space_fluid_type> space_fluid_ptrtype;
    typedef typename space_fluid_type::element_type element_fluid_type;
    typedef boost::shared_ptr<element_fluid_type> element_fluid_ptrtype;
    // subspace velocity
    typedef typename space_fluid_type::template sub_functionspace<0>::type space_fluid_velocity_type;
    typedef typename space_fluid_type::template sub_functionspace<0>::ptrtype space_fluid_velocity_ptrtype;
    typedef typename element_fluid_type::template sub_element<0>::type element_fluid_velocity_type;
    typedef boost::shared_ptr<element_fluid_velocity_type> element_fluid_velocity_ptrtype;
    // subspace pressure
    typedef typename space_fluid_type::template sub_functionspace<1>::type space_fluid_pressure_type;
    typedef typename space_fluid_type::template sub_functionspace<1>::ptrtype space_fluid_pressure_ptrtype;
    typedef typename element_fluid_type::template sub_element<1>::type element_fluid_pressure_type;
    typedef boost::shared_ptr<element_fluid_pressure_type> element_fluid_pressure_ptrtype;
    // function space for lagrange multiplier which impose the mean pressure
    typedef FunctionSpace<mesh_type, bases<basis_l_type> > space_meanpressurelm_type;
    typedef boost::shared_ptr<space_meanpressurelm_type> space_meanpressurelm_ptrtype;
    // function space for Diriclet condition using lagrange multiplier
    typedef FunctionSpace<trace_mesh_type, bases<basis_fluid_u_type> > space_dirichletlm_velocity_type;
    typedef boost::shared_ptr<space_dirichletlm_velocity_type> space_dirichletlm_velocity_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // ALE
#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef MeshALE<convex_type> mesh_ale_type;
    typedef boost::shared_ptr<mesh_ale_type> mesh_ale_ptrtype;
    // ref ALE mesh
    typedef typename mesh_ale_type::mesh_ref_type mesh_ref_type;
    typedef typename mesh_ale_type::mesh_ref_ptrtype mesh_ref_ptrtype;
    // mesh disp
    typedef typename mesh_ale_type::ale_map_functionspace_type space_mesh_disp_type;
    typedef typename mesh_ale_type::ale_map_element_type element_mesh_disp_type;
    typedef boost::shared_ptr<element_mesh_disp_type> element_mesh_disp_ptrtype;
    // mesh velocity (whole domain)
    typedef typename mesh_ale_type::ale_map_element_type element_meshvelocity_type;
    typedef boost::shared_ptr<element_meshvelocity_type> element_meshvelocity_ptrtype;
    // mesh velocity on FSI boundary
    typedef FunctionSpace<mesh_type, bases<basis_fluid_u_type> > space_meshvelocityonboundary_type;
    typedef boost::shared_ptr<space_meshvelocityonboundary_type> space_meshvelocityonboundary_ptrtype;
    typedef typename space_meshvelocityonboundary_type::element_type element_meshvelocityonboundary_type;
    typedef boost::shared_ptr<element_meshvelocityonboundary_type> element_meshvelocityonboundary_ptrtype;
    // save a ALE part of normal stress (usefull semi implicit)
    typedef typename mesh_ale_type::ale_map_functionspacedisc_type space_alemapdisc_type;
    typedef typename mesh_ale_type::ale_map_functionspacedisc_ptrtype space_alemapdisc_ptrtype;
    typedef typename mesh_ale_type::ale_map_elementdisc_type element_alemapdisc_type;
    typedef typename mesh_ale_type::ale_map_elementdisc_ptrtype element_alemapdisc_ptrtype;
    // case where structure displacement is scalar!
    typedef typename space_mesh_disp_type::component_functionspace_type space_mesh_disp_scalar_type;
    typedef boost::shared_ptr<space_mesh_disp_scalar_type> space_mesh_disp_scalar_ptrtype;
    typedef typename space_mesh_disp_scalar_type::element_type element_mesh_disp_scalar_type;
    typedef boost::shared_ptr<element_mesh_disp_scalar_type> element_mesh_disp_scalar_ptrtype;
#endif
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // function space stress
    //typedef bases<Lagrange<nOrderVelocity-1+space_alemapdisc_type::basis_type::nOrder, Vectorial,Discontinuous,PointSetFekete> > basis_stress_type;
    typedef bases<Lagrange<nOrderVelocity-1+mesh_type::nOrder, Vectorial,Discontinuous,PointSetFekete> > basis_stress_type;
    typedef FunctionSpace<mesh_type, basis_stress_type> space_stress_type;
    typedef boost::shared_ptr<space_stress_type> space_stress_ptrtype;
    typedef typename space_stress_type::element_type element_stress_type;
    typedef boost::shared_ptr<element_stress_type> element_stress_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // function space vorticity
    typedef typename mpl::if_< mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,
                               Lagrange<nOrderVelocity, Scalar,Continuous,PointSetFekete>,
                               Lagrange<nOrderVelocity, Vectorial,Continuous,PointSetFekete> >::type basis_vorticity_type;
    typedef FunctionSpace<mesh_type, basis_vorticity_type> space_vorticity_type;
    typedef boost::shared_ptr<space_vorticity_type> space_vorticity_ptrtype;
    typedef typename space_vorticity_type::element_type element_vorticity_type;
    typedef boost::shared_ptr<element_vorticity_type> element_vorticity_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // functionspace for rho, mu, nu
    typedef bases<BasisDVType> basis_densityviscosity_type;
    static const uint16_type nOrderDensityViscosity = BasisDVType::nOrder;
    typedef FunctionSpace<mesh_type, basis_densityviscosity_type> space_densityviscosity_type;
    // viscosity model desc
    typedef DensityViscosityModel<space_densityviscosity_type> densityviscosity_model_type;
    typedef boost::shared_ptr<densityviscosity_model_type> densityviscosity_model_ptrtype;

    typedef bases<Lagrange<nOrderVelocity, Vectorial,Continuous,PointSetFekete> > basis_vectorial_PN_type;
    typedef FunctionSpace<mesh_type, basis_vectorial_PN_type> space_vectorial_PN_type;
    typedef boost::shared_ptr<space_vectorial_PN_type> space_vectorial_PN_ptrtype;
    typedef typename space_vectorial_PN_type::element_type element_vectorial_PN_type;
    typedef boost::shared_ptr<element_vectorial_PN_type> element_vectorial_PN_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // algebraic tools
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;
    typedef typename model_algebraic_factory_type::graph_type graph_type;
    typedef typename model_algebraic_factory_type::graph_ptrtype graph_ptrtype;
    typedef typename model_algebraic_factory_type::indexsplit_type indexsplit_type;
    typedef typename model_algebraic_factory_type::indexsplit_ptrtype indexsplit_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // time
    typedef Bdf<space_fluid_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // windkessel model
    typedef bases<Lagrange<0, Scalar,Continuous>,Lagrange<0, Scalar,Continuous> > basis_fluidoutlet_windkessel_type;
    typedef FunctionSpace<trace_mesh_type, basis_fluidoutlet_windkessel_type > space_fluidoutlet_windkessel_type;
    typedef boost::shared_ptr<space_fluidoutlet_windkessel_type> space_fluidoutlet_windkessel_ptrtype;
    typedef typename space_fluidoutlet_windkessel_type::element_type element_fluidoutlet_windkessel_type;
    typedef boost::shared_ptr<element_fluidoutlet_windkessel_type> element_fluidoutlet_windkessel_ptrtype;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef typename MeshALE<typename trace_mesh_type::shape_type>::ale_map_functionspace_type space_fluidoutlet_windkessel_mesh_disp_type;
    typedef boost::shared_ptr<space_fluidoutlet_windkessel_mesh_disp_type> space_fluidoutlet_windkessel_mesh_disp_ptrtype;
    typedef typename space_fluidoutlet_windkessel_mesh_disp_type::element_type element_fluidoutlet_windkessel_mesh_disp_type;
    typedef boost::shared_ptr<element_fluidoutlet_windkessel_mesh_disp_type> element_fluidoutlet_windkessel_mesh_disp_ptrtype;
    typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<trace_mesh_type>::element_const_iterator,
                         typename MeshTraits<trace_mesh_type>::element_const_iterator> range_fluidoutlet_windkessel_type;
    typedef OperatorInterpolation<space_mesh_disp_type,
                                  space_fluidoutlet_windkessel_mesh_disp_type,
                                  range_fluidoutlet_windkessel_type> op_interpolation_fluidoutlet_windkessel_meshdisp_type;
    typedef boost::shared_ptr<op_interpolation_fluidoutlet_windkessel_meshdisp_type> op_interpolation_fluidoutlet_windkessel_meshdisp_ptrtype;
#endif
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // export
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef Exporter<trace_mesh_type,nOrderGeo> export_trace_type;
    typedef boost::shared_ptr<export_trace_type> export_trace_ptrtype;
    //typedef Exporter<mesh_type,nOrderGeo> gmsh_export_type;
    //typedef boost::shared_ptr<gmsh_export_type> gmsh_export_ptrtype;
    //___________________________________________________________________________________//
    // export ho
#if defined(FEELPP_HAS_VTK)
    //fais comme ca car bug dans opeartorlagrangeP1 pour les champs vectorielles
    typedef FunctionSpace<mesh_type,bases<Lagrange<nOrderVelocity,Scalar,Continuous,PointSetFekete> > > space_create_ho_type;
    // mesh
    typedef Mesh<Simplex<nDim,1,nDim> > mesh_visu_ho_type;
    //function space vectorial
    typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > > space_vectorial_visu_ho_type;
    typedef boost::shared_ptr<space_vectorial_visu_ho_type> space_vectorial_visu_ho_ptrtype;
    typedef typename space_vectorial_visu_ho_type::element_type element_vectorial_visu_ho_type;
    typedef boost::shared_ptr<element_vectorial_visu_ho_type> element_vectorial_visu_ho_ptrtype;
    // function space scalar
    //typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > > space_scalar_visu_ho_type;
    typedef typename space_vectorial_visu_ho_type::component_functionspace_type space_scalar_visu_ho_type;
    typedef boost::shared_ptr<space_scalar_visu_ho_type> space_scalar_visu_ho_ptrtype;
    typedef typename space_scalar_visu_ho_type::element_type element_scalar_visu_ho_type;
    typedef boost::shared_ptr<element_scalar_visu_ho_type> element_scalar_visu_ho_ptrtype;
    // function space vectorial discontinuos
    typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1, Vectorial,Discontinuous,PointSetFekete> > > space_vectorialdisc_visu_ho_type;
    typedef boost::shared_ptr<space_vectorialdisc_visu_ho_type> space_vectorialdisc_visu_ho_ptrtype;
    typedef typename space_vectorialdisc_visu_ho_type::element_type element_vectorialdisc_visu_ho_type;
    typedef boost::shared_ptr<element_vectorialdisc_visu_ho_type> element_vectorialdisc_visu_ho_ptrtype;
    //___________________________________________________________________________________//
    //
    typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<mesh_visu_ho_type>::element_const_iterator,
                         typename MeshTraits<mesh_visu_ho_type>::element_const_iterator> range_visu_ho_type;
    //___________________________________________________________________________________//

    typedef OperatorInterpolation<space_fluid_velocity_type,
                                  space_vectorial_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_vectorial_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_vectorial_type> op_interpolation_visu_ho_vectorial_ptrtype;

    typedef OperatorInterpolation<space_fluid_pressure_type,
                                  space_scalar_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_scalar_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_scalar_type> op_interpolation_visu_ho_scalar_ptrtype;

#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef OperatorInterpolation<space_mesh_disp_type,
                                  space_vectorial_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_meshdisp_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_meshdisp_type> op_interpolation_visu_ho_meshdisp_ptrtype;
#endif

    typedef OperatorInterpolation<space_stress_type,
                                  space_vectorialdisc_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_vectorialdisc_type;
    typedef boost::shared_ptr<op_interpolation_visu_ho_vectorialdisc_type> op_interpolation_visu_ho_vectorialdisc_ptrtype;
    //___________________________________________________________________________________//

    typedef Exporter<mesh_visu_ho_type> export_ho_type;
    typedef boost::shared_ptr<export_ho_type> export_ho_ptrtype;
#endif
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//

    //___________________________________________________________________________________//
    // constructor
    FluidMechanicsBase( //bool __isStationary,
                        std::string prefix,
                        bool __buildMesh=true,
                        WorldComm const& _worldComm=Environment::worldComm(),
                        std::string subPrefix="",
                        std::string appliShortRepository=soption(_name="exporter.directory") );
    FluidMechanicsBase( self_type const & M ) = default;
    //___________________________________________________________________________________//

    static std::string expandStringFromSpec( std::string const& expr );


    void build();
    void init( bool buildMethodNum, typename model_algebraic_factory_type::appli_ptrtype const& app );

    virtual void loadConfigBCFile() = 0;
    virtual void loadConfigMeshFile(std::string const& geofilename) = 0;

    void loadParameterFromOptionsVm();
    void createWorldsComm();
    void createALE();
    void createMesh();
    void createFunctionSpaces();
    void createTimeDiscretisation();
    void createExporters();
    void createOthers();
    void createFunctionSpacesNormalStress();
    void createFunctionSpacesVorticity();
    void createFunctionSpacesSourceAdded();

    void restartExporters();

    void loadMesh(mesh_ptrtype __mesh );

    void updateMarkedZonesInMesh();

    boost::shared_ptr<std::ostringstream> getInfo() const;
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"FluidMechanicsMesh.path"); }
    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }

    space_fluid_ptrtype const& functionSpace() const { return M_Xh; }
    space_fluid_velocity_ptrtype const/*&*/ functionSpaceVelocity() const { return M_Xh->template functionSpace<0>(); }
    space_fluid_pressure_ptrtype const/*&*/ functionSpacePressure() const { return M_Xh->template functionSpace<1>(); }

    element_fluid_ptrtype & fieldVelocityPressurePtr() { return M_Solution; }
    element_fluid_ptrtype const& fieldVelocityPressurePtr() const { return M_Solution; }
    element_fluid_type & fieldVelocityPressure() { return *M_Solution; }
    element_fluid_type const& fieldVelocityPressure() const { return *M_Solution; }
    element_fluid_velocity_type & fieldVelocity() { return M_Solution->template element<0>(); }
    element_fluid_velocity_type const& fieldVelocity() const { return M_Solution->template element<0>(); }
    element_fluid_pressure_type & fieldPressure() { return M_Solution->template element<1>(); }
    element_fluid_pressure_type const& fieldPressure() const { return M_Solution->template element<1>(); }

    element_stress_ptrtype & fieldNormalStressPtr() { return M_fieldNormalStress; }
    element_stress_type const& fieldNormalStress() const { return *M_fieldNormalStress; }
    element_stress_ptrtype & fieldNormalStressRefMeshPtr() { return M_fieldNormalStressRefMesh; }
    element_stress_type const& fieldNormalStressRefMesh() const { return *M_fieldNormalStressRefMesh; }
    element_stress_ptrtype & fieldWallShearStressPtr() { return M_fieldWallShearStress; }
    element_stress_type const& fieldWallShearStress() const { return *M_fieldWallShearStress; }

    bool useExtendedDofTable() const;

    //___________________________________________________________________________________//
    // algebraic data
    backend_ptrtype backend() { return M_backend; }
    backend_ptrtype const& backend() const { return  M_backend; }
    block_pattern_type blockPattern() const;
    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    graph_ptrtype buildMatrixGraph() const;
    int nBlockMatrixGraph() const;
    indexsplit_ptrtype buildIndexSplit() const;
    model_algebraic_factory_ptrtype algebraicFactory() { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    size_type nLocalDof() const;
    std::map<std::string,size_type> const& startDofIndexFieldsInMatrix() const { return M_startDofIndexFieldsInMatrix; }
    BlocksBaseVector<double> blockVectorSolution() { return M_blockVectorSolution; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    void updateBlockVectorSolution();

    //___________________________________________________________________________________//
    // time step scheme
    bdf_ptrtype timeStepBDF() { return M_bdf_fluid; }
    bdf_ptrtype const& timeStepBDF() const { return M_bdf_fluid; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    void updateTimeStepBDF();
    void initTimeStep();
    void updateTimeStep() { this->updateTimeStepBDF(); }


    //___________________________________________________________________________________//
    // export results
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    void setDoExport(bool b);
private :
    void exportResultsImpl( double time );
    void exportResultsImplHO( double time );
public :
    //___________________________________________________________________________________//
    // ale mesh
#if defined( FEELPP_MODELS_HAS_MESHALE )
    mesh_ale_ptrtype meshALE() { return M_meshALE; }
    mesh_ale_ptrtype const& meshALE() const { return M_meshALE; }

    element_mesh_disp_ptrtype meshDisplacementOnInterface() { return M_meshDisplacementOnInterface; }
    element_meshvelocity_type & meshVelocity() { return *M_meshALE->velocity(); }
    element_meshvelocityonboundary_type & meshVelocity2() { return *M_meshVelocityInterface; }
    element_meshvelocityonboundary_ptrtype meshVelocity2Ptr() { return M_meshVelocityInterface; }

    element_meshvelocity_type const & meshVelocity() const { return *M_meshALE->velocity(); }
    element_meshvelocityonboundary_type const & meshVelocity2() const { return *M_meshVelocityInterface; }
    element_meshvelocityonboundary_ptrtype const & meshVelocity2Ptr() const { return M_meshVelocityInterface; }

    element_stress_ptrtype normalStressFromStruct() { return M_normalStressFromStruct; }
    element_stress_ptrtype const& normalStressFromStruct() const { return M_normalStressFromStruct; }
#endif
    //element_fluid_velocity_scalar_type & meshVelocityScalOnInterface() { return *M_meshVelocityScalarOnInterface; }
    //___________________________________________________________________________________//

    bool isMoveDomain() const { return M_isMoveDomain; }

    void pdeType(std::string __type);
    std::string pdeType() const;
    void pdeSolver(std::string __type);
    std::string pdeSolver() const;
    void stressTensorLawType(std::string __type);
    //std::string stressTensorLawType() const;

    bool startBySolveNewtonian() const { return M_startBySolveNewtonian; }
    void startBySolveNewtonian( bool b ) { M_startBySolveNewtonian=b; }
    bool hasSolveNewtonianAtKickOff() const { return M_hasSolveNewtonianAtKickOff; }
    void hasSolveNewtonianAtKickOff( bool b ) { M_hasSolveNewtonianAtKickOff=b; }

    bool startBySolveStokesStationary() const { return M_startBySolveStokesStationary; }
    void startBySolveStokesStationary( bool b ) { M_startBySolveStokesStationary=b; }
    bool hasSolveStokesStationaryAtKickOff() const { return M_hasSolveStokesStationaryAtKickOff; }
    void hasSolveStokesStationaryAtKickOff( bool b ) { M_hasSolveStokesStationaryAtKickOff=b; }

    //___________________________________________________________________________________//
    // fsi parameters

    bool useFSISemiImplicitScheme() const { return M_useFSISemiImplicitScheme; }
    void useFSISemiImplicitScheme(bool b) { M_useFSISemiImplicitScheme=b; }
    std::string couplingFSIcondition() const { return M_couplingFSIcondition; }
    void couplingFSIcondition(std::string s) { M_couplingFSIcondition=s; }

    double couplingFSI_Nitsche_gamma() const { return M_couplingFSI_Nitsche_gamma; }
    void setCouplingFSI_Nitsche_gamma(double d) { M_couplingFSI_Nitsche_gamma=d; }
    double couplingFSI_Nitsche_gamma0() const { return M_couplingFSI_Nitsche_gamma0; }
    void setCouplingFSI_Nitsche_gamma0(double d) { M_couplingFSI_Nitsche_gamma0=d; }
    double couplingFSI_Nitsche_alpha() const { return M_couplingFSI_Nitsche_alpha; }
    void setCouplingFSI_Nitsche_alpha(double d) { M_couplingFSI_Nitsche_alpha=d; }

    bool couplingFSI_RNG_useInterfaceOperator() const { return M_couplingFSI_RNG_useInterfaceOperator; }
    void setCouplingFSI_RNG_useInterfaceOperator(bool b) { M_couplingFSI_RNG_useInterfaceOperator=b; }
    vector_ptrtype const& couplingFSI_RNG_interfaceOperator() const { return M_couplingFSI_RNG_interfaceOperator; }
    void setCouplingFSI_RNG_interfaceOperator(vector_ptrtype const& op ) { M_couplingFSI_RNG_interfaceOperator=op; }
    bool couplingFSI_solidIs1dReduced() const { return M_couplingFSI_solidIs1dReduced; }
    void setCouplingFSI_solidIs1dReduced(bool b) { M_couplingFSI_solidIs1dReduced=b; }
    double couplingFSI_RNG_coeffForm2() const { return M_couplingFSI_RNG_coeffForm2; }
    void setCouplingFSI_RNG_coeffForm2(double d) { M_couplingFSI_RNG_coeffForm2=d; }
    element_meshvelocityonboundary_ptrtype const& couplingFSI_RNG_evalForm1() const { return M_couplingFSI_RNG_evalForm1; }
    void setCouplingFSI_RNG_evalForm1( element_meshvelocityonboundary_ptrtype const& v) { M_couplingFSI_RNG_evalForm1=v; }
    boost::shared_ptr<typename space_fluid_velocity_type::element_type> const& couplingFSI_RNG_evalForm1Bis() const { return M_couplingFSI_RNG_evalForm1Bis; }
    void setCouplingFSI_RNG_evalForm1Bis( boost::shared_ptr<typename space_fluid_velocity_type::element_type> const& v) { M_couplingFSI_RNG_evalForm1Bis=v; }
    sparse_matrix_ptrtype const& couplingFSI_RNG_matrix() const { return M_couplingFSI_RNG_matrix; }
    void couplingFSI_RNG_matrix( sparse_matrix_ptrtype const& mat ) { M_couplingFSI_RNG_matrix=mat; }
    void couplingFSI_RNG_updateForUse();
    void couplingFSI_RNG_updateLinearPDE( vector_ptrtype& F) const;

    //___________________________________________________________________________________//
    // stabilisation parameters
    bool applyCIPStabOnlyOnBoundaryFaces() const { return M_applyCIPStabOnlyOnBoundaryFaces; }
    void applyCIPStabOnlyOnBoundaryFaces(bool b) { M_applyCIPStabOnlyOnBoundaryFaces=b; }
    bool doCIPStabConvection() const { return M_doCIPStabConvection; }
    void doCIPStabConvection(bool b) { M_doCIPStabConvection=b; }
    bool doCIPStabDivergence() const { return M_doCIPStabDivergence; }
    void doCIPStabDivergence(bool b) { M_doCIPStabDivergence=b; }
    bool doCIPStabPressure() const { return M_doCIPStabPressure; }
    void doCIPStabPressure(bool b) { M_doCIPStabPressure=b; }
    double stabCIPConvectionGamma() const { return M_stabCIPConvectionGamma; }
    double stabCIPDivergenceGamma() const { return M_stabCIPDivergenceGamma; }
    double stabCIPPressureGamma() const { return M_stabCIPPressureGamma; }

    bool doStabDivDiv() const { return M_doStabDivDiv; }
    void doStabDivDiv(bool b) { M_doStabDivDiv=b; }

    bool doStabConvectionEnergy() const { return M_doStabConvectionEnergy; }
    void doStabConvectionEnergy(bool b) { M_doStabConvectionEnergy=b; }

    bool definePressureCst() const { return M_definePressureCst; }
    void setDefinePressureCst(bool b) { M_definePressureCst = b; }
    std::string definePressureCstMethod() const { return M_definePressureCstMethod; }
    void setDefinePressureCstMethod(std::string s) { M_definePressureCstMethod = s; }
    double definePressureCstPenalisationBeta() const { return M_definePressureCstPenalisationBeta; }

    //___________________________________________________________________________________//
    // physical parameters rho,mu,nu,...
    densityviscosity_model_ptrtype & densityViscosityModel() { return M_densityViscosityModel; }
    densityviscosity_model_ptrtype const& densityViscosityModel() const { return M_densityViscosityModel; }

    void updateRho(double rho)
    {
        this->densityViscosityModel()->setCstDensity(rho);
    }
    void updateMu(double mu)
    {
        this->densityViscosityModel()->setCstDynamicViscosity(mu);
    }
    template < typename ExprT >
    void updateRho(vf::Expr<ExprT> const& __expr)
    {
        this->densityViscosityModel()->updateDensity( __expr );
    }
    template < typename ExprT >
    void updateMu(vf::Expr<ExprT> const& __expr)
    {
        this->densityViscosityModel()->updateDynamicViscosity( __expr );
    }
    //___________________________________________________________________________________//
    // boundary conditions
    double dirichletBCnitscheGamma() const { return M_dirichletBCnitscheGamma; }
    void setDirichletBCnitscheGamma( double val) { M_dirichletBCnitscheGamma=val; }

    std::list<std::string> const& markersNameMovingBoundary() const { return this->markerALEMeshBC("moving"); }
    //___________________________________________________________________________________//
    // dirichlet with Lagrange multiplier
    trace_mesh_ptrtype const& meshDirichletLM() const { return M_meshDirichletLM; }
    space_dirichletlm_velocity_ptrtype const& XhDirichletLM() const { return M_XhDirichletLM; }
    //___________________________________________________________________________________//
    // impose mean pressure with P0 Lagrange multiplier
    space_meanpressurelm_ptrtype const& XhMeanPressureLM() const { return M_XhMeanPressureLM; }
    //___________________________________________________________________________________//
    // fluid outlets bc
    bool hasFluidOutlet() const { return !M_fluidOutletsBCType.empty(); }

    bool hasFluidOutletFree() const { return this->hasFluidOutlet("free"); }
    bool hasFluidOutletWindkessel() const { return this->hasFluidOutlet("windkessel"); }
    bool hasFluidOutlet(std::string const& type) const
    {
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == type )
                return true;
        return false;
    }
    bool hasFluidOutletWindkesselImplicit() const
    {
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == "windkessel" && std::get<0>( std::get<2>( outletbc ) ) == "implicit" )
                return true;
        return false;
    }
    bool hasFluidOutletWindkesselExplicit() const
    {
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == "windkessel" && std::get<0>( std::get<2>( outletbc ) ) == "explicit" )
                return true;
        return false;
    }
    int nFluidOutlet() const { return M_fluidOutletsBCType.size(); }
    int nFluidOutletWindkesselImplicit() const
    {
        int res=0;
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == "windkessel" && std::get<0>( std::get<2>( outletbc ) ) == "implicit" )
                ++res;
        return res;
    }

    std::map<int,std::vector<double> > const& fluidOutletWindkesselPressureDistalOld() const { return M_fluidOutletWindkesselPressureDistal_old; }
    trace_mesh_ptrtype const& fluidOutletWindkesselMesh() const { return M_fluidOutletWindkesselMesh; }
    space_fluidoutlet_windkessel_ptrtype const& fluidOutletWindkesselSpace() { return M_fluidOutletWindkesselSpace; }

    //___________________________________________________________________________________//

    boost::shared_ptr<typename space_fluid_pressure_type::element_type>/*element_fluid_pressure_ptrtype*/ const& velocityDiv() const { return M_velocityDiv; }
    boost::shared_ptr<typename space_fluid_pressure_type::element_type>/*element_fluid_pressure_ptrtype*/ velocityDiv() { return M_velocityDiv; }
    bool velocityDivIsEqualToZero() const { return M_velocityDivIsEqualToZero; }

    //___________________________________________________________________________________//

    // update normal stress in reference ALE mesh
    void updateNormalStressOnCurrentMesh( std::list<std::string> const& listMarkers = std::list<std::string>() );
    // update normal stress in reference ALE mesh
    void updateNormalStressOnReferenceMesh( std::list<std::string> const& listMarkers = std::list<std::string>() );
    // update normal stress (subfunctions)
    void updateNormalStressOnReferenceMeshStandard( std::list<std::string> const& listMarkers );
    void updateNormalStressOnReferenceMeshOptSI( std::list<std::string> const& listMarkers );
    void updateNormalStressOnReferenceMeshOptPrecompute( std::list<std::string> const& listMarkers );

    void updateWallShearStress();
    void updateVorticity(mpl::int_<2> /***/);
    void updateVorticity(mpl::int_<3> /***/);

    template < typename ExprT >
    void updateVelocity(vf::Expr<ExprT> const& __expr)
    {
        M_Solution->template elementPtr<0>()->on(_range=elements( this->mesh()),_expr=__expr );
    }
    template < typename ExprT >
    void updatePressure(vf::Expr<ExprT> const& __expr)
    {
        M_Solution->template elementPtr<1>()->on(_range=elements( this->mesh()),_expr=__expr );
    }

    template < typename ExprT >
    void updateSourceAdded(vf::Expr<ExprT> const& __expr)
    {
        if (!M_XhSourceAdded) this->createFunctionSpacesSourceAdded();
        M_SourceAdded->on(_range=elements( this->mesh()),_expr=__expr );
        M_haveSourceAdded=true;
    }
    template < typename ExprT >
    void updateVelocityDiv(vf::Expr<ExprT> const& __expr)
    {
        //if (!M_velocityDiv) M_velocityDiv=M_Xh->template functionSpace<1>()->elementPtr();
        if (!M_velocityDiv)
            M_velocityDiv.reset(new typename space_fluid_pressure_type::element_type/*element_fluid_pressure_type*/(M_Xh->template functionSpace<1>(),"velocityDiv") );
        //*M_velocityDiv = vf::project(_space=M_Xh->template functionSpace<1>(),_range=elements( this->mesh()),_expr=__expr);
        M_velocityDiv->on(_range=elements(this->mesh()),_expr=__expr);
        M_velocityDivIsEqualToZero=false;
    }

#if defined( FEELPP_MODELS_HAS_MESHALE )
    template <typename element_mecasol_ptrtype>
    void updateStructureDisplacement(element_mecasol_ptrtype const & structSol);

    //void updateStructureDisplacementBis( typename mesh_ale_type::ale_map_element_ptrtype disp );
    void updateALEmesh();

    template <typename element_vel_mecasol_ptrtype>
    void updateStructureVelocity(element_vel_mecasol_ptrtype velstruct);
#endif
    //___________________________________________________________________________________//

    // save in file value of pressure at point __listPt
    void savePressureAtPoints(const std::list<boost::tuple<std::string,typename mesh_type::node_type> > & __listPt, bool extrapolate=false);

    // compute drag and lift on a markedfaces called markerName
    Eigen::Matrix<value_type,nDim,1> computeForce(std::string markerName);

    double computeFlowRate(std::string marker);
    double computeMeanPressure();
    double computeMeanDivergence();
    double computeNormL2Divergence();

#if 0
    // Averaged Preassure computed on a set of slice (false for compute on actual mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/);
    // Averaged Preassure computed on a set of slice (true for compute on ref mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/);
    // Flow rate computed on a set of slice (false for compute on actual mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/);
    // Flow rate computed on a set of slice (true for compute on ref mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/);
#else
    typedef Mesh<Simplex<1,1,nRealDim> > mesh_slice1d_type;
    typedef boost::shared_ptr<mesh_slice1d_type> mesh_slice1d_ptrtype;
    typedef typename mpl::at_c<typename space_fluid_velocity_type::bases_list,0>::type basis_slice_velocity_type;
    typedef FunctionSpace<mesh_slice1d_type, bases<basis_slice_velocity_type> > space_slice_velocity_type;
    typedef OperatorInterpolation<space_fluid_velocity_type,space_slice_velocity_type> op_interp_velocity_type;
    typedef boost::shared_ptr<op_interp_velocity_type> op_interp_velocity_ptrtype;

    typedef typename mpl::at_c<typename space_fluid_pressure_type::bases_list,0>::type basis_slice_pressure_type;
    typedef FunctionSpace<mesh_slice1d_type, bases<basis_slice_pressure_type> > space_slice_pressure_type;
    typedef OperatorInterpolation<space_fluid_pressure_type,space_slice_pressure_type> op_interp_pressure_type;
    typedef boost::shared_ptr<op_interp_pressure_type> op_interp_pressure_ptrtype;

#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef typename mpl::at_c<typename space_mesh_disp_type::bases_list,0>::type basis_slice_meshdisp_type;
    typedef FunctionSpace<mesh_slice1d_type, bases<basis_slice_meshdisp_type> > space_slice_meshdisp_type;
    typedef OperatorInterpolation<space_mesh_disp_type,space_slice_meshdisp_type> op_interp_meshdisp_type;
    typedef boost::shared_ptr<op_interp_meshdisp_type> op_interp_meshdisp_ptrtype;

    std::vector<double> computeAveragedPreassure( std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
                                                  std::vector<op_interp_pressure_ptrtype> const& opInterp,
                                                  bool computeOnRefMesh=false,
                                                  std::vector<op_interp_meshdisp_ptrtype> const& opInterpMeshDisp = std::vector<op_interp_meshdisp_ptrtype>() );
    std::vector<double> computeFlowRate(std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
                                        std::vector<op_interp_velocity_ptrtype> const& opInterp,
                                        bool computeOnRefMesh=false,
                                        std::vector<op_interp_meshdisp_ptrtype> const& opInterpMeshDisp = std::vector<op_interp_meshdisp_ptrtype>() );
#else
    //TODO?
#endif

#endif
    //___________________________________________________________________________________//

    virtual void solve();

    void updateInHousePreconditioner( sparse_matrix_ptrtype const& mat, vector_ptrtype const& vecSol ) const;
    virtual void updateInHousePreconditionerPCD( sparse_matrix_ptrtype const& mat, vector_ptrtype const& vecSol ) const = 0;

    //___________________________________________________________________________________//

    void updateNewtonInitialGuess(vector_ptrtype& U) const;

    // non linear (newton)
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J , vector_ptrtype& R,
                         bool BuildCstPart,sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                         bool _doClose=true, bool _doBCStrongDirichlet=true ) const;
    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R,
                         bool BuildCstPart, bool UseJacobianLinearTerms,
                         bool _doClose=true, bool _doBCStrongDirichlet=true ) const;
    void updateJacobianModel( element_fluid_type const& U/*const vector_ptrtype& X*/, sparse_matrix_ptrtype& J , vector_ptrtype& R,
                              bool BuildCstPart ) const;
    void updateResidualModel( element_fluid_type const& U/*const vector_ptrtype& X*/, vector_ptrtype& R,
                              bool BuildCstPart, bool UseJacobianLinearTerms ) const;

    virtual void updateInitialNewtonSolutionBCDirichlet(vector_ptrtype& U) const = 0;
    virtual void updateSourceTermResidual( vector_ptrtype& R ) const = 0;
    virtual void updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J,vector_ptrtype& RBis) const = 0;
    virtual void updateBCStrongDirichletResidual(vector_ptrtype& R) const = 0;
    virtual void updateBCDirichletLagMultResidual( vector_ptrtype& R ) const = 0;
    virtual void updateBCDirichletNitscheResidual( vector_ptrtype& R ) const = 0;
    virtual void updateBCNeumannResidual( vector_ptrtype& R ) const = 0;
    virtual void updateBCPressureResidual( vector_ptrtype& R ) const = 0;


    void updateResidualStabilisation(element_fluid_type const& U/*const vector_ptrtype& X*/, vector_ptrtype& R,
                                     bool BuildCstPart, bool UseJacobianLinearTerms) const;
    void updateJacobianStabilisation(element_fluid_type const& U/*const vector_ptrtype& X*/, sparse_matrix_ptrtype& J , vector_ptrtype& R,
                                     bool BuildCstPart ) const;


    // linear
    void updateLinearPDE(const vector_ptrtype& X, sparse_matrix_ptrtype& A, vector_ptrtype& F, bool _buildCstPart,
                         sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                         bool _doClose=true, bool _doBCStrongDirichlet=true ) const;

    void updateOseen( sparse_matrix_ptrtype& A, vector_ptrtype& F, bool _buildCstPart,
                      sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                      bool _doClose=true,
                      bool _doBCStrongDirichlet=true ) const;
    void updateOseenWeakBC( sparse_matrix_ptrtype& A , vector_ptrtype& F, bool _BuildCstPart ) const;

    void updateOseenStabilisation( sparse_matrix_ptrtype& A , vector_ptrtype& F, bool _BuildCstPart,
                                   sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart ) const;

    virtual void updateSourceTermLinearPDE( vector_ptrtype& F, bool BuildCstPart ) const = 0;
    virtual void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const = 0;
    virtual void updateBCDirichletLagMultLinearPDE( vector_ptrtype& F ) const = 0;
    virtual void updateBCDirichletNitscheLinearPDE( vector_ptrtype& F ) const = 0;
    virtual void updateBCNeumannLinearPDE( vector_ptrtype& F ) const = 0;
    virtual void updateBCPressureLinearPDE( vector_ptrtype& F ) const = 0;


#if 0
    void updatePreconditioner(const vector_ptrtype& X,
                              sparse_matrix_ptrtype& A,
                              sparse_matrix_ptrtype& A_extended,
                              sparse_matrix_ptrtype& Prec) const;
#endif
    // non linear (fixed point)
    void updatePtFixe(const vector_ptrtype& Xold, sparse_matrix_ptrtype& A , vector_ptrtype& F,
                      bool _buildCstPart,
                      bool _doClose=true, bool _doBCStrongDirichlet=true ) const;

    double computeDiff(const vector_ptrtype& X1,const vector_ptrtype& X2);

    //___________________________________________________________________________________//


protected:

    bool M_hasBuildFromMesh, M_isUpdatedForUse;
    //----------------------------------------------------
    backend_ptrtype M_backend;
    //----------------------------------------------------
    // mesh
    mesh_ptrtype M_mesh;
    MeshMover<mesh_type> M_mesh_mover;
    // fluid space and solution
    space_fluid_ptrtype M_Xh;
    element_fluid_ptrtype M_Solution;
    // lagrange multiplier space for mean pressure
    space_meanpressurelm_ptrtype M_XhMeanPressureLM;
    // trace mesh and space
    trace_mesh_ptrtype M_meshDirichletLM;
    space_dirichletlm_velocity_ptrtype M_XhDirichletLM;
    // time discrtisation fluid
    bdf_ptrtype M_bdf_fluid;
    //----------------------------------------------------
    // normak boundary stress ans WSS
    space_stress_ptrtype M_XhNormalBoundaryStress;
    element_stress_ptrtype M_fieldNormalStress, M_fieldNormalStressRefMesh;
    element_stress_ptrtype M_fieldWallShearStress;
    // vorticity space
    space_vorticity_ptrtype M_Xh_vorticity;
    element_vorticity_ptrtype M_vorticity;
    //----------------------------------------------------
    // mesh ale tool and space
    bool M_isMoveDomain;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    mesh_ale_ptrtype M_meshALE;
    element_mesh_disp_ptrtype M_meshDisplacementOnInterface;
    space_meshvelocityonboundary_ptrtype M_XhMeshVelocityInterface;
    element_meshvelocityonboundary_ptrtype M_meshVelocityInterface;
    element_stress_ptrtype M_normalStressFromStruct;
    space_alemapdisc_ptrtype M_XhMeshALEmapDisc;
    element_alemapdisc_ptrtype M_saveALEPartNormalStress;
    //----------------------------------------------------
#endif
    // tool solver ( assembly+solver )
    model_algebraic_factory_ptrtype M_algebraicFactory;
    //----------------------------------------------------
    // physical properties/parameters and space
    densityviscosity_model_ptrtype M_densityViscosityModel;
    //----------------------------------------------------
    space_vectorial_PN_ptrtype M_XhSourceAdded;
    element_vectorial_PN_ptrtype M_SourceAdded;
    bool M_haveSourceAdded;
    //----------------------------------------------------
    boost::shared_ptr<typename space_fluid_pressure_type::element_type>/*element_fluid_pressure_ptrtype*/ M_velocityDiv;
    bool M_velocityDivIsEqualToZero;
    //----------------------------------------------------
    std::string M_pdeType;
    std::string M_pdeSolver;
    // fluid outlets bc
    std::vector< std::tuple<std::string,std::string, std::tuple<std::string,double,double,double> > > M_fluidOutletsBCType;

    double M_dirichletBCnitscheGamma;
    bool M_useFSISemiImplicitScheme;
    std::string M_couplingFSIcondition;
    double M_couplingFSI_Nitsche_gamma,M_couplingFSI_Nitsche_gamma0, M_couplingFSI_Nitsche_alpha;
    bool M_couplingFSI_RNG_useInterfaceOperator;
    vector_ptrtype M_couplingFSI_RNG_interfaceOperator;
    bool M_couplingFSI_solidIs1dReduced;
    double M_couplingFSI_RNG_coeffForm2;
    element_meshvelocityonboundary_ptrtype M_couplingFSI_RNG_evalForm1;
    boost::shared_ptr<typename space_fluid_velocity_type::element_type> M_couplingFSI_RNG_evalForm1Bis;
    sparse_matrix_ptrtype M_couplingFSI_RNG_matrix;

    bool M_startBySolveNewtonian, M_hasSolveNewtonianAtKickOff;
    bool M_startBySolveStokesStationary, M_hasSolveStokesStationaryAtKickOff;
    //----------------------------------------------------
    bool M_applyCIPStabOnlyOnBoundaryFaces;
    // stabilisation available
    bool M_doCIPStabConvection,M_doCIPStabDivergence,M_doCIPStabPressure;
    double M_stabCIPConvectionGamma,M_stabCIPDivergenceGamma,M_stabCIPPressureGamma;
    bool M_doStabDivDiv;
    bool M_doStabConvectionEnergy; // see Nobile thesis
    //----------------------------------------------------
    bool M_definePressureCst;
    std::string M_definePressureCstMethod;
    double M_definePressureCstPenalisationBeta;
    //----------------------------------------------------
    // fluid outlet 0d (free, windkessel)
    mutable std::map<int,double> M_fluidOutletWindkesselPressureDistal,M_fluidOutletWindkesselPressureProximal;
    std::map<int,std::vector<double> > M_fluidOutletWindkesselPressureDistal_old;
    trace_mesh_ptrtype M_fluidOutletWindkesselMesh;
    space_fluidoutlet_windkessel_ptrtype M_fluidOutletWindkesselSpace;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    space_fluidoutlet_windkessel_mesh_disp_ptrtype M_fluidOutletWindkesselSpaceMeshDisp;
    element_fluidoutlet_windkessel_mesh_disp_ptrtype M_fluidOutletWindkesselMeshDisp;
    op_interpolation_fluidoutlet_windkessel_meshdisp_ptrtype M_fluidOutletWindkesselOpMeshDisp;
    MeshMover<trace_mesh_type> M_fluidOutletWindkesselMeshMover;
#endif
    //----------------------------------------------------
    // exporter option
    bool M_isHOVisu;
    bool M_doExportVelocity, M_doExportPressure, M_doExportVorticity, M_doExportNormalStress, M_doExportWallShearStress, M_doExportViscosity;
    bool M_doExportMeshDisplacement;
    bool M_doExportMeshALE;
    bool M_doExportMeshDisplacementOnInterface;

    // exporter fluid
    export_ptrtype M_exporter;
    export_trace_ptrtype M_exporterFluidOutlet;
    // exporter fluid ho
#if defined(FEELPP_HAS_VTK)
    export_ho_ptrtype M_exporter_ho;
    space_vectorial_visu_ho_ptrtype M_XhVectorialVisuHO;
    space_scalar_visu_ho_ptrtype M_XhScalarVisuHO;
    space_vectorialdisc_visu_ho_ptrtype M_XhVectorialDiscVisuHO;

    element_vectorial_visu_ho_ptrtype M_velocityVisuHO;
    element_scalar_visu_ho_ptrtype M_pressureVisuHO;
    element_vectorial_visu_ho_ptrtype M_meshdispVisuHO;
    element_vectorialdisc_visu_ho_ptrtype M_normalStressVisuHO;
    element_vectorialdisc_visu_ho_ptrtype M_fieldWallShearStressVisuHO;

    op_interpolation_visu_ho_vectorial_ptrtype M_opIvelocity;
    op_interpolation_visu_ho_scalar_ptrtype M_opIpressure;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    op_interpolation_visu_ho_meshdisp_ptrtype M_opImeshdisp;
    MeshMover<mesh_visu_ho_type> M_meshmover_visu_ho;
#endif
    op_interpolation_visu_ho_vectorialdisc_ptrtype M_opIstress;

#endif
    //----------------------------------------------------
    // start dof index fields in matrix (lm,windkessel,...)
    std::map<std::string,size_type> M_startDofIndexFieldsInMatrix;
    // block vector solution
    BlocksBaseVector<double> M_blockVectorSolution;
    //----------------------------------------------------
    std::set<std::string> M_nameFilesPressureAtPoints;
    //----------------------------------------------------
    // overwrite assembly process : source terms
    typedef boost::function<void ( vector_ptrtype& F, bool buildCstPart )> updateSourceTermLinearPDE_function_type;
    updateSourceTermLinearPDE_function_type M_overwritemethod_updateSourceTermLinearPDE;
    typedef boost::function<void ( vector_ptrtype& R )> updateSourceTermResidual_function_type;
    updateSourceTermResidual_function_type M_overwritemethod_updateSourceTermResidual;

}; // FluidMechanics

//---------------------------------------------------------------------------------------------------------//
#if 0
template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICSBASE_CLASS_NAME::computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/)
{
    using namespace Feel::vf;

    auto solFluid = this->getSolution();

    auto p = solFluid->element<1>();

    auto setMS = setMeshSlices.getContainerOfMeshSlices();

    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];

        double area = integrate(_range=elements(meshSlice),
                                _expr=cst(1.) ).evaluate()(0,0);

        res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                      _expr=idv(p) ).evaluate()(0,0));
    }

    return res;
}

//---------------------------------------------------------------------------------------------------------//

template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICSBASE_CLASS_NAME::computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/)
{
    using namespace Feel::vf;

    this->meshALE()->revertReferenceMesh();

    auto solFluid = this->getSolution();

    auto p = solFluid->element<1>();

    // Identity matrix
    auto Id = BOOST_PP_IF(BOOST_PP_EQUAL(FLUIDMECHANICS_DIM,2),
                          oneX()*trans(oneX()) + oneY()*trans(oneY()),
                          oneX()*trans(oneX()) + oneY()*trans(oneY())+oneZ()*trans(oneZ()) );
    // Deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());
#if (FLUIDMECHANICS_DIM==2)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto detFa = Fa11*Fa22-Fa21*Fa12;
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<2,2>( Fa22,-Fa12,-Fa21,Fa11);
#endif
#if (FLUIDMECHANICS_DIM==3)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa13 = trans(Fa*oneZ())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto Fa23 = trans(Fa*oneZ())*oneY();
    auto Fa31 = trans(Fa*oneX())*oneZ();
    auto Fa32 = trans(Fa*oneY())*oneZ();
    auto Fa33 = trans(Fa*oneZ())*oneZ();
    auto detFa = Fa11*(Fa22*Fa33-Fa23*Fa32) - Fa21*(Fa12*Fa33-Fa13*Fa32) + Fa31*(Fa12*Fa23 - Fa13*Fa22);
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<3,3>( Fa22*Fa33-Fa23*Fa32 , Fa13*Fa32-Fa12*Fa33 , Fa12*Fa23-Fa13*Fa22,
    //                       Fa23*Fa31-Fa21*Fa33 , Fa11*Fa33-Fa13*Fa31 , Fa13*Fa21-Fa11*Fa23,
    //                       Fa21*Fa32-Fa22*Fa31 , Fa12*Fa31-Fa11*Fa32 , Fa11*Fa22-Fa12*Fa21
    //                        );
#endif


    auto setMS = setMeshSlices.getContainerOfMeshSlices();

    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];

        double area = integrate(_range=elements(meshSlice),
                                _expr=cst(1.) ).evaluate()(0,0);

        res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                      _expr=idv(p)*detFa ).evaluate()(0,0));
    }


    this->meshALE()->revertMovingMesh();

    return res;
}

//---------------------------------------------------------------------------------------------------------//

// Flow rate computed on a set of slice
template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICSBASE_CLASS_NAME::computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/)
{
    using namespace Feel::vf;

    auto solFluid = this->getSolution();
    auto u = solFluid->element<0>();

    auto setMS = setMeshSlices.getContainerOfMeshSlices();
    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    auto dirVelocity = vec(cst(1.0),cst(0.));
    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];
        res[i] = integrate(_range=elements(meshSlice),
                           _expr=trans(idv(u))*dirVelocity ).evaluate()(0,0);
    }

    return res;

}

//---------------------------------------------------------------------------------------------------------//

// Flow rate computed on a set of slice
template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICSBASE_CLASS_NAME::computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/)
{
    using namespace Feel::vf;

    auto solFluid = this->getSolution();
    auto u = solFluid->element<0>();

    // Identity matrix
    auto Id = BOOST_PP_IF(BOOST_PP_EQUAL(FLUIDMECHANICS_DIM,2),
                          oneX()*trans(oneX()) + oneY()*trans(oneY()),
                          oneX()*trans(oneX()) + oneY()*trans(oneY())+oneZ()*trans(oneZ()) );
    // Deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());
#if (FLUIDMECHANICS_DIM==2)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto detFa = Fa11*Fa22-Fa21*Fa12;
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<2,2>( Fa22,-Fa12,-Fa21,Fa11);
#endif
#if (FLUIDMECHANICS_DIM==3)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa13 = trans(Fa*oneZ())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto Fa23 = trans(Fa*oneZ())*oneY();
    auto Fa31 = trans(Fa*oneX())*oneZ();
    auto Fa32 = trans(Fa*oneY())*oneZ();
    auto Fa33 = trans(Fa*oneZ())*oneZ();
    auto detFa = Fa11*(Fa22*Fa33-Fa23*Fa32) - Fa21*(Fa12*Fa33-Fa13*Fa32) + Fa31*(Fa12*Fa23 - Fa13*Fa22);
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<3,3>( Fa22*Fa33-Fa23*Fa32 , Fa13*Fa32-Fa12*Fa33 , Fa12*Fa23-Fa13*Fa22,
    //                       Fa23*Fa31-Fa21*Fa33 , Fa11*Fa33-Fa13*Fa31 , Fa13*Fa21-Fa11*Fa23,
    //                       Fa21*Fa32-Fa22*Fa31 , Fa12*Fa31-Fa11*Fa32 , Fa11*Fa22-Fa12*Fa21
    //                        );
#endif



    auto setMS = setMeshSlices.getContainerOfMeshSlices();
    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    auto dirVelocity = vec(cst(1.0),cst(0.));
    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];
        res[i] = integrate(_range=elements(meshSlice),
                           _expr=trans(idv(u))*dirVelocity ).evaluate()(0,0);
    }

    return res;

}
#endif
//---------------------------------------------------------------------------------------------------------//

} // namespace FeelModels
} // namespace Feel


#endif /* FEELPP_FLUIDMECHANICSBASE_HPP */

