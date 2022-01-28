/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
 \file solidmechanics.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-17
 */

#ifndef FEELPP_TOOLBOXES_SOLIDMECHANICS_HPP
#define FEELPP_TOOLBOXES_SOLIDMECHANICS_HPP 1

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/vonmises.hpp>
#include <feel/feelvf/eig.hpp>
#include <feel/feelvf/tresca.hpp>

#include <feel/feelts/bdf.hpp>
#include <feel/feelts/newmark.hpp>

#include <feel/feeldiscr/operatorinterpolation.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

#include <feel/feelmodels/modelcore/options.hpp>

#include <feel/feelmodels/modelvf/solidmecfirstpiolakirchhoff.hpp>

#include <feel/feelmodels/solid/solidmechanics1dreduced.hpp>

namespace Feel
{
namespace FeelModels
{
/**
 * Solid Mechanics Toolbox
 * \ingroup Toolboxes
 */
template< typename ConvexType, typename BasisDisplacementType >
class SolidMechanics : public ModelNumerical,
                       public ModelPhysics<ConvexType::nDim>,
                       public std::enable_shared_from_this< SolidMechanics<ConvexType,BasisDisplacementType> >
{
    typedef ModelPhysics<ConvexType::nDim> super_physics_type;
public:
    typedef ModelNumerical super_type;
    using size_type = typename super_type::size_type;
    typedef SolidMechanics<ConvexType,BasisDisplacementType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

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
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    //___________________________________________________________________________________//
    // basis
    static const uint16_type nOrder = BasisDisplacementType::nOrder;
    static const uint16_type nOrderDisplacement = nOrder;
    static const uint16_type nOrderPressure = (nOrder>1)? nOrder-1:1;
    typedef BasisDisplacementType basis_u_type;
    typedef Lagrange<nOrderPressure, Scalar,Continuous,PointSetFekete> basis_l_type;
    typedef Lagrange<nOrder+1, Vectorial,Discontinuous,PointSetFekete> basis_stress_type;
    typedef Lagrange<nOrder+1, Tensor2,Discontinuous,PointSetFekete> basis_stress_tensor_type;
    typedef Lagrange<0, Vectorial,Continuous> basis_constraint_vec_type;
    //___________________________________________________________________________________//
    // displacement space
    typedef FunctionSpace<mesh_type, bases<basis_u_type>/*,double,NoPeriodicity*/> space_displacement_type;
    typedef std::shared_ptr<space_displacement_type> space_displacement_ptrtype;
    typedef typename space_displacement_type::element_type element_displacement_type;
    typedef std::shared_ptr<element_displacement_type> element_displacement_ptrtype;
    typedef typename space_displacement_type::element_external_storage_type element_displacement_external_storage_type;
    typedef typename space_displacement_type::element_type element_vectorial_type;
    typedef std::shared_ptr<element_vectorial_type> element_vectorial_ptrtype;
    typedef typename space_displacement_type::component_functionspace_type space_displacement_scalar_type;
    typedef typename space_displacement_scalar_type::element_type element_displacement_scalar_type;
    typedef std::shared_ptr<element_displacement_scalar_type> element_displacement_scalar_ptrtype;
    //___________________________________________________________________________________//
    // pressure space
    typedef FunctionSpace<mesh_type, bases<basis_l_type> > space_pressure_type;
    typedef std::shared_ptr<space_pressure_type> space_pressure_ptrtype;
    typedef typename space_pressure_type::element_type element_pressure_type;
    typedef std::shared_ptr<element_pressure_type> element_pressure_ptrtype;
    typedef typename space_pressure_type::element_external_storage_type element_pressure_external_storage_type;
    //___________________________________________________________________________________//
    // vectorial constraint space
    typedef FunctionSpace<mesh_type, bases<basis_constraint_vec_type> > space_constraint_vec_type;
    typedef std::shared_ptr<space_constraint_vec_type> space_constraint_vec_ptrtype;
    typedef typename space_constraint_vec_type::element_type element_constraint_vec_type;
    typedef std::shared_ptr<element_constraint_vec_type> element_constraint_vec_ptrtype;
    //___________________________________________________________________________________//
    // scalar stress space
#if 0
    typedef Lagrange<nOrder, Scalar,Continuous,PointSetFekete> basis_stress_scal_type;
    typedef FunctionSpace<mesh_type, bases<basis_stress_scal_type>/*, double, NoPeriodicity*/> space_stress_scal_type;
    typedef std::shared_ptr<space_stress_scal_type> space_stress_scal_ptrtype;
    typedef typename space_stress_scal_type::element_type element_stress_scal_type;
    typedef std::shared_ptr<element_stress_scal_type> element_stress_scal_ptrtype;
#endif
    // normal stress space
    typedef FunctionSpace<mesh_type, bases<basis_stress_type>/*, double, NoPeriodicity*/> space_normal_stress_type;
    typedef std::shared_ptr<space_normal_stress_type> space_normal_stress_ptrtype;
    typedef typename space_normal_stress_type::element_type element_normal_stress_type;
    typedef std::shared_ptr<element_normal_stress_type> element_normal_stress_ptrtype;
    // stress tensor
    typedef FunctionSpace<mesh_type, bases<basis_stress_tensor_type> > space_stress_tensor_type;
    typedef std::shared_ptr<space_stress_tensor_type> space_stress_tensor_ptrtype;
    typedef typename space_stress_tensor_type::element_type element_stress_tensor_type;
    typedef std::shared_ptr<element_stress_tensor_type> element_stress_tensor_ptrtype;
    // scalar stress space
    typedef typename space_stress_tensor_type::component_functionspace_type space_stress_scal_type;
    typedef typename space_stress_scal_type::element_type element_stress_scal_type;
    typedef std::shared_ptr<element_stress_scal_type> element_stress_scal_ptrtype;
    //___________________________________________________________________________________//
    // methodsnum tool
    // typedef ModelAlgebraicFactory model_algebraic_factory_type;
    // typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;
    // typedef typename model_algebraic_factory_type::graph_type graph_type;
    // typedef typename model_algebraic_factory_type::graph_ptrtype graph_ptrtype;
    //___________________________________________________________________________________//
    // newmark or savets(bdf) class
    typedef Newmark<space_displacement_type> newmark_displacement_type;
    typedef std::shared_ptr<newmark_displacement_type> newmark_displacement_ptrtype;
    typedef Bdf<space_displacement_type> bdf_displacement_type;
    typedef std::shared_ptr<bdf_displacement_type> bdf_displacement_ptrtype;
    typedef Bdf<space_pressure_type> savets_pressure_type;
    typedef std::shared_ptr<savets_pressure_type> savets_pressure_ptrtype;
    //___________________________________________________________________________________//
    // materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;
    //___________________________________________________________________________________//
    // trace mesh
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef std::shared_ptr<trace_mesh_type> trace_mesh_ptrtype;
    typedef Lagrange<nOrderGeo, Vectorial,Continuous> basis_tracemesh_disp_type;
    typedef FunctionSpace<trace_mesh_type, bases<basis_tracemesh_disp_type> > space_tracemesh_disp_type;
    typedef std::shared_ptr<space_tracemesh_disp_type> space_tracemesh_disp_ptrtype;
    typedef typename space_tracemesh_disp_type::element_type element_tracemesh_disp_type;
    typedef std::shared_ptr<element_tracemesh_disp_type> element_tracemesh_disp_ptrtype;
    //___________________________________________________________________________________//
    // exporter
    typedef Exporter<mesh_type,nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;
    //typedef Exporter<mesh_type,nOrderGeo> gmsh_export_type;
    //typedef std::shared_ptr<gmsh_export_type> gmsh_export_ptrtype;
#if 1 //defined(FEELPP_HAS_VTK)
    //fais comme ca car bug dans opeartorlagrangeP1 pour les champs vectorielles
    typedef FunctionSpace<mesh_type,bases<Lagrange<nOrder,Scalar,Continuous,PointSetFekete> > > space_create_ho_type;

    typedef Mesh<Simplex<nDim,1,nDim> > mesh_visu_ho_type;

    typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > > space_vectorial_visu_ho_type;
    typedef std::shared_ptr<space_vectorial_visu_ho_type> space_vectorial_visu_ho_ptrtype;
    typedef typename space_vectorial_visu_ho_type::element_type element_vectorial_visu_ho_type;
    typedef std::shared_ptr<element_vectorial_visu_ho_type> element_vectorial_visu_ho_ptrtype;


    typedef elements_reference_wrapper_t<mesh_visu_ho_type> range_visu_ho_type;
    typedef faces_reference_wrapper_t<mesh_visu_ho_type> range_visu_boundaryfaces_ho_type;

    typedef OperatorInterpolation<space_displacement_type,
                                  space_vectorial_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_disp_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_disp_type> op_interpolation_visu_ho_disp_ptrtype;

    typedef OperatorInterpolation<space_normal_stress_type,
                                  space_vectorial_visu_ho_type,
                                  range_visu_boundaryfaces_ho_type/*range_visu_ho_type*/> op_interpolation_visu_ho_normalstress_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_normalstress_type> op_interpolation_visu_ho_normalstress_ptrtype;

    //typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > > space_scalar_visu_ho_type;
    typedef typename space_vectorial_visu_ho_type::component_functionspace_type space_scalar_visu_ho_type;
    typedef std::shared_ptr<space_scalar_visu_ho_type> space_scalar_visu_ho_ptrtype;
    typedef typename space_scalar_visu_ho_type::element_type element_scalar_visu_ho_type;
    typedef std::shared_ptr<element_scalar_visu_ho_type> element_scalar_visu_ho_ptrtype;

    typedef OperatorInterpolation<space_pressure_type,
                                  space_scalar_visu_ho_type,
                                  range_visu_ho_type> op_interpolation_visu_ho_pressure_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_pressure_type> op_interpolation_visu_ho_pressure_ptrtype;

    typedef Exporter<mesh_visu_ho_type> export_ho_type;
    typedef std::shared_ptr<export_ho_type> export_ho_ptrtype;
#endif


    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //1d reduced_part
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//

    using solid_1dreduced_type = SolidMechanics1dReduced<Simplex<1,nOrderGeo,nRealDim>,basis_u_type>;
    using solid_1dreduced_ptrtype = std::shared_ptr<solid_1dreduced_type>;


    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // interpolation between 1d and 2d/3d
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//

    typedef faces_reference_wrapper_t<mesh_type> range_face_type;
    typedef OperatorInterpolation<typename solid_1dreduced_type::space_displacement_type/* space_vect_1dreduced_type*/, space_displacement_type, range_face_type> op_interpolation1dTo2d_disp_type;
    typedef std::shared_ptr<op_interpolation1dTo2d_disp_type> op_interpolation1dTo2d_disp_ptrtype;


    //___________________________________________________________________________________//

     struct FieldTag
     {
         static auto displacement( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
         static auto pressure( self_type const* t ) { return ModelFieldTag<self_type,1>( t ); }
     };

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//


    SolidMechanics( std::string const& prefix,
                    std::string const& keyword = "solid",
                    worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                    std::string const& subPrefix = "",
                    ModelBaseRepository const& modelRep = ModelBaseRepository() );
    SolidMechanics( self_type const& ) = default;

    static self_ptrtype New( std::string const& prefix,
                             std::string const& keyword = "solid",
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    static std::string expandStringFromSpec( std::string const& expr );

private :

    void loadParameterFromOptionsVm();
    void initMaterialProperties();
    void initMesh();
    void initFunctionSpaces();
    void initBoundaryConditions();

    void createExporters();


    void createAdditionalFunctionSpacesNormalStress();
    //    void createAdditionalFunctionSpacesStressTensor();


public :
    void createsSolid1dReduced();

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"SolidMechanicsMesh.path"); }

    //-----------------------------------------------------------------------------------//

    void init( bool buildAlgebraicFactory = true );
    void solve( bool upVelAcc=true );

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;


    bool hasSolidEquation1dReduced() const { return M_solid1dReduced.use_count() > 0; }
    bool hasSolidEquationStandard() const { return !this->hasSolidEquation1dReduced(); } // TOOD

    bool isStandardModel() const { return this->hasSolidEquationStandard(); } // TODO DEPRECATED
    bool is1dReducedModel() const { return this->hasSolidEquation1dReduced(); }

    bool hasDisplacementPressureFormulation() const
        {
            bool res = false;
            for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
            {
                if ( !this->materialsProperties()->hasPhysic( physicName ) )
                    continue;

                auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
                if ( physicSolidData->useDisplacementPressureFormulation() )
                {
                    res = true;
                    break;
                }
            }
            return res;
        }

    void setSolver( std::string const& type );
    //-----------------------------------------------------------------------------------//
    // all models
    //-----------------------------------------------------------------------------------//

    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { CHECK( !this->isUpdatedForUse() ) << "setMaterialsProperties can be called only before called isUpdatedForUse";  M_materialsProperties = mp; }


    bool hasDirichletBC() const
        {
            return ( !M_bcDirichlet.empty() ||
                     !M_bcDirichletComponents.find(Component::X)->second.empty() ||
                     !M_bcDirichletComponents.find(Component::Y)->second.empty() ||
                     !M_bcDirichletComponents.find(Component::Z)->second.empty() );
        }

    std::string const& timeStepping() const { return M_timeStepping; }

#if 0
    std::shared_ptr<TSBase> timeStepBase()
    {
        if ( this->hasSolidEquationStandard() )
        {
            if ( M_timeStepping == "Newmark" )
                return this->timeStepNewmark();
            else
                return this->timeStepBdfDisplacement();
        }
        else// if (this->hasSolidEquation1dReduced())
            return M_solid1dReduced->timeStepBase();
    }
#endif
    std::shared_ptr<TSBase> timeStepBase() const
    {
        if ( this->hasSolidEquationStandard() )
        {
            if ( M_timeStepping == "Newmark" )
                return this->timeStepNewmark();
            else
                return this->timeStepBdfDisplacement();
        }
        else// if (this->hasSolidEquation1dReduced())
            return M_solid1dReduced->timeStepBase();
    }
    void initTimeStep();
    void startTimeStep();
    void updateTimeStep();
    void updateVelocity();

    // post process
    void initPostProcess() override;

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr )
        {
            if ( this->hasSolidEquationStandard() )
                this->exportResults( time, this->modelFields(), symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );

            if ( this->hasSolidEquation1dReduced() )
                M_solid1dReduced->exportResults( time, symbolsExpr );
        }

    template <typename ModelFieldsType,typename SymbolsExpr>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr );


    template <typename SymbExprType>
    auto exprPostProcessExportsToolbox( SymbExprType const& se, std::string const& prefix ) const
        {
            using _expr_firstPiolaKirchhof_type = std::decay_t<decltype(Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(this->fieldDisplacement(),M_fieldPressure,std::shared_ptr<ModelPhysicSolid<nDim>>{}, this->materialsProperties()->materialProperties(""),se))>;
            using _expr_vonmises_type = std::decay_t<decltype( vonmises(_expr_firstPiolaKirchhof_type{}) )>;
            using _expr_tresca_type = std::decay_t<decltype( tresca(_expr_firstPiolaKirchhof_type{}) )>;
            using _expr_princial_stress_type = std::decay_t<decltype( eig(_expr_firstPiolaKirchhof_type{}) )>;
            std::map<std::string,std::vector<std::tuple< _expr_vonmises_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprVonMisses;
            std::map<std::string,std::vector<std::tuple< _expr_tresca_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprTresca;
            std::map<std::string,std::vector<std::tuple< _expr_princial_stress_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprPrincipalStress;
            if ( this->hasSolidEquationStandard() )
            {
                auto const& u = this->fieldDisplacement();

                for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
                {
                    auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
                    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
                    {
                        auto const& matProperties = this->materialsProperties()->materialProperties( matName );
                        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );

                        auto fpk = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(u,M_fieldPressure,physicSolidData,matProperties,se);
                        auto vonmisesExpr = vonmises( fpk );
                        mapExprVonMisses[prefixvm(prefix,"von-mises-criterion")].push_back( std::make_tuple( vonmisesExpr, range, "element" ) );

                        auto trescaExpr = tresca( fpk );
                        mapExprTresca[prefixvm(prefix,"tresca-criterion")].push_back( std::make_tuple( trescaExpr, range, "element" ) );

                        auto principalStressExpr = eig( fpk );
                        mapExprPrincipalStress[prefixvm(prefix,"principal-stresses")].push_back( std::make_tuple( principalStressExpr, range, "element" ) );
                    }
                }
            }

            return hana::make_tuple( mapExprVonMisses,mapExprTresca,mapExprPrincipalStress );
        }
    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            return hana::concat( this->materialsProperties()->exprPostProcessExports( this->mesh(), this->physicsAvailable(),se ),
                                 this->exprPostProcessExportsToolbox( se, prefix ) );
        }


#if 0
    void exportFields( double time );
    bool updateExportedFields( exporter_ptrtype exporter, std::set<std::string> const& fields, double time );
    bool updateExportedFields1dReduced( exporter_1dreduced_ptrtype exporter, std::set<std::string> const& fields, double time );
    void exportMeasures( double time );
#endif
    void restartExporters() { this->restartExporters( this->timeInitial() ); }
    void restartExporters( double time );
private :
    //void exportFieldsImpl( double time );
    // void exportFieldsImplHO( double time );

    void updateTimeStepCurrentResidual();
public :

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    void predictorDispl();

    block_pattern_type blockPattern() const override;


    //-----------------------------------------------------------------------------------//
    // standart model
    //-----------------------------------------------------------------------------------//

    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_displacement_ptrtype const& functionSpace() const { return this->functionSpaceDisplacement(); }
    space_displacement_ptrtype const& functionSpaceDisplacement() const { return M_XhDisplacement; }
    space_pressure_ptrtype const& functionSpacePressure() const { CHECK( M_XhPressure ) << "space pressure not define";return M_XhPressure; }

    element_displacement_type & fieldDisplacement() { return *M_fieldDisplacement; }
    element_displacement_type const& fieldDisplacement() const { return *M_fieldDisplacement; }
    element_displacement_ptrtype & fieldDisplacementPtr() { return M_fieldDisplacement; }
    element_displacement_ptrtype const& fieldDisplacementPtr() const { return M_fieldDisplacement; }
    element_pressure_type & fieldPressure() { CHECK( M_fieldPressure ) << "field pressure not define"; return *M_fieldPressure; }
    element_pressure_type const& fieldPressure() const { CHECK( M_fieldPressure ) << "field pressure not define"; return *M_fieldPressure; }
    element_pressure_ptrtype & fieldPressurePtr() { return M_fieldPressure; }
    element_pressure_ptrtype const& fieldPressurePtr() const { return M_fieldPressure; }

    newmark_displacement_ptrtype & timeStepNewmark() { return M_timeStepNewmark; }
    newmark_displacement_ptrtype const& timeStepNewmark() const { return M_timeStepNewmark; }
    savets_pressure_ptrtype const& timeStepSavetsPressure() const { CHECK( M_savetsPressure ) << "savets pressure not define"; return M_savetsPressure; }
    bdf_displacement_ptrtype timeStepBdfDisplacement() const { return M_timeStepBdfDisplacement; }
    bdf_displacement_ptrtype timeStepBdfVelocity() const { return M_timeStepBdfVelocity; }
    double timeStepThetaValue() const { return M_timeStepThetaValue; }

    element_displacement_ptrtype & fieldVelocityPtr() { return M_fieldVelocity; }
    element_displacement_ptrtype const& fieldVelocityPtr() const { return M_fieldVelocity; }
    element_displacement_type & fieldVelocity() { return *M_fieldVelocity;; }
    element_displacement_type const& fieldVelocity() const { return *M_fieldVelocity; }

    element_displacement_type & fieldAcceleration() { return *M_fieldAcceleration; }
    element_displacement_type const& fieldAcceleration() const { return *M_fieldAcceleration; }
    element_displacement_ptrtype const& fieldAccelerationPtr() const { return M_fieldAcceleration; }

    element_normal_stress_ptrtype & fieldNormalStressFromStructPtr() { return M_fieldNormalStressFromStruct; }
    element_normal_stress_ptrtype const& fieldNormalStressFromStructPtr() const { return M_fieldNormalStressFromStruct; }

    // stress tensor ( tensor2 )
    element_stress_tensor_ptrtype const& fieldStressTensorPtr() const { return M_fieldStressTensor; }
    element_stress_tensor_type const& fieldStressTensor() const { return *M_fieldStressTensor; }



    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return this->modelFields( this->fieldDisplacementPtr(), this->fieldVelocityPtr(), this->fieldPressurePtr(), prefix );
        }

    auto modelFields( vector_ptrtype sol, size_type startBlockSpaceIndex = 0, std::string const& prefix = "" ) const
        {
            std::shared_ptr<element_displacement_external_storage_type> field_d;
            std::shared_ptr<element_displacement_external_storage_type> field_v;
            std::shared_ptr<element_pressure_external_storage_type> field_p;
            if ( this->hasSolidEquationStandard() )
            {
                field_d = this->functionSpaceDisplacement()->elementPtr( *sol, startBlockSpaceIndex + this->startSubBlockSpaceIndex( "displacement" ) );
                if ( M_timeSteppingUseMixedFormulation )
                    field_v = this->functionSpaceDisplacement()->elementPtr( *sol, startBlockSpaceIndex + this->startSubBlockSpaceIndex( "velocity" ) );
                if ( this->hasDisplacementPressureFormulation() )
                    field_p = this->functionSpacePressure()->elementPtr( *sol, startBlockSpaceIndex + this->startSubBlockSpaceIndex( "pressure" ) );
            }
            return this->modelFields( field_d, field_v, field_p, prefix );
        }

    template <typename DisplacementFieldType, typename VelocityFieldType,typename PressureFieldType>
    auto modelFields( DisplacementFieldType const& field_s, VelocityFieldType const& field_v, PressureFieldType const& field_p, std::string const& prefix = "" ) const
        {
            //auto mfield_disp = modelField<FieldCtx::ID|FieldCtx::MAGNITUDE,element_displacement_ptrtype>( FieldTag::displacement(this) );
            auto mfield_disp = modelField<FieldCtx::FULL,DisplacementFieldType>( FieldTag::displacement(this) );
            mfield_disp.add( FieldTag::displacement(this), prefix,"displacement", field_s, "s", this->keyword() );
            if constexpr ( std::is_same_v<DisplacementFieldType,VelocityFieldType> )
                mfield_disp.add( FieldTag::displacement(this), prefix,"velocity", field_v, "v", this->keyword() );
            //mfield_disp.add( FieldTag::displacement(this), prefix,"acceleration", this->fieldAccelerationPtr(), "a", this->keyword() );
            auto mfield_pressure = modelField<FieldCtx::ID>( FieldTag::pressure(this), prefix,"pressure", field_p, "p", this->keyword() );

            return Feel::FeelModels::modelFields( mfield_disp, mfield_pressure );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            // TODO add velocity/pressure
            return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::displacement(this), "displacement", startBlockSpaceIndex + this->startSubBlockSpaceIndex( "displacement" ) ) );
        }


    //___________________________________________________________________________________//
    // symbols expression
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
#if 1
            auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            return Feel::vf::symbolsExpr( seToolbox, seParam, seMat, seFields );
#else
            return symbols_expression_empty_t{};
#endif
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType>
    auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
        {
            auto const& u = mfields.field( FieldTag::displacement(this), "displacement" );

            using _expr_fpkt_type = std::decay_t<decltype(this->firstPiolaKirchhoffTensorExpr( u ))>;
            symbol_expression_t<_expr_fpkt_type> se_fpkt;
            std::string symbol_fpkt = fmt::format("{}_stress_P",this->keyword());
            se_fpkt.add( symbol_fpkt, this->firstPiolaKirchhoffTensorExpr( u ), SymbolExprComponentSuffix(nDim,nDim) );

            return Feel::vf::symbolsExpr( se_fpkt );
        }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
        {
            return mfields.trialSymbolsExpr( tsmf );
        }


    
    template <typename DispFieldType,typename SymbolsExprType = symbols_expression_empty_t>
    auto firstPiolaKirchhoffTensorExpr( DispFieldType const& u, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            // TODO pressure
            auto p = this->fieldPressurePtr();
            using _stesstensor_expr_type = std::decay_t<decltype( Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(u,p,
                                                                                                                      std::static_pointer_cast<ModelPhysicSolid<nDim>>( this->physicsFromCurrentType().begin()->second ),
                                                                                                                      this->materialsProperties()->materialProperties(""),se) )>;
            std::vector<std::pair<std::string,_stesstensor_expr_type>> theExprs;
            for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
            {
                auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
                for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
                {
                    auto const& matProperties = this->materialsProperties()->materialProperties( matName );
                    auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );

                    auto stressTensorExpr = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(u,p,physicSolidData,matProperties,se);
                    theExprs.push_back( std::make_pair( matName, std::move( stressTensorExpr ) ) );
                }
            }
            return expr<typename mesh_type::index_type>( this->materialsProperties()->exprSelectorByMeshElementMapping(), theExprs );

        }

    //___________________________________________________________________________________//
    // model context helper
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
        {
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type,typename solid_1dreduced_type::mesh_type>();
            return Feel::FeelModels::modelContext( mfields, std::move( se ) );
        }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type,typename solid_1dreduced_type::mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
    auto modelContext( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type,typename solid_1dreduced_type::mesh_type>();
            auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( rowStartInVector ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }

    //----------------------------------//
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;

    void updateMassMatrixLumped();
    bool useMassMatrixLumped() const { return M_useMassMatrixLumped; }
    void setUseMassMatrixLumped( bool b ) { M_useMassMatrixLumped = b; }
    sparse_matrix_ptrtype const& massMatrixLumped() const { return M_massMatrixLumped; }
    vector_ptrtype const& vecDiagMassMatrixLumped() const { return M_vecDiagMassMatrixLumped; }


    //-----------------------------------------------------------------------------------//

    solid_1dreduced_ptrtype solid1dReduced() const { return M_solid1dReduced; }

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


    // fsi coupling
    bool useFSISemiImplicitScheme() const { return M_useFSISemiImplicitScheme; }
    void useFSISemiImplicitScheme(bool b) { M_useFSISemiImplicitScheme=b; }
    std::string couplingFSIcondition() const { return M_couplingFSIcondition; }
    void couplingFSIcondition(std::string s) { M_couplingFSIcondition=s; }

    std::set<std::string> const& markerNameFSI() const { return M_bcFSIMarkerManagement.markerFluidStructureInterfaceBC(); }

    void updateNormalStressFromStruct();
    //void updateStressCriterions();

    void updatePreStress() { *U_displ_struct_prestress=*M_fieldDisplacement; }

#if 0 // TODO
    //usefull for 1d reduced model
    void updateInterfaceDispFrom1dDisp();
    void updateInterfaceVelocityFrom1dVelocity();

    element_vect_1dreduced_ptrtype
    extendVelocity1dReducedVectorial( element_1dreduced_type const& vel1d ) const;
#endif
    //-----------------------------------------------------------------------------------//
    // post processing computation
    //-----------------------------------------------------------------------------------//

    double computeExtremumValue( std::string const& field, std::set<std::string> const& markers, std::string const& type ) const;
    double computeVolumeVariation( elements_reference_wrapper_t<mesh_type> const& rangeElt ) const;


    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//

    // assembly methods for linear system
    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const;

    // assembly methods for nonlinear system
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    template <typename ModelContextType>
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    template <typename ModelContextType>
    void updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    template <typename ModelContextType>
    void updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;

private :


    //! solver
    std::string M_solverName;

    // physical parameter
    materialsproperties_ptrtype M_materialsProperties;


    // boundary conditions
    map_vector_field<nDim,1,2> M_bcDirichlet;
    std::map<ComponentType,map_scalar_field<2> > M_bcDirichletComponents;
    map_scalar_field<2> M_bcNeumannScalar,M_bcInterfaceFSI;
    map_vector_field<nDim,1,2> M_bcNeumannVectorial;
    map_matrix_field<nDim,nDim,2> M_bcNeumannTensor2;
    map_vector_fields<nDim,1,2> M_bcRobin;
    map_scalar_field<2> M_bcNeumannEulerianFrameScalar;
    map_vector_field<nDim,1,2> M_bcNeumannEulerianFrameVectorial;
    map_matrix_field<nDim,nDim,2> M_bcNeumannEulerianFrameTensor2;
    map_vector_field<nDim,1,2> M_volumicForcesProperties;
    MarkerManagementDirichletBC M_bcDirichletMarkerManagement;
    MarkerManagementNeumannBC M_bcNeumannMarkerManagement;
    MarkerManagementNeumannEulerianFrameBC M_bcNeumannEulerianFrameMarkerManagement;
    MarkerManagementRobinBC M_bcRobinMarkerManagement;
    MarkerManagementFluidStructureInterfaceBC M_bcFSIMarkerManagement;

    //-------------------------------------------//
    // standard model
    //-------------------------------------------//
    // mesh
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;
    MeshMover<mesh_type> M_meshMover;
    MeshMover<typename mesh_type::trace_mesh_type> M_meshMoverTrace;
    // function space
    space_displacement_ptrtype M_XhDisplacement;
    element_displacement_ptrtype M_fieldDisplacement, M_fieldVelocity, M_fieldAcceleration;
    space_pressure_ptrtype M_XhPressure;
    element_pressure_ptrtype M_fieldPressure;
    space_constraint_vec_ptrtype M_XhConstraintVec;
    element_constraint_vec_ptrtype M_fieldConstraintVec;
    element_vectorial_ptrtype U_displ_struct_prestress;
    // normal stress space
    space_normal_stress_ptrtype M_XhNormalStress;
    //element_normal_stress_ptrtype M_fieldNormalStressFromFluid;
    element_normal_stress_ptrtype M_fieldNormalStressFromStruct;
    // stress tensor space
    space_stress_tensor_ptrtype M_XhStressTensor;
    element_stress_tensor_ptrtype M_fieldStressTensor;
    // time discretisation
    newmark_displacement_ptrtype M_timeStepNewmark;
    savets_pressure_ptrtype M_savetsPressure;

    // mass matrix lumped
    bool M_useMassMatrixLumped;
    sparse_matrix_ptrtype M_massMatrixLumped;
    vector_ptrtype M_vecDiagMassMatrixLumped;

    // exporter
    exporter_ptrtype M_exporter;
    // ho exporter
    bool M_isHOVisu;
#if 1 //defined(FEELPP_HAS_VTK)
    export_ho_ptrtype M_exporter_ho;
    space_vectorial_visu_ho_ptrtype M_XhVectorialVisuHO;
    element_vectorial_visu_ho_ptrtype M_displacementVisuHO;
    op_interpolation_visu_ho_disp_ptrtype M_opIdisplacement;
    op_interpolation_visu_ho_normalstress_ptrtype M_opInormalstress;
    space_scalar_visu_ho_ptrtype M_XhScalarVisuHO;
    element_scalar_visu_ho_ptrtype M_pressureVisuHO;
    op_interpolation_visu_ho_pressure_ptrtype M_opIpressure;
#endif

    std::map<std::string,std::set<std::string> > M_postProcessVolumeVariation; // (name,list of markers)
    //-------------------------------------------//
    // 1d_reduced model
    //-------------------------------------------//
    solid_1dreduced_ptrtype M_solid1dReduced;

    //-------------------------------------------//
    // others
    //-------------------------------------------//

    bool M_timeSteppingUseMixedFormulation;
    std::string M_timeStepping;
    //generalised-alpha parameters
    double M_genAlpha_rho;
    double M_genAlpha_alpha_m;
    double M_genAlpha_alpha_f;
    double M_genAlpha_gamma;
    double M_genAlpha_beta;


    bdf_displacement_ptrtype M_timeStepBdfDisplacement, M_timeStepBdfVelocity;
    double M_timeStepThetaValue;
    vector_ptrtype M_timeStepThetaSchemePreviousContrib;

    // fsi
    bool M_useFSISemiImplicitScheme;
    std::string M_couplingFSIcondition;


}; // SolidMechanics

template< typename ConvexType, typename BasisDisplacementType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
SolidMechanics<ConvexType,BasisDisplacementType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    if ( !this->hasSolidEquationStandard() )
        return;

    this->log("SolidMechanics","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if constexpr ( nOrderGeo <= 2 )
    {
        this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr, exportsExpr );
    }

    this->executePostProcessMeasures( time, mfields, symbolsExpr );
    // TODO
    //this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : M_bdfTemperature->iteration(), mfields );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("SolidMechanics","exportResults", "finish");
}

template< typename ConvexType, typename BasisDisplacementType>
template <typename ModelFieldsType, typename SymbolsExpr>
void
SolidMechanics<ConvexType,BasisDisplacementType>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr )
{
    if ( !this->hasSolidEquationStandard() )
        return;

    //auto const& u = mfields.field( FieldTag::displacement(this), "displacement" );

    // volume variation
    for ( auto const& ppvv : M_postProcessVolumeVariation )
    {
        std::string const& vvname = ppvv.first;
        auto const& vvmarkers = ppvv.second;
        elements_reference_wrapper_t<mesh_type> vvrange = ( vvmarkers.size() == 1 && vvmarkers.begin()->empty() )?
            M_rangeMeshElements : markedelements( this->mesh(),vvmarkers );
        double volVar = this->computeVolumeVariation( vvrange );
        this->postProcessMeasures().setValue( vvname, volVar );
    }

    model_measures_quantities_empty_t mquantities;

    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), M_rangeMeshElements, symbolsExpr, mfields, mquantities );
}


} // FeelModels
} // Feel

#include <feel/feelmodels/solid/solidmechanicsassembly.hpp>


#endif /* FEELPP_TOOLBOXES_SOLIDMECHANICS_HPP */
