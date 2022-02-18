/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-12-12

  Copyright (C) 2016 Feel++ Consortium

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
   \file thermoelectric.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-12-12
 */

#ifndef FEELPP_TOOLBOXES_ELECTRIC_HPP
#define FEELPP_TOOLBOXES_ELECTRIC_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/electric/electricboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{
/** 
 * Toolbox Electric 
 * @ingroup Toolboxes
 */
template< typename ConvexType, typename BasisPotentialType>
class Electric : public ModelNumerical,
                 public ModelPhysics<ConvexType::nDim>,
                 public std::enable_shared_from_this< Electric<ConvexType,BasisPotentialType> >
{
    typedef ModelPhysics<ConvexType::nDim> super_physics_type;
public:
    typedef ModelNumerical super_type;
    typedef Electric<ConvexType,BasisPotentialType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type mesh_electric_type;

    // function space electric-potential
    typedef BasisPotentialType basis_electricpotential_type;
    static const uint16_type nOrderPolyElectricPotential = basis_electricpotential_type::nOrder;
    typedef FunctionSpace<mesh_type, bases<basis_electricpotential_type> > space_electricpotential_type;
    typedef std::shared_ptr<space_electricpotential_type> space_electricpotential_ptrtype;
    typedef typename space_electricpotential_type::element_type element_electricpotential_type;
    typedef std::shared_ptr<element_electricpotential_type> element_electricpotential_ptrtype;
    typedef typename space_electricpotential_type::element_external_storage_type element_electricpotential_external_storage_type;

    // materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    struct FieldTag
    {
        static auto potential( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };

    //___________________________________________________________________________________//
    // constructor
    Electric( std::string const& prefix,
              std::string const& keyword = "electric",
              worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_electricpotential_ptrtype const& spaceElectricPotential() const { return M_XhElectricPotential; }
    element_electricpotential_ptrtype const& fieldElectricPotentialPtr() const { return M_fieldElectricPotential; }
    element_electricpotential_type const& fieldElectricPotential() const { return *M_fieldElectricPotential; }

    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initInitialConditions();
    void initBoundaryConditions();
    void initPostProcess() override;

public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;
    void initAlgebraicFactory();

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    //___________________________________________________________________________________//
    // execute post-processing
    //___________________________________________________________________________________//

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr )
        {
            return this->exportResults( time, this->modelFields(), symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );
        }

    template <typename ModelFieldsType, typename SymbolsExpr>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr );

    //___________________________________________________________________________________//
    // export expressions
    //___________________________________________________________________________________//

    template <typename SymbExprType>
    auto exprPostProcessExportsToolbox( SymbExprType const& se, std::string const& prefix ) const
        {
            auto const& v = this->fieldElectricPotential();
            auto electricFieldExpr = -trans(gradv( v ) );
            std::map<std::string,std::vector<std::tuple<std::decay_t<decltype(electricFieldExpr)>, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprElectricField;
            mapExprElectricField[prefixvm(prefix,"electric-field")].push_back( std::make_tuple( electricFieldExpr, M_rangeMeshElements, "element" ) );

            typedef decltype(expr(typename ModelExpression::expr_scalar_type{},se)) _expr_scalar_type;
            std::map<std::string,std::vector<std::tuple<_expr_scalar_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprScalar;

            typedef std::decay_t<decltype(-_expr_scalar_type{}*trans(gradv(v)))> expr_current_density_type;
            std::map<std::string,std::vector<std::tuple< expr_current_density_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprCurrentDensity;

            typedef std::decay_t<decltype(_expr_scalar_type{}*inner(gradv(v)))> expr_joules_losses_type;
            std::map<std::string,std::vector<std::tuple< expr_joules_losses_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprJoulesLosses;

            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                if ( this->materialsProperties()->hasElectricConductivity( matName ) )
                {
                    auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
                    auto electricConductivityExpr = expr( electricConductivity.expr(), se );

                    auto currentDensityExpr = -electricConductivityExpr*trans(gradv(v)) ;
                    mapExprCurrentDensity[prefixvm(prefix,"current-density")].push_back( std::make_tuple( currentDensityExpr, range, "element" ) );

                    auto joulesLossesExpr = electricConductivityExpr*inner(gradv(v));
                    mapExprJoulesLosses[prefixvm(prefix,"joules-losses")].push_back( std::make_tuple( joulesLossesExpr, range, "element" ) );
                }
            }
            return hana::make_tuple( mapExprElectricField, mapExprCurrentDensity, mapExprJoulesLosses );
        }
        template <typename SymbExprType>
        auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
            {
                return hana::concat( this->materialsProperties()->exprPostProcessExports( this->mesh(),this->physicsAvailable(),se ),
                                     this->exprPostProcessExportsToolbox( se, prefix ) );
            }

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return this->modelFields( this->fieldElectricPotentialPtr(), prefix );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto field_p = this->spaceElectricPotential()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( "potential-electric" ) );
            return this->modelFields( field_p, prefix );
        }
    template <typename PotentialFieldType>
    auto modelFields( PotentialFieldType const& field_p, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( modelField<FieldCtx::FULL>( FieldTag::potential(this), prefix, "electric-potential", field_p, "P", this->keyword() ) );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::potential(this), "electric-potential", startBlockSpaceIndex ) );
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr(); // generate symbols electric_P, electric_grad_P(_x,_y,_z), electric_dn_P
            return Feel::vf::symbolsExpr( seToolbox, seParam, seMat, seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType>
    auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
        {
            auto const& v = mfields.field( FieldTag::potential(this), "electric-potential" );

            // generate symbol electric_matName_current_density
            typedef decltype( this->currentDensityExpr(v,"") ) _expr_currentdensity_type;
            std::vector<std::tuple<std::string,_expr_currentdensity_type,SymbolExprComponentSuffix>> currentDensitySymbs;
            symbol_expression_t<_expr_currentdensity_type> se_currentdensity;
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
            {
                std::string symbolcurrentDensityStr = prefixvm( this->keyword(), (boost::format("%1%_current_density") %matName).str(), "_");
                auto _currentDensityExpr = this->currentDensityExpr( v, matName );
                se_currentdensity.add( symbolcurrentDensityStr, _currentDensityExpr, SymbolExprComponentSuffix( nDim,1 ) );
            }

            return Feel::vf::symbolsExpr( se_currentdensity );
        }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
        {
            return mfields.trialSymbolsExpr( tsmf );
        }

    //___________________________________________________________________________________//
    // model context helper
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
        {
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( mfields, std::move( se ) );
        }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
    auto modelContext( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( rowStartInVector ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }

    //___________________________________________________________________________________//
    // apply assembly and solver
    //___________________________________________________________________________________//

    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mfields ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mfields ) const;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    template <typename ModelContextType>
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mfields ) const;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    template <typename ModelContextType>
    void updateJacobian( DataUpdateJacobian & data, ModelContextType const& mfields ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    template <typename ModelContextType>
    void updateResidual( DataUpdateResidual & data, ModelContextType const& mfields ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;


    //___________________________________________________________________________________//
    // toolbox expressions
    //___________________________________________________________________________________//

    auto electricFieldExpr() const { return this->electricFieldExpr( this->fieldElectricPotential() ); }

    template <typename FieldElectricPotentialType>
    auto electricFieldExpr( FieldElectricPotentialType const& v ) const
        {
            return -trans(gradv( v ) );
        }

    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto currentDensityExpr( std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            return this->currentDensityExpr( this->fieldElectricPotential(), matName, symbolsExpr );
        }
    template <typename FieldElectricPotentialType, typename SymbolsExpr = symbols_expression_empty_t>
    auto currentDensityExpr( FieldElectricPotentialType const& v, std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
            auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
            return -sigmaExpr*trans(gradv(v));
        }

    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto joulesLossesExpr(std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            return this->joulesLossesExpr( this->fieldElectricPotential(), matName, symbolsExpr );
        }
    template <typename FieldElectricPotentialType, typename SymbolsExpr = symbols_expression_empty_t>
    auto joulesLossesExpr( FieldElectricPotentialType const& v, std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
            auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
            return sigmaExpr*inner(gradv(v));
        }

private :

    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_electricpotential_ptrtype M_XhElectricPotential;
    element_electricpotential_ptrtype M_fieldElectricPotential;

    // physical parameter
    materialsproperties_ptrtype M_materialsProperties;

    // boundary conditions
    ElectricBoundaryConditions M_boundaryConditions;
#if 0
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_scalar_field<2> M_volumicForcesProperties;
    MarkerManagementDirichletBC M_bcDirichletMarkerManagement;
    MarkerManagementNeumannBC M_bcNeumannMarkerManagement;
    MarkerManagementRobinBC M_bcRobinMarkerManagement;
#endif
    // post-process
    export_ptrtype M_exporter;
};



template< typename ConvexType, typename BasisPotentialType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
Electric<ConvexType,BasisPotentialType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("Electric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, mfields, symbolsExpr );
    this->executePostProcessSave( invalid_uint32_type_value, mfields );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Electric","exportResults", "finish");
}


template< typename ConvexType, typename BasisPotentialType>
template <typename ModelFieldsType,typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr )
{
    model_measures_quantities_empty_t mquantities;

    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), M_rangeMeshElements, symbolsExpr, mfields, mquantities );
}

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/electric/electricassembly.hpp>

#endif // FEELPP_TOOLBOXES_ELECTRIC_HPP
