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

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisPotentialType>
class Electric : public ModelNumerical,
                 public ModelPhysics<ConvexType::nDim>,
                 public MarkerManagementDirichletBC,
                 public MarkerManagementNeumannBC,
                 public MarkerManagementRobinBC,
                 public std::enable_shared_from_this< Electric<ConvexType,BasisPotentialType> >
{

public:
    typedef ModelNumerical super_type;
    typedef Electric<ConvexType,BasisPotentialType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
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
    typedef MaterialsProperties<mesh_type> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    // measure tools for points evaluation
    typedef MeasurePointsEvaluation<space_electricpotential_type> measure_points_evaluation_type;
    typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    Electric( std::string const& prefix,
              std::string const& keyword = "electric",
              worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"ElectricMesh.path"); }
    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( pt::ptree & p ) override;

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initInitialConditions();
    void initBoundaryConditions();
    void initPostProcess() override;

    template <typename SymbExprType>
    auto symbolsExprFit( SymbExprType const& se ) const { return super_type::symbolsExprFit( se ); }


public :
    void setMesh(mesh_ptrtype const& mesh) { M_mesh = mesh; }
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;
    void initAlgebraicFactory();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr )
        {
            this->exportResults( time, symbolsExpr, hana::concat( this->materialsProperties()->exprPostProcessExports( this->physic(),symbolsExpr ),
                                                                  this->exprPostProcessExports( symbolsExpr ) ) );
        }
    template <typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename TupleFieldsType, typename SymbolsExpr>
    void executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr );

    void updateParameterValues();

    auto allFields( std::string const& prefix = "" ) const
        {
            return hana::make_tuple( std::make_pair( prefixvm(prefix,"electric-potential"),this->fieldElectricPotentialPtr() ) );
        }

    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
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

            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physic() ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( matName );
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


    template <typename FieldElectricPotentialType>
    auto symbolsExpr( FieldElectricPotentialType const& v, std::string const& prefix_symbol = "electric_" ) const
        {
            auto seField = this->symbolsExprField( v, prefix_symbol );
            auto seFit = this->symbolsExprFit( seField );
            auto seMat = this->symbolsExprMaterial( Feel::vf::symbolsExpr( seField, seFit ), prefix_symbol );
            return Feel::vf::symbolsExpr( seField, seFit, seMat );
        }
    auto symbolsExpr( std::string const& prefix_symbol = "electric_" ) const { return this->symbolsExpr( this->fieldElectricPotential(), prefix_symbol ); }

    constexpr auto symbolsExprField( std::string const& prefix_symbol = "electric_" ) const { return this->symbolsExprField( this->fieldElectricPotential(), prefix_symbol ); }
    template <typename FieldElectricPotentialType>
    constexpr auto symbolsExprField( FieldElectricPotentialType const& v,  std::string const& prefix_symbol = "electric_" ) const
        {
            // generate symbols electric_P, electric_grad_P(_x,_y,_z), electric_dn_P
            return Feel::vf::symbolsExpr( symbolExpr( (boost::format("%1%P")%prefix_symbol).str(),idv(v) ),
                                          symbolExpr( (boost::format("%1%grad_P")%prefix_symbol).str(),gradv(v), SymbolExprComponentSuffix( 1,nDim,true ) ),
                                          symbolExpr( (boost::format("%1%dn_P")%prefix_symbol).str(),dnv(v) )
                                          );
        }

    template <typename SymbExprType>
    auto symbolsExprMaterial( SymbExprType const& se, std::string const& prefix_symbol = "electric_" ) const
        {
            return this->materialsProperties()->symbolsExpr( se, prefix_symbol );
        }

    template <typename FieldElectricPotentialType, typename SymbExprType>
    auto symbolsExprPostProcess( FieldElectricPotentialType const& v, SymbExprType const& se, std::string const& prefix_symbol = "electric_" ) const
        {
            typedef decltype( this->currentDensityExpr(v,"",se) ) _expr_currentdensity_type;
            std::vector<std::tuple<std::string,_expr_currentdensity_type,SymbolExprComponentSuffix>> currentDensitySymbs;

            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physic() ) )
            {
                std::string symbolcurrentDensityStr = (boost::format("%1%%2%_current_density")%prefix_symbol %matName).str();
                auto _currentDensityExpr = this->currentDensityExpr( v, matName, se );
                currentDensitySymbs.push_back( std::make_tuple( symbolcurrentDensityStr, _currentDensityExpr, SymbolExprComponentSuffix( nDim,1,true ) ) );
            }
            return Feel::vf::symbolsExpr( symbolExpr( currentDensitySymbs ) );
        }

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_electricpotential_ptrtype const& spaceElectricPotential() const { return M_XhElectricPotential; }
    element_electricpotential_ptrtype const& fieldElectricPotentialPtr() const { return M_fieldElectricPotential; }
    element_electricpotential_type const& fieldElectricPotential() const { return *M_fieldElectricPotential; }

    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename SymbolsExpr>
    void updateLinearPDE( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    template <typename SymbolsExpr>
    void updateJacobian( DataUpdateJacobian & data, SymbolsExpr const& symbolsExpr ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    template <typename SymbolsExpr>
    void updateResidual( DataUpdateResidual & data, SymbolsExpr const& symbolsExpr ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;


    //___________________________________________________________________________________//

    auto electricFieldExpr() const { return this->electricFieldExpr( this->fieldElectricPotential() ); }

    template <typename FieldElectricPotentialType>
    auto electricFieldExpr( FieldElectricPotentialType const& v ) const
        {
            return -trans(gradv( v ) );
        }

    template <typename SymbolsExpr>
    auto currentDensityExpr( std::string const& matName, SymbolsExpr const& symbolsExpr ) const { return this->currentDensityExpr( this->fieldElectricPotential(), matName, symbolsExpr ); }

    template <typename FieldElectricPotentialType, typename SymbolsExpr>
    auto currentDensityExpr( FieldElectricPotentialType const& v, std::string const& matName, SymbolsExpr const& symbolsExpr ) const
        {
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
            auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
            return -sigmaExpr*trans(gradv(v));
        }

    template <typename SymbolsExpr>
    auto joulesLossesExpr(std::string const& matName, SymbolsExpr const& symbolsExpr ) const { return this->joulesLossesExpr( this->fieldElectricPotential(), matName, symbolsExpr ); }

    template <typename FieldElectricPotentialType, typename SymbolsExpr>
    auto joulesLossesExpr( FieldElectricPotentialType const& v, std::string const& matName, SymbolsExpr const& symbolsExpr ) const
        {
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
            auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
            return sigmaExpr*inner(gradv(v));
        }

private :

    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_electricpotential_ptrtype M_XhElectricPotential;
    element_electricpotential_ptrtype M_fieldElectricPotential;

    // physical parameter
    materialsproperties_ptrtype M_materialsProperties;
    // boundary conditions
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_scalar_field<2> M_volumicForcesProperties;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;

    // post-process
    export_ptrtype M_exporter;
    measure_points_evaluation_ptrtype M_measurePointsEvaluation;
};



template< typename ConvexType, typename BasisPotentialType>
template <typename SymbolsExpr, typename ExportsExprType>
void
Electric<ConvexType,BasisPotentialType>::exportResults( double time, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("Electric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    auto fields = this->allFields();
    this->executePostProcessExports( M_exporter, time, fields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, fields, symbolsExpr );
    this->executePostProcessSave( invalid_uint32_type_value, fields );

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
template <typename TupleFieldsType,typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;
    bool hasMeasureNorm = this->executePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, tupleFields, symbolsExpr );
    bool hasMeasureStatistics = this->executePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, tupleFields, symbolsExpr );
    bool hasMeasurePoint = this->executePostProcessMeasuresPoint( M_measurePointsEvaluation, tupleFields );
    if ( hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint )
        hasMeasure = true;

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }
}

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/electric/electricassembly.hpp>

#endif // FEELPP_TOOLBOXES_ELECTRIC_HPP
