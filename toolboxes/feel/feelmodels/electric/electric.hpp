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
#include <feel/feelmodels/electric/electricpropertiesdescription.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisPotentialType>
class Electric : public ModelNumerical,
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
    // function space electric-field
    typedef Lagrange<nOrderPolyElectricPotential, Vectorial,Discontinuous/*Continuous*/,PointSetFekete> basis_electricfield_type;
    typedef FunctionSpace<mesh_electric_type, bases<basis_electricfield_type> > space_electricfield_type;
    typedef std::shared_ptr<space_electricfield_type> space_electricfield_ptrtype;
    typedef typename space_electricfield_type::element_type element_electricfield_type;
    typedef std::shared_ptr<element_electricfield_type> element_electricfield_ptrtype;

    typedef typename space_electricfield_type::component_functionspace_type space_component_electricfield_type;
    typedef std::shared_ptr<space_component_electricfield_type> space_component_electricfield_ptrtype;
    typedef typename space_component_electricfield_type::element_type element_component_electricfield_type;
    typedef std::shared_ptr<element_component_electricfield_type> element_component_electricfield_ptrtype;

    // mechanical properties desc
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef std::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef ElectricPropertiesDescription<space_scalar_P0_type> electricproperties_type;
    typedef std::shared_ptr<electricproperties_type> electricproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    // measure tools for points evaluation
    typedef MeasurePointsEvaluation<space_electricpotential_type,space_electricfield_type> measure_points_evaluation_type;
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

    template <typename FieldElectricPotentialType>
    constexpr auto symbolsExprField( FieldElectricPotentialType const& v, hana::int_<2> /**/ ) const
        {
            return Feel::vf::symbolsExpr( symbolExpr("electric_P",idv(v) ),
                                          symbolExpr("electric_dxP",dxv(v) ),
                                          symbolExpr("electric_dyP",dyv(v) ),
                                          symbolExpr("electric_dnP",dnv(v) )
                                          );
        }
    template <typename FieldElectricPotentialType>
    constexpr auto symbolsExprField( FieldElectricPotentialType const& v, hana::int_<3> /**/ ) const
        {
            return Feel::vf::symbolsExpr( symbolExpr("electric_P",idv(v) ),
                                          symbolExpr("electric_dxP",dxv(v) ),
                                          symbolExpr("electric_dyP",dyv(v) ),
                                          symbolExpr("electric_dzP",dzv(v) ),
                                          symbolExpr("electric_dnP",dnv(v) )
                                          );
        }
    //auto symbolsExprFit() const { return symbolsExprFit( this->symbolsExprField() ); }

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
    void exportResults( double time, SymbolsExpr const& symbolsExpr );
    template <typename TupleFieldsType, typename SymbolsExpr>
    void executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr );

    void updateParameterValues();

    auto allFields( std::string const& prefix = "" ) const
        {
            return hana::make_tuple( std::make_pair( prefixvm(prefix,"electric-potential"),this->fieldElectricPotentialPtr() ),
                                     std::make_pair( prefixvm(prefix,"electric-field"),this->fieldElectricFieldPtr() ),
                                     //std::make_pair( prefixvm(prefix,"electric-conductivity"),M_electricProperties->fieldElectricConductivityPtr() ),
                                     //std::make_pair( prefixvm(prefix,"current-density"),this->fieldCurrentDensityPtr() ),
                                     std::make_pair( prefixvm(prefix,"joules-losses"),this->fieldJoulesLossesPtr() )
                                     );
        }

    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            typedef decltype(expr(typename ModelExpression::expr_scalar_type{},se)) _expr_scalar_type;
            std::map<std::string,std::vector<std::tuple<_expr_scalar_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprScalar;

            auto const& v = this->fieldElectricPotential();
            typedef std::decay_t<decltype(-_expr_scalar_type{}*trans(gradv(v)))> expr_current_density_type;
            std::map<std::string,std::vector<std::tuple< expr_current_density_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprCurrentDensity;
            for ( auto const& rangeData : this->electricProperties()->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& range = rangeData.second;
                if ( this->electricProperties()->hasElectricConductivity( _matName ) )
                {
                    auto const& electricConductivity = this->electricProperties()->electricConductivity( _matName );
                    auto electricConductivityExpr = expr( electricConductivity.expr(), se );
                    mapExprScalar[prefixvm(prefix,"electric-conductivity")].push_back( std::make_tuple( electricConductivityExpr, range, "element" ) );

                    auto currentDensityExpr = -electricConductivityExpr*trans(gradv(v)) ;
                    mapExprCurrentDensity[prefixvm(prefix,"current-density")].push_back( std::make_tuple( currentDensityExpr, range, "element" ) );
                }
            }
            return hana::make_tuple( mapExprScalar, mapExprCurrentDensity );
        }


    template <typename FieldElectricPotentialType>
    /*constexpr*/auto symbolsExpr( FieldElectricPotentialType const& v ) const
        {
            auto seField = this->symbolsExprField( v );
            auto seFit = this->symbolsExprFit( seField );
            auto seMat = this->symbolsExprMaterial( Feel::vf::symbolsExpr( seField, seFit ) );
            return Feel::vf::symbolsExpr( seField, seFit, seMat );
        }
    auto symbolsExpr() const { return this->symbolsExpr( this->fieldElectricPotential() ); }

    constexpr auto symbolsExprField() const { return this->symbolsExprField( this->fieldElectricPotential() ); }
    template <typename FieldElectricPotentialType>
    constexpr auto symbolsExprField( FieldElectricPotentialType const& v ) const { return this->symbolsExprField( v, hana::int_<nDim>() ); }

    template <typename SymbExprType>
    auto symbolsExprMaterial( SymbExprType const& se ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},se)) _expr_type;
            std::vector<std::pair<std::string,_expr_type>> matPropSymbs;
            for ( auto const& [_matname, _expr] : this->electricProperties()->electricConductivityByMaterial() )
            {
                matPropSymbs.push_back( std::make_pair( (boost::format("electric_%1%_sigma")%_matname).str(), expr( _expr.expr(), se ) ) );
            }
            return Feel::vf::symbolsExpr( symbolExpr( matPropSymbs ) );
        }

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_electricpotential_ptrtype const& spaceElectricPotential() const { return M_XhElectricPotential; }
    element_electricpotential_ptrtype const& fieldElectricPotentialPtr() const { return M_fieldElectricPotential; }
    element_electricpotential_type const& fieldElectricPotential() const { return *M_fieldElectricPotential; }

    space_electricfield_ptrtype const& spaceElectricField() const { return M_XhElectricField; }
    element_electricfield_ptrtype const& fieldElectricFieldPtr() const { return M_fieldElectricField; }
    element_electricfield_type const& fieldElectricField() const { return *M_fieldElectricField; }
    element_electricfield_ptrtype const& fieldCurrentDensityPtr() const { return M_fieldCurrentDensity; }
    element_electricfield_type const& fieldCurrentDensity() const { return *M_fieldCurrentDensity; }
    element_component_electricfield_type const& fieldJoulesLosses() const { return *M_fieldJoulesLosses; }
    element_component_electricfield_ptrtype const& fieldJoulesLossesPtr() const { return M_fieldJoulesLosses; }

    electricproperties_ptrtype const& electricProperties() const { return M_electricProperties; }

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

    template <typename SymbolsExpr>
    void updateFields( SymbolsExpr const& symbolsExpr )
        {
            //this->electricProperties()->updateFields( symbolsExpr );
            this->updateElectricField();
            this->updateCurrentDensity( symbolsExpr );
            this->updateJoulesLosses( symbolsExpr );
        }
private :
    void updateElectricField();
    void updateCurrentDensity();
    template<typename SymbolsExpr>
    void updateCurrentDensity( SymbolsExpr const& symbolsExpr )
        {
            auto const& v = this->fieldElectricPotential();
            for ( auto const& rangeData : this->electricProperties()->rangeMeshElementsByMaterial() )
            {
                std::string const& matName = rangeData.first;
                auto const& range = rangeData.second;
                auto const& electricConductivity = this->electricProperties()->electricConductivity( matName );
                auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
                M_fieldCurrentDensity->on(_range=range, _expr=-sigmaExpr*trans(gradv(v)) );
            }
        }
    template<typename SymbolsExpr>
    void updateJoulesLosses( SymbolsExpr const& symbolsExpr )
        {
            auto const& v = this->fieldElectricPotential();
            for ( auto const& rangeData : this->electricProperties()->rangeMeshElementsByMaterial() )
            {
                std::string const& matName = rangeData.first;
                auto const& range = rangeData.second;
                auto const& electricConductivity = this->electricProperties()->electricConductivity( matName );
                auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
                M_fieldJoulesLosses->on(_range=range, _expr=sigmaExpr*inner(gradv(v)) );
            }
        }

private :
    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_electricpotential_ptrtype M_XhElectricPotential;
    element_electricpotential_ptrtype M_fieldElectricPotential;
    space_electricfield_ptrtype M_XhElectricField;
    element_electricfield_ptrtype M_fieldElectricField;
    element_electricfield_ptrtype M_fieldCurrentDensity;
    element_component_electricfield_ptrtype M_fieldJoulesLosses;

    // physical parameter
    electricproperties_ptrtype M_electricProperties;
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
template <typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::exportResults( double time, SymbolsExpr const& symbolsExpr )
{
    this->log("Electric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    auto fields = this->allFields();
    this->executePostProcessExports( M_exporter, time, fields, symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );
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
