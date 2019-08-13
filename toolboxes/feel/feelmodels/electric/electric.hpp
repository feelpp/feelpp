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
    void initBoundaryConditions();
    void initPostProcess();

    constexpr auto symbolsExprField( hana::int_<2> /**/ ) const
        {
            return Feel::vf::symbolsExpr( symbolExpr("electric_P",idv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dxP",dxv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dyP",dyv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dnP",dnv(this->fieldElectricPotential()) )
                                          );
        }
    constexpr auto symbolsExprField( hana::int_<3> /**/ ) const
        {
            return Feel::vf::symbolsExpr( symbolExpr("electric_P",idv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dxP",dxv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dyP",dyv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dzP",dzv(this->fieldElectricPotential()) ),
                                          symbolExpr("electric_dnP",dnv(this->fieldElectricPotential()) )
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

    void exportFields( double time );
    std::set<std::string> postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix = "" ) const;
    bool updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time );
    bool hasPostProcessFieldExported( std::string const& key ) const { return M_postProcessFieldExported.find( key ) != M_postProcessFieldExported.end(); }

    void exportMeasures( double time );
    template <typename SymbolsExpr>
    void exportMeasures( double time, SymbolsExpr const& symbolsExpr );

    void updateParameterValues();

    template <typename SymbExprType>
    /*constexpr*/auto symbolsExpr( SymbExprType const& se ) const
        {
            auto seFit = this->symbolsExprFit( se );
            auto seMat = this->symbolsExprMaterial( Feel::vf::symbolsExpr( se, seFit ) );
            return Feel::vf::symbolsExpr( se, seFit, seMat );
        }
    auto symbolsExpr() const { return this->symbolsExpr( this->symbolsExprField() ); }

    constexpr auto symbolsExprField() const { return this->symbolsExprField( hana::int_<nDim>() ); }

    template <typename SymbExprType>
    auto symbolsExprMaterial( SymbExprType const& se ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},se)) _expr_type;
            std::vector<std::pair<std::string,_expr_type>> matPropSymbs;
            for ( auto const& [_matname, _expr] : this->electricProperties()->electricConductivityByMaterial() )
            {
                matPropSymbs.push_back( std::make_pair( std::string("electric_sigma_"+ _matname), expr( _expr.expr(), se ) ) );
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
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;


    //___________________________________________________________________________________//
    void updateElectricField();
    void updateCurrentDensity();

    template<typename ExprT>
    void updateCurrentDensity( Expr<ExprT> const& expr, elements_reference_wrapper_t<mesh_type> range )
        {
            M_fieldCurrentDensity->on(_range=range, _expr=expr );
        }
private :
    void updateLinearPDEWeakBC( sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart ) const;
    void updateJacobianWeakBC( element_electricpotential_external_storage_type const& v, sparse_matrix_ptrtype& J, bool buildCstPart ) const;
    void updateResidualWeakBC( element_electricpotential_external_storage_type const& v, vector_ptrtype& R, bool buildCstPart ) const;

private :
    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_electricpotential_ptrtype M_XhElectricPotential;
    element_electricpotential_ptrtype M_fieldElectricPotential;
    space_electricfield_ptrtype M_XhElectricField;
    element_electricfield_ptrtype M_fieldElectricField;
    element_electricfield_ptrtype M_fieldCurrentDensity;

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
    std::set<std::string> M_postProcessFieldExported;
    std::set<std::string> M_postProcessUserFieldExported;
    measure_points_evaluation_ptrtype M_measurePointsEvaluation;
};



template< typename ConvexType, typename BasisPotentialType>
template <typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::exportResults( double time, SymbolsExpr const& symbolsExpr )
{
    this->log("Electric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->exportFields( time );

    this->exportMeasures( time, symbolsExpr );

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
template <typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::exportMeasures( double time, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    auto fieldTuple = hana::make_tuple( std::make_pair( "electric-potential",this->fieldElectricPotentialPtr() ),
                                        std::make_pair( "electric-field",this->fieldElectricFieldPtr() ) );
    for ( auto const& ppNorm : this->modelProperties().postProcess().measuresNorm( this->keyword() ) )
    {
        std::map<std::string,double> resPpNorms;
        measureNormEvaluation( this->mesh(), M_rangeMeshElements, ppNorm, resPpNorms, symbolsExpr, fieldTuple );
        for ( auto const& resPpNorm : resPpNorms )
        {
            this->postProcessMeasuresIO().setMeasure( resPpNorm.first, resPpNorm.second );
            hasMeasure = true;
        }
    }
    for ( auto const& ppStat : this->modelProperties().postProcess().measuresStatistics( this->keyword() ) )
    {
        std::map<std::string,double> resPpStats;
        measureStatisticsEvaluation( this->mesh(), M_rangeMeshElements, ppStat, resPpStats, symbolsExpr, fieldTuple );
        for ( auto const& resPpStat : resPpStats )
        {
            this->postProcessMeasuresIO().setMeasure( resPpStat.first, resPpStat.second );
            hasMeasure = true;
        }
    }

    std::map<std::string,double> resPpPoints;
    M_measurePointsEvaluation->eval( this->modelProperties().postProcess().measuresPoint( this->keyword() ), resPpPoints, fieldTuple );
    for ( auto const& resPpPoint : resPpPoints )
    {
        this->postProcessMeasuresIO().setMeasure( resPpPoint.first, resPpPoint.second );
        hasMeasure = true;
    }

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

#endif // FEELPP_TOOLBOXES_ELECTRIC_HPP
