/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_THERMALPROPERTIES_DESCRIPTION_H
#define FEELPP_TOOLBOXES_THERMALPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelvf/exprselectorbymeshelement.hpp>

namespace Feel
{
namespace FeelModels
{

template<class MeshType>
class ThermalPropertiesDescription
{
    typedef ThermalPropertiesDescription<MeshType> self_type;
public :
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    static const uint16_type nDim = mesh_type::nDim;

    ThermalPropertiesDescription( std::string const& prefix, std::string const& exprRepository )
        :
        M_exprRepository( exprRepository ),
        M_thermalConductivityDefaultValue( doption(_name="thermal-conductivity",_prefix=prefix) ),// [ W/(m*K) ]
        M_heatCapacityDefaultValue( doption(_name="heat-capacity",_prefix=prefix) ), // [ J/(kg*K) ]
        M_rhoDefaultValue( doption(_name="rho",_prefix=prefix) ),
        M_thermalExpansionDefaultValue( doption(_name="thermal-expansion",_prefix=prefix) ) // [ 1/K ]
        {}

    ThermalPropertiesDescription( ThermalPropertiesDescription const& ) = default;

    void updateForUse( mesh_ptrtype const& mesh , ModelMaterials const& mats )
        {
            std::set<std::string> eltMarkersInMesh;
            for (auto const& markPair : mesh->markerNames() )
            {
                std::string meshMarker = markPair.first;
                if ( mesh->hasElementMarker( meshMarker ) )
                    eltMarkersInMesh.insert( meshMarker );
            }

            std::map<std::string,std::set<std::string>> markersByMaterial;
            M_markers.clear();
            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
                auto const& mat = m.second;
                if ( mat.hasPhysics() && !mat.hasPhysics( { "heat","aerothermal","thermo-electric" } ) )
                    continue;

                for ( std::string const& matmarker : mat.meshMarkers() )
                {
                    if ( eltMarkersInMesh.find( matmarker ) == eltMarkersInMesh.end() )
                        continue;
                    M_markers.insert( matmarker );
                    markersByMaterial[matName].insert( matmarker );
                }
            }

            M_isDefinedOnWholeMesh = ( M_markers.size() == eltMarkersInMesh.size() );

            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
                auto const& mat = m.second;
                auto itFindMat = markersByMaterial.find( matName );
                if ( itFindMat == markersByMaterial.end() )
                    continue;
                if ( itFindMat->second.empty() )
                    continue;
                auto const& matmarkers = itFindMat->second;
                auto range = markedelements( mesh,matmarkers );
                M_rangeMeshElementsByMaterial[matName] = range;

                M_thermalConductivityByMaterial[matName];
                if ( mat.hasProperty("k") )
                {
                    M_thermalConductivityByMaterial[matName] = mat.property( "k");
                    auto & prop = M_thermalConductivityByMaterial[matName];
                    if ( prop.hasAtLeastOneExpr() && ( !prop.hasExprScalar() && !prop.template hasExprMatrix<nDim,nDim>() )  )
                        CHECK( false ) << "thermal conductivty should be a scalar or a matrix";
                }

                M_rhoByMaterial[matName];
                if ( mat.hasPropertyExprScalar("rho") )
                {
                    auto const& expr = mat.propertyExprScalar("rho");
                    M_rhoByMaterial[matName].setExpr( expr );
                }
                else M_rhoByMaterial[matName].setExpr( expr("0") );

                M_heatCapacityByMaterial[matName];
                if ( mat.hasPropertyExprScalar("Cp") )
                {
                    auto const& expr = mat.propertyExprScalar("Cp");
                    M_heatCapacityByMaterial[matName].setExpr( expr );
                }
                else M_heatCapacityByMaterial[matName].setExpr( expr("0") );

                M_thermalExpansionByMaterial[matName];
                if ( mat.hasPropertyExprScalar("beta") )
                {
                    auto const& expr = mat.propertyExprScalar("beta");
                    M_thermalExpansionByMaterial[matName].setExpr( expr );
                }
                else M_thermalExpansionByMaterial[matName].setExpr( expr("0") );

                // rho * Cp
                M_rhoHeatCapacityByMaterial[matName];
                if ( this->rho( matName ).hasExpr()  && this->heatCapacity( matName ).hasExpr() )
                {
                    auto rhoExpr = M_rhoByMaterial[matName].expr();
                    auto expr = expr_mult<2>( rhoExpr,M_heatCapacityByMaterial[matName].expr(),"",mesh->worldComm(),M_exprRepository );
                    M_rhoHeatCapacityByMaterial[matName].setExpr( expr );
                }
#if 0
                if ( M_rhoByMaterial[matName].isConstant() )
                {
                    double rhoValue = M_rhoByMaterial[matName].value();
                    if ( M_heatCapacityByMaterial[matName].isConstant() )
                        M_rhoHeatCapacityByMaterial[matName].setValue( rhoValue*M_heatCapacityByMaterial[matName].value() );
                    else
                    {
                        auto expr = expr_mult( M_heatCapacityByMaterial[matName].expr(),rhoValue,"",mesh->worldComm(),M_exprRepository );
                        M_rhoHeatCapacityByMaterial[matName].setExpr( expr );
                    }
                }
                else
                {
                    auto rhoExpr = M_rhoByMaterial[matName].expr();
                    if ( M_heatCapacityByMaterial[matName].isConstant() )
                    {
                        auto expr = expr_mult( rhoExpr,M_heatCapacityByMaterial[matName].value(),"",mesh->worldComm(),M_exprRepository );
                        M_rhoHeatCapacityByMaterial[matName].setExpr( expr );
                    }
                    else
                    {
                        auto expr = expr_mult<2>( rhoExpr,M_heatCapacityByMaterial[matName].expr(),"",mesh->worldComm(),M_exprRepository );
                        M_rhoHeatCapacityByMaterial[matName].setExpr( expr );
                    }
                }
#endif
            }

            M_exprSelectorByMeshElementMapping = std::make_shared<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>>();
            M_exprSelectorByMeshElementMapping->template updateForUse<mesh_type>( this->rangeMeshElementsByMaterial() );
        }

    std::set<std::string> const& markers() const { return M_markers; }

    bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElementsByMaterial( std::string const& matName ) const
        {
            CHECK( this->hasMaterial(matName) ) << "no material with name " << matName;
            return M_rangeMeshElementsByMaterial.find( matName )->second;
        }

    bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

    std::shared_ptr<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>> exprSelectorByMeshElementMapping() const { return M_exprSelectorByMeshElementMapping; }

    // thermal conductivity
    bool hasThermalConductivity( std::string const& matName ) const
        {
            return M_thermalConductivityByMaterial.find( matName ) != M_thermalConductivityByMaterial.end();
        }
    ModelExpression const& thermalConductivity( std::string const& matName ) const
        {
            CHECK( this->hasThermalConductivity( matName ) ) << "material name not registered : " << matName;
            return M_thermalConductivityByMaterial.find( matName )->second;
        }

    bool allThermalConductivitiesAreScalar() const
        {
            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& thermalConductivity = this->thermalConductivity( _matName );
                if ( !thermalConductivity.isScalar() )
                    return false;
            }
            return true;
        }

    template <typename SymbolsExpr>
    auto thermalConductivityExpr( SymbolsExpr const& symbolsExpr ) const
        {
            auto Id = eye<nDim,nDim>();
            typedef decltype(expr(scalar_field_expression<2>{},symbolsExpr)*Id) _expr_scalar_type;
            std::vector<std::pair<std::string,_expr_scalar_type>> exprs_k;
            typedef decltype(expr(matrix_field_expression<nDim,nDim,2>{},symbolsExpr)) _expr_matrix_type;
            std::vector<std::pair<std::string,_expr_matrix_type>> exprs_k_matrix;

            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& thermalConductivity = this->thermalConductivity( _matName );
                if ( thermalConductivity.isMatrix() )
                {
                    auto thermalConductivityExpr = expr( thermalConductivity.template expr<nDim,nDim>(), symbolsExpr );
                    exprs_k_matrix.push_back( std::make_pair( _matName, thermalConductivityExpr ) );
                }
                else
                {
                    auto thermalConductivityExpr = expr( thermalConductivity.exprScalar(), symbolsExpr );
                    exprs_k.push_back( std::make_pair( _matName, thermalConductivityExpr*Id ) );
                }
            }
            auto kappa = expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, exprs_k, exprs_k_matrix );
            return kappa;
        }

    template <typename SymbolsExpr>
    auto thermalConductivityScalarExpr( SymbolsExpr const& symbolsExpr ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},symbolsExpr)) _expr_scalar_type;
            std::vector<std::pair<std::string,_expr_scalar_type>> exprs_k;
            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& thermalConductivity = this->thermalConductivity( _matName );
                if ( thermalConductivity.isScalar() )
                {
                    auto thermalConductivityExpr = expr( thermalConductivity.exprScalar(), symbolsExpr );
                    exprs_k.push_back( std::make_pair( _matName, thermalConductivityExpr ) );
                }
            }
            auto kappa = expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, exprs_k );
            return kappa;
        }



    // rho
    bool hasRho( std::string const& matName ) const
        {
            return M_rhoByMaterial.find( matName ) != M_rhoByMaterial.end();
        }
    ModelExpressionScalar const& rho( std::string const& matName ) const
        {
            CHECK( this->hasRho( matName ) ) << "material name not registered : " << matName;
            return M_rhoByMaterial.find( matName )->second;
        }
    // heat capacity
    bool hasHeatCapacity( std::string const& matName ) const
        {
            return M_heatCapacityByMaterial.find( matName ) != M_heatCapacityByMaterial.end();
        }
    ModelExpressionScalar const& heatCapacity( std::string const& matName ) const
        {
            CHECK( this->hasHeatCapacity( matName ) ) << "material name not registered : " << matName;
            return M_heatCapacityByMaterial.find( matName )->second;
        }
    // thermal expansion
    bool hasThermalExpansion( std::string const& matName ) const
        {
            return M_thermalExpansionByMaterial.find( matName ) != M_thermalExpansionByMaterial.end();
        }
    ModelExpressionScalar const& thermalExpansion( std::string const& matName ) const
        {
            CHECK( this->hasThermalExpansion( matName ) ) << "material name not registered : " << matName;
            return M_thermalExpansionByMaterial.find( matName )->second;
        }
    // rho * Cp
    ModelExpressionScalar const& rhoHeatCapacity( std::string const& matName ) const
        {
            CHECK( this->hasMaterial( matName ) ) << "material name not registered : " << matName;
            return M_rhoHeatCapacityByMaterial.find( matName )->second;
        }


    std::shared_ptr<std::ostringstream>
    getInfoMaterialParameters() const
        {
            std::shared_ptr<std::ostringstream> ostr( new std::ostringstream() );
            *ostr << "\n   Materials parameters";
            *ostr << "\n     -- defined on whole mesh : " << this->isDefinedOnWholeMesh();
            *ostr << "\n     -- number of materials : " << M_rangeMeshElementsByMaterial.size();
            for ( auto const& matRange : M_rangeMeshElementsByMaterial)
            {
                std::string const& matName = matRange.first;
                if ( this->rho( matName ).hasExpr() )
                {
                    *ostr << "\n     -- [" << matName << "] rho : ";
                    if ( this->rho(matName).isConstant() )
                        *ostr << this->rho(matName).value();
                    else
                        *ostr << str( this->rho(matName).expr().expression() );
                }
                if ( this->thermalConductivity( matName ).hasExprScalar() || this->thermalConductivity( matName ).template hasExpr<nDim,nDim>() )
                {
                    *ostr << "\n     -- [" << matName << "] thermal conductivity : ";
                    auto const& thermalConductivity = this->thermalConductivity(matName);
                    if ( thermalConductivity.isMatrix() )
                        *ostr << str( thermalConductivity.template exprMatrix<nDim,nDim>().expression() );
                    else if ( thermalConductivity.isConstant() )
                        *ostr << thermalConductivity.value();
                    else
                        *ostr << str( thermalConductivity.expr().expression() );
                }
                if ( this->heatCapacity( matName ).hasExpr() )
                {
                    *ostr << "\n     -- [" << matName << "] heat capacity : ";
                    if ( this->heatCapacity(matName).isConstant() )
                        *ostr << this->heatCapacity(matName).value();
                    else
                        *ostr << str( this->heatCapacity(matName).expr().expression() );
                }
                if ( this->thermalExpansion( matName ).hasExpr() )
                {
                    *ostr << "\n     -- [" << matName << "] thermal expansion : ";
                    if ( this->thermalExpansion(matName).isConstant() )
                        *ostr << this->thermalExpansion(matName).value();
                    else
                        *ostr << str( this->thermalExpansion(matName).expr().expression() );
                }
            }
            return ostr;
        }

    void updateInformationObject( pt::ptree & p )
        {
            p.put( "number of materials", M_rangeMeshElementsByMaterial.size() );
            for ( auto const& matRange : M_rangeMeshElementsByMaterial)
            {
                pt::ptree matPt;
                std::string const& matName = matRange.first;
                if ( this->rho( matName ).hasExpr() )
                {
                    if ( this->rho(matName).isConstant() )
                        matPt.put( "rho", this->rho(matName).value() );
                    else
                        matPt.put( "rho", str( this->rho(matName).expr().expression() ) );
                }
                if ( this->thermalConductivity( matName ).hasExprScalar() || this->thermalConductivity( matName ).template hasExpr<nDim,nDim>() )
                {
                    auto const& thermalConductivity = this->thermalConductivity(matName);
                    if ( thermalConductivity.isMatrix() )
                        matPt.put( "thermal conductivity", str( thermalConductivity.template exprMatrix<nDim,nDim>().expression() ) );
                    else if ( thermalConductivity.isConstant() )
                        matPt.put( "thermal conductivity", thermalConductivity.value() );
                    else
                        matPt.put( "thermal conductivity", str( thermalConductivity.expr().expression() ) );
                }
                if ( this->heatCapacity( matName ).hasExpr() )
                {
                    if ( this->heatCapacity(matName).isConstant() )
                        matPt.put( "heat capacity", this->heatCapacity(matName).value() );
                    else
                        matPt.put( "heat capacity", str( this->heatCapacity(matName).expr().expression() ) );
                }
                if ( this->thermalExpansion( matName ).hasExpr() )
                {
                    if ( this->thermalExpansion(matName).isConstant() )
                        matPt.put( "thermal expansion", this->thermalExpansion(matName).value() );
                    else
                        matPt.put( "thermal expansion", str( this->thermalExpansion(matName).expr().expression() ) );
                }
                p.add_child( matName, matPt );
            }
        }

    bool hasThermalConductivityDependingOnSymbol( std::string const& symbolStr ) const
        {
            for ( auto const& conductivityData : M_thermalConductivityByMaterial )
            {
                auto const& thermalConductivity = conductivityData.second;
                if ( thermalConductivity.isMatrix() )
                {
                    if ( thermalConductivity.template exprMatrix<nDim,nDim>().expression().hasSymbol( symbolStr ) )
                        return true;
                    continue;
                }
                if ( thermalConductivity.isConstant() )
                    continue;
                if ( thermalConductivity.expr().expression().hasSymbol( symbolStr ) )
                    return true;
            }
            return false;
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & prop : M_thermalConductivityByMaterial )
                prop.second.setParameterValues( mp );
            for ( auto & prop : M_heatCapacityByMaterial )
                prop.second.setParameterValues( mp );
            for ( auto & prop : M_rhoByMaterial )
                prop.second.setParameterValues( mp );
            for ( auto & prop : M_thermalExpansionByMaterial )
                prop.second.setParameterValues( mp );
        }

    template <typename SymbExprType>
    auto symbolsExpr( SymbExprType const& se, std::string const& prefix_symbol ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},se)) _expr_scalar_type;
            std::vector<std::pair<std::string,_expr_scalar_type>> matPropSymbsScalar;
            typedef decltype(expr(matrix_field_expression<nDim,nDim,2>{},se)) _expr_matrix_type;
            std::vector<std::tuple<std::string,_expr_matrix_type,SymbolExprComponentSuffix>> matPropSymbsMatrix;

            // generate symbols heat_matName_k or heat_matName_k(_xx,_xy,...,_zz)
            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& thermalConductivity = this->thermalConductivity( _matName );
                if ( thermalConductivity.isMatrix() )
                {
                    matPropSymbsMatrix.push_back( std::make_tuple( (boost::format("%1%%2%_k")%prefix_symbol %_matName).str(), expr( thermalConductivity.template exprMatrix<nDim,nDim>(), se ), SymbolExprComponentSuffix( nDim,nDim,true )  ) );
                }
                else
                    matPropSymbsScalar.push_back( std::make_pair( (boost::format("%1%%2%_k")%prefix_symbol %_matName).str(), expr( thermalConductivity.exprScalar(), se ) ) );
            }

            typedef decltype( this->thermalConductivityScalarExpr( se ) ) _expr_scalar_selector_type;
            std::vector<std::pair<std::string,_expr_scalar_selector_type>> matPropSymbsScalarSelector;
            typedef decltype( this->thermalConductivityExpr( se ) ) _expr_conductivity_matrix_selector_type;
            std::vector<std::tuple<std::string,_expr_conductivity_matrix_selector_type,SymbolExprComponentSuffix>> matPropSymbsConductivityMatrixSelector;
            // generate the symbol heat_k if the all conductivities are scalar or heat_k(_xx,_xy,...) if at least one conductivity is a matrix
            if ( this->allThermalConductivitiesAreScalar() )
            {
                auto expr_k = this->thermalConductivityScalarExpr( se );
                matPropSymbsScalarSelector.push_back( std::make_pair( (boost::format("%1%k")% prefix_symbol).str(), expr_k ) );
            }
            else
            {
                auto expr_k = this->thermalConductivityExpr( se );
                matPropSymbsConductivityMatrixSelector.push_back( std::make_tuple( (boost::format("%1%k")% prefix_symbol).str(), expr_k,  SymbolExprComponentSuffix( nDim,nDim,true ) ) );
            }
            return Feel::vf::symbolsExpr( symbolExpr( matPropSymbsScalar ), symbolExpr( matPropSymbsMatrix ), symbolExpr(matPropSymbsScalarSelector), symbolExpr( matPropSymbsConductivityMatrixSelector ) );
        }

private :
    std::string M_exprRepository;
    std::set<std::string> M_markers;
    bool M_isDefinedOnWholeMesh;
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;
    std::shared_ptr<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>> M_exprSelectorByMeshElementMapping;

    std::map<std::string, ModelExpressionScalar> /*M_thermalConductivityByMaterial,*/ M_heatCapacityByMaterial, M_rhoByMaterial, M_thermalExpansionByMaterial, M_rhoHeatCapacityByMaterial;
    std::map<std::string, ModelExpression> M_thermalConductivityByMaterial;
    double M_thermalConductivityDefaultValue, M_heatCapacityDefaultValue, M_rhoDefaultValue, M_thermalExpansionDefaultValue;
};


} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
