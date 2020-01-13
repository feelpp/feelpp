/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELMATERIALS_H
#define FEELPP_TOOLBOXES_MODELMATERIALS_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelvf/exprselectorbymeshelement.hpp>
#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

class MaterialProperty : public ModelExpression
{
public :
    MaterialProperty( std::string const& symbol, ModelExpression const& expr )
        :
        ModelExpression( expr ),
        M_symbol( symbol )
        //M_expr( expr )
        {}
    MaterialProperty( MaterialProperty const& ) = default;
    MaterialProperty( MaterialProperty && ) = default;

    //ModelExpression const& expr() const { return M_expr; }
private :
    std::string M_name;
    std::string M_symbol;
    //ModelExpression M_expr;
    //std::set<std::string> M_shapes;
};


template<class MeshType>
class MaterialProperties : public std::map<std::string, MaterialProperty>
{
public :
    MaterialProperties( std::string const& name ) : M_materialName( name ) {}
    MaterialProperties( MaterialProperties const& ) = default;
    MaterialProperties( MaterialProperties && ) = default;

    void add( std::string const& propName, ModelExpression const& expr )
        {
            this->/*M_propertyToExpression.*/emplace(propName, MaterialProperty(propName, expr) );
        }

    bool has( std::string const& propName ) const
        {
            if ( this->/*M_propertyToExpression.*/find( propName ) != this->/*M_propertyToExpression.*/end() )
                return true;
            return false;
        }

    MaterialProperty const&
    property( std::string const& propName ) const
        {
            auto itFind = this->/*M_propertyToExpression.*/find( propName );
            CHECK( itFind != this->/*M_propertyToExpression.*/end() ) << "material property" << propName << " not found";
            return itFind->second;
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & [propName,matProp] : *this )
                matProp.setParameterValues( mp );
        }

    // ModelExpression const&
    // expr( std::string const& propName ) const
    //     {
    //         return this->property( propName ).expr();
    //     }

private :
    std::string M_materialName;
    //std::map<std::string, MaterialProperty> M_propertyToExpression;
};

template<class MeshType>
class MaterialsProperties
{
    typedef MaterialsProperties<MeshType> self_type;
public :
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    static const uint16_type nDim = mesh_type::nDim;

    MaterialsProperties( std::string const& prefix, std::string const& exprRepository )
        :
        M_exprRepository( exprRepository )
        //M_thermalConductivityDefaultValue( doption(_name="thermal-conductivity",_prefix=prefix) ),// [ W/(m*K) ]
        //M_heatCapacityDefaultValue( doption(_name="heat-capacity",_prefix=prefix) ), // [ J/(kg*K) ]
        //M_rhoDefaultValue( doption(_name="rho",_prefix=prefix) ),
        //M_thermalExpansionDefaultValue( doption(_name="thermal-expansion",_prefix=prefix) ) // [ 1/K ]
        {}

    MaterialsProperties( MaterialsProperties const& ) = default;

    void updateForUse( mesh_ptrtype const& mesh, ModelMaterials const& mats, ModelPhysics const& modelphysics )
        {
            std::set<std::string> eltMarkersInMesh;
            for (auto const& markPair : mesh->markerNames() )
            {
                std::string meshMarker = markPair.first;
                if ( mesh->hasElementMarker( meshMarker ) )
                    eltMarkersInMesh.insert( meshMarker );
            }

            auto const& mapPhysicsToSubphysics = modelphysics.mapPhysicsToSubphysics();
            auto const& defaultPhysics = modelphysics.physic();
            auto const& physicsAvailable = modelphysics.physics();

            std::map<std::string,std::set<std::string>> markersByMaterial;
            M_markers.clear();
            std::set<std::string> currentPhysics;

            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
                auto const& mat = m.second;
#if 0
                if ( mat.hasPhysics() && !mat.hasPhysics( { "heat","aerothermal","thermo-electric" } ) )
                    continue;

                if (  mat.hasPhysics() && !mat.hasPhysics( physicsAvailable ) )
                    continue;
#endif
                std::set<std::string> currentPhysicsReaded;
                if ( mat.hasPhysics() )
                {
                    if ( !mat.hasPhysics( physicsAvailable ) )
                        continue;
                    for ( std::string const& p : mat.physics() )
                        if ( physicsAvailable.find( p ) != physicsAvailable.end() )
                            currentPhysicsReaded.insert( p );
                }
                else
                {
                    if ( physicsAvailable.find( defaultPhysics ) != physicsAvailable.end() )
                        currentPhysicsReaded.insert( defaultPhysics );
                }


                for ( std::string const& p : currentPhysicsReaded )
                {
                    currentPhysics.insert( p );
                    auto itFindPhysics = mapPhysicsToSubphysics.find( p );
                    if ( itFindPhysics != mapPhysicsToSubphysics.end() )
                        currentPhysics.insert( itFindPhysics->second.begin(), itFindPhysics->second.end() );
                }

                if ( currentPhysics.empty() )
                    continue;

                for ( std::string const& matmarker : mat.meshMarkers() )
                {
                    if ( eltMarkersInMesh.find( matmarker ) == eltMarkersInMesh.end() )
                        continue;
                    for ( std::string const& p : currentPhysics )
                        M_markers[p].insert( matmarker );
                    markersByMaterial[matName].insert( matmarker );
                }
            }

            for ( std::string const& p : currentPhysics )
                M_isDefinedOnWholeMesh[p] = ( this->markers( p ).size() == eltMarkersInMesh.size() );

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

                auto [itProp,isAdded] = M_materialNameToProperties.emplace( matName, MaterialProperties<mesh_type>( matName ) );
                auto & matProperties = itProp->second;

#if 0
                M_thermalConductivityByMaterial[matName];
                if ( mat.hasProperty("k") )
                {
                    M_thermalConductivityByMaterial[matName] = mat.property( "k");
                    auto & prop = M_thermalConductivityByMaterial[matName];
                    if ( prop.hasAtLeastOneExpr() && ( !prop.hasExprScalar() && !prop.template hasExprMatrix<nDim,nDim>() )  )
                        CHECK( false ) << "thermal conductivty should be a scalar or a matrix";
                }
#endif
                using shape_dim_type = std::pair<uint16_type,uint16_type>;
                shape_dim_type scalarShape = std::make_pair(1,1);
                shape_dim_type matrixShape = std::make_pair(nDim,nDim);
                std::map<std::string,std::tuple<std::string,std::vector<shape_dim_type>>> mapPropNameToPropDesc;
                mapPropNameToPropDesc["density"] = std::make_tuple("rho", std::vector<shape_dim_type>( { scalarShape } ) );
                mapPropNameToPropDesc["specific-heat-capacity"] = std::make_tuple("Cp", std::vector<shape_dim_type>( { scalarShape } ) );
                mapPropNameToPropDesc["thermal-expansion"] = std::make_tuple("beta", std::vector<shape_dim_type>( { scalarShape } ) );
                mapPropNameToPropDesc["thermal-conductivity"] = std::make_tuple("k", std::vector<shape_dim_type>( { scalarShape, matrixShape } /*"scalar_or_matrix"*/ ) );

                for ( auto const& [propName,desc] : mapPropNameToPropDesc )
                {
                    std::string const& symbol = std::get<0>( desc );
                    if ( !mat.hasProperty( symbol ) )
                        continue;
                    auto const& possibleExprShapes = std::get<1>( desc );
                    auto const& matExprProperty = mat.property( symbol );
                    for ( auto [nComp1,nComp2] : std::get<1>( desc ) )
                    {
                        if ( matExprProperty.hasExpr(nComp1,nComp2) )
                        {
                            matProperties.add( propName/*symbol*/, matExprProperty );
                            break;
                        }
                    }
                }

                // rho * Cp
                M_rhoHeatCapacityByMaterial[matName];
                if ( this->hasDensity( matName ) && this->density( matName ).hasExprScalar()  && this->hasHeatCapacity( matName ) && this->heatCapacity( matName ).hasExprScalar() )
                {
                    auto rhoExpr = this->density( matName ).expr();
                    auto CpExpr = this->heatCapacity( matName ).expr();
                    auto expr = expr_mult<2>( rhoExpr,CpExpr,"",mesh->worldComm(),M_exprRepository );
                    M_rhoHeatCapacityByMaterial[matName].setExpr( expr );
                }
            }

            M_exprSelectorByMeshElementMapping = std::make_shared<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>>();
            M_exprSelectorByMeshElementMapping->template updateForUse<mesh_type>( this->rangeMeshElementsByMaterial() );
        }

    std::set<std::string> /*const&*/ markers( std::string const& p ) const
        {
            auto itFindMarkers = M_markers.find( p );
            if ( itFindMarkers == M_markers.end() )
                return std::set<std::string>{};
            else
                return itFindMarkers->second;
        }

    bool isDefinedOnWholeMesh( std::string const& p ) const
        {
            auto itFindMarkers = M_isDefinedOnWholeMesh.find( p );
            if ( itFindMarkers == M_isDefinedOnWholeMesh.end() )
                return false;
            else
                return itFindMarkers->second;
        }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElementsByMaterial( std::string const& matName ) const
        {
            CHECK( this->hasMaterial(matName) ) << "no material with name " << matName;
            return M_rangeMeshElementsByMaterial.find( matName )->second;
        }

    bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

    std::shared_ptr<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>> exprSelectorByMeshElementMapping() const { return M_exprSelectorByMeshElementMapping; }

    MaterialProperties<mesh_type> const&
    materialProperties( std::string const& matName ) const
        {
            auto itFindMatProp = M_materialNameToProperties.find( matName );
            CHECK( itFindMatProp != M_materialNameToProperties.end() ) << "material name not registered : " << matName;
            return itFindMatProp->second;
        }

    //! return true is the property \propName is defined in material \matName
    bool hasProperty( std::string const& matName, std::string const& propName ) const
        {
            if ( !this->hasMaterial( matName ) )
                return false;
            return this->materialProperties( matName ).has( propName );
        }

    //! return the property \propName which is defined in material \matName
    MaterialProperty const&
    materialProperty( std::string const& matName, std::string const& propName ) const
        {
            CHECK( this->hasProperty( matName, propName ) ) << "material property not registered";
            return this->materialProperties( matName ).property( propName );
        }

    //! return true is the property \propName is defined as a scalar expression in all materials
    bool materialPropertyIsScalarInAllMaterials( std::string const& propName ) const
        {
            for ( auto & [matName,matProps] : M_materialNameToProperties )
            {
                if ( !matProps.has( propName ) )
                    continue;
                if ( !matProps.property( propName ).isScalar() )
                    return false;
            }
            return true;
        }

    bool hasMaterialPropertyDependingOnSymbol( std::string const& propName, std::string const& symbolStr ) const
        {
            for ( auto & [matName,matProps] : M_materialNameToProperties )
            {
                if ( !matProps.has( propName ) )
                    continue;
                auto const& matProp = matProps.property( propName );

                if ( matProp.hasSymbolDependency( symbolStr ) )
                    return true;
            }
            return false;
        }


    // thermal conductivity
    bool hasThermalConductivity( std::string const& matName ) const
        {
            return this->hasProperty( matName, "thermal-conductivity" );
        }
    MaterialProperty const& thermalConductivity( std::string const& matName ) const
        {
            return this->materialProperty( matName, "thermal-conductivity" );
        }
    bool allThermalConductivitiesAreScalar() const
        {
            return this->materialPropertyIsScalarInAllMaterials( "thermal-conductivity" );
#if 0
            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& thermalConductivity = this->thermalConductivity( _matName );
                if ( !thermalConductivity.isScalar() )
                    return false;
            }
            return true;
#endif
        }

        bool hasThermalConductivityDependingOnSymbol( std::string const& symbolStr ) const
        {
            return this->hasMaterialPropertyDependingOnSymbol( "thermal-conductivity", symbolStr );
#if 0
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
#endif
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


    // density
    bool hasDensity( std::string const& matName ) const
        {
            return this->hasProperty( matName, "density" );
        }
    MaterialProperty/*ModelExpression*/ const& density( std::string const& matName ) const
        {
            return this->materialProperty( matName, "density" );
        }
    MaterialProperty/*ModelExpression*/ const& rho( std::string const& matName ) const { return this->density( matName ); } // DEPRECATED
    bool hasRho( std::string const& matName ) const { return this->hasDensity( matName ); } // DEPRECATED

        // heat capacity
    bool hasHeatCapacity( std::string const& matName ) const
        {
            return this->hasProperty( matName, "specific-heat-capacity" );
        }
    MaterialProperty const& heatCapacity( std::string const& matName ) const
        {
            return this->materialProperty( matName, "specific-heat-capacity" );
        }
    // thermal expansion
    bool hasThermalExpansion( std::string const& matName ) const
        {
            return this->hasProperty( matName, "thermal-expansion" );
        }
    MaterialProperty const& thermalExpansion( std::string const& matName ) const
        {
            return this->materialProperty( matName, "thermal-expansion" );
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
            for ( auto const&[physic,b] : M_isDefinedOnWholeMesh )
                *ostr << "\n     -- defined on whole mesh [" << physic << "] : " << b;
            *ostr << "\n     -- number of materials : " << M_rangeMeshElementsByMaterial.size();
            for ( auto const& matRange : M_rangeMeshElementsByMaterial)
            {
                std::string const& matName = matRange.first;
                auto const& matProps = this->materialProperties( matName );
                for ( auto const& [propName,matProp] : matProps )
                {
                    *ostr << "\n     -- [" << matName << "] " << propName << " : ";
                    *ostr << matProp.exprToString();
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
                for ( auto const& [propName,matProp] : this->materialProperties( matName ) )
                {
                    matPt.put( propName, matProp.exprToString() );
                }
                p.add_child( matName, matPt );
            }
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & [matName,matProps] : M_materialNameToProperties )
                matProps.setParameterValues( mp );
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
    std::map<std::string,std::set<std::string>> M_markers; // physic -> markers
    std::map<std::string,bool> M_isDefinedOnWholeMesh; // physics -> bool
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial; // matName -> range
    std::shared_ptr<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>> M_exprSelectorByMeshElementMapping;

    std::map<std::string, ModelExpressionScalar> /*M_thermalConductivityByMaterial,*/ /*M_heatCapacityByMaterial, M_rhoByMaterial, M_thermalExpansionByMaterial,*/ M_rhoHeatCapacityByMaterial;
    //std::map<std::string, ModelExpression> M_thermalConductivityByMaterial;
    //double M_thermalConductivityDefaultValue, M_heatCapacityDefaultValue, M_rhoDefaultValue, M_thermalExpansionDefaultValue;

    std::map<std::string,MaterialProperties<mesh_type>> M_materialNameToProperties;
};


} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
