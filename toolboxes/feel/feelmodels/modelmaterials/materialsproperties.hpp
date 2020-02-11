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
    MaterialProperty( std::string const& name, ModelExpression const& expr )
        :
        ModelExpression( expr ),
        M_name( name )
        {}
    MaterialProperty( MaterialProperty const& ) = default;
    MaterialProperty( MaterialProperty && ) = default;

    std::string const& name() const { return M_name; }
private :
    std::string M_name;
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
            this->emplace(propName, MaterialProperty(propName, expr) );
        }

    bool has( std::string const& propName ) const
        {
            if ( this->find( propName ) != this->end() )
                return true;
            return false;
        }

    MaterialProperty const&
    property( std::string const& propName ) const
        {
            auto itFind = this->find( propName );
            CHECK( itFind != this->end() ) << "material property" << propName << " not found";
            return itFind->second;
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & [propName,matProp] : *this )
                matProp.setParameterValues( mp );
        }

private :
    std::string M_materialName;
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

    void updateForUse( mesh_ptrtype const& mesh, ModelMaterials const& mats, ModelPhysics<nDim> const& modelphysics,
                       std::set<std::string> const& onlyTheseMaterialNames = std::set<std::string>{},
                       std::set<std::string> const& onlyTheseMarkers = std::set<std::string>{} )
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
                if ( !onlyTheseMaterialNames.empty() && (onlyTheseMaterialNames.find( matName ) == onlyTheseMaterialNames.end()) )
                    continue;
                auto const& mat = m.second;

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

                for ( std::string const& p : currentPhysics )
                    M_materialsNames[p].insert( matName );

                for ( std::string const& matmarker : mat.meshMarkers() )
                {
                    if ( !onlyTheseMarkers.empty() && (onlyTheseMarkers.find( matName ) == onlyTheseMarkers.end()) )
                        continue;
                    if ( eltMarkersInMesh.find( matmarker ) == eltMarkersInMesh.end() )
                        continue;
                    for ( std::string const& p : currentPhysics )
                        M_markers[p].insert( matmarker );
                    markersByMaterial[matName].insert( matmarker );
                }
            }

            for ( std::string const& p : currentPhysics )
                M_isDefinedOnWholeMesh[p] = ( this->markers( p ).size() == eltMarkersInMesh.size() );


            for ( auto const& [propName, propDesc] : modelphysics.matertialPropertyDescription() )
                M_matertialPropertyDescription[propName] = propDesc;

            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
                if ( !this->hasMaterial( matName ) )
                    continue;
                //if ( !onlyTheseMaterialNames.empty() && (onlyTheseMaterialNames.find( matName ) == onlyTheseMaterialNames.end()) )
                //continue;

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

                for ( auto const& [propName,desc] : M_matertialPropertyDescription  )
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

    //! return true if the material matName has been defined
    bool hasMaterial( std::string const& matName ) const
        {
            for ( auto const& [physic, matNames] : M_materialsNames )
                if ( matNames.find( matName ) != matNames.end() )
                    return true;
            return false;
        }

    //! return the materials name from a physic
    std::set<std::string> physicToMaterials( std::string const& physic ) const
        {
            auto itFindMat = M_materialsNames.find( physic );
            if( itFindMat != M_materialsNames.end() )
                return itFindMat->second;
            return std::set<std::string>{};
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

    //bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

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

    //! return true is the property \propName is defined as a scalar or a matrix TheDim x TheDim in all materials
    template <int TheDim>
    bool materialPropertyIsScalarOrMatrixInAllMaterials( std::string const& propName ) const
        {
            auto itFindPropDesc = M_matertialPropertyDescription.find( propName );
            if ( itFindPropDesc == M_matertialPropertyDescription.end() )
                return false;
            auto const& theShapes = std::get<1>( itFindPropDesc->second );
            if ( theShapes.size() != 2 )
                return false;
            if ( theShapes[0].first == 1 && theShapes[0].second == 1 && theShapes[1].first == TheDim && theShapes[1].second == TheDim )
                return true;
            if ( theShapes[1].first == 1 && theShapes[1].second == 1 && theShapes[0].first == TheDim && theShapes[0].second == TheDim )
                return true;
            return false;
        }

    //! return true if the  property \propName depends on symbol named \symbolStr
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

    template <int TheDim, typename SymbolsExpr>
    auto materialPropertyExprScalarOrMatrix( std::string const& propName, SymbolsExpr const& symbolsExpr ) const
        {
            auto Id = eye<TheDim,TheDim>();
            typedef decltype(expr(scalar_field_expression<2>{},symbolsExpr)*Id) _expr_scalar_type;
            std::vector<std::pair<std::string,_expr_scalar_type>> exprs_scalar_as_matrix;
            typedef decltype(expr(matrix_field_expression<TheDim,TheDim,2>{},symbolsExpr)) _expr_matrix_type;
            std::vector<std::pair<std::string,_expr_matrix_type>> exprs_matrix;

            for ( auto & [matName,matProps] : M_materialNameToProperties )
            {
                if ( !matProps.has( propName ) )
                    continue;
                auto const& matProp = matProps.property( propName );

                if ( matProp.template hasExpr<TheDim,TheDim>() )
                {
                    auto matPropExpr = expr( matProp.template expr<TheDim,TheDim>(), symbolsExpr );
                    exprs_matrix.push_back( std::make_pair( matName, matPropExpr ) );
                }
                else if ( matProp.hasExprScalar() )
                {
                    auto matPropExpr = expr( matProp.exprScalar(), symbolsExpr );
                    exprs_scalar_as_matrix.push_back( std::make_pair( matName, matPropExpr*Id ) );
                }
            }

            return expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, exprs_scalar_as_matrix, exprs_matrix );
        }

    template <typename SymbolsExpr>
    auto materialPropertyExprScalar( std::string const& propName, SymbolsExpr const& symbolsExpr ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},symbolsExpr)) _expr_scalar_type;
            std::vector<std::pair<std::string,_expr_scalar_type>> exprs_scalar;
            for ( auto & [matName,matProps] : M_materialNameToProperties )
            {
                if ( !matProps.has( propName ) )
                    continue;
                auto const& matProp = matProps.property( propName );

                if ( matProp.isScalar() )
                {
                    auto matPropExpr = expr( matProp.exprScalar(), symbolsExpr );
                    exprs_scalar.push_back( std::make_pair( matName, matPropExpr ) );
                }
            }

            return expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, exprs_scalar );
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
        }

    bool hasThermalConductivityDependingOnSymbol( std::string const& symbolStr ) const
        {
            return this->hasMaterialPropertyDependingOnSymbol( "thermal-conductivity", symbolStr );
        }

    template <typename SymbolsExpr>
    auto thermalConductivityExpr( SymbolsExpr const& symbolsExpr ) const
        {
            return this->materialPropertyExprScalarOrMatrix<nDim>( "thermal-conductivity", symbolsExpr );
        }

    template <typename SymbolsExpr>
    auto thermalConductivityScalarExpr( SymbolsExpr const& symbolsExpr ) const
        {
            return this->materialPropertyExprScalar( "thermal-conductivity", symbolsExpr );
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


    // electric conductivity
    bool hasElectricConductivity( std::string const& matName ) const
        {
            return this->hasProperty( matName, "electric-conductivity" );
        }
    MaterialProperty const& electricConductivity( std::string const& matName ) const
        {
            return this->materialProperty( matName, "electric-conductivity" );
        }
    bool hasElectricConductivityDependingOnSymbol( std::string const& symbolStr ) const
        {
            return this->hasMaterialPropertyDependingOnSymbol( "electric-conductivity", symbolStr );
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
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                for ( auto const& [propName,matProp] : matProps )
                {
                    auto itFindPropDesc = M_matertialPropertyDescription.find( propName );
                    if ( itFindPropDesc == M_matertialPropertyDescription.end() )
                        continue;
                    std::string const& symbolProp = std::get<0>( itFindPropDesc->second );
                    std::string symbolPrefixMatProp = (boost::format("%1%%2%_%3%")%prefix_symbol %matName %symbolProp).str();
                    if ( matProp.template hasExpr<nDim,nDim>() )
                    {
                        auto matPropExpr = expr( matProp.template expr<nDim,nDim>(), se );
                        matPropSymbsMatrix.push_back( std::make_tuple( symbolPrefixMatProp, matPropExpr, SymbolExprComponentSuffix( nDim,nDim,true ) ) );

                    }
                    else if ( matProp.hasExprScalar() )
                    {
                        auto matPropExpr = expr( matProp.exprScalar(), se );
                        matPropSymbsScalar.push_back( std::make_pair( symbolPrefixMatProp, matPropExpr ) );
                    }
                }
            }

            typedef decltype( this->materialPropertyExprScalar( "", se ) ) _expr_scalar_selector_type;
            std::vector<std::pair<std::string,_expr_scalar_selector_type>> matPropSymbsScalarSelector;
            typedef decltype( this->materialPropertyExprScalarOrMatrix<nDim>( "", se ) ) _expr_scalar_or_matrix_selector_type;
            std::vector<std::tuple<std::string,_expr_scalar_or_matrix_selector_type,SymbolExprComponentSuffix>> matPropSymbsScalarOrMatrixSelector;

            // generate the symbol heat_k if the all conductivities are scalar or heat_k(_xx,_xy,...) if at least one conductivity is a matrix
            for ( auto const& [propName,propDesc] : M_matertialPropertyDescription )
            {
                std::string const& symbolProp = std::get<0>( propDesc );
                if ( this->materialPropertyIsScalarInAllMaterials( propName ) )
                {
                    auto propExpr = this->materialPropertyExprScalar( propName,se );
                    matPropSymbsScalarSelector.push_back( std::make_pair( (boost::format("%1%%2%")% prefix_symbol %symbolProp).str(), propExpr ) );
                }
                else if ( this->materialPropertyIsScalarOrMatrixInAllMaterials<nDim>( propName ) )
                {
                    auto propExpr = materialPropertyExprScalarOrMatrix<nDim>( propName,se );
                    matPropSymbsScalarOrMatrixSelector.push_back( std::make_tuple( (boost::format("%1%%2%")% prefix_symbol %symbolProp).str(), propExpr,  SymbolExprComponentSuffix( nDim,nDim,true ) ) );
                }
            }

            return Feel::vf::symbolsExpr( symbolExpr( matPropSymbsScalar ), symbolExpr( matPropSymbsMatrix ), symbolExpr(matPropSymbsScalarSelector), symbolExpr( matPropSymbsScalarOrMatrixSelector ) );
        }

private :
    std::string M_exprRepository;
    std::map<std::string,std::set<std::string>> M_markers; // physic -> markers
    std::map<std::string,bool> M_isDefinedOnWholeMesh; // physics -> bool
    std::map<std::string,std::set<std::string>> M_materialsNames; // physic -> matNames
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial; // matName -> range
    std::shared_ptr<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>> M_exprSelectorByMeshElementMapping;

    std::map<std::string, ModelExpressionScalar> M_rhoHeatCapacityByMaterial;

    std::map<std::string,MaterialProperties<mesh_type>> M_materialNameToProperties;
    std::map<std::string,typename ModelPDE<nDim>::material_property_description_type> M_matertialPropertyDescription; // name -> (symbol, shapes.. )
};


} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
