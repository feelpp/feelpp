/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MODELMATERIALS_H
#define FEELPP_TOOLBOXES_MODELMATERIALS_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/ones.hpp>
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


class MaterialProperties : public std::map<std::string, MaterialProperty>
{
public :
    explicit MaterialProperties( std::string const& name ) : M_materialName( name ) {}
    MaterialProperties( MaterialProperties const& ) = default;
    MaterialProperties( MaterialProperties && ) = default;

    void add( std::string const& propName, ModelExpression const& expr )
        {
            auto itFind = this->find( propName );
            if ( itFind != this->end() )
                this->erase( itFind );
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

    std::set<std::string> markers() const { return M_markers; }
    void setMarkers( std::set<std::string> const& m ) { M_markers = m; }

    std::string const& materialName() const { return M_materialName; }

private :
    std::string M_materialName;
    std::set<std::string> M_markers;
};

struct MaterialsOnMeshBase;
template<class MeshType>
class MaterialsOnMesh;

template<uint16_type Dim>
class MaterialsProperties
{
    typedef MaterialsProperties<Dim> self_type;
    using modelphysics_type = ModelPhysics<Dim>;
    using modelphysics_ptrtype = std::shared_ptr<modelphysics_type>;
    using modelphysics_weakptrtype = std::weak_ptr<modelphysics_type>;
public :
    using physic_id_type = typename modelphysics_type::physic_id_type;
    static const uint16_type nDim = Dim;

    MaterialsProperties( modelphysics_ptrtype const& mphysics )
        :
        M_modelPhysics( mphysics )
        {}

    MaterialsProperties( MaterialsProperties const& ) = default;

    void updateForUse( ModelMaterials const& mats,
                       std::set<std::string> const& onlyTheseMaterialNames = std::set<std::string>{},
                       std::set<std::string> const& onlyTheseMarkers = std::set<std::string>{} )
        {
            modelphysics_ptrtype mphysics = M_modelPhysics.lock();
            worldcomm_t const& worldComm = mphysics->worldComm();
            std::string const& exprRepository = mphysics->repository().expr();

            std::map<std::string,std::set<std::string>> markersByMaterial;

            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
                if ( !onlyTheseMaterialNames.empty() && (onlyTheseMaterialNames.find( matName ) == onlyTheseMaterialNames.end()) )
                    continue;
                auto const& mat = m.second;

                std::set<physic_id_type> currentPhysics;
                for ( auto const& [physicId,physicObj] : mphysics->physics() )
                {
                    auto const& matOnPhysic = physicObj->materialNames();
                    if ( matOnPhysic.empty() || matOnPhysic.find( matName ) != matOnPhysic.end() ) // empty say all physics applied on mat
                        currentPhysics.insert( physicId );
                }

                if ( currentPhysics.empty() )
                    continue;

                for ( physic_id_type const& p : currentPhysics )
                    M_materialsNames[p].insert( matName );

                for ( std::string const& matmarker : mat.meshMarkers() )
                {
                    if ( !onlyTheseMarkers.empty() && (onlyTheseMarkers.find( matmarker ) == onlyTheseMarkers.end()) )
                        continue;
                    // if ( M_eltMarkersInMesh.find( matmarker ) == M_eltMarkersInMesh.end() )
                    //     continue;
                    markersByMaterial[matName].insert( matmarker );
                }
            }



            std::map<std::string,std::string> propSymbolToPropNameInDescription;
            for ( auto const& [physicName,physicData] : mphysics->physics() )
            {
                for ( auto const& [propName, propDesc] : physicData->materialPropertyDescription() )
                {
                    M_materialPropertyPhysicDescription[propName] = propDesc;
                    propSymbolToPropNameInDescription[std::get<0>( propDesc )] = propName;
                }
            }

            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
                bool _hasThisMat = false;
                for ( auto const& [physic, matNames] : M_materialsNames )
                    if ( matNames.find( matName ) != matNames.end() )
                        _hasThisMat = true;
                if ( !_hasThisMat )
                    continue;

                auto const& mat = m.second;
                auto itFindMat = markersByMaterial.find( matName );
                if ( itFindMat == markersByMaterial.end() )
                    continue;
                if ( itFindMat->second.empty() )
                    continue;

                auto [itProp,isAdded] = M_materialNameToProperties.emplace( matName, MaterialProperties( matName ) );
                auto & matProperties = itProp->second;
                // attach markers
                matProperties.setMarkers( itFindMat->second );

                for ( auto const& [propSymbol,propExpr] : mat.properties() )
                {
                    std::string propName = propSymbol;
                    auto itFindSymbolInDesc = propSymbolToPropNameInDescription.find( propSymbol );
                    if ( itFindSymbolInDesc != propSymbolToPropNameInDescription.end() )
                        propName = itFindSymbolInDesc->second;
                    this->addProperty( matProperties, propName, propExpr.mexpr() );
                }

                // TODO : move in heat toolbox
                // rho * Cp
                M_rhoHeatCapacityByMaterial[matName];
                if ( this->hasDensity( matName ) && this->density( matName ).hasExprScalar()  && this->hasHeatCapacity( matName ) && this->heatCapacity( matName ).hasExprScalar() )
                {
                    auto rhoExpr = this->density( matName ).expr();
                    auto CpExpr = this->heatCapacity( matName ).expr();
                    auto expr = expr_mult<2>( rhoExpr,CpExpr,"",worldComm,exprRepository );
                    M_rhoHeatCapacityByMaterial[matName].setExpr( expr );
                }

                // update Lame's parameters and bulk modulus if required
                bool hasLameFirstParameterInPhysicDesc = M_materialPropertyPhysicDescription.find( "Lame-first-parameter" ) != M_materialPropertyPhysicDescription.end();
                bool hasLameSecondParameterInPhysicDesc = M_materialPropertyPhysicDescription.find( "Lame-second-parameter" ) != M_materialPropertyPhysicDescription.end();
                bool hasBulkModulusInPhysicDesc = M_materialPropertyPhysicDescription.find( "bulk-modulus" ) != M_materialPropertyPhysicDescription.end();
                auto itFindYoungModulusInPhysicDesc = M_materialPropertyPhysicDescription.find( "Young-modulus" );
                auto itFindPoissonRatioInPhysicDesc = M_materialPropertyPhysicDescription.find( "Poisson-ratio" );
                bool hasYoungModulusInPhysicDesc = itFindYoungModulusInPhysicDesc != M_materialPropertyPhysicDescription.end();
                bool hasPoissonRatioInPhysicDesc = itFindPoissonRatioInPhysicDesc != M_materialPropertyPhysicDescription.end();
                if ( hasLameFirstParameterInPhysicDesc && hasYoungModulusInPhysicDesc && hasPoissonRatioInPhysicDesc &&
                     !matProperties.has( "Lame-first-parameter" ) && matProperties.has( "Young-modulus" ) && matProperties.has( "Poisson-ratio" ) )
                {
                    std::string const& YoungModulusSymb = itFindYoungModulusInPhysicDesc->second.symbol();
                    std::string const& PoissonRatioSymb = itFindPoissonRatioInPhysicDesc->second.symbol();
                    std::string lame1ExprStr = (boost::format("%1%*%2%/((1+%2%)*(1-2*%2%)):%1%:%2%") %YoungModulusSymb %PoissonRatioSymb).str();
                    ModelExpression lame1Expr;
                    lame1Expr.setExpr( lame1ExprStr,worldComm,exprRepository);
                    this->addProperty( matProperties, "Lame-first-parameter", lame1Expr, true );
                }
                if ( hasLameSecondParameterInPhysicDesc && hasYoungModulusInPhysicDesc && hasPoissonRatioInPhysicDesc &&
                     !matProperties.has( "Lame-second-parameter" ) && matProperties.has( "Young-modulus" ) && matProperties.has( "Poisson-ratio" ) )
                {
                    std::string const& YoungModulusSymb = itFindYoungModulusInPhysicDesc->second.symbol();
                    std::string const& PoissonRatioSymb = itFindPoissonRatioInPhysicDesc->second.symbol();
                    std::string lame2ExprStr = (boost::format("%1%/(2*(1+%2%)):%1%:%2%") %YoungModulusSymb %PoissonRatioSymb).str();
                    ModelExpression lame2Expr;
                    lame2Expr.setExpr( lame2ExprStr,worldComm,exprRepository);
                    this->addProperty( matProperties, "Lame-second-parameter", lame2Expr, true );
                }
                if ( hasBulkModulusInPhysicDesc && hasYoungModulusInPhysicDesc && hasPoissonRatioInPhysicDesc &&
                     !matProperties.has( "bulk-modulus" ) && matProperties.has( "Young-modulus" ) && matProperties.has( "Poisson-ratio" ) )
                {
                    std::string const& YoungModulusSymb = itFindYoungModulusInPhysicDesc->second.symbol();
                    std::string const& PoissonRatioSymb = itFindPoissonRatioInPhysicDesc->second.symbol();
                    std::string bulkModulusExprStr = (boost::format("%1%/(3*(1-2*%2%)):%1%:%2%") %YoungModulusSymb %PoissonRatioSymb).str();
                    ModelExpression bulkModulusExpr;
                    bulkModulusExpr.setExpr( bulkModulusExprStr,worldComm,exprRepository);
                    this->addProperty( matProperties, "bulk-modulus", bulkModulusExpr, true );
                }

                // rename symbols in current material : <propName> -> materials_<matName>_<propName>
                std::string prefix_symbol = "materials_";
                std::map<std::string,std::string> old2newSymbols;
                for ( auto const& [propName,matProp] : matProperties )
                {
                    auto itFindPropDesc = M_materialPropertyDescription.find( propName );
                    if ( itFindPropDesc == M_materialPropertyDescription.end() )
                        continue;

                    std::string const& symbolProp = std::get<0>( itFindPropDesc->second );
                    std::set<std::string> allSymbNames;
                    matProp.updateSymbolNames( symbolProp, allSymbNames );
                    for ( std::string const& _symbName : allSymbNames )
                        old2newSymbols[_symbName] = (boost::format("%1%%2%_%3%")%prefix_symbol %matName %_symbName).str();
                }
                for ( auto & [propName,matProp] : matProperties )
                    matProp.renameSymbols( old2newSymbols );
            }

            M_exprSelectorByMeshElementMapping = std::make_shared<ExprSelectorByMeshElementMapping<typename MeshBase<>/*mesh_type*/::index_type>>();
        }


    //! return the number of materials
    int numberOfMaterials() const { return M_materialNameToProperties.size(); }

    //! return true if the material matName has been defined
    bool hasMaterial( std::string const& matName ) const
        {
            return (M_materialNameToProperties.find( matName ) != M_materialNameToProperties.end());
        }

    //! return a mapping between physics and materials names
    std::map<physic_id_type,std::set<std::string>> const& physicToMaterials() const { return M_materialsNames; }

    //! return the materials names used with a set of physic
    std::set<std::string> physicToMaterials( std::set<physic_id_type> const& physicIds ) const
        {
            std::set<std::string> res;
            for ( auto const& physicId : physicIds )
            {
                auto itFindMat = M_materialsNames.find( physicId );
                if( itFindMat != M_materialsNames.end() )
                    res.insert( itFindMat->second.begin(), itFindMat->second.end() );
            }
            return res;
        }

    //! return the materials names used with a physic
    std::set<std::string> physicToMaterials( physic_id_type const& physicId ) const
        {
            return this->physicToMaterials( std::set<physic_id_type>({ physicId }) );
        }

    //! return true if the physic is defined in a material
    bool hasPhysic( physic_id_type const& physicId ) const
        {
            return M_materialsNames.find( physicId ) != M_materialsNames.end();
        }


    //! return the physics that are used in the material \matName from \physicsCollection
    std::map<physic_id_type,typename modelphysics_type::model_physic_ptrtype> physicsFromMaterial( std::string const& matName, std::map<physic_id_type,typename modelphysics_type::model_physic_ptrtype> const& physicsCollection ) const
        {
            std::map<physic_id_type,typename modelphysics_type::model_physic_ptrtype> res;
            for ( auto const& [physicId,physicData] : physicsCollection )
            {
                auto const& matNames = this->physicToMaterials( physicId );
                if ( matNames.find( matName ) == matNames.end() )
                    continue;
                res[physicId] = physicData;
            }
            return res;
        }


    void addProperty( MaterialProperties & matProperties, std::string const& propName, ModelExpression const& propExpr, bool onlyInPhysicDescription = false )
        {
            auto itFindPropPhysicDesc = M_materialPropertyPhysicDescription.find( propName );
            if ( itFindPropPhysicDesc != M_materialPropertyPhysicDescription.end() )
            {
                auto const& desc = itFindPropPhysicDesc->second;
                bool findProp = false;
                for ( auto [nComp1,nComp2] : desc.shapes() )
                {
                    if ( propExpr.hasExpr(nComp1,nComp2) )
                    {
                        matProperties.add( propName/*symbol*/, propExpr );
                        findProp = true;
                        break;
                    }
                }
                CHECK( findProp ) << "shape of the material property " << propName << " is not compatible with the physic";
                M_materialPropertyDescription.try_emplace( propName,desc );
            }
            else if ( !onlyInPhysicDescription )
            {
                std::string const& propSymbol = propName;
                hana::for_each( ModelExpression::expr_shapes, [this,&propName,&propExpr,&propSymbol,&matProperties]( auto const& e_ij )
                                {
                                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                                    if ( propExpr.template hasExpr<ni,nj>() )
                                    {
                                        matProperties.add( propName, propExpr );
                                        auto propShape = MaterialPropertyDescription::shape(ni,nj);
                                        auto itFindPropDesc = M_materialPropertyDescription.find( propName);
                                        if( itFindPropDesc == M_materialPropertyDescription.end() )
                                            M_materialPropertyDescription.emplace( propName, MaterialPropertyDescription( propSymbol, { propShape } ) );
                                        else
                                            itFindPropDesc->second.add( propShape );
                                    }
                                });
            }
        }

    std::map<std::string,MaterialProperties> const& materialNameToProperties() const { return M_materialNameToProperties; }

    std::map<std::string,MaterialPropertyDescription> const& materialPropertyDescription() const { return M_materialPropertyDescription; }
    std::map<std::string,MaterialPropertyDescription> const& materialPropertyPhysicDescription() const { return M_materialPropertyPhysicDescription; }

    std::shared_ptr<ExprSelectorByMeshElementMapping<typename MeshBase<>/*mesh_type*/::index_type>> exprSelectorByMeshElementMapping() const { return M_exprSelectorByMeshElementMapping; }


    MaterialProperties const&
    materialProperties( std::string const& matName ) const
        {
            auto itFindMatProp = M_materialNameToProperties.find( matName );
            CHECK( itFindMatProp != M_materialNameToProperties.end() ) << "material name not registered : " << matName;
            return itFindMatProp->second;
        }

    MaterialProperties &
    materialProperties( std::string const& matName )
        {
            auto itFindMatProp = M_materialNameToProperties.find( matName );
            CHECK( itFindMatProp != M_materialNameToProperties.end() ) << "material name not registered : " << matName;
            return itFindMatProp->second;
        }

    //! return true is the property \propName is defined in at leat one material
    bool hasMaterialProperty( std::string const& propName ) const
        {
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                if ( matProps.has( propName ) )
                    return true;
            }
            return false;
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
            CHECK( this->hasProperty( matName, propName ) ) << "material property " << propName << " is not registered in material " << matName;
            return this->materialProperties( matName ).property( propName );
        }

    //! return true is the property \propName has the same expression shape (M \times N) in all materials
    template <int M,int N>
    bool materialPropertyHasSameExprShapeInAllMaterials( std::string const& propName ) const
        {
            auto itFindPropDesc = M_materialPropertyDescription.find( propName );
            if ( itFindPropDesc == M_materialPropertyDescription.end() )
                return false;

            // look if the prop has only one shape available and exists
            auto const& theShapes = std::get<1>( itFindPropDesc->second );
            if ( theShapes.size() == 1 )
                if ( theShapes[0].first == M && theShapes[0].second == N && this->hasMaterialProperty( propName ) )
                    return true;

            // else loop on prop for each material
            bool findProp = false ;
            for ( auto & [matName,matProps] : M_materialNameToProperties )
            {
                if ( !matProps.has( propName ) )
                    continue;
                findProp = true;
                if ( !matProps.property( propName ).template hasExpr<M,N>() )
                    return false;
            }
            return findProp;
        }

    //! return true is the property \propName is defined as a scalar expression in all materials
    bool materialPropertyIsScalarInAllMaterials( std::string const& propName ) const
        {
            return this->materialPropertyHasSameExprShapeInAllMaterials<1,1>( propName );
        }


    //! return true is the property \propName is defined as a scalar or a matrix TheDim x TheDim in all materials
    template <int TheDim>
    bool materialPropertyIsScalarOrMatrixInAllMaterials( std::string const& propName ) const
        {
            auto itFindPropDesc = M_materialPropertyDescription.find( propName );
            if ( itFindPropDesc == M_materialPropertyDescription.end() )
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

    template <int TheDim, typename SymbolsExprType = symbols_expression_empty_t>
    auto materialPropertyExprScalarOrMatrix( std::string const& propName, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            auto Id = eye<TheDim,TheDim>();
            using _expr_scalar_type = std::decay_t< decltype( ModelExpression{}.template expr<1,1>() ) >;
            using _expr_matrix_type = std::decay_t< decltype( ModelExpression{}.template expr<TheDim,TheDim>() ) >;
            using _expr_scalar_as_matrix_type = std::decay_t< decltype(_expr_scalar_type{}*Id) >;
            std::vector<std::pair<std::string,_expr_scalar_as_matrix_type>> exprs_scalar_as_matrix;
            std::vector<std::pair<std::string,_expr_matrix_type>> exprs_matrix;

            for ( auto const& [matName,matProps] : this->materialNameToProperties() )
            {
                if ( !matProps.has( propName ) )
                    continue;
                auto const& matProp = matProps.property( propName );

                if ( matProp.template hasExpr<TheDim,TheDim>() )
                {
                    auto const& matPropExpr = matProp.template expr<TheDim,TheDim>();
                    exprs_matrix.push_back( std::make_pair( matName, matPropExpr ) );
                }
                else if ( matProp.hasExprScalar() )
                {
                    auto const& matPropExpr = matProp.exprScalar();
                    exprs_scalar_as_matrix.push_back( std::make_pair( matName, matPropExpr*Id ) );
                }
            }

            if constexpr ( std::is_same_v<SymbolsExprType,symbols_expression_empty_t> )
                return expr<typename MeshBase<>/*mesh_type*/::index_type>( M_exprSelectorByMeshElementMapping, exprs_scalar_as_matrix, exprs_matrix );
            else
                return expr<typename MeshBase<>/*mesh_type*/::index_type>( M_exprSelectorByMeshElementMapping, exprs_scalar_as_matrix, exprs_matrix ).applySymbolsExpr( se );
        }

    template <int M,int N, typename SymbolsExprType = symbols_expression_empty_t>
    auto materialPropertyExpr( std::string const& propName, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            using _expr_type = std::decay_t< decltype( ModelExpression{}.template expr<M,N>() ) >;
            std::vector<std::pair<std::string,_expr_type>> theExprs;
            for ( auto const& [matName,matProps] : this->materialNameToProperties() )
            {
                if ( !matProps.has( propName ) )
                    continue;
                auto const& matProp = matProps.property( propName );

                if ( matProp.template hasExpr<M,N>() )
                {
                    auto const& matPropExpr = matProp.template expr<M,N>();
                    theExprs.push_back( std::make_pair( matName, matPropExpr ) );
                }
            }

            if constexpr ( std::is_same_v<SymbolsExprType,symbols_expression_empty_t> )
                return expr<typename MeshBase<>/*mesh_type*/::index_type>( M_exprSelectorByMeshElementMapping, theExprs );
            else
                return expr<typename MeshBase<>/*mesh_type*/::index_type>( M_exprSelectorByMeshElementMapping, theExprs ).applySymbolsExpr( se );
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

    template <typename SymbolsExprType = symbols_expression_empty_t>
    auto thermalConductivityExpr( SymbolsExprType const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            return this->materialPropertyExprScalarOrMatrix<nDim>( "thermal-conductivity", symbolsExpr );
        }

    template <typename SymbolsExprType = symbols_expression_empty_t>
    auto thermalConductivityScalarExpr( SymbolsExprType const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            return this->materialPropertyExpr<1,1>( "thermal-conductivity", symbolsExpr );
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
#if 0 // TODO
            for ( auto const&[physic,b] : M_isDefinedOnWholeMesh )
                *ostr << "\n     -- defined on whole mesh [" << physic << "] : " << b;
#endif
            *ostr << "\n     -- number of materials : " << this->numberOfMaterials();
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                for ( auto const& [propName,matProp] : matProps )
                {
                    *ostr << "\n     -- [" << matName << "] " << propName << " : ";
                    *ostr << matProp.exprToString();
                }
            }
            return ostr;
        }

    void updateInformationObject( nl::json & p ) const
        {
            std::string prefix_symbol = "materials_";
            int nMat = this->numberOfMaterials();
            p.emplace( "number of materials", nMat );

            std::map<std::string,std::string> propertyNamesToSymbol;
            auto & jByMat = p["by_material"];
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                nl::json::array_t jaOnMat;
                for ( auto const& [propName,matProp] : matProps )
                {
                    //p[matName].emplace( propName, matProp.exprToString() );
                    if ( !matProp.hasAtLeastOneExpr() )
                        continue;

                    auto itFindPropDesc = M_materialPropertyDescription.find( propName );
                    if ( itFindPropDesc == M_materialPropertyDescription.end() )
                        continue;
                    std::string const& symbolProp = std::get<0>( itFindPropDesc->second );
                    std::string symbolPrefixMatProp = (boost::format("%1%%2%_%3%")%prefix_symbol %matName %symbolProp).str();
                    auto [exprStr,compInfo] =  matProp.exprInformations();
                    jaOnMat.push_back( symbolExprInformations( symbolPrefixMatProp, exprStr, compInfo, propName ) );

                    if ( nMat>1 )
                        propertyNamesToSymbol[propName] = symbolProp; // need for define multi-material symbols
                }
                jByMat.emplace( matName, jaOnMat );
            }


            // multi-material symbols (same shape)
            nl::json::array_t jaMultiMat;
            hana::for_each( ModelExpression::expr_shapes, [this,&prefix_symbol,&jaMultiMat,&propertyNamesToSymbol](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;

                    std::vector<typename std::map<std::string,std::string>::iterator> propItOfMultiMatSameShape;
                    for ( auto it = propertyNamesToSymbol.begin(), en = propertyNamesToSymbol.end(); it != en; ++it )
                    {
                        auto const& [propName,propSymbol] = *it;
                        if ( this->materialPropertyHasSameExprShapeInAllMaterials<ni,nj>( propName ) )
                        {
                            std::string symbolGlobalMatProp = (boost::format("%1%%2%")% prefix_symbol %propSymbol).str();
                            auto compInfo =  SymbolExprComponentSuffix( ni,nj );
                            jaMultiMat.push_back( symbolExprInformations( symbolGlobalMatProp, "concat(...)", compInfo, propName ) );
                            propItOfMultiMatSameShape.push_back( it );
                        }
                    }
                    for ( auto const& it : propItOfMultiMatSameShape )
                        propertyNamesToSymbol.erase( it );
                });

            // multi-material symbols (get symbols heat_k_xx, heat_k_xy,... if a prop is scalar in one mat and matrix in another mat)
            for ( auto const& [propName,propSymbol] : propertyNamesToSymbol )
            {
                if ( this->materialPropertyIsScalarOrMatrixInAllMaterials<nDim>( propName ) )
                {
                    std::string symbolGlobalMatProp = (boost::format("%1%%2%")% prefix_symbol %propSymbol).str();
                    auto compInfo =  SymbolExprComponentSuffix( nDim,nDim );
                    jaMultiMat.push_back( symbolExprInformations( symbolGlobalMatProp, "concat(...)", compInfo, propName ) );
                }
            }

            if ( !jaMultiMat.empty() )
                p["multi-material"] = jaMultiMat;
        }

    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
        {
            auto tabInfo = TabulateInformationsSections::New();

            Feel::Table tabInfoOthers;
            TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoOthers, jsonInfo, tabInfoProp );
            tabInfoOthers.format()
                .setShowAllBorders( false )
                .setColumnSeparator(":")
                .setHasRowSeparator( false );
            tabInfo->add( "", TabulateInformations::New( tabInfoOthers ) );

            if ( jsonInfo.contains("by_material") )
            {
                for ( auto const& el : jsonInfo.at("by_material").items() )
                {
                    std::string const& matName = el.key();
                    tabInfo->add( "Material : "+matName, TabulateInformationTools::FromJSON::tabulateInformationsSymbolsExpr( el.value(),tabInfoProp,true) );
                }
            }

            if ( jsonInfo.contains( "multi-material") )
                tabInfo->add( "Multi-Materials", TabulateInformationTools::FromJSON::tabulateInformationsSymbolsExpr( jsonInfo.at("multi-material"),tabInfoProp.newByIncreasingVerboseLevel(),true) );

            return tabInfo;
        }


    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & [matName,matProps] : M_materialNameToProperties )
                matProps.setParameterValues( mp );
        }

    std::map<std::string,double>
    toParameterValues( std::string const& prefix_symbol = "materials_" ) const
        {
            std::map<std::string,double> pv;
            int nMat = this->numberOfMaterials();
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                for ( auto const& [propName,matProp] : matProps )
                {
                    auto itFindPropDesc = M_materialPropertyDescription.find( propName );
                    if ( itFindPropDesc == M_materialPropertyDescription.end() )
                        continue;

                    std::string const& symbolProp = std::get<0>( itFindPropDesc->second );
                    std::set<std::string> symbolMatProp( { (boost::format("%1%%2%_%3%")%prefix_symbol %matName %symbolProp).str() } );
                    if ( nMat == 1 )
                        symbolMatProp.insert( (boost::format("%1%%2%")%prefix_symbol %symbolProp).str() );

                    matProp.updateParameterValues( symbolMatProp, pv );
                }
            }
            return pv;
        }


    //! update constant (and scalar) material properties into the mapping of values \mp
    void updateParameterValues( std::map<std::string,double> & mp, std::string const& prefix_symbol = "materials_" )// const
        {
            // start by updated the expression from mp
            this->setParameterValues( mp );

            // get all symbol names from material properties
            int nMat = this->numberOfMaterials();
            std::set<std::string> allSymbNames;
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                for ( auto const& [propName,matProp] : matProps )
                {
                    auto itFindPropDesc = M_materialPropertyDescription.find( propName );
                    if ( itFindPropDesc == M_materialPropertyDescription.end() )
                        continue;

                    std::string const& symbolProp = std::get<0>( itFindPropDesc->second );
                    std::set<std::string> symbolMatProp( { (boost::format("%1%%2%_%3%")%prefix_symbol %matName %symbolProp).str() } );
                    if ( nMat == 1 )
                        symbolMatProp.insert( (boost::format("%1%%2%")%prefix_symbol %symbolProp).str() );;
                    matProp.updateSymbolNames( symbolMatProp, allSymbNames );
                }
            }

            // erase parameters values from material properties
            for ( auto & [matName,matProps] : M_materialNameToProperties )
                for ( auto & [propName,matProp] : matProps )
                    matProp.eraseParameterValues( allSymbNames );

            // get value of symbols evaluables
            int previousParam = 0;
            std::map<std::string,double> curmp;
            while ( true )
            {
                curmp = this->toParameterValues( prefix_symbol );
                if ( curmp.size() == previousParam )
                    break;
                previousParam = curmp.size();
                this->setParameterValues( curmp );
            }

            // update output parameter values
            for ( auto const& [_name,_value] : curmp )
                mp[_name] = _value;
        }

    auto symbolsExpr( std::string const& prefix_symbol = "materials_" ) const
        {
            int nMat = this->numberOfMaterials();
            std::map<std::string,std::string> propertyNamesToSymbol; //usefull for expression on whole mesh
            // generate symbols heat_matName_k or heat_matName_k(_xx,_xy,...,_zz)
            auto tupleSymbolExprsByMaterial = hana::transform( ModelExpression::expr_shapes, [this,&prefix_symbol,&nMat,&propertyNamesToSymbol](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;

                    using _expr_type = std::decay_t< decltype( ModelExpression{}.template expr<ni,nj>() ) >;
                    symbol_expression_t<_expr_type> se;
                    for ( auto const& [matName,matProps] : M_materialNameToProperties )
                    {
                        for ( auto const& [propName,matProp] : matProps )
                        {
                            auto itFindPropDesc = M_materialPropertyDescription.find( propName );
                            if ( itFindPropDesc == M_materialPropertyDescription.end() )
                                continue;

                            if ( !matProp.template hasExpr<ni,nj>() )
                                continue;

                            std::string const& symbolProp = std::get<0>( itFindPropDesc->second );

                            if ( nMat>1 )
                                propertyNamesToSymbol[propName] = symbolProp; // need for define the global expr

                            auto const& matPropExpr = matProp.template expr<ni,nj>();
                            if ( matPropExpr.expression().isEvaluable() )
                                continue;

                            std::string symbolPrefixMatProp = (boost::format("%1%%2%_%3%")%prefix_symbol %matName %symbolProp).str();
                            std::string symbolGlobalMatProp = (boost::format("%1%%2%")%prefix_symbol %symbolProp).str();

                            se.add( symbolPrefixMatProp, matPropExpr, SymbolExprComponentSuffix( ni,nj ) );
                            if ( nMat == 1 )
                                se.add( symbolGlobalMatProp, matPropExpr, SymbolExprComponentSuffix( ni,nj ) );
                        }
                    }
                    return se;
                });

            // generate symbol heat_k(_xx,_xy,...) if the all property k has same shape
            auto tupleSymbolExprsOnAllMaterials = hana::transform( ModelExpression::expr_shapes, [this,&prefix_symbol,&propertyNamesToSymbol](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                    using _expr_type =  std::decay_t< decltype( this->materialPropertyExpr<ni,nj>( "" ) ) >;
                    symbol_expression_t<_expr_type> se;
                    std::vector<typename std::map<std::string,std::string>::iterator> propItOfMultiMatSameShape;
                    for ( auto it = propertyNamesToSymbol.begin(), en = propertyNamesToSymbol.end(); it != en; ++it )
                    {
                        auto const& [propName,propSymbol] = *it;
                        if ( this->materialPropertyHasSameExprShapeInAllMaterials<ni,nj>( propName ) )
                        {
                            std::string symbolGlobalMatProp = (boost::format("%1%%2%")% prefix_symbol %propSymbol).str();
                            auto globalMatPropExpr = this->materialPropertyExpr<ni,nj>( propName );
                            se.add( symbolGlobalMatProp, globalMatPropExpr, SymbolExprComponentSuffix( ni,nj ) );
                            propItOfMultiMatSameShape.push_back( it );
                        }
                    }
                    for ( auto const& it : propItOfMultiMatSameShape )
                        propertyNamesToSymbol.erase( it );
                    return se;
                });

            // generate symbols heat_k_xx, heat_k_xy,... if a prop is scalar in one mat and matrix in another mat
            typedef decltype( this->materialPropertyExprScalarOrMatrix<nDim>( "" ) ) _expr_scalar_or_matrix_selector_type;
            symbol_expression_t<_expr_scalar_or_matrix_selector_type> seScalarOrMatrixSelector;
            for ( auto const& [propName,propSymbol] : propertyNamesToSymbol )
            {
                if ( this->materialPropertyIsScalarOrMatrixInAllMaterials<nDim>( propName ) )
                {
                    auto propExpr = this->materialPropertyExprScalarOrMatrix<nDim>( propName );
                    seScalarOrMatrixSelector.add( (boost::format("%1%%2%")% prefix_symbol %propSymbol).str(), propExpr,  SymbolExprComponentSuffix( nDim,nDim ) );
                }
            }

            return Feel::vf::symbolsExpr( Feel::vf::SymbolsExpr( tupleSymbolExprsByMaterial ),
                                          Feel::vf::SymbolsExpr( tupleSymbolExprsOnAllMaterials ),
                                          seScalarOrMatrixSelector
                                          );
        }

    template <typename MeshType,typename SymbExprType>
    auto exprPostProcessExports( std::shared_ptr<MeshType> mesh, std::set<physic_id_type> const& physics, SymbExprType const& se, std::string const& prefix = "materials" ) const
        {
            auto mom = this->materialsOnMesh( mesh );
            auto setOfMatNameUsed = this->physicToMaterials( physics );
            std::map<std::string,std::map<std::string,const MaterialProperty* >> propertiesNotHaveSameShape;

            auto tupleExportsExprs = hana::transform( ModelExpression::expr_shapes, [this,&setOfMatNameUsed,&se,&prefix,&propertiesNotHaveSameShape,&mom](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;

                    using _expr_type = std::decay_t< decltype( expr( ModelExpression{}.template expr<ni,nj>(), se ) ) >;
                    std::map<std::string,std::vector<std::tuple<_expr_type, elements_reference_wrapper_t<MeshType>, std::string > > > mapExprExported;

                    if ( !mom )
                        return mapExprExported;

                    for ( auto const& [matName,matProps] : M_materialNameToProperties )
                    {
                        if ( setOfMatNameUsed.find( matName ) == setOfMatNameUsed.end() )
                            continue;
                        auto const& range = mom->rangeMeshElementsByMaterial( matName );
                        for ( auto const& [propName,matProp] : matProps )
                        {
                            // only properties in physic descroption
                            auto itFindPropDesc = M_materialPropertyPhysicDescription.find( propName );
                            if ( itFindPropDesc == M_materialPropertyPhysicDescription.end() )
                                continue;

                            if ( !matProp.template hasExpr<ni,nj>() )
                                continue;
                            if ( this->materialPropertyHasSameExprShapeInAllMaterials<ni,nj>( propName ) )
                            {
                                auto matPropExpr = expr( matProp.template expr<ni,nj>(), se );
                                mapExprExported[prefixvm(prefix,propName)].push_back( std::make_tuple( matPropExpr, range, "element" ) );
                            }
                            else
                            {
                                propertiesNotHaveSameShape[matName][propName] = &matProp;
                            }
                        }
                    }
                    return mapExprExported;
                });

            typedef decltype(expr(ModelExpression{}.template expr<nDim,nDim>(),se)) _expr_tensor2_type;
            std::map<std::string,std::vector<std::tuple<_expr_tensor2_type, elements_reference_wrapper_t<MeshType>, std::string > > > mapExprTensor2;

            auto Id = eye<nDim,nDim>();
            typedef decltype(expr(typename ModelExpression::expr_scalar_type{},se)*Id) _expr_tensor2_from_scalar_type;
            std::map<std::string,std::vector<std::tuple<_expr_tensor2_from_scalar_type, elements_reference_wrapper_t<MeshType>, std::string > > > mapExprTensor2FromScalar;

            for ( auto const& [matName,matProps] : propertiesNotHaveSameShape )
            {
                auto const& range = mom->rangeMeshElementsByMaterial( matName );
                for ( auto const& [propName,matPropPtr] : matProps )
                {
                    if ( this->materialPropertyIsScalarOrMatrixInAllMaterials<nDim>( propName ) )
                    {
                        auto const& matProp = *matPropPtr;
                        if ( matProp.template hasExpr<nDim,nDim>() )
                        {
                            auto matPropExpr = expr( matProp.template expr<nDim,nDim>(), se );
                            // TODO : insert mapExprTensor2 into tupleExportsExprs
                            mapExprTensor2[prefixvm(prefix,propName)].push_back( std::make_tuple( matPropExpr, range, "element" ) );
                        }
                        else
                        {
                            auto matPropExpr = expr( matProp.exprScalar(), se );
                            mapExprTensor2FromScalar[prefixvm(prefix,propName)].push_back( std::make_tuple( matPropExpr*Id, range, "element" ) );
                        }
                    }
                }

            }

            return  hana::concat( tupleExportsExprs, hana::make_tuple( mapExprTensor2,mapExprTensor2FromScalar ) );
        }

    template <typename MeshType>
    std::set<std::string> postProcessExportsAllFieldsAvailable( std::shared_ptr<MeshType> mesh, std::set<physic_id_type> const& physics, std::string const& prefix = "materials" ) const
        {
            std::set<std::string> res;
            auto mom = this->materialsOnMesh( mesh );
            if ( !mom )
                return res;
            auto setOfMatNameUsed = this->physicToMaterials( physics );
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                if ( setOfMatNameUsed.find( matName ) == setOfMatNameUsed.end() )
                    continue;
                auto range = mom->rangeMeshElementsByMaterial( matName );
                for ( auto const& [propName,matProp] : matProps )
                {
                    auto itFindPropDesc =  M_materialPropertyPhysicDescription.find( propName );
                    if ( itFindPropDesc ==  M_materialPropertyPhysicDescription.end() )
                        continue;
                    res.insert(prefixvm(prefix,propName));
                }
            }
            return res;
        }


    template <typename MeshType>
    void addMesh( std::shared_ptr<MeshType> mesh );

    template <typename MeshType>
    void removeMesh( std::shared_ptr<MeshType> mesh );

    template <typename MeshType>
    std::shared_ptr<MaterialsOnMesh<MeshType>> materialsOnMesh( std::shared_ptr<MeshType> mesh ) const
        {
            auto itFindMatMesh = M_materialsOnMesh.find( mesh );
            if ( itFindMatMesh == M_materialsOnMesh.end() )
                return std::shared_ptr<MaterialsOnMesh<MeshType>>{};
            return std::dynamic_pointer_cast<MaterialsOnMesh<MeshType>>( itFindMatMesh->second );
        }
    template <typename MeshType>
    elements_reference_wrapper_t<MeshType> const& rangeMeshElementsByMaterial( std::shared_ptr<MeshType> const& mesh, std::string const& matName ) const
        {
            auto mom = this->materialsOnMesh(mesh);
            CHECK( mom ) << "no materialmesh";
            return mom->rangeMeshElementsByMaterial( matName );
        }


private :
    modelphysics_weakptrtype M_modelPhysics;

    std::map<physic_id_type,std::set<std::string>> M_materialsNames; // physic -> matNames

    std::map<std::string, ModelExpressionScalar> M_rhoHeatCapacityByMaterial;

    std::map<std::string,MaterialProperties> M_materialNameToProperties;
    std::map<std::string,MaterialPropertyDescription> M_materialPropertyPhysicDescription; // name -> (symbol, shapes.. )   only defined in phycis
    std::map<std::string,MaterialPropertyDescription> M_materialPropertyDescription; // name -> (symbol, shapes.. )   all props read

    std::shared_ptr<ExprSelectorByMeshElementMapping<typename MeshBase<>/*mesh_type*/::index_type>> M_exprSelectorByMeshElementMapping;

    std::map<std::shared_ptr<MeshBase<>>, std::shared_ptr<MaterialsOnMeshBase>> M_materialsOnMesh;
};

struct MaterialsOnMeshBase
{
    virtual ~MaterialsOnMeshBase() {}
};


template<class MeshType>
class MaterialsOnMesh : public MaterialsOnMeshBase
{
    typedef MaterialsOnMesh<MeshType> self_type;
public :
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    static const uint16_type nDim = mesh_type::nDim;
    static const uint16_type nRealDim = mesh_type::nRealDim;

    using materials_properties_type = MaterialsProperties<nRealDim>;
    using materials_properties_ptrtype = std::shared_ptr<materials_properties_type>;
    using physic_id_type = typename materials_properties_type::physic_id_type;

    MaterialsOnMesh( materials_properties_type const& materialsProperties, mesh_ptrtype mesh )
        :
        M_materialsProperties( materialsProperties )
        {
            M_eltMarkersInMesh.clear();
            for (auto const& markPair : mesh->markerNames() )
            {
                std::string meshMarker = markPair.first;
                if ( mesh->hasElementMarker( meshMarker ) )
                    M_eltMarkersInMesh.insert( meshMarker );
            }

            std::map<std::string,std::set<std::string>> materialToMarkers;
            for ( auto const& [p,matNames] : M_materialsProperties.physicToMaterials() )
            {
                for ( std::string const& matName : matNames )
                {
                    auto const& matProps = M_materialsProperties.materialProperties( matName );
                    for ( std::string const& m : matProps.markers() )
                        if ( M_eltMarkersInMesh.find( m ) != M_eltMarkersInMesh.end() )
                        {
                            M_markers[p].insert( m );
                            materialToMarkers[matName].insert( m );
                        }
                }
                M_isDefinedOnWholeMesh[p] = ( this->markers( p ).size() == M_eltMarkersInMesh.size() );
            }

            for ( auto const& [matName,matmarkers] : materialToMarkers )
            {
                auto range = markedelements( mesh,matmarkers );
                M_rangeMeshElementsByMaterial[matName] = range;
            }

        }

    materials_properties_type const& materialsProperties() const { return M_materialsProperties; }

    void markers( physic_id_type const& p, std::set<std::string> & res ) const
        {
            auto itFindMarkers = M_markers.find( p );
            if ( itFindMarkers != M_markers.end() )
                res.insert( itFindMarkers->second.begin(), itFindMarkers->second.end() );
        }

    std::set<std::string> markers( physic_id_type const& p ) const
        {
            std::set<std::string> res;
            this->markers( p, res );
            return res;
        }

    std::set<std::string> markers( std::set<physic_id_type> const& setOfPhysics ) const
        {
            std::set<std::string> res;
            for ( physic_id_type const& p : setOfPhysics )
                this->markers( p, res );
            return res;
        }

    bool isDefinedOnWholeMesh( physic_id_type const& p ) const
        {
            auto itFindMarkers = M_isDefinedOnWholeMesh.find( p );
            if ( itFindMarkers == M_isDefinedOnWholeMesh.end() )
                return false;
            else
                return itFindMarkers->second;
        }


    bool isDefinedOnWholeMesh( std::set<physic_id_type> const& setOfPhysics ) const
        {
            return this->markers( setOfPhysics ).size() == M_eltMarkersInMesh.size();
        }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElementsByMaterial( std::string const& matName ) const
        {
            // CHECK( this->materialsProperties()->hasMaterial(matName) ) << "no material with name " << matName;
            auto itFindMat = M_rangeMeshElementsByMaterial.find( matName );
            CHECK( itFindMat != M_rangeMeshElementsByMaterial.end() ) << "no material with name " << matName;
            return itFindMat->second;
        }


private :
    materials_properties_type const& M_materialsProperties;

    std::set<std::string> M_eltMarkersInMesh;
    std::map<physic_id_type,bool> M_isDefinedOnWholeMesh; // physics -> bool
    std::map<physic_id_type,std::set<std::string>> M_markers; // physic -> markers
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial; // matName -> range

};


template<uint16_type Dim>
template <typename MeshType>
void
MaterialsProperties<Dim>::addMesh( std::shared_ptr<MeshType> mesh )
{
    if ( M_materialsOnMesh.find( mesh->shared_from_this_meshbase() ) != M_materialsOnMesh.end() )
        return;
    auto mom = std::make_shared<MaterialsOnMesh<MeshType>>( *this, mesh );
    M_materialsOnMesh[ mesh->shared_from_this_meshbase() ] = mom;
    M_exprSelectorByMeshElementMapping->template updateForUse<MeshType>( mom->rangeMeshElementsByMaterial() );
}

template<uint16_type Dim>
template <typename MeshType>
void
MaterialsProperties<Dim>::removeMesh( std::shared_ptr<MeshType> mesh )
{
    auto itFindMesh = M_materialsOnMesh.find( mesh->shared_from_this_meshbase() );
    if ( itFindMesh == M_materialsOnMesh.end() )
        return;
    M_materialsOnMesh.erase( itFindMesh );

    // clear elt mapping (TODO : only for this mesh)
    M_exprSelectorByMeshElementMapping->clear();

}

} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
