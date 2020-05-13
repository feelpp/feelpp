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
        {}

    MaterialsProperties( MaterialsProperties const& ) = default;

    void updateForUse( mesh_ptrtype const& mesh, ModelMaterials const& mats, ModelPhysics<nDim> const& modelphysics,
                       std::set<std::string> const& onlyTheseMaterialNames = std::set<std::string>{},
                       std::set<std::string> const& onlyTheseMarkers = std::set<std::string>{} )
        {
            M_eltMarkersInMesh.clear();
            for (auto const& markPair : mesh->markerNames() )
            {
                std::string meshMarker = markPair.first;
                if ( mesh->hasElementMarker( meshMarker ) )
                    M_eltMarkersInMesh.insert( meshMarker );
            }

            //auto const& mapPhysicsToSubphysics = modelphysics.mapPhysicsToSubphysics();
            auto const& defaultPhysics = modelphysics.physicDefault();
            auto const& physicsAvailable = modelphysics.physicsAvailable();

            std::map<std::string,std::set<std::string>> markersByMaterial;
            M_markers.clear();

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


                std::set<std::string> currentPhysics;
                for ( std::string const& p : currentPhysicsReaded )
                {
                    auto thePhysicsShared = modelphysics.physicsShared( p );
                    currentPhysics.insert( thePhysicsShared.begin(), thePhysicsShared.end() );
                    // auto itFindPhysics = mapPhysicsToSubphysics.find( p );
                    // if ( itFindPhysics != mapPhysicsToSubphysics.end() )
                    //     currentPhysics.insert( itFindPhysics->second.begin(), itFindPhysics->second.end() );
                }

                if ( currentPhysics.empty() )
                    continue;

                for ( std::string const& p : currentPhysics )
                    M_materialsNames[p].insert( matName );

                for ( std::string const& matmarker : mat.meshMarkers() )
                {
                    if ( !onlyTheseMarkers.empty() && (onlyTheseMarkers.find( matmarker ) == onlyTheseMarkers.end()) )
                        continue;
                    if ( M_eltMarkersInMesh.find( matmarker ) == M_eltMarkersInMesh.end() )
                        continue;
                    for ( std::string const& p : currentPhysics )
                        M_markers[p].insert( matmarker );
                    markersByMaterial[matName].insert( matmarker );
                }
            }

            for ( std::string const& p : physicsAvailable /*currentPhysics*/ )
                M_isDefinedOnWholeMesh[p] = ( this->markers( p ).size() == M_eltMarkersInMesh.size() );


            std::map<std::string,std::string> propSymbolToPropNameInDescription;
#if 0
            for ( auto const& [propName, propDesc] : modelphysics.materialPropertyDescription() )
            {
                 M_materialPropertyPhysicDescription[propName] = propDesc;
                propSymbolToPropNameInDescription[std::get<0>( propDesc )] = propName;
            }
#else
            for ( auto const& [physicName,physicData] : modelphysics.physics() )
            {
                for ( auto const& [propName, propDesc] : physicData->materialPropertyDescription() )
                {
                    M_materialPropertyPhysicDescription[propName] = propDesc;
                    propSymbolToPropNameInDescription[std::get<0>( propDesc )] = propName;
                }
            }
#endif

            for( auto const& m : mats )
            {
                std::string const& matName = m.first;
#if 0
                if ( !this->hasMaterial( matName ) )
                    continue;
#else

                bool _hasThisMat = false;
                for ( auto const& [physic, matNames] : M_materialsNames )
                    if ( matNames.find( matName ) != matNames.end() )
                        _hasThisMat = true;
                if ( !_hasThisMat )
                    continue;
#endif


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

                for ( auto const& [propSymbol,propExpr] : mat.properties() )
                {
                    auto itFindSymbolInDesc = propSymbolToPropNameInDescription.find( propSymbol );
                    if ( itFindSymbolInDesc != propSymbolToPropNameInDescription.end() )
                    {
                        std::string propName = itFindSymbolInDesc->second;
                        auto const& desc = M_materialPropertyPhysicDescription.find( propName )->second;
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
                    else
                    {
                        // we need to copy these variable elses the lambda don't capture, I don't know why
                        auto propSymbol2 = propSymbol;
                        auto propExpr2 = propExpr;
                        hana::for_each( ModelExpression::expr_shapes, [this,&propSymbol2,&propExpr2,&matProperties]( auto const& e_ij )
                                        {
                                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                                            if ( propExpr2.template hasExpr<ni,nj>() )
                                            {
                                                std::string propName = propSymbol2;
                                                matProperties.add( propName, propExpr2 );
                                                auto propShape = MaterialPropertyDescription::shape(ni,nj);
                                                auto itFindPropDesc = M_materialPropertyDescription.find( propName);
                                                if( itFindPropDesc == M_materialPropertyDescription.end() )
                                                    M_materialPropertyDescription.emplace( propName, MaterialPropertyDescription( propSymbol2, { propShape } ) );
                                                else
                                                    itFindPropDesc->second.add( propShape );
                                            }
                                        });
                    }
                }

                // TODO : move in heat toolbox
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


    //! return the number of materials
    int numberOfMaterials() const { return M_materialNameToProperties.size(); }

    //! return true if the material matName has been defined
    bool hasMaterial( std::string const& matName ) const
        {
            return (M_materialNameToProperties.find( matName ) != M_materialNameToProperties.end());
#if 0
            for ( auto const& [physic, matNames] : M_materialsNames )
                if ( matNames.find( matName ) != matNames.end() )
                    return true;
            return false;
#endif
        }

    //! return the materials names used with a set of physic
    std::set<std::string> physicToMaterials( std::set<std::string> const& physics ) const
        {
            std::set<std::string> res;
            for ( std::string const& physic : physics )
            {
                auto itFindMat = M_materialsNames.find( physic );
                if( itFindMat != M_materialsNames.end() )
                    res.insert( itFindMat->second.begin(), itFindMat->second.end() );
            }
            return res;
        }

    //! return the materials names used with a physic
    std::set<std::string> physicToMaterials( std::string const& physic ) const
        {
            return this->physicToMaterials( std::set<std::string>({ physic }) );
        }

    void markers( std::string const& p, std::set<std::string> & res ) const
        {
            auto itFindMarkers = M_markers.find( p );
            if ( itFindMarkers != M_markers.end() )
                res.insert( itFindMarkers->second.begin(), itFindMarkers->second.end() );
        }

    std::set<std::string> markers( std::string const& p ) const
        {
            std::set<std::string> res;
            this->markers( p, res );
            return res;
        }

    std::set<std::string> markers( std::set<std::string> const& setOfPhysics ) const
        {
            std::set<std::string> res;
            for ( std::string p : setOfPhysics )
                this->markers( p, res );
            return res;
        }

    bool isDefinedOnWholeMesh( std::string const& p ) const
        {
            auto itFindMarkers = M_isDefinedOnWholeMesh.find( p );
            if ( itFindMarkers == M_isDefinedOnWholeMesh.end() )
                return false;
            else
                return itFindMarkers->second;
        }


    bool isDefinedOnWholeMesh( std::set<std::string> const& setOfPhysics ) const
        {
            return this->markers( setOfPhysics ).size() == M_eltMarkersInMesh.size();
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

    template <int TheDim, typename SymbolsExpr = symbols_expression_empty_t>
    auto materialPropertyExprScalarOrMatrix( std::string const& propName, SymbolsExpr const& se = symbols_expression_empty_t{} ) const
        {
            auto Id = eye<TheDim,TheDim>();
            using _expr_scalar_type = std::decay_t< decltype( ModelExpression{}.template expr<1,1>() ) >;
            using _expr_matrix_type = std::decay_t< decltype( ModelExpression{}.template expr<TheDim,TheDim>() ) >;
            using _expr_scalar_as_matrix_type = std::decay_t< decltype(_expr_scalar_type{}*Id) >;
            std::vector<std::pair<std::string,_expr_scalar_as_matrix_type>> exprs_scalar_as_matrix;
            std::vector<std::pair<std::string,_expr_matrix_type>> exprs_matrix;

            for ( auto & [matName,matProps] : M_materialNameToProperties )
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

            if constexpr ( std::is_same_v<SymbolsExpr,symbols_expression_empty_t> )
                return expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, exprs_scalar_as_matrix, exprs_matrix );
            else
                return expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, exprs_scalar_as_matrix, exprs_matrix ).applySymbolsExpr( se );
        }

    template <int M,int N, typename SymbolsExpr = symbols_expression_empty_t>
    auto materialPropertyExpr( std::string const& propName, SymbolsExpr const& se = symbols_expression_empty_t{} ) const
        {
            using _expr_type = std::decay_t< decltype( ModelExpression{}.template expr<M,N>() ) >;
            std::vector<std::pair<std::string,_expr_type>> theExprs;
            for ( auto & [matName,matProps] : M_materialNameToProperties )
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

            if constexpr ( std::is_same_v<SymbolsExpr,symbols_expression_empty_t> )
                return expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, theExprs );
            else
                return expr<typename mesh_type::index_type>( M_exprSelectorByMeshElementMapping, theExprs ).applySymbolsExpr( se );
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

    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto thermalConductivityExpr( SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            return this->materialPropertyExprScalarOrMatrix<nDim>( "thermal-conductivity", symbolsExpr );
        }

    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto thermalConductivityScalarExpr( SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
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

    std::map<std::string,double>
    toParameterValues( std::string const& prefix_symbol ) const
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
            std::map<std::string,std::string> propertyNamesToSymbol2;
            auto tupleSymbolExprsOnAllMaterials = hana::transform( ModelExpression::expr_shapes, [this,&prefix_symbol/*,&nMat*/,&propertyNamesToSymbol,&propertyNamesToSymbol2](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                    using _expr_type =  std::decay_t< decltype( this->materialPropertyExpr<ni,nj>( "" ) ) >;
                    symbol_expression_t<_expr_type> se;
                    for ( auto const& [propName,propSymbol] : propertyNamesToSymbol )
                    {
                        if ( this->materialPropertyHasSameExprShapeInAllMaterials<ni,nj>( propName ) )
                        {
                            std::string symbolGlobalMatProp = (boost::format("%1%%2%")% prefix_symbol %propSymbol).str();
                            auto globalMatPropExpr = this->materialPropertyExpr<ni,nj>( propName );
                            se.add( symbolGlobalMatProp, globalMatPropExpr, SymbolExprComponentSuffix( ni,nj ) );
                        }
                        else
                            propertyNamesToSymbol2[propName] = propSymbol;
                    }
                    return se;
                });

            // generate symbols heat_k_xx, heat_k_xy,... if a prop is scalar in one mat and matrix in another mat
            typedef decltype( this->materialPropertyExprScalarOrMatrix<nDim>( "" ) ) _expr_scalar_or_matrix_selector_type;
            symbol_expression_t<_expr_scalar_or_matrix_selector_type> seScalarOrMatrixSelector;
            for ( auto const& [propName,propSymbol] : propertyNamesToSymbol2 )
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

        template <typename SymbExprType>
        auto exprPostProcessExports( std::set<std::string> const& physics, SymbExprType const& se, std::string const& prefix = "materials" ) const
        {
            auto setOfMatNameUsed = this->physicToMaterials( physics );
            std::map<std::string,std::map<std::string,const MaterialProperty* >> propertiesNotHaveSameShape;

            auto tupleExportsExprs = hana::transform( ModelExpression::expr_shapes, [this,&setOfMatNameUsed,&se,&prefix,&propertiesNotHaveSameShape](auto const& e_ij) {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;

                    using _expr_type = std::decay_t< decltype( expr( ModelExpression{}.template expr<ni,nj>(), se ) ) >;
                    std::map<std::string,std::vector<std::tuple<_expr_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprExported;

                    for ( auto const& [matName,matProps] : M_materialNameToProperties )
                    {
                        if ( setOfMatNameUsed.find( matName ) == setOfMatNameUsed.end() )
                            continue;
                        auto const& range = this->rangeMeshElementsByMaterial( matName );
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
            std::map<std::string,std::vector<std::tuple<_expr_tensor2_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprTensor2;

            auto Id = eye<nDim,nDim>();
            typedef decltype(expr(typename ModelExpression::expr_scalar_type{},se)*Id) _expr_tensor2_from_scalar_type;
            std::map<std::string,std::vector<std::tuple<_expr_tensor2_from_scalar_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprTensor2FromScalar;

            for ( auto const& [matName,matProps] : propertiesNotHaveSameShape )
            {
                auto const& range = this->rangeMeshElementsByMaterial( matName );
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

    std::set<std::string> postProcessExportsAllFieldsAvailable( std::set<std::string> const& physics, std::string const& prefix = "materials" ) const
        {
            std::set<std::string> res;
            auto setOfMatNameUsed = this->physicToMaterials( physics );
            for ( auto const& [matName,matProps] : M_materialNameToProperties )
            {
                if ( setOfMatNameUsed.find( matName ) == setOfMatNameUsed.end() )
                    continue;
                auto range = this->rangeMeshElementsByMaterial( matName );
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

private :
    std::string M_exprRepository;
    std::set<std::string> M_eltMarkersInMesh;
    std::map<std::string,std::set<std::string>> M_markers; // physic -> markers
    std::map<std::string,bool> M_isDefinedOnWholeMesh; // physics -> bool
    std::map<std::string,std::set<std::string>> M_materialsNames; // physic -> matNames
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial; // matName -> range
    std::shared_ptr<ExprSelectorByMeshElementMapping<typename mesh_type::index_type>> M_exprSelectorByMeshElementMapping;

    std::map<std::string, ModelExpressionScalar> M_rhoHeatCapacityByMaterial;

    std::map<std::string,MaterialProperties<mesh_type>> M_materialNameToProperties;
    std::map<std::string,MaterialPropertyDescription> M_materialPropertyPhysicDescription; // name -> (symbol, shapes.. )   only defined in phycis
    std::map<std::string,MaterialPropertyDescription> M_materialPropertyDescription; // name -> (symbol, shapes.. )   all props read
};


} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
