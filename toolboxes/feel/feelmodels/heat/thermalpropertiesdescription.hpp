/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_THERMALPROPERTIES_DESCRIPTION_H
#define FEELPP_TOOLBOXES_THERMALPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelmodels/modelexpression.hpp>

namespace Feel
{
namespace FeelModels
{

template<class SpaceType>
class ThermalPropertiesDescription
{
    typedef ThermalPropertiesDescription<SpaceType> self_type;
public :
    typedef SpaceType space_type;
    typedef std::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
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
            if ( M_isDefinedOnWholeMesh )
                M_space = space_type::New(_mesh=mesh );
            else
                M_space = space_type::New(_mesh=mesh,_range=markedelements(mesh,M_markers) );
            M_fieldThermalConductivity = M_space->elementPtr( vf::cst( this->cstThermalConductivity() ) );
            M_fieldHeatCapacity = M_space->elementPtr( vf::cst( this->cstHeatCapacity() ) );
            M_fieldRho = M_space->elementPtr( vf::cst( this->cstRho() ) );
            M_fieldThermalExpansion = M_space->elementPtr( vf::cst( this->cstThermalExpansion() ) );

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
                    if ( prop.template hasExprMatrix<nDim,nDim>() )
                    {
                        //M_fieldThermalConductivity->on(_range=range,_expr=expr);
                    }
                    else if ( prop.hasExprScalar() )
                    {
                        auto const& expr = prop.expr();
                        M_fieldThermalConductivity->on(_range=range,_expr=expr);
                    }
                }

                M_rhoByMaterial[matName];
                if ( mat.hasPropertyExprScalar("rho") )
                {
                    auto const& expr = mat.propertyExprScalar("rho");
                    M_rhoByMaterial[matName].setExpr( expr );
                    M_fieldRho->on(_range=range,_expr=expr);
                }
                else M_rhoByMaterial[matName].setExpr( expr("0") );

                M_heatCapacityByMaterial[matName];
                if ( mat.hasPropertyExprScalar("Cp") )
                {
                    auto const& expr = mat.propertyExprScalar("Cp");
                    M_heatCapacityByMaterial[matName].setExpr( expr );
                    M_fieldHeatCapacity->on(_range=range,_expr=expr);
                }
                else M_heatCapacityByMaterial[matName].setExpr( expr("0") );

                M_thermalExpansionByMaterial[matName];
                if ( mat.hasPropertyExprScalar("beta") )
                {
                    auto const& expr = mat.propertyExprScalar("beta");
                    M_thermalExpansionByMaterial[matName].setExpr( expr );
                    M_fieldThermalExpansion->on(_range=range,_expr=expr);
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
        }

    std::set<std::string> const& markers() const { return M_markers; }

    bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }

    bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

    element_type const& fieldThermalConductivity() const { return *M_fieldThermalConductivity; }
    element_type const& fieldHeatCapacity() const { return *M_fieldHeatCapacity; }
    element_type const& fieldRho() const { return *M_fieldRho; }
    element_type const& fieldThermalExpansion() const { return *M_fieldThermalExpansion; }
    element_ptrtype const& fieldThermalConductivityPtr() const { return M_fieldThermalConductivity; }
    element_ptrtype const& fieldHeatCapacityPtr() const { return M_fieldHeatCapacity; }
    element_ptrtype const& fieldRhoPtr() const { return M_fieldRho; }
    element_ptrtype const& fieldThermalExpansionPtr() const { return M_fieldThermalExpansion; }

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
    double cstThermalConductivity( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_thermalConductivityByMaterial.empty() )
                    return M_thermalConductivityDefaultValue;
                else
                    return M_thermalConductivityByMaterial.begin()->second.value();
            }
            auto itFindMat = M_thermalConductivityByMaterial.find( matName );
            CHECK( itFindMat != M_thermalConductivityByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
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
    double cstRho( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_rhoByMaterial.empty() )
                    return M_rhoDefaultValue;
                else
                    return M_rhoByMaterial.begin()->second.value();
            }
            auto itFindMat = M_rhoByMaterial.find( matName );
            CHECK( itFindMat != M_rhoByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
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
    double cstHeatCapacity( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_heatCapacityByMaterial.empty() )
                    return M_heatCapacityDefaultValue;
                else
                    return M_heatCapacityByMaterial.begin()->second.value();
            }
            auto itFindMat = M_heatCapacityByMaterial.find( matName );
            CHECK( itFindMat != M_heatCapacityByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
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
    double cstThermalExpansion( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_thermalExpansionByMaterial.empty() )
                    return M_thermalExpansionDefaultValue;
                else
                    return M_thermalExpansionByMaterial.begin()->second.value();
            }
            auto itFindMat = M_thermalExpansionByMaterial.find( matName );
            CHECK( itFindMat != M_thermalExpansionByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
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

private :
    std::string M_exprRepository;
    std::set<std::string> M_markers;
    bool M_isDefinedOnWholeMesh;
    space_ptrtype M_space;
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;

    std::map<std::string, ModelExpressionScalar> /*M_thermalConductivityByMaterial,*/ M_heatCapacityByMaterial, M_rhoByMaterial, M_thermalExpansionByMaterial, M_rhoHeatCapacityByMaterial;
    std::map<std::string, ModelExpression> M_thermalConductivityByMaterial;
    element_ptrtype M_fieldThermalConductivity, M_fieldHeatCapacity, M_fieldRho, M_fieldThermalExpansion;
    double M_thermalConductivityDefaultValue, M_heatCapacityDefaultValue, M_rhoDefaultValue, M_thermalExpansionDefaultValue;
};


} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
