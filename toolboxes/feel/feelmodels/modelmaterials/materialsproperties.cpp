/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

namespace Feel
{
namespace FeelModels
{

    MaterialProperties::MaterialProperties( MaterialProperties const& matProp ):
        super_type( matProp ),
        M_materialName( matProp.M_materialName ),
        M_markers( matProp.M_markers ),
        M_materialLaws( matProp.M_materialLaws )
    {
        for ( auto const& [subMatName, subMatProps]: matProp.M_subMaterialProperties )
        {
            M_subMaterialProperties.emplace( subMatName, std::make_unique<MaterialProperties>( *subMatProps ) );
        }
    }

    void MaterialProperties::add( std::string const& propName, ModelExpression const& expr )
    {
        auto itFind = this->find( propName );
        if ( itFind != this->end() )
            this->erase( itFind );
        this->emplace(propName, MaterialProperty(propName, expr) );
    }

    bool MaterialProperties::has( std::string const& propName ) const
    {
        if ( this->find( propName ) != this->end() )
            return true;
        return false;
    }

    MaterialProperty const&
    MaterialProperties::property( std::string const& propName ) const
    {
        auto itFind = this->find( propName );
        CHECK( itFind != this->end() ) << "material property" << propName << " not found";
        return itFind->second;
    }

    void MaterialProperties::setParameterValues( std::map<std::string,double> const& mp )
    {
        for ( auto & [propName,matProp] : *this )
            matProp.setParameterValues( mp );
    }

    std::shared_ptr<ModelMaterialLaw> MaterialProperties::law( std::string const& l ) const
    {
        auto itFind = M_materialLaws.find( l );
        CHECK( itFind != M_materialLaws.end() ) << "material law " << l << " not found";
        return itFind->second;
    }

    MaterialProperties const& MaterialProperties::subMaterialProperties( std::string const& subMatName ) const
    {
        auto itFind = M_subMaterialProperties.find( subMatName );
        CHECK( itFind != M_subMaterialProperties.end() ) << "sub material " << subMatName << " not found";
        return *(itFind->second);
    }
    MaterialProperties & MaterialProperties::subMaterialProperties( std::string const& subMatName )
    {
        auto [it, inserted] = M_subMaterialProperties.try_emplace( subMatName, nullptr );
        if( inserted )
            it->second = std::make_unique<MaterialProperties>( subMatName );
        return *it->second;
    }
    void MaterialProperties::addSubMaterialProperties( std::string const& subMatName, MaterialProperties const& subMatProps ) 
    {
        auto [it, inserted] = M_subMaterialProperties.try_emplace( subMatName, nullptr );
        CHECK( inserted ) << "sub material " << subMatName << " already exists";
        if( inserted )
            it->second = std::make_unique<MaterialProperties>( subMatProps );
    }
    void MaterialProperties::addSubMaterialProperty( std::string const& subMatName, std::string const& propName, ModelExpression const& expr )
    {
        M_subMaterialProperties[subMatName]->add( propName, expr );
    }

} // ns FeelModels
} // ns Feel
