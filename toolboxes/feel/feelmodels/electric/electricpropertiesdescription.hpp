/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_MODELS_ELECTRICPROPERTIES_DESCRIPTION_H
#define FEELPP_MODELS_ELECTRICPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelmodels/modelexpression.hpp>

namespace Feel
{
namespace FeelModels
{

template<class SpaceType>
class ElectricPropertiesDescription
{
    typedef ElectricPropertiesDescription<SpaceType> self_type;
public :
    typedef SpaceType space_type;
    typedef std::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    ElectricPropertiesDescription( std::string const& prefix )
        :
        M_isDefinedOnWholeMesh( true ),
        M_electricConductivityDefaultValue( doption(_name="electric-conductivity",_prefix=prefix) )
        {}

    ElectricPropertiesDescription( ElectricPropertiesDescription const& ) = default;

    void updateForUse( mesh_ptrtype const& mesh , ModelMaterials const& mats, worldscomm_ptr_t const& worldsComm )
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
                if ( mat.hasPhysics() && !mat.hasPhysics( { "electric","thermo-electric" } ) )
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
                M_space = space_type::New(_mesh=mesh, _worldscomm=worldsComm );
            else
                M_space = space_type::New(_mesh=mesh, _worldscomm=worldsComm,_range=markedelements(mesh,M_markers) );
            M_fieldElectricConductivity = M_space->elementPtr( vf::cst( this->cstElectricConductivity() ) );

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

                M_electricConductivityByMaterial[matName];
                if ( mat.hasPropertyExprScalar("sigma") )
                {
                    auto const& expr = mat.propertyExprScalar("sigma");
                    M_electricConductivityByMaterial[matName].setExpr( expr );
                    M_fieldElectricConductivity->on(_range=range,_expr=expr);
                }
            }
        }

    std::set<std::string> const& markers() const { return M_markers; }
    bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }

    bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

    std::map<std::string, ModelExpressionScalar> const& electricConductivityByMaterial() const { return M_electricConductivityByMaterial; }

    double cstElectricConductivity( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_electricConductivityByMaterial.empty() )
                    return M_electricConductivityDefaultValue;
                else
                    return M_electricConductivityByMaterial.begin()->second.value();
            }
            auto itFindMat = M_electricConductivityByMaterial.find( matName );
            CHECK( itFindMat != M_electricConductivityByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
        }

    element_type const& fieldElectricConductivity() const { return *M_fieldElectricConductivity; }
    element_ptrtype const& fieldElectricConductivityPtr() const { return M_fieldElectricConductivity; }

    bool hasElectricConductivity( std::string const& matName ) const
        {
            return M_electricConductivityByMaterial.find( matName ) != M_electricConductivityByMaterial.end();
        }
    ModelExpressionScalar const& electricConductivity( std::string const& matName ) const
        {
            CHECK( this->hasElectricConductivity( matName ) ) << "material name not registered : " << matName;
            return M_electricConductivityByMaterial.find( matName )->second;
        }

    std::shared_ptr<std::ostringstream>
    getInfoMaterialParameters() const
        {
            std::shared_ptr<std::ostringstream> ostr( new std::ostringstream() );
            *ostr << "\n   Materials parameters";
            std::set<std::string> matNames;
            for ( auto const& matRange : M_rangeMeshElementsByMaterial)
                matNames.insert( matRange.first );
            // if ( matNames.empty() )
            //     matNames.insert( std::string("") );
            *ostr << "\n     -- number of materials : " << matNames.size();
            for ( std::string const& matName : matNames)
            {
                if ( this->electricConductivity( matName ).hasExpr() )
                {
                    *ostr << "\n     -- [" << matName << "] electric conductivity : ";
                    if ( this->electricConductivity(matName).isConstant() )
                        *ostr << this->electricConductivity(matName).value();
                    else
                        *ostr << str( this->electricConductivity(matName).expr().expression() );
                }
            }
            return ostr;
        }

    bool hasElectricConductivityDependingOnSymbol( std::string const& symbolStr ) const
        {
            for ( auto const& conductivityData : M_electricConductivityByMaterial )
            {
                auto const& electricConductivity = conductivityData.second;
                if ( electricConductivity.isConstant() )
                    continue;
                if ( electricConductivity.expr().expression().hasSymbol( symbolStr ) )
                    return true;
            }
            return false;
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & prop : M_electricConductivityByMaterial )
                prop.second.setParameterValues( mp );
        }


private :
    std::set<std::string> M_markers;
    bool M_isDefinedOnWholeMesh;
    space_ptrtype M_space;
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;
    std::map<std::string, ModelExpressionScalar> M_electricConductivityByMaterial;
    element_ptrtype M_fieldElectricConductivity;
    double M_electricConductivityDefaultValue;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_ELECTRICPROPERTIES_DESCRIPTION_H
