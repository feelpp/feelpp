#ifndef FEELPP_MODELS_MAXWELLPROPERTIES_DESCRIPTION_H
#define FEELPP_MODELS_MAXWELLPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelmodels/modelexpression.hpp>

namespace Feel
{
namespace FeelModels
{

template<class SpaceType>
class MaxwellPropertiesDescription
{
    typedef MaxwellPropertiesDescription<SpaceType> self_type;
public :
    typedef SpaceType space_type;
    typedef std::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    explicit MaxwellPropertiesDescription( std::string const& prefix )
        :
        M_isDefinedOnWholeMesh( true ),
        M_magneticPermeabilityDefaultValue( doption(_name="magnetic-permeability",_prefix=prefix) )
        {}

    MaxwellPropertiesDescription( MaxwellPropertiesDescription const& ) = default;

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
                if ( mat.hasPhysics() && !mat.hasPhysics( { "magneto","magnetostatic" } ) )
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
            M_fieldMagneticPermeability = M_space->elementPtr( vf::cst( this->cstMagneticPermeability() ) );

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

                M_magneticPermeabilityByMaterial[matName];
                if ( mat.hasProperty("mu") )
                {
                    auto const& expr = mat.property("mu").template expr<1,1>();
                    M_magneticPermeabilityByMaterial[matName].setExpr( expr );
                    M_fieldMagneticPermeability->on(_range=range,_expr=expr);
                }
            }
        }

    std::set<std::string> const& markers() const { return M_markers; }
    bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

    std::map<std::string, Range<mesh_type,MESH_ELEMENTS> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }

    bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

    std::map<std::string, ModelExpressionScalar> const& magneticPermeabilityByMaterial() const { return M_magneticPermeabilityByMaterial; }

    double cstMagneticPermeability( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_magneticPermeabilityByMaterial.empty() )
                    return M_magneticPermeabilityDefaultValue;
                else
                    return M_magneticPermeabilityByMaterial.begin()->second.value();
            }
            auto itFindMat = M_magneticPermeabilityByMaterial.find( matName );
            CHECK( itFindMat != M_magneticPermeabilityByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
        }

    element_type const& fieldMagneticPermeability() const { return *M_fieldMagneticPermeability; }
    element_ptrtype const& fieldMagneticPermeabilityPtr() const { return M_fieldMagneticPermeability; }

    bool hasMagneticPermeability( std::string const& matName ) const
        {
            return M_magneticPermeabilityByMaterial.find( matName ) != M_magneticPermeabilityByMaterial.end();
        }
    ModelExpressionScalar const& magneticPermeability( std::string const& matName ) const
        {
            CHECK( this->hasMagneticPermeability( matName ) ) << "material name not registered : " << matName;
            return M_magneticPermeabilityByMaterial.find( matName )->second;
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
                *ostr << "\n     -- [" << matName << "] magnetic permeability : ";
                if ( this->magneticPermeability(matName).isConstant() )
                    *ostr << this->magneticPermeability(matName).value();
                else
                    *ostr << str( this->magneticPermeability(matName).expr().expression() );
            }
            return ostr;
        }

    bool hasMagneticPermeabilityDependingOnSymbol( std::string const& symbolStr ) const
        {
            for ( auto const& conductivityData : M_magneticPermeabilityByMaterial )
            {
                auto const& magneticPermeability = conductivityData.second;
                if ( magneticPermeability.isConstant() )
                    continue;
                if ( magneticPermeability.expr().expression().hasSymbol( symbolStr ) )
                    return true;
            }
            return false;
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            for ( auto & prop : M_magneticPermeabilityByMaterial )
                prop.second.setParameterValues( mp );
        }


private :
    std::set<std::string> M_markers;
    bool M_isDefinedOnWholeMesh;
    space_ptrtype M_space;
    std::map<std::string, Range<mesh_type,MESH_ELEMENTS> > M_rangeMeshElementsByMaterial;
    std::map<std::string, ModelExpressionScalar> M_magneticPermeabilityByMaterial;
    element_ptrtype M_fieldMagneticPermeability;
    double M_magneticPermeabilityDefaultValue;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_MAXWELLPROPERTIES_DESCRIPTION_H
