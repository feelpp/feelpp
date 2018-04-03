/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_MODELS_ELECTRICPROPERTIES_DESCRIPTION_H
#define FEELPP_MODELS_ELECTRICPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>

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
        typedef boost::shared_ptr<SpaceType> space_ptrtype;
        typedef typename SpaceType::element_type element_type;
        typedef boost::shared_ptr<element_type> element_ptrtype;
        typedef typename space_type::mesh_type mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        static std::string defaultMaterialName() { return std::string("FEELPP_DEFAULT_MATERIAL_NAME"); }

        ElectricPropertiesDescription( std::string const& prefix )
            :
            M_isDefinedOnWholeMesh( true )
            {
                M_cstElectricConductivity[self_type::defaultMaterialName()] = doption(_name="electric-conductivity",_prefix=prefix);
            }

        ElectricPropertiesDescription( ElectricPropertiesDescription const& app ) = default;

        void updateForUse( mesh_ptrtype const& mesh , ModelMaterials const& mat, std::vector<WorldComm> const& worldsComm )
            {
                std::set<std::string> eltMarkersInMesh;
                for (auto const& markPair : mesh->markerNames() )
                {
                    std::string meshMarker = markPair.first;
                    if ( mesh->hasElementMarker( meshMarker ) )
                        eltMarkersInMesh.insert( meshMarker );
                }

                M_markers.clear();
                for( auto const& m : mat )
                {
                    auto const& matmarker = m.first;
                    if ( eltMarkersInMesh.find( matmarker ) == eltMarkersInMesh.end() )
                        continue;
                    auto const& mat = m.second;
                    std::string matphysics = ( mat.physics().empty() )? "electric" : mat.physics();
                    if ( ( matphysics != "electric" ) && ( matphysics != "thermo-electric" ) )
                        continue;
                    M_markers.insert( matmarker );
                }

                M_isDefinedOnWholeMesh = ( M_markers.size() == eltMarkersInMesh.size() );
                if ( M_isDefinedOnWholeMesh )
                    M_space = space_type::New(_mesh=mesh, _worldscomm=worldsComm );
                else
                    M_space = space_type::New(_mesh=mesh, _worldscomm=worldsComm,_range=markedelements(mesh,M_markers) );
                M_fieldElectricConductivity = M_space->elementPtr( vf::cst( this->cstElectricConductivity() ) );

                for( auto const& m : mat )
                {
                    auto const& mat = m.second;
                    auto const& matmarker = m.first;
                    if ( M_markers.find( matmarker ) == M_markers.end() )
                        continue;

                    if ( mat.hasPropertyExprScalar("sigma") )
                        this->setElectricConductivity( mat.propertyExprScalar("sigma"),matmarker );
                    else
                        this->setCstElectricConductivity( mat.propertyConstant("sigma"), matmarker );
                }

            }

        double cstElectricConductivity( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstElectricConductivity.find( markerUsed );
            CHECK( itFindMarker != M_cstElectricConductivity.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }

        std::set<std::string> const& markers() const { return M_markers; }

        bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

        element_type const& fieldElectricConductivity() const { return *M_fieldElectricConductivity; }
        element_ptrtype const& fieldElectricConductivityPtr() const { return M_fieldElectricConductivity; }

        void setCstElectricConductivity( double val, std::string const& marker = "", bool update = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstElectricConductivity[markerUsed]=val;
            if ( update )
                this->updateElectricConductivity( cst(val), marker );
        }
        template < typename ExprT >
        void setElectricConductivity( vf::Expr<ExprT> const& vfexpr, std::string const& marker = "" )
            {
                this->updateElectricConductivity( vfexpr,marker);
                if ( M_fieldElectricConductivity )
                    this->setCstElectricConductivity( M_fieldElectricConductivity->min(), marker, false );
            }

        template < typename ExprT >
        void updateElectricConductivity(vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldElectricConductivity ) return;
            auto rangeEltUsed = (marker.empty())? markedelements( M_space->mesh(), M_markers ) : markedelements( M_space->mesh(),marker );
            M_fieldElectricConductivity->on(_range=rangeEltUsed,_expr=__expr);
        }

        boost::shared_ptr<std::ostringstream>
        getInfoMaterialParameters() const
        {
            boost::shared_ptr<std::ostringstream> ostr( new std::ostringstream() );
            *ostr << "\n   Materials parameters";
            std::set<std::string> matmarkerset = this->markers();
            if (  matmarkerset.empty() ) matmarkerset.insert(std::string(""));
            *ostr << "\n     -- number of materials : " << matmarkerset.size();
            for ( std::string const& matmarker : matmarkerset)
            {
                std::string matmarkertag = matmarker.empty()? std::string("") : (boost::format("[%1%] ")%matmarker).str();
                *ostr << "\n     -- " << matmarkertag << "electric conductivity : " << this->cstElectricConductivity(matmarker);
            }
            return ostr;
        }

    private :
        std::map<std::string,double> M_cstElectricConductivity;// [ W/(m*K) ]
        std::set<std::string> M_markers;
        bool M_isDefinedOnWholeMesh;
        space_ptrtype M_space;
        element_ptrtype M_fieldElectricConductivity;
    };


    template<class SpaceType>
    ElectricPropertiesDescription<SpaceType>
    electricPropertiesDesc( boost::shared_ptr<SpaceType> const& space, std::string const& prefix )
    {
        ElectricPropertiesDescription<SpaceType> res(space,prefix);
        return res;
    }

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_ELECTRICPROPERTIES_DESCRIPTION_H
