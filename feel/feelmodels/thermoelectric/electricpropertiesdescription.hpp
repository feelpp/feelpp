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
        static std::string defaultMaterialName() { return std::string("FEELPP_DEFAULT_MATERIAL_NAME"); }
        ElectricPropertiesDescription( space_ptrtype const& space, std::string const& prefix )
            {
                M_cstElectricConductivity[self_type::defaultMaterialName()] = doption(_name="electric-conductivity",_prefix=prefix);
                this->initFromSpace( space );
            }

        ElectricPropertiesDescription( std::string const& prefix )
            {
                M_cstElectricConductivity[self_type::defaultMaterialName()] = doption(_name="electric-conductivity",_prefix=prefix);
            }

        ElectricPropertiesDescription( ElectricPropertiesDescription const& app  ) = default;

        void initFromSpace( space_ptrtype const& space )
        {
            M_space = space;
            M_fieldElectricConductivity = space->elementPtr( vf::cst( this->cstElectricConductivity() ) );
        }

        double cstElectricConductivity( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstElectricConductivity.find( markerUsed );
            CHECK( itFindMarker != M_cstElectricConductivity.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }

        std::set<std::string> const& markers() const { return M_markers; }

        element_type const& fieldElectricConductivity() const { return *M_fieldElectricConductivity; }
        element_ptrtype const& fieldElectricConductivityPtr() const { return M_fieldElectricConductivity; }

        void setCstElectricConductivity( double val, std::string const& marker = "" )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstElectricConductivity[markerUsed]=val;
            this->updateElectricConductivity( cst(val), marker );
        }

        template < typename ExprT >
        void updateElectricConductivity(vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldElectricConductivity ) return;
            if ( marker.empty() )
                M_fieldElectricConductivity->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldElectricConductivity->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
        }

        void updateFromModelMaterials( ModelMaterials const& mat )
        {
            if ( mat.empty() ) return;

            for( auto const& m : mat )
            {
                auto const& mat = m.second;
                auto const& matmarker = m.first;
                //LOG(INFO) << "set material " << mat.name() << " associated to marker : " << matmarker<< "\n";
                M_markers.insert( matmarker );
                this->setCstElectricConductivity( mat.sigma(), matmarker );
            }
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
