/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef __THERMALPROPERTIES_DESCRIPTION_H
#define __THERMALPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>

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
        typedef boost::shared_ptr<SpaceType> space_ptrtype;
        typedef typename SpaceType::element_type element_type;
        typedef boost::shared_ptr<element_type> element_ptrtype;
        static std::string defaultMaterialName() { return std::string("FEELPP_DEFAULT_MATERIAL_NAME"); }
        ThermalPropertiesDescription( space_ptrtype const& space, std::string const& prefix )
            {
                M_cstThermalConductivity[self_type::defaultMaterialName()] = doption(_name="thermal-conductivity",_prefix=prefix);
                M_cstHeatCapacity[self_type::defaultMaterialName()] = doption(_name="heat-capacity",_prefix=prefix);
                M_cstRho[self_type::defaultMaterialName()] = doption(_name="rho",_prefix=prefix);
                M_cstThermalExpansion[self_type::defaultMaterialName()] = doption(_name="thermal-expansion",_prefix=prefix);
                this->initFromSpace( space );
            }

        ThermalPropertiesDescription( std::string const& prefix )
            {
                M_cstThermalConductivity[self_type::defaultMaterialName()] = doption(_name="thermal-conductivity",_prefix=prefix);
                M_cstHeatCapacity[self_type::defaultMaterialName()] = doption(_name="heat-capacity",_prefix=prefix);
                M_cstRho[self_type::defaultMaterialName()] = doption(_name="rho",_prefix=prefix);
                M_cstThermalExpansion[self_type::defaultMaterialName()] = doption(_name="thermal-expansion",_prefix=prefix);
            }

        ThermalPropertiesDescription( ThermalPropertiesDescription const& app  ) = default;

        void initFromSpace( space_ptrtype const& space )
        {
            M_space = space;
            M_fieldThermalConductivity = space->elementPtr( vf::cst( this->cstThermalConductivity() ) );
            M_fieldHeatCapacity = space->elementPtr( vf::cst( this->cstHeatCapacity() ) );
            M_fieldRho = space->elementPtr( vf::cst( this->cstRho() ) );
            M_fieldThermalExpansion = space->elementPtr( vf::cst( this->cstThermalExpansion() ) );
        }

        double cstThermalConductivity( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstThermalConductivity.find( markerUsed );
            CHECK( itFindMarker != M_cstThermalConductivity.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }
        double cstHeatCapacity( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstHeatCapacity.find( markerUsed );
            CHECK( itFindMarker != M_cstHeatCapacity.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }
        double cstRho( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstRho.find( markerUsed );
            CHECK( itFindMarker != M_cstRho.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }
        double cstThermalExpansion( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstThermalExpansion.find( markerUsed );
            CHECK( itFindMarker != M_cstThermalExpansion.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }

        std::set<std::string> const& markers() const { return M_markers; }

        element_type const& fieldThermalConductivity() const { return *M_fieldThermalConductivity; }
        element_type const& fieldHeatCapacity() const { return *M_fieldHeatCapacity; }
        element_type const& fieldRho() const { return *M_fieldRho; }
        element_type const& fieldThermalExpansion() const { return *M_fieldThermalExpansion; }
        element_ptrtype const& fieldThermalConductivityPtr() const { return M_fieldThermalConductivity; }
        element_ptrtype const& fieldHeatCapacityPtr() const { return M_fieldHeatCapacity; }
        element_ptrtype const& fieldRhoPtr() const { return M_fieldRho; }
        element_ptrtype const& fieldThermalExpansionPtr() const { return M_fieldThermalExpansion; }

        void setCstThermalConductivity( double val, std::string const& marker = "" )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstThermalConductivity[markerUsed]=val;
            this->updateThermalConductivity( cst(val), marker );
        }
        void setCstHeatCapacity( double val, std::string const& marker = "" )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstHeatCapacity[markerUsed]=val;
            this->updateHeatCapacity( cst(val), marker );
        }
        void setCstRho( double val, std::string const& marker = "" )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstRho[markerUsed]=val;
            this->updateRho( vf::cst(val), marker );
        }
        void setCstThermalExpansion( double val, std::string const& marker = "" )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstThermalExpansion[markerUsed]=val;
            this->updateThermalExpansion( vf::cst(val), marker );
        }

        template < typename ExprT >
        void updateThermalConductivity(vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldThermalConductivity ) return;
            if ( marker.empty() )
                M_fieldThermalConductivity->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldThermalConductivity->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
        }
        template < typename ExprT >
        void updateHeatCapacity(vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldHeatCapacity ) return;
            if ( marker.empty() )
                M_fieldHeatCapacity->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldHeatCapacity->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
        }

        template < typename ExprT >
        void updateRho(vf::Expr<ExprT> const& __expr, std::string const& marker = "")
        {
            if ( !M_fieldRho ) return;
            if ( marker.empty() )
                M_fieldRho->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldRho->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
        }

        template < typename ExprT >
        void updateThermalExpansion(vf::Expr<ExprT> const& __expr, std::string const& marker = "")
        {
            if ( !M_fieldThermalExpansion ) return;
            if ( marker.empty() )
                M_fieldThermalExpansion->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldThermalExpansion->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
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
                this->setCstRho( mat.rho(), matmarker );
                this->setCstThermalConductivity( mat.k11(), matmarker );
                this->setCstHeatCapacity( mat.Cp(), matmarker );
                this->setCstThermalExpansion( mat.beta(), matmarker );
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
                *ostr << "\n     -- " << matmarkertag << "rho : " << this->cstRho(matmarker)
                      << "\n     -- " << matmarkertag << "thermal conductivity : " << this->cstThermalConductivity(matmarker)
                      << "\n     -- " << matmarkertag << "heat capacity : " << this->cstHeatCapacity(matmarker)
                      << "\n     -- " << matmarkertag << "thermal expansion : " << this->cstThermalExpansion(matmarker);
            }
            return ostr;
        }

    private :
        std::map<std::string,double> M_cstThermalConductivity;// [ W/(m*K) ]
        std::map<std::string,double> M_cstHeatCapacity;// [ J/(kg*K) ]
        std::map<std::string,double> M_cstRho;
        std::map<std::string,double> M_cstThermalExpansion;// [ 1/K ]
        std::set<std::string> M_markers;
        space_ptrtype M_space;
        element_ptrtype M_fieldThermalConductivity, M_fieldHeatCapacity, M_fieldRho;
        element_ptrtype M_fieldThermalExpansion;
    };


    template<class SpaceType>
    ThermalPropertiesDescription<SpaceType>
    thermalPropertiesDesc( boost::shared_ptr<SpaceType> const& space, std::string const& prefix )
    {
        ThermalPropertiesDescription<SpaceType> res(space,prefix);
        return res;
    }

} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
