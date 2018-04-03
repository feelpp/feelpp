/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_THERMALPROPERTIES_DESCRIPTION_H
#define FEELPP_TOOLBOXES_THERMALPROPERTIES_DESCRIPTION_H 1

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
        typedef typename space_type::mesh_type mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        static std::string defaultMaterialName() { return std::string("FEELPP_DEFAULT_MATERIAL_NAME"); }

        ThermalPropertiesDescription( std::string const& prefix )
            {
                M_cstThermalConductivity[self_type::defaultMaterialName()] = doption(_name="thermal-conductivity",_prefix=prefix);
                M_cstHeatCapacity[self_type::defaultMaterialName()] = doption(_name="heat-capacity",_prefix=prefix);
                M_cstRho[self_type::defaultMaterialName()] = doption(_name="rho",_prefix=prefix);
                M_cstThermalExpansion[self_type::defaultMaterialName()] = doption(_name="thermal-expansion",_prefix=prefix);
            }

        ThermalPropertiesDescription( ThermalPropertiesDescription const& app  ) = default;

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
                    std::string matphysics = ( mat.physics().empty() )? "heat-transfert" : mat.physics();
                    if ( ( matphysics != "heat-transfert" ) && ( matphysics != "aerothermal" ) && ( matphysics != "thermo-electric" ) )
                        continue;
                    M_markers.insert( matmarker );
                }

                M_isDefinedOnWholeMesh = ( M_markers.size() == eltMarkersInMesh.size() );
                if ( M_isDefinedOnWholeMesh )
                    M_space = space_type::New(_mesh=mesh, _worldscomm=worldsComm );
                else
                    M_space = space_type::New(_mesh=mesh, _worldscomm=worldsComm,_range=markedelements(mesh,M_markers) );
                M_fieldThermalConductivity = M_space->elementPtr( vf::cst( this->cstThermalConductivity() ) );
                M_fieldHeatCapacity = M_space->elementPtr( vf::cst( this->cstHeatCapacity() ) );
                M_fieldRho = M_space->elementPtr( vf::cst( this->cstRho() ) );
                M_fieldThermalExpansion = M_space->elementPtr( vf::cst( this->cstThermalExpansion() ) );

                for( auto const& m : mat )
                {
                    auto const& mat = m.second;
                    auto const& matmarker = m.first;
                    if ( M_markers.find( matmarker ) == M_markers.end() )
                        continue;

                    if ( mat.hasPropertyExprScalar("rho") )
                        this->setRho( mat.propertyExprScalar("rho"),matmarker );
                    else
                        this->setCstRho( mat.propertyConstant("rho"), matmarker );

                    if ( mat.hasPropertyExprScalar("k11") )
                        this->setRho( mat.propertyExprScalar("k11"),matmarker );
                    else
                        this->setCstThermalConductivity( mat.propertyConstant("k11"), matmarker );

                    if ( mat.hasPropertyExprScalar("Cp") )
                        this->setRho( mat.propertyExprScalar("Cp"),matmarker );
                    else
                        this->setCstHeatCapacity( mat.propertyConstant("Cp"), matmarker );

                    if ( mat.hasPropertyExprScalar("beta") )
                        this->setRho( mat.propertyExprScalar("beta"),matmarker );
                    else
                        this->setCstThermalExpansion( mat.propertyConstant("beta"), matmarker );
                }
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

        bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

        element_type const& fieldThermalConductivity() const { return *M_fieldThermalConductivity; }
        element_type const& fieldHeatCapacity() const { return *M_fieldHeatCapacity; }
        element_type const& fieldRho() const { return *M_fieldRho; }
        element_type const& fieldThermalExpansion() const { return *M_fieldThermalExpansion; }
        element_ptrtype const& fieldThermalConductivityPtr() const { return M_fieldThermalConductivity; }
        element_ptrtype const& fieldHeatCapacityPtr() const { return M_fieldHeatCapacity; }
        element_ptrtype const& fieldRhoPtr() const { return M_fieldRho; }
        element_ptrtype const& fieldThermalExpansionPtr() const { return M_fieldThermalExpansion; }

        void setCstThermalConductivity( double val, std::string const& marker = "", bool update = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstThermalConductivity[markerUsed]=val;
            if ( update )
                this->updateThermalConductivity( cst(val), marker );
        }
        void setCstHeatCapacity( double val, std::string const& marker = "", bool update = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstHeatCapacity[markerUsed]=val;
            if ( update )
                this->updateHeatCapacity( cst(val), marker );
        }
        void setCstRho( double val, std::string const& marker = "", bool update = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstRho[markerUsed]=val;
            if ( update )
                this->updateRho( vf::cst(val), marker );
        }
        void setCstThermalExpansion( double val, std::string const& marker = "", bool update = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstThermalExpansion[markerUsed]=val;
            if ( update )
                this->updateThermalExpansion( vf::cst(val), marker );
        }

        template < typename ExprT >
        void setThermalConductivity( vf::Expr<ExprT> const& vfexpr, std::string const& marker = "" )
        {
            this->updateThermalConductivity( vfexpr,marker);
            if ( M_fieldThermalConductivity )
                this->setCstThermalConductivity( M_fieldThermalConductivity->min(), marker, false );
        }
        template < typename ExprT >
        void updateThermalConductivity(vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldThermalConductivity ) return;
            auto rangeEltUsed = ( marker.empty() )? elements( M_space->mesh()) : markedelements( M_space->mesh(),marker );
            M_fieldThermalConductivity->on(_range=rangeEltUsed,_expr=__expr);
        }

        template < typename ExprT >
        void setHeatCapacity( vf::Expr<ExprT> const& vfexpr, std::string const& marker = "" )
        {
            this->updateHeatCapacity( vfexpr,marker);
            if ( M_fieldHeatCapacity )
                this->setCstHeatCapacity( M_fieldHeatCapacity->min(), marker, false );
        }
        template < typename ExprT >
        void updateHeatCapacity(vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldHeatCapacity ) return;
            auto rangeEltUsed = ( marker.empty() )? elements( M_space->mesh()) : markedelements( M_space->mesh(),marker );
            M_fieldHeatCapacity->on(_range=rangeEltUsed,_expr=__expr);
        }

        template < typename ExprT >
        void setRho( vf::Expr<ExprT> const& vfexpr, std::string const& marker = "" )
        {
            this->updateRho( vfexpr,marker);
            if ( M_fieldRho )
                this->setCstRho( M_fieldRho->min(), marker, false );
        }
        template < typename ExprT >
        void updateRho(vf::Expr<ExprT> const& __expr, std::string const& marker = "")
        {
            if ( !M_fieldRho ) return;
            auto rangeEltUsed = ( marker.empty() )? elements( M_space->mesh()) : markedelements( M_space->mesh(),marker );
            M_fieldRho->on(_range=rangeEltUsed,_expr=__expr);
        }

        template < typename ExprT >
        void setThermalExpansion( vf::Expr<ExprT> const& vfexpr, std::string const& marker = "" )
        {
            this->updateThermalExpansion( vfexpr,marker);
            if ( M_fieldThermalExpansion )
                this->setCstThermalExpansion( M_fieldRho->min(), marker, false );
        }
        template < typename ExprT >
        void updateThermalExpansion(vf::Expr<ExprT> const& __expr, std::string const& marker = "")
        {
            if ( !M_fieldThermalExpansion ) return;
            auto rangeEltUsed = ( marker.empty() )? elements( M_space->mesh()) : markedelements( M_space->mesh(),marker );
            M_fieldThermalExpansion->on(_range=rangeEltUsed,_expr=__expr);
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

                if ( mat.hasPropertyExprScalar("rho") )
                    this->setRho( mat.propertyExprScalar("rho"),matmarker );
                else
                    this->setCstRho( mat.propertyConstant("rho"), matmarker );

                if ( mat.hasPropertyExprScalar("k11") )
                    this->setRho( mat.propertyExprScalar("k11"),matmarker );
                else
                    this->setCstThermalConductivity( mat.propertyConstant("k11"), matmarker );

                if ( mat.hasPropertyExprScalar("Cp") )
                    this->setRho( mat.propertyExprScalar("Cp"),matmarker );
                else
                    this->setCstHeatCapacity( mat.propertyConstant("Cp"), matmarker );

                if ( mat.hasPropertyExprScalar("beta") )
                    this->setRho( mat.propertyExprScalar("beta"),matmarker );
                else
                    this->setCstThermalExpansion( mat.propertyConstant("beta"), matmarker );
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
        bool M_isDefinedOnWholeMesh;
        space_ptrtype M_space;
        element_ptrtype M_fieldThermalConductivity, M_fieldHeatCapacity, M_fieldRho;
        element_ptrtype M_fieldThermalExpansion;
    };


} // namespace FeelModels
} // namespace Feel

#endif // __THERMALPROPERTIES_DESCRIPTION_H
