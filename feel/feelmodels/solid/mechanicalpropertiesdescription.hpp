/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef __MECHANICALPROPERTIES_DESCRIPTION_H
#define __MECHANICALPROPERTIES_DESCRIPTION_H 1

#include <feel/feelvf/cst.hpp>

namespace Feel
{
namespace FeelModels
{

    template<class SpaceType>
    class MechanicalPropertiesDescription
    {
        typedef MechanicalPropertiesDescription<SpaceType> self_type;
    public :
        typedef SpaceType space_type;
        typedef boost::shared_ptr<SpaceType> space_ptrtype;
        typedef typename SpaceType::element_type element_type;
        typedef boost::shared_ptr<element_type> element_ptrtype;
        static std::string defaultMaterialName() { return std::string("FEELPP_DEFAULT_MATERIAL_NAME"); }
        MechanicalPropertiesDescription( space_ptrtype const& space, bool useDisplacementPressureFormulation, std::string const& prefix )
            :
            M_materialLaw( soption(_name="material_law",_prefix=prefix) ),
            M_decouplingEnergyVolumicLaw( soption(_name="mechanicalproperties.compressible.volumic-law",_prefix=prefix) ),
            M_compressibleNeoHookeanVariantName( soption(_name="mechanicalproperties.compressible.neohookean.variant",_prefix=prefix) ),
            M_useDisplacementPressureFormulation( useDisplacementPressureFormulation )
            {
                M_cstYoungModulus[self_type::defaultMaterialName()] = doption(_name="youngmodulus",_prefix=prefix);// E
                M_cstCoeffPoisson[self_type::defaultMaterialName()] = doption(_name="coeffpoisson",_prefix=prefix);// nu
                M_cstRho[self_type::defaultMaterialName()] = doption(_name="rho",_prefix=prefix);// rho
                this->initFromSpace( space );
                CHECK( M_materialLaw == "StVenantKirchhoff" || M_materialLaw == "NeoHookean" ) << "invalid material law : " << M_materialLaw;
                CHECK( M_decouplingEnergyVolumicLaw == "classic" || M_decouplingEnergyVolumicLaw == "simo1985" ) << "invalid decouplingEnergyVolumicLaw : " << M_decouplingEnergyVolumicLaw;
                CHECK( M_compressibleNeoHookeanVariantName == "default" || M_compressibleNeoHookeanVariantName == "molecular-theory" ||
                       M_compressibleNeoHookeanVariantName == "molecular-theory-simo1985" ) << "invalid compressibleNeoHookeanVariantName : " <<  M_compressibleNeoHookeanVariantName;
            }

        MechanicalPropertiesDescription( std::string const& prefix )
            :
            M_materialLaw( soption(_name="material_law",_prefix=prefix) ),
            M_decouplingEnergyVolumicLaw( soption(_name="mechanicalproperties.compressible.volumic-law",_prefix=prefix) ),
            M_compressibleNeoHookeanVariantName( soption(_name="mechanicalproperties.compressible.neohookean.variant",_prefix=prefix) ),
            M_useDisplacementPressureFormulation( false )
            {
                M_cstYoungModulus[self_type::defaultMaterialName()] = doption(_name="youngmodulus",_prefix=prefix);// E
                M_cstCoeffPoisson[self_type::defaultMaterialName()] = doption(_name="coeffpoisson",_prefix=prefix);// nu
                M_cstRho[self_type::defaultMaterialName()] = doption(_name="rho",_prefix=prefix);// rho
                CHECK( M_materialLaw == "StVenantKirchhoff" || M_materialLaw == "NeoHookean" ) << "invalid material law : " << M_materialLaw;
                CHECK( M_decouplingEnergyVolumicLaw == "classic" || M_decouplingEnergyVolumicLaw == "simo1985" ) << "invalid decouplingEnergyVolumicLaw : " << M_decouplingEnergyVolumicLaw;
                CHECK( M_compressibleNeoHookeanVariantName == "default" || M_compressibleNeoHookeanVariantName == "molecular-theory" ||
                       M_compressibleNeoHookeanVariantName == "molecular-theory-simo1985" ) << "invalid compressibleNeoHookeanVariantName : " <<  M_compressibleNeoHookeanVariantName;
            }

        MechanicalPropertiesDescription( MechanicalPropertiesDescription const& app  ) = default;

        void initFromSpace( space_ptrtype const& space )
        {
            M_space = space;
            M_fieldYoungModulus = space->elementPtr( vf::cst( this->cstYoungModulus() ) );
            M_fieldCoeffPoisson = space->elementPtr( vf::cst( this->cstCoeffPoisson() ) );
            M_fieldCoeffLame1 = space->elementPtr( vf::cst( this->cstCoeffLame1() ) );
            M_fieldCoeffLame2 = space->elementPtr( vf::cst( this->cstCoeffLame2() ) );
            M_fieldBulkModulus = space->elementPtr( vf::cst( this->cstBulkModulus() ) );
            M_fieldRho = space->elementPtr( vf::cst( this->cstRho() ) );
        }
        std::string const& materialLaw() const { return M_materialLaw; }
        std::string const& decouplingEnergyVolumicLaw() const { return M_decouplingEnergyVolumicLaw; }
        std::string const& compressibleNeoHookeanVariantName() const { return M_compressibleNeoHookeanVariantName; }

        bool useDisplacementPressureFormulation() const { return M_useDisplacementPressureFormulation; }
        void setUseDisplacementPressureFormulation( bool b ) { M_useDisplacementPressureFormulation = b; }

        double cstYoungModulus( std::string const& marker = "" ) const  // E
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstYoungModulus.find( markerUsed );
            CHECK( itFindMarker != M_cstYoungModulus.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }
        double cstCoeffPoisson( std::string const& marker = "" ) const // nu
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstCoeffPoisson.find( markerUsed );
            CHECK( itFindMarker != M_cstCoeffPoisson.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }
        double cstRho( std::string const& marker = "" ) const
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            auto itFindMarker = M_cstRho.find( markerUsed );
            CHECK( itFindMarker != M_cstRho.end() ) << "invalid marker not registered " << markerUsed;
            return itFindMarker->second;
        }

        double cstCoeffLame2( std::string const& marker = "" ) const // mu
        {
            return this->cstYoungModulus(marker)/(2*(1+this->cstCoeffPoisson(marker)));
        }
        double cstCoeffLame1( std::string const& marker = "" ) const // lambda
        {
            double mycoeffPoisson =this->cstCoeffPoisson(marker);
            return this->cstYoungModulus(marker)*mycoeffPoisson/((1+mycoeffPoisson)*(1-2*mycoeffPoisson));
        }
        double cstBulkModulus( std::string const& marker = "" ) const
        {
            return this->cstCoeffLame1(marker) + (2./3.)*this->cstCoeffLame2(marker);
        }
        std::set<std::string> const& markers() const { return M_markers; }

        element_type const& fieldYoungModulus() const { return *M_fieldYoungModulus; }
        element_type const& fieldCoeffPoisson() const { return *M_fieldCoeffPoisson; }
        element_type const& fieldCoeffLame1() const { return *M_fieldCoeffLame1; }
        element_type const& fieldCoeffLame2() const { return *M_fieldCoeffLame2; }
        element_type const& fieldBulkModulus() const { return *M_fieldBulkModulus; }
        element_type const& fieldRho() const { return *M_fieldRho; }
        element_ptrtype const& fieldYoungModulusPtr() const { return M_fieldYoungModulus; }
        element_ptrtype const& fieldCoeffPoissonPtr() const { return M_fieldCoeffPoisson; }
        element_ptrtype const& fieldCoeffLame1Ptr() const { return M_fieldCoeffLame1; }
        element_ptrtype const& fieldCoeffLame2Ptr() const { return M_fieldCoeffLame2; }
        element_ptrtype const& fieldBulkModulusPtr() const { return M_fieldBulkModulus; }
        element_ptrtype const& fieldRhoPtr() const { return M_fieldRho; }

        void setCstYoungModulus( double val, std::string const& marker = "", bool updateOthersFields = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstYoungModulus[markerUsed]=val;
            this->updateYoungModulus( cst(val), marker, updateOthersFields );
        }
        void setCstCoeffPoisson( double val, std::string const& marker = "", bool updateOthersFields = true )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstCoeffPoisson[markerUsed]=val;
            this->updateCoeffPoisson( cst(val), marker, updateOthersFields );
        }
        void setCstRho( double val, std::string const& marker = "" )
        {
            std::string markerUsed = ( marker.empty() )? self_type::defaultMaterialName() : marker;
            M_cstRho[markerUsed]=val;
            this->updateRho( vf::cst(val), marker );
        }

        template < typename ExprT >
        void updateYoungModulus(vf::Expr<ExprT> const& __expr, std::string const& marker = "", bool updateOthersFields = true )
        {
            if ( !M_fieldYoungModulus ) return;
            if ( marker.empty() )
                M_fieldYoungModulus->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldYoungModulus->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
            if ( updateOthersFields )
                this->updateOtherMaterialFields( marker );
        }
        template < typename ExprT >
        void updateCoeffPoisson(vf::Expr<ExprT> const& __expr, std::string const& marker = "", bool updateOthersFields = true )
        {
            if ( !M_fieldCoeffPoisson ) return;
            if ( marker.empty() )
                M_fieldCoeffPoisson->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldCoeffPoisson->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
            if ( updateOthersFields )
                this->updateOtherMaterialFields( marker );
        }

        template < typename ExprT >
        void updateCoeffLame1(vf::Expr<ExprT> const& __expr, std::string const& marker = "", bool updateOthersFields = true )
        {
            if ( marker.empty() )
                M_fieldCoeffLame1->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldCoeffLame1->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
            if ( updateOthersFields )
                this->updateBulkModulus( marker );
        }
        template < typename ExprT >
        void updateCoeffLame2(vf::Expr<ExprT> const& __expr, std::string const& marker = "", bool updateOthersFields = true)
        {
            if ( marker.empty() )
                M_fieldCoeffLame2->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldCoeffLame2->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
            if ( updateOthersFields )
                this->updateBulkModulus( marker );
        }
        void updateBulkModulus( std::string const& marker = "" )
        {
            auto bulkModulusExpr = idv(M_fieldCoeffLame1) + (2/3.)*idv(M_fieldCoeffLame2);
            if ( marker.empty() )
                M_fieldBulkModulus->on(_range=elements( M_space->mesh()),_expr=bulkModulusExpr );
            else
                M_fieldBulkModulus->on(_range=markedelements( M_space->mesh(),marker ),_expr=bulkModulusExpr );
        }
        template < typename ExprT >
        void updateRho(vf::Expr<ExprT> const& __expr, std::string const& marker = "")
        {
            if ( marker.empty() )
                M_fieldRho->on(_range=elements( M_space->mesh()),_expr=__expr);
            else
                M_fieldRho->on(_range=markedelements( M_space->mesh(),marker ),_expr=__expr);
        }

        void updateOtherMaterialFields( std::string const& marker = "" )
        {
            if ( M_fieldCoeffLame1 && M_fieldCoeffLame2 )
            {
                //this->updateCoeffLame1( vf::cst( this->cstCoeffLame1() ), marker, false );
                //this->updateCoeffLame2( vf::cst( this->cstCoeffLame2() ), marker, true );
                auto coeffLame1Expr = idv(M_fieldYoungModulus)*idv(M_fieldCoeffPoisson)/((1+idv(M_fieldCoeffPoisson))*(1-2*idv(M_fieldCoeffPoisson)));
                auto coeffLame2Expr = idv(M_fieldYoungModulus)/(2*(1+idv(M_fieldCoeffPoisson)));
                this->updateCoeffLame1( coeffLame1Expr, marker, false );
                this->updateCoeffLame2( coeffLame2Expr, marker, true );
            }
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
                this->setCstYoungModulus( mat.E(), matmarker, false );
                this->setCstCoeffPoisson( mat.nu(), matmarker, true );
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
                     << "\n     -- " << matmarkertag << "young modulus : " << this->cstYoungModulus(matmarker)
                     << "\n     -- " << matmarkertag << "coeff poisson : " << this->cstCoeffPoisson(matmarker)
                     << "\n     -- " << matmarkertag << "coeff Lame 1 : " << this->cstCoeffLame1(matmarker)
                     << "\n     -- " << matmarkertag << "coeff Lame 2 : " << this->cstCoeffLame2(matmarker);
            }
            return ostr;
        }

    private :
        std::string M_materialLaw, M_decouplingEnergyVolumicLaw, M_compressibleNeoHookeanVariantName;
        bool M_useDisplacementPressureFormulation;
        std::map<std::string,double> M_cstYoungModulus;// E;
        std::map<std::string,double> M_cstCoeffPoisson;// nu;
        std::map<std::string,double> M_cstRho;
        std::set<std::string> M_markers;
        //boost::reference_wrapper<const element_muP0_type> M_coefflame1P0, M_coefflame2P0;
        space_ptrtype M_space;
        element_ptrtype M_fieldYoungModulus, M_fieldCoeffPoisson, M_fieldCoeffLame1, M_fieldCoeffLame2, M_fieldBulkModulus;
        element_ptrtype M_fieldRho;
    };


    template<class SpaceType>
    MechanicalPropertiesDescription<SpaceType>
    mechanicalPropertiesDesc( boost::shared_ptr<SpaceType> const& space, bool useDisplacementPressureFormulation, std::string const& prefix )
    {
        MechanicalPropertiesDescription<SpaceType> res(space,useDisplacementPressureFormulation,prefix);
        return res;
    }

} // namespace FeelModels
} // namespace Feel

#endif // __MECHANICALPROPERTIES_DESCRIPTION_H
