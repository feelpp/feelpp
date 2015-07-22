/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

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
    public :
        typedef SpaceType space_type;
        typedef boost::shared_ptr<SpaceType> space_ptrtype;
        typedef typename SpaceType::element_type element_type;
        typedef boost::shared_ptr<element_type> element_ptrtype;

        MechanicalPropertiesDescription( space_ptrtype const& space, bool useDisplacementPressureFormulation, std::string const& prefix )
            :
            M_materialLaw( soption(_name="material_law",_prefix=prefix) ),
            M_decouplingEnergyVolumicLaw( soption(_name="mechanicalproperties.compressible.volumic-law",_prefix=prefix) ),
            M_compressibleNeoHookeanVariantName( soption(_name="mechanicalproperties.compressible.neohookean.variant",_prefix=prefix) ),
            M_useDisplacementPressureFormulation( useDisplacementPressureFormulation ),
            M_youngModulus( doption(_name="youngmodulus",_prefix=prefix) ),// E
            M_coeffPoisson( doption(_name="coeffpoisson",_prefix=prefix) ),// sigma
            M_rho( doption(_name="rho",_prefix=prefix) ),// rho
            M_space( space ),
            M_fieldCoeffLame1( space->elementPtr( vf::cst( this->cstCoeffLame1() ) ) ),
            M_fieldCoeffLame2( space->elementPtr( vf::cst( this->cstCoeffLame2() ) ) ),
            M_fieldBulkModulus( space->elementPtr( vf::cst( this->cstBulkModulus() ) ) ),
            M_fieldRho( space->elementPtr( vf::cst(M_rho) ) )
            {
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
            M_useDisplacementPressureFormulation( false ),
            M_youngModulus( doption(_name="youngmodulus",_prefix=prefix) ),// E
            M_coeffPoisson( doption(_name="coeffpoisson",_prefix=prefix) ),// sigma
            M_rho( doption(_name="rho",_prefix=prefix) )// rho
            {
                CHECK( M_materialLaw == "StVenantKirchhoff" || M_materialLaw == "NeoHookean" ) << "invalid material law : " << M_materialLaw;
                CHECK( M_decouplingEnergyVolumicLaw == "classic" || M_decouplingEnergyVolumicLaw == "simo1985" ) << "invalid decouplingEnergyVolumicLaw : " << M_decouplingEnergyVolumicLaw;
                CHECK( M_compressibleNeoHookeanVariantName == "default" || M_compressibleNeoHookeanVariantName == "molecular-theory" ||
                       M_compressibleNeoHookeanVariantName == "molecular-theory-simo1985" ) << "invalid compressibleNeoHookeanVariantName : " <<  M_compressibleNeoHookeanVariantName;
            }

        MechanicalPropertiesDescription( MechanicalPropertiesDescription const& app  ) = default;

        void initFromSpace( space_ptrtype const& space )
        {
            M_space = space;
            M_fieldCoeffLame1 = space->elementPtr( vf::cst( this->cstCoeffLame1() ) );
            M_fieldCoeffLame2 = space->elementPtr( vf::cst( this->cstCoeffLame2() ) );
            M_fieldBulkModulus = space->elementPtr( vf::cst( this->cstBulkModulus() ) );
            M_fieldRho = space->elementPtr( vf::cst(M_rho) );
        }
        std::string const& materialLaw() const { return M_materialLaw; }
        std::string const& decouplingEnergyVolumicLaw() const { return M_decouplingEnergyVolumicLaw; }
        std::string const& compressibleNeoHookeanVariantName() const { return M_compressibleNeoHookeanVariantName; }

        bool useDisplacementPressureFormulation() const { return M_useDisplacementPressureFormulation; }
        void setUseDisplacementPressureFormulation( bool b ) { M_useDisplacementPressureFormulation = b; }

        double cstYoungModulus() const { return M_youngModulus; } // E
        double cstCoeffPoisson() const { return M_coeffPoisson; } // sigma;
        double cstCoeffLame2() const { return M_youngModulus/(2*(1+M_coeffPoisson)); }// mu
        double cstCoeffLame1() const { return M_youngModulus*M_coeffPoisson/((1+M_coeffPoisson)*(1-2*M_coeffPoisson)); }// lambda
        double cstBulkModulus() const { return this->cstCoeffLame1() + (2./3.)*this->cstCoeffLame2(); }
        double cstRho() const { return M_rho; }

        element_type const& fieldCoeffLame1() const { return *M_fieldCoeffLame1; }
        element_type const& fieldCoeffLame2() const { return *M_fieldCoeffLame2; }
        element_type const& fieldBulkModulus() const { return *M_fieldBulkModulus; }
        element_type const& fieldRho() const { return *M_fieldRho; }
        element_ptrtype const& fieldCoeffLame1Ptr() const { return M_fieldCoeffLame1; }
        element_ptrtype const& fieldCoeffLame2Ptr() const { return M_fieldCoeffLame2; }
        element_ptrtype const& fieldBulkModulusPtr() const { return M_fieldBulkModulus; }
        element_ptrtype const& fieldRhoPtr() const { return M_fieldRho; }

        void setCstYoungModulus( double val ) { M_youngModulus=val;this->updateLameCoeffFromYoungPoisson(); }
        void setCstCoeffPoisson( double val ) { M_coeffPoisson=val;this->updateLameCoeffFromYoungPoisson(); }
        void setCstRho( double val ) { M_rho=val; this->updateRho( vf::cst(M_rho) ); }

        template < typename ExprT >
        void updateCoeffLame1(vf::Expr<ExprT> const& __expr, bool updateOthersFields = true )
        {
            M_fieldCoeffLame1->on(_range=elements( M_space->mesh()),_expr=__expr);
            if ( updateOthersFields )
                this->updateBulkModulus();
        }
        template < typename ExprT >
        void updateCoeffLame2(vf::Expr<ExprT> const& __expr, bool updateOthersFields = true)
        {
            M_fieldCoeffLame2->on(_range=elements( M_space->mesh()),_expr=__expr);
            if ( updateOthersFields )
                this->updateBulkModulus();
        }
        void updateBulkModulus()
        {
            M_fieldBulkModulus->on(_range=elements( M_space->mesh()),_expr=idv(M_fieldCoeffLame1) + (2/3.)*idv(M_fieldCoeffLame2) );
        }
        template < typename ExprT >
        void updateRho(vf::Expr<ExprT> const& __expr)
        {
            M_fieldRho->on(_range=elements( M_space->mesh()),_expr=__expr);
        }

        void updateLameCoeffFromYoungPoisson()
        {
            if ( M_fieldCoeffLame1 && M_fieldCoeffLame2 )
            {
                this->updateCoeffLame1( vf::cst( this->cstCoeffLame1() ), false );
                this->updateCoeffLame2( vf::cst( this->cstCoeffLame2() ), true );
            }
        }

    private :
        std::string M_materialLaw, M_decouplingEnergyVolumicLaw, M_compressibleNeoHookeanVariantName;
        bool M_useDisplacementPressureFormulation;
        double M_youngModulus;// E;
        double M_coeffPoisson;// sigma;
        double M_rho;
        //boost::reference_wrapper<const element_muP0_type> M_coefflame1P0, M_coefflame2P0;
        space_ptrtype M_space;
        element_ptrtype M_fieldCoeffLame1, M_fieldCoeffLame2, M_fieldBulkModulus;
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
