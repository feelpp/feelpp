/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2014-03-21

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_MATERIALPROPERTIES_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_MATERIALPROPERTIES_HPP 1

//#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>
//#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelmaterials.hpp>

namespace Feel
{
namespace FeelModels
{

struct DynamicViscosityLaw: public ModelMaterialLaw
{
    using super_type = ModelMaterialLaw;

    struct Law { 
        virtual ~Law() = default;
        virtual std::string name() const = 0;
    };
    struct NewtonianLaw: Law { virtual std::string name() const override { return "newtonian"; } };
    struct PowerLaw: Law { virtual std::string name() const override { return "power_law"; } };
    struct CarreauLaw: Law { virtual std::string name() const override { return "carreau_law"; } };
    struct CarreauYasudaLaw: Law { virtual std::string name() const override { return "carreau-yasuda_law"; } };
    struct MultifluidLaw: Law { 
        MultifluidLaw() = default;
        MultifluidLaw( const std::string & oF, const std::set<std::string> & iFs ):
            outerFluid{ oF }, innerFluids{ iFs }
        {}
        virtual std::string name() const override { return "multifluid"; }
        uint32_t nFluids() const { return 1 + innerFluids.size(); }
        std::string outerFluid;
        std::set<std::string> innerFluids;
    };

    DynamicViscosityLaw( std::string const& law ):
        super_type( ),
        M_law( createLaw( law ) )
    {}
    DynamicViscosityLaw( std::unique_ptr<Law> && law ):
        super_type( ),
        M_law( std::move(law) )
    {}
    DynamicViscosityLaw( DynamicViscosityLaw const& ) = default;
    DynamicViscosityLaw( DynamicViscosityLaw && ) = default;

    void setup( nl::json const& jarg ) override;

    void setLaw( std::string const& law ) { M_law.reset( createLaw( law ) ); }
    void setLaw( std::unique_ptr<Law> && law ) { M_law = std::move(law); }
    template< typename DVLawT >
    DVLawT const & law() const { return static_cast<DVLawT const &>( *M_law ); }

    std::string lawName() const { return M_law->name(); }

    bool isNewtonianLaw() const { return this->lawName() == "newtonian"; }
    bool isPowerLaw() const { return this->lawName() == "power_law"; }
    bool isCarreauLaw() const { return this->lawName() == "carreau_law"; }
    bool isCarreauYasudaLaw() const { return this->lawName() == "carreau-yasuda_law"; }
    bool isMultifluidLaw() const { return this->lawName() == "multifluid"; }

    bool isConstant() const { return this->isNewtonianLaw(); }

private:
    static Law* createLaw( std::string const& law );

private:
    std::unique_ptr<Law> M_law;
};

#if 0
/**
 * Power Law parameters
 */
struct DynamicViscosityPowerLaw
{
    DynamicViscosityPowerLaw( std::string const& prefix )
        :
        M_k( doption(_name="power_law.k",_prefix=prefix) ),
        M_n( doption(_name="power_law.n",_prefix=prefix) ),
        M_muMin( doption(_name="viscosity.min",_prefix=prefix) ),
        M_muMax( doption(_name="viscosity.max",_prefix=prefix) )
        {}
    DynamicViscosityPowerLaw( DynamicViscosityPowerLaw const& ) = default;
    DynamicViscosityPowerLaw( DynamicViscosityPowerLaw && ) = default;
    DynamicViscosityPowerLaw& operator=( DynamicViscosityPowerLaw const& ) = default;
    DynamicViscosityPowerLaw& operator=( DynamicViscosityPowerLaw && ) = default;

    double k() const { return M_k; }
    double n() const { return M_n; }
    double muMin() const { return M_muMin; }
    double muMax() const { return M_muMax; }

    void setK( double k ) { M_k = k; }
    void setN( double n ) { M_n = n; }
    void setMuMin( double v ) { M_muMin = v; }
    void setMuMax( double v ) { M_muMax = v; }

    void setup( pt::ptree const& pt )
        {
            if ( auto k  = pt.get_optional<double>("k") )
                M_k = *k;
            if ( auto n  = pt.get_optional<double>("n") )
                M_n = *n;
            if ( auto muMin  = pt.get_optional<double>("muMin") )
                M_muMin = *muMin;
            if ( auto muMax  = pt.get_optional<double>("muMax") )
                M_muMax = *muMax;
        }

private :
    double M_k, M_n, M_muMin, M_muMax;
};

/**
 * Carreau Law parameters
 */
struct DynamicViscosityCarreauLaw
{
    DynamicViscosityCarreauLaw( std::string const& prefix )
        :
        M_mu0( doption(_name="viscosity.zero_shear",_prefix=prefix) ),
        M_muInf( doption(_name="viscosity.infinite_shear",_prefix=prefix) ),
        M_lambda( doption(_name="carreau_law.lambda",_prefix=prefix) ),
        M_n( doption(_name="carreau_law.n",_prefix=prefix) )
        {}
    DynamicViscosityCarreauLaw( DynamicViscosityCarreauLaw const& ) = default;
    DynamicViscosityCarreauLaw( DynamicViscosityCarreauLaw && ) = default;
    DynamicViscosityCarreauLaw& operator=( DynamicViscosityCarreauLaw const& ) = default;
    DynamicViscosityCarreauLaw& operator=( DynamicViscosityCarreauLaw && ) = default;

    double mu0() const { return M_mu0; }
    double muInf() const { return M_muInf; }
    double lambda() const { return M_lambda; }
    double n() const { return M_n; }

    void setup( pt::ptree const& pt )
        {
            if ( auto mu0  = pt.get_optional<double>("mu0") )
                M_mu0 = *mu0;
            if ( auto muInf  = pt.get_optional<double>("muInf") )
                M_muInf = *muInf;
            if ( auto l  = pt.get_optional<double>("lambda") )
                M_lambda = *l;
            if ( auto n  = pt.get_optional<double>("n") )
                M_n = *n;
        }

private :
    double M_mu0, M_muInf;
    double M_lambda, M_n;
};

/**
 * Carreau-Yasuda Law parameters
 */
struct DynamicViscosityCarreauYasudaLaw
{
    DynamicViscosityCarreauYasudaLaw( std::string const& prefix )
        :
        M_mu0( doption(_name="viscosity.zero_shear",_prefix=prefix) ),
        M_muInf( doption(_name="viscosity.infinite_shear",_prefix=prefix) ),
        M_lambda( doption(_name="carreau-yasuda_law.lambda",_prefix=prefix) ),
        M_n( doption(_name="carreau-yasuda_law.n",_prefix=prefix) ),
        M_a( doption(_name="carreau-yasuda_law.a",_prefix=prefix) )
        {}
    DynamicViscosityCarreauYasudaLaw( DynamicViscosityCarreauYasudaLaw const& ) = default;
    DynamicViscosityCarreauYasudaLaw( DynamicViscosityCarreauYasudaLaw && ) = default;
    DynamicViscosityCarreauYasudaLaw& operator=( DynamicViscosityCarreauYasudaLaw const& ) = default;
    DynamicViscosityCarreauYasudaLaw& operator=( DynamicViscosityCarreauYasudaLaw && ) = default;

    double mu0() const { return M_mu0; }
    double muInf() const { return M_muInf; }
    double lambda() const { return M_lambda; }
    double n() const { return M_n; }
    double a() const { return M_a; }

    void setup( pt::ptree const& pt )
        {
            if ( auto mu0  = pt.get_optional<double>("mu0") )
                M_mu0 = *mu0;
            if ( auto muInf  = pt.get_optional<double>("muInf") )
                M_muInf = *muInf;
            if ( auto l  = pt.get_optional<double>("lambda") )
                M_lambda = *l;
            if ( auto n  = pt.get_optional<double>("n") )
                M_n = *n;
            if ( auto a  = pt.get_optional<double>("a") )
                M_a = *a;
        }

private :
    double M_mu0, M_muInf;
    double M_lambda, M_n, M_a;
};

/**
 * Walburn-Schneck Law parameters
 */
struct DynamicViscosityWalburnSchneckLaw
{
    DynamicViscosityWalburnSchneckLaw( std::string const& prefix )
        :
        M_C1( doption(_name="walburn-schneck_law.C1",_prefix=prefix) ),
        M_C2( doption(_name="walburn-schneck_law.C2",_prefix=prefix) ),
        M_C3( doption(_name="walburn-schneck_law.C3",_prefix=prefix) ),
        M_C4( doption(_name="walburn-schneck_law.C4",_prefix=prefix) ),
        M_hematocrit( doption(_name="hematocrit",_prefix=prefix) ),
        M_TPMA( doption(_name="TPMA",_prefix=prefix) ),
        M_powerLaw( prefix )
        {
            this->updateForUse();
        }
    DynamicViscosityWalburnSchneckLaw( DynamicViscosityWalburnSchneckLaw const& ) = default;
    DynamicViscosityWalburnSchneckLaw( DynamicViscosityWalburnSchneckLaw && ) = default;
    DynamicViscosityWalburnSchneckLaw& operator=( DynamicViscosityWalburnSchneckLaw const& ) = default;
    DynamicViscosityWalburnSchneckLaw& operator=( DynamicViscosityWalburnSchneckLaw && ) = default;

    double C1() const { return M_C1; }
    double C2() const { return M_C2; }
    double C3() const { return M_C3; }
    double C4() const { return M_C4; }
    double hematocrit() const { return M_hematocrit; }
    double TPMA() const { return M_TPMA; }

    DynamicViscosityPowerLaw const& powerLaw() const { return M_powerLaw; }

    void setup( pt::ptree const& pt )
        {
            if ( auto C1  = pt.get_optional<double>("C1") )
                M_C1 = *C1;
            if ( auto C2  = pt.get_optional<double>("C2") )
                M_C2 = *C2;
            if ( auto C3  = pt.get_optional<double>("C3") )
                M_C3 = *C3;
            if ( auto C4  = pt.get_optional<double>("C4") )
                M_C4 = *C4;
            if ( auto hematocrit  = pt.get_optional<double>("hematocrit") )
                M_hematocrit = *hematocrit;
            if ( auto TPMA  = pt.get_optional<double>("TPMA") )
                M_TPMA = *TPMA;
            if ( auto muMin = pt.get_optional<double>("muMin") )
                M_powerLaw.setMuMin( *muMin );
            if ( auto muMax = pt.get_optional<double>("muMax") )
                M_powerLaw.setMuMax( *muMax );
            this->updateForUse();
        }
private :
    void updateForUse()
    {
        M_powerLaw.setK( M_C1*math::exp(M_hematocrit*M_C2)*math::exp( M_C4*M_TPMA/math::pow(M_hematocrit,2) ) );
        M_powerLaw.setN( 1.-M_C3*M_hematocrit );
    }
    double M_C1, M_C2, M_C3, M_C4;
    double M_hematocrit;
    double M_TPMA; //Total Proteins Minus Albumin (TPMA)
    DynamicViscosityPowerLaw M_powerLaw;
};

/**
 * Manage dynamic viscosity for a material
 */
class DynamicViscosityMaterialProperties
{
public :
    DynamicViscosityMaterialProperties( std::string const& prefix, ModelMaterial const& mat )
        :
        M_lawName( soption(_name="viscosity.law",_prefix=prefix ) )
        {
            if ( mat.hasPropertyExprScalar("mu") )
            {
                auto const& expr = mat.propertyExprScalar("mu");
                M_newtonian.setExpr( expr );
            }

            pt::ptree const& pt = mat.pTree();
            if ( auto viscosityLawOpt = pt.template get_optional<std::string>("mu_law") )
                M_lawName = *viscosityLawOpt;
            CHECK( this->checkLawName() ) << "invalid law name : " << this->lawName();

            if ( this->lawName() == "power_law" )
            {
                M_powerLaw = boost::optional<DynamicViscosityPowerLaw>( prefix );
                if ( auto ptPowerLawOpt = pt.get_child_optional("mu_power-law") )
                    M_powerLaw->setup( *ptPowerLawOpt );
            }
            else if ( this->lawName() == "carreau_law" )
            {
                M_carreauLaw = boost::optional<DynamicViscosityCarreauLaw>( prefix );
                if ( auto ptCarreauLawOpt = pt.get_child_optional("mu_carreau-law") )
                    M_carreauLaw->setup( *ptCarreauLawOpt );
            }
            else if ( this->lawName() == "carreau-yasuda_law" )
            {
                M_carreauYasudaLaw = boost::optional<DynamicViscosityCarreauYasudaLaw>( prefix );
                if ( auto ptCarreauYasudaLawOpt = pt.get_child_optional("mu_carreau-yasuda-law") )
                    M_carreauYasudaLaw->setup( *ptCarreauYasudaLawOpt );
            }
            else if ( this->lawName() == "walburn-schneck_law" )
            {
                M_walburnSchneckLaw = boost::optional<DynamicViscosityWalburnSchneckLaw>( prefix );
                if ( auto ptWalburnSchneckLawLawOpt = pt.get_child_optional("mu_walburn-schneck-law") )
                    M_walburnSchneckLaw->setup( *ptWalburnSchneckLawLawOpt );
            }

        }

    DynamicViscosityMaterialProperties( DynamicViscosityMaterialProperties const& ) = default;
    DynamicViscosityMaterialProperties( DynamicViscosityMaterialProperties && ) = default;

    std::string const& lawName() const { return M_lawName; }

    bool isNewtonianLaw() const { return (this->lawName() == "newtonian"); }
    bool isPowerLaw() const { return (this->lawName() == "power_law"); }
    bool isCarreauLaw() const { return (this->lawName() == "carreau_law"); }
    bool isCarreauYasudaLaw() const { return (this->lawName() == "carreau-yasuda_law"); }
    bool isWalburnSchneckLaw() const { return (this->lawName() == "walburn-schneck_law"); }

    bool checkLawName() const
        {
            return ( this->isNewtonianLaw() || this->isPowerLaw() || this->isCarreauLaw() || this->isCarreauYasudaLaw() || this->isWalburnSchneckLaw() );
        }

    ModelExpressionScalar const& newtonian() const { return M_newtonian; }
    ModelExpressionScalar & newtonian() { return M_newtonian; }

    bool hasPowerLaw() const { return (M_powerLaw)? true : false; }
    bool hasCarreauLaw() const { return (M_carreauLaw)? true : false; }
    bool hasCarreauYasudaLaw() const { return (M_carreauYasudaLaw)? true : false; }
    bool hasWalburnSchneckLaw() const { return (M_walburnSchneckLaw)? true : false; }

    DynamicViscosityPowerLaw const& powerLaw() const { CHECK( this->hasPowerLaw() ) << "no PowerLaw";return *M_powerLaw; }
    DynamicViscosityCarreauLaw const& carreauLaw() const { CHECK( this->hasCarreauLaw() ) << "no CarreauLaw";return *M_carreauLaw; }
    DynamicViscosityCarreauYasudaLaw const& carreauYasudaLaw() const { CHECK( this->hasCarreauYasudaLaw() ) << "no CarreauYasudaLaw";return *M_carreauYasudaLaw; }
    DynamicViscosityWalburnSchneckLaw const& walburnSchneckLaw() const { CHECK( this->hasWalburnSchneckLaw() ) << "no WalburnSchneckLaw";return *M_walburnSchneckLaw; }

    void getInfo( std::shared_ptr<std::ostringstream> ostr ) const
        {
            *ostr << "\n          + law : " << this->lawName();
            if ( this->isNewtonianLaw() )
            {
                *ostr << "\n          + mu : ";
                if ( this->newtonian().isConstant() )
                    *ostr << this->newtonian().value();
                else
                    *ostr << this->newtonian().expr().expression();
            }
            else if ( this->isPowerLaw() )
            {
                *ostr << "\n          + k : " << this->powerLaw().k()
                      << "\n          + n : " << this->powerLaw().n()
                      << "\n          + muMin : " << this->powerLaw().muMin()
                      << "\n          + muMax : " << this->powerLaw().muMax();
            }
            else if ( this->isCarreauLaw() )
            {
                *ostr << "\n          + lambda : " << this->carreauLaw().lambda()
                      << "\n          + n : " << this->carreauLaw().n()
                      << "\n          + mu0 : " << this->carreauLaw().mu0()
                      << "\n          + muInf : " << this->carreauLaw().muInf();
            }
            else if ( this->isCarreauYasudaLaw() )
            {
                *ostr << "\n          + lambda : " << this->carreauYasudaLaw().lambda()
                      << "\n          + n : " << this->carreauYasudaLaw().n()
                      << "\n          + a : " << this->carreauYasudaLaw().a()
                      << "\n          + mu0 : " << this->carreauYasudaLaw().mu0()
                      << "\n          + muInf : " << this->carreauYasudaLaw().muInf();
            }
            else if ( this->isWalburnSchneckLaw() )
            {
                *ostr << "\n          + C1 : " << this->walburnSchneckLaw().C1()
                      << "\n          + C2 : " << this->walburnSchneckLaw().C2()
                      << "\n          + C3 : " << this->walburnSchneckLaw().C3()
                      << "\n          + C4 : " << this->walburnSchneckLaw().C4()
                      << "\n          + hematocrit : " << this->walburnSchneckLaw().hematocrit()
                      << "\n          + TPMA : " << this->walburnSchneckLaw().TPMA()
                      << "\n          + muMin : " << this->walburnSchneckLaw().powerLaw().muMin()
                      << "\n          + muMax : " << this->walburnSchneckLaw().powerLaw().muMax();
            }
        }
private :

    std::string M_lawName;

    ModelExpressionScalar M_newtonian;

    boost::optional<DynamicViscosityPowerLaw> M_powerLaw;
    boost::optional<DynamicViscosityCarreauLaw> M_carreauLaw;
    boost::optional<DynamicViscosityCarreauYasudaLaw> M_carreauYasudaLaw;
    boost::optional<DynamicViscosityWalburnSchneckLaw> M_walburnSchneckLaw;

    double M_mu0, M_muInf;
    double M_carreau_lambda, M_carreau_n;
    double M_carreauYasuda_lambda, M_carreauYasuda_n, M_carreauYasuda_a;

    double M_walburnSchneck_C1,M_walburnSchneck_C2,M_walburnSchneck_C3,M_walburnSchneck_C4;
    double M_non_newtonian_hematocrit;
    double M_non_newtonian_TPMA; //Total Proteins Minus Albumin (TPMA)
};

/**
 * Manage Fluid Mechanics material properties
 */
template<class SpaceType>
class FluidMechanicsMaterialProperties
{
    typedef FluidMechanicsMaterialProperties<SpaceType> self_type;
public :
    typedef SpaceType space_type;
    typedef std::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    FluidMechanicsMaterialProperties( std::string const& prefix )
        :
        M_prefix( prefix ),
        M_isDefinedOnWholeMesh( true ),
        M_dynamicViscosityDefaultValue( doption(_name="mu",_prefix=prefix) ),
        M_densityDefaultValue( doption(_name="rho",_prefix=prefix) )
        {}
    FluidMechanicsMaterialProperties( FluidMechanicsMaterialProperties const& ) = default;

    void updateForUse( mesh_ptrtype const& mesh , ModelMaterials const& mats, bool useExtendedDofTable )
        {
            this->updateForUseImpl( mesh,mats,useExtendedDofTable );
        }
    void updateForUse( space_ptrtype const& space, ModelMaterials const& mats )
        {
            this->updateForUseImpl( space->mesh(),mats,false,space );
        }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements( std::string const& matName ) const
        {
            CHECK( this->hasMaterial( matName ) ) << "no material " << matName;
            return M_rangeMeshElementsByMaterial.find( matName )->second;
        }
    bool hasMaterial( std::string const& matName ) const { return M_rangeMeshElementsByMaterial.find( matName ) != M_rangeMeshElementsByMaterial.end(); }

    space_ptrtype const& dynamicViscositySpace() const { return M_space; }

    std::set<std::string> const& markers() const { return M_markers; }

    bool isDefinedOnWholeMesh() const { return M_isDefinedOnWholeMesh; }

    // return true one of fluid material use a non newtonian law
    bool hasNonNewtonianLaw() const
        {
            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                std::string const& matName = rangeData.first;
                auto const& dynamicViscosity = this->dynamicViscosity( matName );
                if ( !dynamicViscosity.isNewtonianLaw() )
                    return true;
            }
            return false;
        }


    //! return dynamic viscosity projection field
    element_type const& fieldMu() const { return this->fieldDynamicViscosity(); }
    element_type const& fieldDynamicViscosity() const { return *M_fieldDynamicViscosity; }
    element_ptrtype const& fieldDynamicViscosityPtr() const { return M_fieldDynamicViscosity; }
    //! update dynamic viscosity projection field
    template < typename ExprT >
    void updateDynamicViscosityField( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldDynamicViscosity ) return;
            auto rangeEltUsed = ( marker.empty() )? elements(M_fieldDynamicViscosity->mesh()) : markedelements(M_fieldDynamicViscosity->mesh(),marker);
            M_fieldDynamicViscosity->on(_range=rangeEltUsed,_expr=__expr );
        }

    //! return true if has DynamicViscosity for this mat
    bool hasDynamicViscosity( std::string const& matName ) const
        {
            return M_dynamicViscosityModelByMaterial.find( matName ) != M_dynamicViscosityModelByMaterial.end();
        }
    //! return the expression manager of dynamicViscosity
    DynamicViscosityMaterialProperties/*ModelExpression*/ const& dynamicViscosity( std::string const& matName ) const
        {
            CHECK( this->hasDynamicViscosity( matName ) ) << "material name not registered : " << matName;
            return M_dynamicViscosityModelByMaterial.find( matName )->second;
        }
    //! return constant DynamicViscosity
    double cstMu( std::string const& matName = "" ) const { return this->cstDynamicViscosity( matName ); }
    //! return constant DynamicViscosity
    double cstDynamicViscosity( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_dynamicViscosityModelByMaterial.empty() )
                    return M_dynamicViscosityDefaultValue;
                else
                    return M_dynamicViscosityModelByMaterial.begin()->second.newtonian().value();
            }
            auto itFindMat = M_dynamicViscosityModelByMaterial.find( matName );
            CHECK( itFindMat != M_dynamicViscosityModelByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.newtonian().value();
        }

    void setCstDynamicViscosity( double d, std::string const& matName = "", bool updateField = true )
        {
            CHECK( false ) << "TODO";
        }


    //! return density projection field
    element_type const& fieldRho() const { return this->fieldDensity(); }
    element_type const& fieldDensity() const { return *M_fieldDensity; }
    element_ptrtype const& fieldDensityPtr() const { return M_fieldDensity; }
    //! update density projection field
    template < typename ExprT >
    void updateDensityField( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
        {
            if ( !M_fieldDensity ) return;
            auto rangeEltUsed = ( marker.empty() )? elements(M_fieldDensity->mesh()) : markedelements(M_fieldDensity->mesh(),marker);
            M_fieldDensity->on(_range=rangeEltUsed,_expr=__expr );
        }

    //! return true if has density description for this material name
    bool hasDensity( std::string const& matName ) const
        {
            return M_densityByMaterial.find( matName ) != M_densityByMaterial.end();
        }
    //! return the density description for a material
    ModelExpressionScalar const& density( std::string const& matName ) const
        {
            CHECK( this->hasDensity( matName ) ) << "material name not registered : " << matName;
            return M_densityByMaterial.find( matName )->second;
        }
    //! return constant density for a material
    double cstDensity( std::string const& matName = "" ) const
        {
            if ( matName.empty() )
            {
                if ( M_densityByMaterial.empty() )
                    return M_densityDefaultValue;
                else
                    return M_densityByMaterial.begin()->second.value();
            }
            auto itFindMat = M_densityByMaterial.find( matName );
            CHECK( itFindMat != M_densityByMaterial.end() ) << "material name not registered : " << matName;
            return itFindMat->second.value();
        }

    void setCstDensity( double d, std::string const& matName = "", bool updateField = true )
        {
            CHECK( false ) << "TODO";
        }



    std::shared_ptr<std::ostringstream>
    getInfo() const
    {
        std::shared_ptr<std::ostringstream> ostr( new std::ostringstream() );

        *ostr << "\n   Materials parameters";
        *ostr << "\n     -- defined on whole mesh : " << this->isDefinedOnWholeMesh();
        *ostr << "\n     -- number of materials : " << M_rangeMeshElementsByMaterial.size();
        for ( auto const& matRange : M_rangeMeshElementsByMaterial)
        {
            std::string const& matName = matRange.first;
            *ostr << "\n     -- [" << matName << "] rho : ";
            *ostr << "\n       - density : ";
            if ( this->density(matName).isConstant() )
                *ostr << this->density(matName).value();
            else
                *ostr << str( this->density(matName).expr().expression() );

            *ostr << "\n       - dynamic viscosity : ";
            this->dynamicViscosity( matName ).getInfo( ostr );
        }
        return ostr;
    }

    private :
    void updateForUseImpl( mesh_ptrtype const& mesh , ModelMaterials const& mats, bool useExtendedDofTable, space_ptrtype space=space_ptrtype() )
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
            if ( mat.hasPhysics() && !mat.hasPhysics( { "fluid","aerothermal","heat-fluid" } ) )
                continue;
            for ( std::string const& matmarker : mat.meshMarkers() )
            {
                if ( eltMarkersInMesh.find( matmarker ) == eltMarkersInMesh.end() )
                    continue;
                M_markers.insert( matmarker );
                markersByMaterial[matName].insert( matmarker );
            }
        }

        if ( !space )
        {
            if( M_markers.size() > 0 )
                M_isDefinedOnWholeMesh = ( M_markers.size() == eltMarkersInMesh.size() );
            else
                M_isDefinedOnWholeMesh = true;
            if ( M_isDefinedOnWholeMesh )
                M_space = space_type::New(_mesh=mesh, _extended_doftable=useExtendedDofTable );
            else
                M_space = space_type::New(_mesh=mesh,_range=markedelements(mesh,M_markers), _extended_doftable=useExtendedDofTable );
        }
        else
        {
            M_isDefinedOnWholeMesh = true;
            M_space = space;
        }
        M_fieldDynamicViscosity = M_space->elementPtr( cst( this->cstDynamicViscosity() ) );
        M_fieldDensity = M_space->elementPtr( cst( this->cstDensity() ) );

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

            M_dynamicViscosityModelByMaterial.insert( std::make_pair( matName, DynamicViscosityMaterialProperties(M_prefix,mat ) ) );
            auto const& dynamicViscosity = M_dynamicViscosityModelByMaterial.find( matName )->second;
            // update dynamic viscosity field (for newtonian)
            if ( dynamicViscosity.newtonian().isConstant() )
            {
                double value = dynamicViscosity.newtonian().value();
                M_fieldDynamicViscosity->on(_range=range,_expr=cst(value));
            }
            else
            {
                auto const& expr = dynamicViscosity.newtonian().expr();
                M_fieldDynamicViscosity->on(_range=range,_expr=expr);
            }

            M_densityByMaterial[matName];
            if ( mat.hasPropertyExprScalar("rho") )
            {
                auto const& expr = mat.propertyExprScalar("rho");
                M_densityByMaterial[matName].setExpr( expr );
                M_fieldDensity->on(_range=range,_expr=expr);
            }
            else M_densityByMaterial[matName].setExpr( expr("0") );
        }
    }


    template <typename SymbolsExpr>
    void updateFields( SymbolsExpr const& symbolsExpr )
        {
            // TODO
        }

private :
    std::string M_prefix;
    space_ptrtype M_space;
    std::set<std::string> M_markers;
    bool M_isDefinedOnWholeMesh;
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;

    element_ptrtype M_fieldDynamicViscosity;
    double M_dynamicViscosityDefaultValue;
    std::map<std::string, DynamicViscosityMaterialProperties> M_dynamicViscosityModelByMaterial;

    element_ptrtype M_fieldDensity;
    double M_densityDefaultValue;
    std::map<std::string, ModelExpressionScalar> M_densityByMaterial;
};
#endif

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_FLUIDMECHANICS_MATERIALPROPERTIES_HPP
