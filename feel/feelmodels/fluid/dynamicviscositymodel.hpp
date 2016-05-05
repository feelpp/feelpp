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
/**
 \file dynamicviscositymodel.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2014-03-21
 */

#ifndef FEELPP_FLUIDMECHANICS_DYNAMICVISCOSITYMODEL_H
#define FEELPP_FLUIDMECHANICS_DYNAMICVISCOSITYMODEL_H 1

//#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{

template <class SpaceType>
class DynamicViscosityModel
{
    typedef DynamicViscosityModel<SpaceType> self_type;

  public:
    typedef SpaceType space_type;
    typedef boost::shared_ptr<SpaceType> space_ptrtype;
    typedef typename SpaceType::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    static std::string defaultMaterialName() { return std::string( "FEELPP_DEFAULT_MATERIAL_NAME" ); }

    DynamicViscosityModel( std::string prefix )
        : M_viscosityModelName( soption( _name = "viscosity.law", _prefix = prefix ) ),
          M_powerLaw_n( doption( _name = "power_law.n", _prefix = prefix ) ),
          M_powerLaw_k( doption( _name = "power_law.k", _prefix = prefix ) ),
          M_powerLaw_n_generic( ( M_powerLaw_n - 1. ) / 2. ), M_powerLaw_k_generic( M_powerLaw_k ),
          M_mu_0( doption( _name = "viscosity.zero_shear", _prefix = prefix ) ),
          M_mu_inf( doption( _name = "viscosity.infinite_shear", _prefix = prefix ) ),
          M_carreau_lambda( doption( _name = "carreau_law.lambda", _prefix = prefix ) ),
          M_carreau_n( doption( _name = "carreau_law.n", _prefix = prefix ) ),
          M_carreauYasuda_lambda( doption( _name = "carreau-yasuda_law.lambda", _prefix = prefix ) ),
          M_carreauYasuda_n( doption( _name = "carreau-yasuda_law.n", _prefix = prefix ) ),
          M_carreauYasuda_a( doption( _name = "carreau-yasuda_law.a", _prefix = prefix ) ),
          M_walburnSchneck_C1( doption( _name = "walburn-schneck_law.C1", _prefix = prefix ) ),
          M_walburnSchneck_C2( doption( _name = "walburn-schneck_law.C2", _prefix = prefix ) ),
          M_walburnSchneck_C3( doption( _name = "walburn-schneck_law.C3", _prefix = prefix ) ),
          M_walburnSchneck_C4( doption( _name = "walburn-schneck_law.C4", _prefix = prefix ) ),
          M_non_newtonian_hematocrit( doption( _name = "hematocrit", _prefix = prefix ) ),
          M_non_newtonian_TPMA( doption( _name = "TPMA", _prefix = prefix ) )
    {
        M_cstDynamicViscosity[self_type::defaultMaterialName()] = doption( _name = "mu", _prefix = prefix );

        CHECK( this->checkDynamicViscosityLaw() ) << "invalid viscosity model : " << this->dynamicViscosityLaw();
        if ( this->dynamicViscosityLaw() == "walburn-schneck_law" )
        {
            double hematocrit = this->nonNewtonianHematocrit();
            M_powerLaw_k_generic = this->walburnSchneck_C1() * math::exp( hematocrit * this->walburnSchneck_C2() ) *
                                   math::exp( this->walburnSchneck_C4() * this->nonNewtonianTPMA() / math::pow( hematocrit, 2 ) );
            M_powerLaw_n_generic = -this->walburnSchneck_C3() * hematocrit;
        }
    }
    DynamicViscosityModel( DynamicViscosityModel const& app ) = default;

    void initFromMesh( mesh_ptrtype const& mesh, bool useExtendedDofTable )
    {
        M_space = space_type::New( _mesh = mesh, _worldscomm = std::vector<WorldComm>( 1, mesh->worldComm() ),
                                   _extended_doftable = std::vector<bool>( 1, useExtendedDofTable ) );
        M_fieldDynamicViscosity = M_space->elementPtr( cst( this->cstDynamicViscosity() ) );
    }
    void initFromSpace( space_ptrtype const& space )
    {
        M_space = space;
        M_fieldDynamicViscosity = M_space->elementPtr( cst( this->cstDynamicViscosity() ) );
    }

    bool checkDynamicViscosityLaw() const
    {
        return ( this->dynamicViscosityLaw() == "newtonian" ||
                 this->dynamicViscosityLaw() == "power_law" || this->dynamicViscosityLaw() == "walburn-schneck_law" ||
                 this->dynamicViscosityLaw() == "carreau_law" || this->dynamicViscosityLaw() == "carreau-yasuda_law" );
    }
    std::string const& dynamicViscosityLaw() const { return M_viscosityModelName; }
    void setDynamicViscosityLaw( std::string s )
    {
        M_viscosityModelName = s;
        CHECK( this->checkDynamicViscosityLaw() ) << "invalid viscosity model : " << s;
    }

    space_ptrtype const& dynamicViscositySpace() const { return M_space; }

    std::set<std::string> const& markers() const { return M_markers; }

    double cstMu( std::string const& marker = "" ) const { return this->cstDynamicViscosity( marker ); }
    double cstDynamicViscosity( std::string const& marker = "" ) const
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        auto itFindMarker = M_cstDynamicViscosity.find( markerUsed );
        CHECK( itFindMarker != M_cstDynamicViscosity.end() ) << "invalid marker not registered " << markerUsed;
        return itFindMarker->second;
    }

    element_type const& fieldMu() const { return this->fieldDynamicViscosity(); }
    element_type const& fieldDynamicViscosity() const { return *M_fieldDynamicViscosity; }
    element_ptrtype const& fieldDynamicViscosityPtr() const { return M_fieldDynamicViscosity; }

    void setCstDynamicViscosity( double d, std::string const& marker = "" )
    {
        std::string markerUsed = ( marker.empty() ) ? self_type::defaultMaterialName() : marker;
        M_cstDynamicViscosity[markerUsed] = d;
        this->updateDynamicViscosity( cst( d ), marker );
    }

    template <typename ExprT>
    void updateDynamicViscosity( vf::Expr<ExprT> const& __expr, std::string const& marker = "" )
    {
        if ( !M_fieldDynamicViscosity ) return;
        if ( marker.empty() )
            M_fieldDynamicViscosity->on( _range = elements( M_fieldDynamicViscosity->mesh() ), _expr = __expr );
        else
            M_fieldDynamicViscosity->on( _range = markedelements( M_fieldDynamicViscosity->mesh(), marker ), _expr = __expr );
    }

    double powerLaw_n() const { return M_powerLaw_n; }
    void powerLaw_n( double d ) { M_powerLaw_n = d; }
    double powerLaw_k() const { return M_powerLaw_k; }
    void powerLaw_k( double d ) { M_powerLaw_k = d; }
    double powerLaw_n_generic() const { return M_powerLaw_n_generic; }
    double powerLaw_k_generic() const { return M_powerLaw_k_generic; }

    double mu_0() const { return M_mu_0; }
    void mu_0( double d ) { M_mu_0 = d; }
    double mu_inf() const { return M_mu_inf; }
    void mu_inf( double d ) { M_mu_inf = d; }
    double carreau_lambda() const { return M_carreau_lambda; }
    void carreau_lambda( double d ) { M_carreau_lambda = d; }
    double carreau_n() const { return M_carreau_n; }
    void carreau_n( double d ) { M_carreau_n = d; }
    double carreauYasuda_lambda() const { return M_carreauYasuda_lambda; }
    void carreauYasuda_lambda( double d ) { M_carreauYasuda_lambda = d; }
    double carreauYasuda_n() const { return M_carreauYasuda_n; }
    void carreauYasuda_n( double d ) { M_carreauYasuda_n = d; }
    double carreauYasuda_a() const { return M_carreauYasuda_a; }
    void carreauYasuda_a( double d ) { M_carreauYasuda_a = d; }

    double nonNewtonianHematocrit() const { return M_non_newtonian_hematocrit; }
    void nonNewtonianHematocrit( double d ) { M_non_newtonian_hematocrit = d; }
    double nonNewtonianTPMA() const { return M_non_newtonian_TPMA; }
    void nonNewtonianTPMA( double d ) { M_non_newtonian_TPMA = d; }

    double walburnSchneck_C1() const { return M_walburnSchneck_C1; }
    void walburnSchneck_C1( double d ) { M_walburnSchneck_C1 = d; }
    double walburnSchneck_C2() const { return M_walburnSchneck_C2; }
    void walburnSchneck_C2( double d ) { M_walburnSchneck_C2 = d; }
    double walburnSchneck_C3() const { return M_walburnSchneck_C3; }
    void walburnSchneck_C3( double d ) { M_walburnSchneck_C3 = d; }
    double walburnSchneck_C4() const { return M_walburnSchneck_C4; }
    void walburnSchneck_C4( double d ) { M_walburnSchneck_C4 = d; }

    void updateFromModelMaterials( ModelMaterials const& mat )
    {
        if ( mat.empty() ) return;
        for ( auto const& m : mat )
        {
            auto const& mat = m.second;
            auto const& matmarker = m.first;
            M_markers.insert( matmarker );
            this->setCstDynamicViscosity( mat.mu(), matmarker );
        }
    }

    boost::shared_ptr<std::ostringstream>
    getInfo() const
    {
        boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
        if ( this->dynamicViscosityLaw() == "newtonian" )
        {
        }
        else if ( this->dynamicViscosityLaw() == "power_law" )
        {
            *_ostr << "\n   Power Law "
                   << "\n     -- n   : " << this->powerLaw_n()
                   << "\n     -- k   : " << this->powerLaw_k();
        }
        else if ( this->dynamicViscosityLaw() == "walburn-schneck_law" )
        {
            *_ostr << "\n   Walburn-Schneck "
                   << "\n     -- C1  : " << this->walburnSchneck_C1()
                   << "\n     -- C2  : " << this->walburnSchneck_C2()
                   << "\n     -- C3  : " << this->walburnSchneck_C3()
                   << "\n     -- C4  : " << this->walburnSchneck_C4()
                   << "\n     -- hematocrit  : " << this->nonNewtonianHematocrit()
                   << "\n     -- TPMA  : " << this->nonNewtonianTPMA();
        }
        else if ( this->dynamicViscosityLaw() == "carreau_law" )
        {
            *_ostr << "\n   Carreau "
                   << "\n     -- lambda  : " << this->carreau_lambda()
                   << "\n     -- n       : " << this->carreau_n()
                   << "\n     -- mu_0    : " << this->mu_0()
                   << "\n     -- mu_inf  : " << this->mu_inf();
        }
        else if ( this->dynamicViscosityLaw() == "carreau-yasuda_law" )
        {
            *_ostr << "\n   Carreau-Yasuda "
                   << "\n     -- lambda  : " << this->carreauYasuda_lambda()
                   << "\n     -- n       : " << this->carreauYasuda_n()
                   << "\n     -- a       : " << this->carreauYasuda_a()
                   << "\n     -- mu_0    : " << this->mu_0()
                   << "\n     -- mu_inf  : " << this->mu_inf();
        }

        return _ostr;
    }

  private:
    std::string M_viscosityModelName;
    space_ptrtype M_space;
    std::set<std::string> M_markers;

    element_ptrtype M_fieldDynamicViscosity;
    std::map<std::string, double> M_cstDynamicViscosity;

    double M_powerLaw_n, M_powerLaw_k, M_powerLaw_n_generic, M_powerLaw_k_generic;
    double M_mu_0, M_mu_inf;
    double M_carreau_lambda, M_carreau_n;
    double M_carreauYasuda_lambda, M_carreauYasuda_n, M_carreauYasuda_a;

    double M_walburnSchneck_C1, M_walburnSchneck_C2, M_walburnSchneck_C3, M_walburnSchneck_C4;
    double M_non_newtonian_hematocrit;
    double M_non_newtonian_TPMA; //Total Proteins Minus Albumin (TPMA)
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_FLUIDMECHANICS_DYNAMICVISCOSITYMODEL_H
