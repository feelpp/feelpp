/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
 \file viscositymodeldescription.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2014-03-21
 */

#ifndef FEELPP_VISCOSITY_MODEL_DESCRIPTION_H
#define FEELPP_VISCOSITY_MODEL_DESCRIPTION_H 1

namespace Feel
{
namespace FeelModels
{

template<class SpaceType>
class ViscosityModelDescription
{
public :
    typedef SpaceType space_type;
    typedef typename space_type::element_type element_type;

    ViscosityModelDescription( std::string const& _viscosityModel,
                               element_type const& muP0,
                               std::string prefix )
        :
        M_viscosityModelName( _viscosityModel ), M_muP0( muP0 ),
        M_powerLaw_n( doption(_name="power_law.n",_prefix=prefix) ),
        M_powerLaw_k( doption(_name="power_law.k",_prefix=prefix) ),
        M_powerLaw_n_generic( (M_powerLaw_n-1.)/2. ), M_powerLaw_k_generic( M_powerLaw_k ),
        M_mu_0( doption(_name="viscosity.zero_shear",_prefix=prefix) ),
        M_mu_inf( doption(_name="viscosity.infinite_shear",_prefix=prefix) ),
        M_carreau_lambda( doption(_name="carreau_law.lambda",_prefix=prefix) ),
        M_carreau_n( doption(_name="carreau_law.n",_prefix=prefix) ),
        M_carreauYasuda_lambda( doption(_name="carreau-yasuda_law.lambda",_prefix=prefix) ),
        M_carreauYasuda_n( doption(_name="carreau-yasuda_law.n",_prefix=prefix) ),
        M_carreauYasuda_a( doption(_name="carreau-yasuda_law.a",_prefix=prefix) ),
        M_walburnSchneck_C1( doption(_name="walburn-schneck_law.C1",_prefix=prefix) ),
        M_walburnSchneck_C2( doption(_name="walburn-schneck_law.C2",_prefix=prefix) ),
        M_walburnSchneck_C3( doption(_name="walburn-schneck_law.C3",_prefix=prefix) ),
        M_walburnSchneck_C4( doption(_name="walburn-schneck_law.C4",_prefix=prefix) ),
        M_non_newtonian_hematocrit( doption(_name="hematocrit",_prefix=prefix) ),
        M_non_newtonian_TPMA( doption(_name="TPMA",_prefix=prefix) )
        {
            if (this->name() == "walburn-schneck_law")
            {
                double hematocrit = this->nonNewtonianHematocrit();
                M_powerLaw_k_generic=this->walburnSchneck_C1()*math::exp(hematocrit*this->walburnSchneck_C2())*math::exp( this->walburnSchneck_C4()*this->nonNewtonianTPMA()/math::pow(hematocrit,2) );
                M_powerLaw_n_generic=-this->walburnSchneck_C3()*hematocrit;
            }

        }

    ViscosityModelDescription( std::string const& _viscosityModel,
                               element_type const& muP0,
                               double _powerLaw_n, double _powerLaw_k, double _mu_0, double _mu_inf,
                               double _carreau_lambda, double _carreau_n,
                               double _carreauYasuda_lambda, double _carreauYasuda_n, double _carreauYasuda_a,
                               double _walburnSchneck_C1, double _walburnSchneck_C2, double _walburnSchneck_C3, double _walburnSchneck_C4,
                               double _non_newtonian_hematocrit, double _non_newtonian_TPMA )
    :
    M_viscosityModelName( _viscosityModel ), M_muP0( muP0 ),
        M_powerLaw_n( _powerLaw_n ), M_powerLaw_k( _powerLaw_k ),
        M_powerLaw_n_generic( (M_powerLaw_n-1.)/2. ), M_powerLaw_k_generic( M_powerLaw_k ),
        M_mu_0( _mu_0 ), M_mu_inf( _mu_inf ),
        M_carreau_lambda( _carreau_lambda ), M_carreau_n( _carreau_n ),
        M_carreauYasuda_lambda( _carreauYasuda_lambda ), M_carreauYasuda_n( _carreauYasuda_n ), M_carreauYasuda_a( _carreauYasuda_a ),
        M_walburnSchneck_C1( _walburnSchneck_C1 ),M_walburnSchneck_C2( _walburnSchneck_C2 ),
        M_walburnSchneck_C3( _walburnSchneck_C3 ),M_walburnSchneck_C4( _walburnSchneck_C4 ),
        M_non_newtonian_hematocrit( _non_newtonian_hematocrit ), M_non_newtonian_TPMA( _non_newtonian_TPMA )
        {
            if (this->name() == "walburn-schneck_law")
            {
                double hematocrit = this->nonNewtonianHematocrit();
                M_powerLaw_k_generic=this->walburnSchneck_C1()*math::exp(hematocrit*this->walburnSchneck_C2())*math::exp( this->walburnSchneck_C4()*this->nonNewtonianTPMA()/math::pow(hematocrit,2) );
                M_powerLaw_n_generic=-this->walburnSchneck_C3()*hematocrit;
            }

        }

    ViscosityModelDescription( ViscosityModelDescription const& app  )
        :
        M_viscosityModelName( app.M_viscosityModelName ), M_muP0( app.M_muP0 ),
        M_powerLaw_n( app.M_powerLaw_n ), M_powerLaw_k( app.M_powerLaw_k ),
        M_powerLaw_n_generic( app.M_powerLaw_n_generic ), M_powerLaw_k_generic( app.M_powerLaw_k_generic ),
        M_mu_0( app.M_mu_0 ), M_mu_inf( app.M_mu_inf ),
        M_carreau_lambda( app.M_carreau_lambda ), M_carreau_n( app.M_carreau_n ),
        M_carreauYasuda_lambda( app.M_carreauYasuda_lambda ), M_carreauYasuda_n( app.M_carreauYasuda_n ), M_carreauYasuda_a( app.M_carreauYasuda_a ),
        M_walburnSchneck_C1( app.M_walburnSchneck_C1 ), M_walburnSchneck_C2( app.M_walburnSchneck_C2 ),
        M_walburnSchneck_C3( app.M_walburnSchneck_C3 ), M_walburnSchneck_C4( app.M_walburnSchneck_C4 ),
        M_non_newtonian_hematocrit( app.M_non_newtonian_hematocrit ), M_non_newtonian_TPMA( app.M_non_newtonian_TPMA )
        {}

    std::string const& name() const { return M_viscosityModelName; }
    void name( std::string s ) { M_viscosityModelName=s; }

    element_type const& muP0() const { return M_muP0; }

    double powerLaw_n() const { return M_powerLaw_n; }
    void powerLaw_n( double d ) { M_powerLaw_n=d; }
    double powerLaw_k() const { return M_powerLaw_k; }
    void powerLaw_k( double d ) { M_powerLaw_k=d; }
    double powerLaw_n_generic() const { return M_powerLaw_n_generic; }
    double powerLaw_k_generic() const { return M_powerLaw_k_generic; }

    double mu_0() const { return M_mu_0; }
    void mu_0( double d ) { M_mu_0=d; }
    double mu_inf() const { return M_mu_inf; }
    void mu_inf( double d ) { M_mu_inf=d; }
    double carreau_lambda() const { return M_carreau_lambda; }
    void carreau_lambda( double d ) { M_carreau_lambda=d; }
    double carreau_n() const { return M_carreau_n; }
    void carreau_n( double d ) { M_carreau_n=d; }
    double carreauYasuda_lambda() const { return M_carreauYasuda_lambda; }
    void carreauYasuda_lambda( double d ){ M_carreauYasuda_lambda=d; }
    double carreauYasuda_n() const { return M_carreauYasuda_n; }
    void carreauYasuda_n( double d ) { M_carreauYasuda_n=d; }
    double carreauYasuda_a() const { return M_carreauYasuda_a; }
    void carreauYasuda_a( double d ) { M_carreauYasuda_a=d; }

    double nonNewtonianHematocrit() const { return M_non_newtonian_hematocrit; }
    void nonNewtonianHematocrit(double d) { M_non_newtonian_hematocrit=d; }
    double nonNewtonianTPMA() const { return M_non_newtonian_TPMA; }
    void nonNewtonianTPMA(double d) { M_non_newtonian_TPMA=d; }

    double walburnSchneck_C1() const { return M_walburnSchneck_C1; }
    void walburnSchneck_C1(double d) { M_walburnSchneck_C1=d; }
    double walburnSchneck_C2() const { return M_walburnSchneck_C2; }
    void walburnSchneck_C2(double d) { M_walburnSchneck_C2=d; }
    double walburnSchneck_C3() const { return M_walburnSchneck_C3; }
    void walburnSchneck_C3(double d) { M_walburnSchneck_C3=d; }
    double walburnSchneck_C4() const { return M_walburnSchneck_C4; }
    void walburnSchneck_C4(double d) { M_walburnSchneck_C4=d; }


    boost::shared_ptr<std::ostringstream>
    getInfo() const
        {
#if 0

            if ( this->name() == "newtonian" )
                M_stressTensorLaw  = "newtonian";
            else if ( this->name() == "power_law" )
                M_stressTensorLaw  = "power_law";
            else if ( this->name() == "walburn-schneck_law" )
                M_stressTensorLaw  = "walburn-schneck_law";
            else if ( this->name() == "carreau_law")
                M_stressTensorLaw  = "carreau_law";
            else if (  this->name() == "carreau-yasuda_law")
                M_stressTensorLaw  = "carreau-yasuda_law";
            else
                CHECK( false ) << "invalid stressTensorLawType\n";
#endif
            boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
            if ( this->name() ==  "newtonian" )
            {
            }
            else if ( this->name() == "power_law" )
            {
                *_ostr << "\n   Power Law "
                       << "\n     -- n   : " << this->powerLaw_n()
                       << "\n     -- k   : " << this->powerLaw_k();
            }
            else if ( this->name() == "walburn-schneck_law" )
            {
                *_ostr << "\n   Walburn-Schneck "
                       << "\n     -- C1  : " << this->walburnSchneck_C1()
                       << "\n     -- C2  : " << this->walburnSchneck_C2()
                       << "\n     -- C3  : " << this->walburnSchneck_C3()
                       << "\n     -- C4  : " << this->walburnSchneck_C4()
                       << "\n     -- hematocrit  : " << this->nonNewtonianHematocrit()
                       << "\n     -- TPMA  : " << this->nonNewtonianTPMA();
            }
            else if ( this->name() == "carreau_law")
            {
                *_ostr << "\n   Carreau "
                       << "\n     -- lambda  : " << this->carreau_lambda()
                       << "\n     -- n       : " << this->carreau_n()
                       << "\n     -- mu_0    : " << this->mu_0()
                       << "\n     -- mu_inf  : " << this->mu_inf();
            }
            else if (  this->name() == "carreau-yasuda_law")
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

private :
    std::string M_viscosityModelName;
    boost::reference_wrapper<const element_type> M_muP0;
    double M_powerLaw_n, M_powerLaw_k, M_powerLaw_n_generic, M_powerLaw_k_generic;
    double M_mu_0, M_mu_inf;
    double M_carreau_lambda, M_carreau_n;
    double M_carreauYasuda_lambda, M_carreauYasuda_n, M_carreauYasuda_a;

    double M_walburnSchneck_C1,M_walburnSchneck_C2,M_walburnSchneck_C3,M_walburnSchneck_C4;
    double M_non_newtonian_hematocrit;
    double M_non_newtonian_TPMA; //Total Proteins Minus Albumin (TPMA)

};



template<class SpaceType>
ViscosityModelDescription<SpaceType>
viscosityModelDesc( std::string const& _viscosityModel,
                    typename SpaceType::element_type const& muP0,
                    double _powerLaw_n, double _powerLaw_k, double _mu_0, double _mu_inf,
                    double _carreau_lambda, double _carreau_n,
                    double _carreauYasuda_lambda, double _carreauYasuda_n, double _carreauYasuda_a,
                    double _walburnSchneck_C1, double _walburnSchneck_C2, double _walburnSchneck_C3, double _walburnSchneck_C4,
                    double _non_newtonian_hematocrit, double _non_newtonian_TPMA )
{
    ViscosityModelDescription<SpaceType> res(_viscosityModel, muP0,
                                             _powerLaw_n,_powerLaw_k,_mu_0,_mu_inf,
                                             _carreau_lambda,_carreau_n,
                                             _carreauYasuda_lambda,_carreauYasuda_n,_carreauYasuda_a,
                                             _walburnSchneck_C1, _walburnSchneck_C2, _walburnSchneck_C3, _walburnSchneck_C4,
                                             _non_newtonian_hematocrit, _non_newtonian_TPMA);
    return res;
}

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_VISCOSITY_MODEL_DESCRIPTION_H
