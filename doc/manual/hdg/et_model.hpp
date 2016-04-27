/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Daniele Prada <daniele.prada85@gmail.com>
       Date: 2016-02-10

  Copyright (C) 2016 Feel++ Consortium

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
//#include <feel/feel.hpp>
#include "mixedpoisson.hpp"

namespace Feel {

inline
po::options_description
makeETHDGOptions()
{
    po::options_description options ( "Electro-Thermal options");
    options.add_options()
        ( "picard.itol", po::value<double>()->default_value( 1e-14 ), "tolerance" )
        ( "picard.itmax", po::value<int>()->default_value( 10 ), "iterations max" )
        ;
    options.add( makeMixedPoissonOptions("E"));
    options.add( makeMixedPoissonOptions("T"));
    return options;
}

inline
po::options_description
makeETHDGLibOptions()
{
    po::options_description options ( "Electro-Thermal lib options");
    options.add( makeMixedPoissonLibOptions("E"));
    options.add( makeMixedPoissonLibOptions("T"));
    return options;
}

template<int Dim, int Order>
class ElectroThermal
{
public:
    using convex_type = Simplex<Dim,1>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    using electro_type = MixedPoisson<Dim,Order>;
    using thermal_type = MixedPoisson<Dim,Order>;

    using electro_flux_ptrtype = typename electro_type::Vh_element_ptr_t;
    using electro_potential_ptrtype = typename electro_type::Wh_element_ptr_t;
    using thermal_flux_ptrtype = typename thermal_type::Vh_element_ptr_t;
    using thermal_potential_ptrtype = typename thermal_type::Wh_element_ptr_t;

private:
    mesh_ptrtype M_mesh;

    electro_type M_ElectroModel;
    thermal_type M_ThermalModel;

    electro_flux_ptrtype M_j;
    electro_potential_ptrtype M_V;
    thermal_flux_ptrtype M_GT;
    thermal_potential_ptrtype M_T;

public:
    ElectroThermal();
    void run();
};

template<int Dim, int Order>
ElectroThermal<Dim, Order>::ElectroThermal()
{
    M_mesh = loadMesh( new mesh_type );

    M_ElectroModel = electro_type( "E" );
    M_ThermalModel = thermal_type( "T" );
    M_ElectroModel.init( M_mesh);
    M_ThermalModel.init( M_mesh);
}

template<int Dim, int Order>
void
ElectroThermal<Dim, Order>::run()
{
    auto electroModelProp = M_ElectroModel.modelProperties();
    auto thermalModelProp = M_ThermalModel.modelProperties();

    auto Vo = M_ElectroModel.potentialSpace()->element();
    auto To = M_ThermalModel.potentialSpace()->element();
    auto tol = doption("picard.itol");
    auto itmax = ioption("picard.itmax");
    auto i = 0;
    double incrV = 0, incrT = 0;

    M_ElectroModel.solve();
    M_j = M_ElectroModel.fluxField();
    M_V = M_ElectroModel.potentialField();

    M_ThermalModel.assembleLinear();
    for( auto const& pairMat : electroModelProp.materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        auto expr = inner(idv(*M_j))/material.getScalar("sigma0");
        M_ThermalModel.updatePotentialRHS(expr, marker);
    }
    for( auto const& pairMat : thermalModelProp.materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        auto expr = material.getScalar("k0");
        M_ThermalModel.updateConductivity(expr, marker);
    }

    M_ThermalModel.solveNL();
    M_T = M_ThermalModel.potentialField();
    M_GT = M_ThermalModel.fluxField();

    cout << " #iteration incrV incrT";
    for ( auto const& marker : M_ElectroModel.integralMarkersList())
        cout << " current_" << marker << std::endl;

    incrV = normL2( _range=elements(M_mesh), _expr=idv(*M_V)-idv(Vo) );
    incrT = normL2( _range=elements(M_mesh), _expr=idv(*M_T)-idv(To) );
    Vo = *M_V;
    To = *M_T;

    cout << " picard #" << i << " " << incrV << " " << incrT;

    for ( auto const& marker : M_ElectroModel.integralMarkersList())
    {
        double I = integrate( _range=markedfaces(M_mesh, marker ),
                              _expr=inner(idv(*M_j),N()) ).evaluate()(0,0);
        cout << " " << I << std::endl;
    }

    while ( ( incrV > tol || incrT > tol ) && ( ++i < itmax ) )
    {
        M_ElectroModel.assembleLinear();

        for( auto const& pairMat : electroModelProp.materials() )
        {
            auto marker = pairMat.first;
            auto material = pairMat.second;
            auto alpha = material.getDouble("alpha");
            auto T0 = material.getDouble("T0");
            auto sigma0 = material.getDouble("sigma0");
            auto sigma = material.getScalar("sigma", "T", idv(*M_T));
            sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            M_ElectroModel.updateConductivity( sigma, marker);
        }

        M_ElectroModel.solveNL();
        M_j = M_ElectroModel.fluxField();
        M_V = M_ElectroModel.potentialField();

        M_ThermalModel.assembleLinear();
        for( auto const& pairMat : electroModelProp.materials() )
        {
            auto marker = pairMat.first;
            auto material = pairMat.second;
            auto alpha = material.getDouble("alpha");
            auto T0 = material.getDouble("T0");
            auto sigma0 = material.getDouble("sigma0");
            auto sigma = material.getScalar("sigma", "T", idv(*M_T));
            sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            auto expr = inner(idv(*M_j))/sigma;
            M_ThermalModel.updatePotentialRHS(expr, marker);
        }
        for( auto const& pairMat : thermalModelProp.materials() )
        {
            auto marker = pairMat.first;
            auto material = pairMat.second;
            auto alpha = material.getDouble("alpha");
            auto T0 = material.getDouble("T0");
            auto k0 = material.getDouble("k0");
            auto k = material.getScalar("k", "T", idv(*M_T));
            k.setParameterValues({{"k0",k0},{"T0",T0},{"alpha",alpha}});
            M_ThermalModel.updateConductivity(k,marker);
        }

        M_ThermalModel.solveNL();
        M_T = M_ThermalModel.potentialField();
        M_GT = M_ThermalModel.fluxField();

        for ( auto const& marker : M_ElectroModel.integralMarkersList())
            cout << " current_" << marker << std::endl;

        incrV = normL2( _range=elements(M_mesh), _expr=idv(*M_V)-idv(Vo) );
        incrT = normL2( _range=elements(M_mesh), _expr=idv(*M_T)-idv(To) );
        Vo = *M_V;
        To = *M_T;

        cout << " picard #" << i << " " << incrV << " " << incrT;

        for ( auto const& marker : M_ElectroModel.integralMarkersList())
        {
            double I = integrate( _range=markedfaces(M_mesh, marker ),
                                  _expr=inner(idv(*M_j),N()) ).evaluate()(0,0);
            cout << " " << I << std::endl;
        }
    }
}


} // Feel
