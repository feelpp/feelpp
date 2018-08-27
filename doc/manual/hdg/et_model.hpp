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
namespace FeelModels{

inline
po::options_description
makeETHDGOptions()
{
    po::options_description options ( "Electro-Thermal options");
    options.add_options()
        ( "picard.itol", po::value<double>()->default_value( 1e-14 ), "tolerance" )
        ( "picard.itmax", po::value<int>()->default_value( 10 ), "iterations max" )
        ( "linear", po::value<bool>()->default_value( false ), "update conductivities or not")
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
    return options.add( feel_options());
}

template<int Dim, int Order>
class ElectroThermal
{
public:
    using convex_type = Simplex<Dim,1>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

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

    void updateElectroAssembly( sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
    void updateThermalAssembly( sparse_matrix_ptrtype& A, vector_ptrtype& F) const;

public:
    ElectroThermal() : M_ElectroModel("E"), M_ThermalModel("T") {
        M_mesh = loadMesh( new mesh_type );
        M_ElectroModel.init( M_mesh);
        M_ThermalModel.init( M_mesh);
    }
    void run();
};

template<int Dim, int Order>
void
ElectroThermal<Dim, Order>::updateElectroAssembly( sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    auto Vh = M_ElectroModel.fluxSpace();
    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );

    auto a11 = form2( _trial=Vh, _test=Vh,_matrix=A );

    auto electroModelProp = M_ElectroModel.modelProperties();

    for( auto const& pairMat : electroModelProp.materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        auto alpha = material.getDouble("alpha");
        auto T0 = material.getDouble("T0");
        auto sigma0 = material.getDouble("sigma0");
        auto sigma = material.getScalar("sigma", "T", idv(*M_T));
        sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
        a11 += integrate( _range=markedelements(M_mesh, marker), _expr=inner(idt(u),id(v))/sigma);
    }
}

template<int Dim, int Order>
void
ElectroThermal<Dim, Order>::updateThermalAssembly( sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    auto Vh = M_ThermalModel.fluxSpace();
    auto Wh = M_ThermalModel.potentialSpace();
    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto w = Wh->element();

    auto a11 = form2( _trial=Vh, _test=Vh,_matrix=A );
    auto rhs = form1( _test=Wh, _vector=F, _rowstart=1);

    auto thermalModelProp = M_ThermalModel.modelProperties();
    auto electroModelProp = M_ElectroModel.modelProperties();

    for( auto const& pairMat : electroModelProp.materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        auto alpha = material.getDouble("alpha");
        auto T0 = material.getDouble("T0");
        auto sigma0 = material.getDouble("sigma0");
        auto sigma = material.getScalar("sigma", {"T"}, {idv(*M_T)}, {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
        auto expr = inner(idv(*M_j))/sigma;
        rhs += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(w)) );
    }
    for( auto const& pairMat : thermalModelProp.materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        auto alpha = material.getDouble("alpha");
        auto T0 = material.getDouble("T0");
        auto k0 = material.getDouble("k0");
        auto k = material.getScalar("k", {"T"}, {idv(*M_T)}, {{"k0",k0},{"T0",T0},{"alpha",alpha}});
        a11 += integrate( _range=markedelements(M_mesh, marker), _expr=inner(idt(u),id(v))/k);
    }
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

    auto e = exporter( M_mesh);

    M_ElectroModel.solve();
    M_j = M_ElectroModel.fluxField();
    M_V = M_ElectroModel.potentialField();

    M_ThermalModel.assembleF();
    M_ThermalModel.updateConductivityTerm();
    for( auto const& pairMat : electroModelProp.materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        auto expr = inner(idv(*M_j))/material.getScalar("sigma0");
        M_ThermalModel.updatePotentialRHS(expr, marker);
    }

    M_ThermalModel.solveNL();
    M_T = M_ThermalModel.potentialField();
    M_GT = M_ThermalModel.fluxField();

    e->step(i)->add("potential", *M_V);
    e->step(i)->add("current", *M_j);
    e->step(i)->add("temperature", *M_T);
    e->step(i)->add("thermal-flux", *M_GT);
    e->save();

    cout << " #iteration incrV incrT";
    for ( auto const& marker : M_ElectroModel.integralMarkersList())
        cout << " current_" << marker << std::endl;

    incrV = normL2( _range=elements(M_mesh), _expr=idv(*M_V)-idv(Vo) );
    incrT = normL2( _range=elements(M_mesh), _expr=idv(*M_T)-idv(To) );
    Vo = *M_V;
    To = *M_T;

    cout << " picard #" << i << " " << incrV << " " << incrT;

    auto itField = electroModelProp.boundaryConditions().find( "flux");
    if ( itField != electroModelProp.boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                double meas = integrate( _range=markedfaces(M_mesh,marker), _expr=cst(1.0)).evaluate()(0,0);
                double I_target = boost::lexical_cast<double>(exAtMarker.expression());
                double I = integrate( _range=markedfaces(M_mesh, marker ),
                                      _expr=inner(idv(*M_j),N()) ).evaluate()(0,0);
                cout << " " << std::abs(I-I_target) << std::endl;
            }
        }
    }

    M_ElectroModel.M_updateAssembly = boost::bind( &ElectroThermal<Dim, Order>::updateElectroAssembly,
                                                   this, _1, _2 );
    M_ThermalModel.M_updateAssembly = boost::bind( &ElectroThermal<Dim, Order>::updateThermalAssembly,
                                                   this, _1, _2 );

    while ( ( incrV > tol || incrT > tol ) && ( ++i < itmax ) )
    {
        M_ElectroModel.assembleF();
        M_ElectroModel.solveNL();
        M_j = M_ElectroModel.fluxField();
        M_V = M_ElectroModel.potentialField();

        M_ThermalModel.assembleF();
        M_ThermalModel.solveNL();
        M_T = M_ThermalModel.potentialField();
        M_GT = M_ThermalModel.fluxField();

        e->step(i)->add("potential", *M_V);
        e->step(i)->add("current", *M_j);
        e->step(i)->add("temperature", *M_T);
        e->step(i)->add("thermal-flux", *M_GT);
        e->save();

        incrV = normL2( _range=elements(M_mesh), _expr=idv(*M_V)-idv(Vo) );
        incrT = normL2( _range=elements(M_mesh), _expr=idv(*M_T)-idv(To) );
        Vo = *M_V;
        To = *M_T;

        cout << " picard #" << i << " " << incrV << " " << incrT;

        auto itField = electroModelProp.boundaryConditions().find( "flux");
        if ( itField != electroModelProp.boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    double meas = integrate( _range=markedfaces(M_mesh,marker), _expr=cst(1.0)).evaluate()(0,0);
                    double I_target = boost::lexical_cast<double>(exAtMarker.expression());
                    double I = integrate( _range=markedfaces(M_mesh, marker ),
                                          _expr=inner(idv(*M_j),N()) ).evaluate()(0,0);
                    cout << " " << std::abs(I-I_target) << std::endl;
                }
            }
        }
    } // Picard Loop

}

} // FeelModels
} // Feel
