//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/feelcrb/crbplugin.hpp>
#include "poissonCRB.hpp"

namespace Feel {

Poisson::Poisson() : super_type("poisson") {}

int Poisson::Qa()
{
    return 2;
}

int Poisson::Nl()
{
    return 1;
}

int Poisson::Ql( int l)
{
    switch(l) {
    case 0:
        return 1;
    default:
        return 0;
    }
}

void Poisson::initModel()
{
    M_modelProps = std::make_shared<prop_type>(Environment::expand( soption("poisson.filename")));
    auto parameters = M_modelProps->parameters();
    Dmu->setDimension(parameters.size());
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        mu_min(i) = parameterPair.second.min();
        mu_max(i) = parameterPair.second.max();
        Dmu->setParameterName(i++, parameterPair.first );
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);

    M_mesh = loadMesh( new mesh_type );
    this->setFunctionSpaces(functionspace_type::New( M_mesh ) );

    M_Aq.resize(Qa(), backend()->newMatrix(Xh, Xh ) );
    M_betaAq.resize(Qa());
    M_Fq.resize(Nl());
    M_betaFq.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        M_Fq[l].resize(Ql(l), backend()->newVector(Xh) );
        M_betaFq[l].resize(Ql(l));
    }

    auto u = Xh->element();
    auto v = Xh->element();
    auto gamma = doption("poisson.gamma");

    // Omega 1 (ext)
    auto a1 = form2(_test=Xh, _trial=Xh);
    a1 = integrate( markedelements(M_mesh, "omega1"),
                    inner(gradt(u),grad(v)) );
    // Dirichlet condition
    a1 += integrate( markedfaces(M_mesh, "top"),
                     gamma/hFace()*inner(idt(u),id(v))
                     -inner(gradt(u)*N(),id(v))
                     -inner(grad(v)*N(),idt(u)) );
    M_Aq[0] = a1.matrixPtr();

    // Omega 0 (int)
    auto a2 = form2( _test=Xh, _trial=Xh);
    a2 = integrate( markedelements(M_mesh, "omega0"),
                    inner(gradt(u),grad(v)) );
    M_Aq[1] = a2.matrixPtr();

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(u),grad(v)) );
    M_energy_matrix = m.matrixPtr();

    // Neumann condition
    auto f1 = form1(_test=Xh);
    f1 = integrate( markedfaces(M_mesh, "base"),
                    id(v) );
    M_Fq[0][0] = f1.vectorPtr();
}

double Poisson::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index == 0 )
    {
        for ( int q = 0; q < Ql(0); q++ )
        {
            element_ptrtype eltF( new element_type( Xh ) );
            *eltF = *M_Fq[output_index][q];
            output += M_betaFq[output_index][q]*dot( *eltF, u );
            //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

auto Poisson::computeBetaQ( parameter_type const& mu ) ->
    boost::tuple<beta_vector_light_type, std::vector<beta_vector_light_type> >
{
    M_betaAq[0] = 1;
    M_betaAq[1] = mu.parameterNamed("kappa");
    M_betaFq[0][0] = 1;
    M_betaFq[0][1] = mu.parameterNamed("flux");
    return boost::make_tuple( M_betaAq, M_betaFq );
}

FEELPP_CRB_PLUGIN( Poisson, poisson )
}
