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
//! @author Thomas Saigre @thomas-saigre
//! @date 20 Fev 2024
//! @copyright 2024 Feel++ Consortium
//!
//!

#include <feel/feelmor/crbplugin.hpp>
#include "stokes.hpp"


namespace Feel
{

FEELPP_EXPORT po::options_description
makeEye2Brain_StokesOptions()
{
    po::options_description options( "Eye2BrainStokes" );
    return options;
}

FEELPP_EXPORT AboutData
makeEye2Brain_StokesAbout( std::string const& str = "eye2brain-stokes" )
{
    AboutData about( /*AppName  */ str.c_str(),
                     /*ProgName */ str.c_str(),
                     /*Version  */ "0.1",
                     /*ShortDesc*/ "Eye2Brain 3D Heat Application",
                     /*Licence  */ AboutData::License_GPL,
                     /*Copyright*/ "Copyright (c) 2024 Feel++ Consortium" );
    return about;
}

template<int Order, int Dim>
Eye2BrainStokes<Order, Dim>::Eye2BrainStokes(): super_type( fmt::format("Eye2BrainStokes_{}D_P{}", Dim, Order) ) {}

template<int Order, int Dim>
Eye2BrainStokes<Order, Dim>::Eye2BrainStokes(mesh_ptrtype mesh): super_type( fmt::format("Eye2BrainStokes_{}D_P{}", Dim, Order) )
{
    this->M_mesh = mesh;
}


template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::Qa() const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::Nl() const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::Ql(int l) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::mQA(int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::mLQF(int l, int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::mCompliantQ(int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::mIntensityQ(int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2BrainStokes<Order, Dim>::mAverageTempQ(int q) const
{
    // TODO
    return 1;
}


template<int Order, int Dim>
void Eye2BrainStokes<Order, Dim>::resizeQm( bool resizeMatrix )
{
    // TODO
}


template<int Order, int Dim>
void Eye2BrainStokes<Order, Dim>::initModel()
{
    Feel::cout << "initModel" << std::endl; LOG(INFO) << "initModel" << std::endl;

    M_mesh = loadMesh( _mesh = new mesh_type,
                       _update = MESH_UPDATE_FACES | MESH_UPDATE_EDGES | MESH_NO_UPDATE_MEASURES,
                       _savehdf5 = 0 );

    this->setFunctionSpaces( functionspace_type::New( _mesh=M_mesh ));
    this->setSymbolicExpressionBuildDir("$repository/crb/Eye2BrainStokes/symbolicexpr/");

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    // Trial function
    auto UP = Xh->element();
    auto u = UP.template element<0>();
    auto p = UP.template element<1>();

    // Test function
    auto PHI = Xh->element();
    auto v = PHI.template element<0>();
    auto q = PHI.template element<1>();

    double muAH_ref = 1080, muAH_min = 1000, muAH_max = 2000;     // random values
    std::vector<double> muRef = {muAH_ref, 1};
    auto energy = backend()->newMatrix( _test = this->Xh, _trial = this->Xh );


    this->Dmu->setDimension( 2 );
    auto mu_min = this->Dmu->element();
    mu_min << muAH_min, 1;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << muAH_max, 1;
    this->Dmu->setMax( mu_max );
    LOG(INFO) << "[Eye2brainStokes::initModel] Set parameter space : mu_min " << mu_min << "\n" << std::endl;
    LOG(INFO) << "[Eye2brainStokes::initModel] Set parameter space : mu_max " << mu_max << "\n" << std::endl;


    // Left hand side
    auto a0 = form2( _test  =Xh, _trial = Xh);
    a0 = integrate( _range = elements(M_mesh),
                    _expr = inner( 2 * sym(gradt(u)), grad(v) ) );
    a0.matrixPtr()->close();
    this->addLhs( {a0, "mu0"} );
    energy->addMatrix( muRef[0], a0.matrixPtr() );

    auto a1 = form2( _test = Xh, _trial = Xh);
    a1 = integrate( _range = elements(M_mesh),
                    _expr = -idt(p) * div(v) );
    a1 += integrate( _range = elements(M_mesh),
                     _expr = id(q) * divt(u) );
    a1.matrixPtr()->close();
    this->addLhs( {a1, "mu1"} );
    energy->addMatrix( muRef[1], a1.matrixPtr() );

    // Right hand side
    auto f = expr<Dim,1>( soption(_name="functions.f") );
    auto l = form1( _test = Xh );
    l = integrate( _range = elements(M_mesh),
                   _expr = inner( expr<Dim,1>( f ), id(v) ) );
    l.vectorPtr()->close();
    this->addRhs( {l, "mu3"} );

    // Energy
    energy->close();
    this->addEnerfyMatrix( energy );

    // Output: for now avrage pressure
    auto avgP = form1( _test = Xh );
    double meas = integrate( _range = elements(M_mesh), _expr = cst(1.) ).evaluate()(0,0);
    avgP = integrate( _range = elements(M_mesh), _expr = id(p) / cst(meas) );
    this->addOutput( {avgP, "1"} );

}



}  // namespace Feel