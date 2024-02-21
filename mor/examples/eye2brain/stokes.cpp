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
    po::options_description options( "Eye2Brain_Stokes" );
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
Eye2Brain_Stokes<Order, Dim>::Eye2Brain_Stokes(): super_type( fmt::format("Eye2Brain_Stokes_{}D_P{}", Dim, Order) ) {}

template<int Order, int Dim>
Eye2Brain_Stokes<Order, Dim>::Eye2Brain_Stokes(mesh_ptrtype mesh): super_type( fmt::format("Eye2Brain_Stokes_{}D_P{}", Dim, Order) )
{
    this->M_mesh = mesh;
}


template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::Qa() const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::Nl() const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::Ql(int l) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::mQA(int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::mLQF(int l, int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::mCompliantQ(int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::mIntensityQ(int q) const
{
    // TODO
    return 1;
}

template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::mAverageTempQ(int q) const
{
    // TODO
    return 1;
}


template<int Order, int Dim>
int Eye2Brain_Stokes<Order, Dim>::resizeQm(int q) const
{
    // TODO
    return 1;
}


template<int Order, int Dim>
void Eye2Brain_Stokes<Order, Dim>::resizeQm( bool resizeMatrix )
{
    // TODO
}


template<int Order, int Dim>
void Eye2Brain_Stokes<Order, Dim>::initModel()
{
    Feel::cout << "initModel" << std::endl; LOG(INFO) << "initModel" << std::endl;


}


template<int Order, int Dim>
void Eye2Brain_Stokes<Order, Dim>::decomposition()
{
    auto UP = Xh->element();
    auto u = up.template element<0>();
    auto p = up.template element<1>();

    auto PHI = Xh->element();
    auto v = phi.template element<0>();
    auto q = phi.template element<1>();


    // Left hand side

}


}  // namespace Feel