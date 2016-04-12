/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 Mar 2015

 Copyright (C) 2015 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <feel/feelcore/environment.hpp>
#include <feel/feelmodels/modelproperties.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv );

    ModelProperties model_props( "test.feelpp" );

    auto mats = model_props.materials();
    for ( auto matPair : mats )
    {
        std::cout << "properties for " << matPair.first << std::endl;
        auto mat = matPair.second;
        auto rho = mat.getScalar( "rho" );
        auto mu = mat.getVector<2>( "mu" );
        auto nu = mat.getVector<3>( "nu" );
        auto curlnu = curl(nu);
        auto chi = mat.getMatrix<2>( "chi" );
        auto xhi = mat.getMatrix<3>( "xhi" );

        std::cout << "\t" << rho << std::endl;
        std::cout << "\t" << mu << std::endl;
        std::cout << "\t" << nu << std::endl;
        std::cout << "\t" << chi << std::endl;
        std::cout << "\t" << xhi << std::endl;
        std::cout << "\t" << curlnu << std::endl;
    }
}
