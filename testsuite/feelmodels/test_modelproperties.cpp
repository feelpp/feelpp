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
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv );

    auto mesh = loadMesh(new Mesh<Simplex<3> >);
    auto Xh = Pch<1>(mesh);
    auto g = Xh->element();
    auto d = Xh->element();

    ModelProperties model_props( "test.feelpp" );

    auto mats = model_props.materials();
    for ( auto matPair : mats )
    {
        cout << "properties for " << matPair.first << std::endl;
        auto mat = matPair.second;
        auto name = mat.getString("name");
        auto rhoInt = mat.getInt("rho");
        auto etaDouble = mat.getDouble("eta");

        auto rho = mat.getScalar( "rho" );
        auto f = mat.getScalar("f","g",idv(g));
        auto fParam = mat.getScalar("f",{{"g",1.}});
        auto hList = mat.getScalar("h",{"g","t"},{idv(g),idv(d)});
        auto fVec = mat.getScalar("f",std::vector<std::string>(1,"g"), std::vector<decltype(idv(g))>(1,idv(g)));
        auto hParams = mat.getScalar("h",{"g"},{idv(g)}, {{"t",1.}});

        auto nu = mat.getVector<3>( "nu" );
        auto curlnu = curl(nu);
        auto muPair = mat.getVector<2>("mu",{"t",1.});
        auto muMap = mat.getVector<2>("mu", {{"t",1.}});

        auto chi = mat.getMatrix<2>( "chi" );
        auto xhiPair = mat.getMatrix<3>( "xhi", {"t",2.} );
        auto xhiMap = mat.getMatrix<3,3>( "xhi", {{"t",3.}});

        cout << "\t" << name << std::endl;
        cout << "\t" << rhoInt << std::endl;
        cout << "\t" << etaDouble << std::endl;
        cout << "\t" << rho << std::endl;
        cout << "\t" << nu << std::endl;
        cout << "\t" << curlnu << std::endl;
        cout << "\t" << muPair << std::endl;
        cout << "\t" << muMap << std::endl;
        cout << "\t" << chi << std::endl;
        cout << "\t" << xhiPair << std::endl;
        cout << "\t" << xhiMap << std::endl;
    }
}
