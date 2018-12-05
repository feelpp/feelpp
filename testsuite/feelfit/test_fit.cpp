/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil;
   c-basic-offset: 4; show-trailing-whitespace: t -*-
vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2014-01-14

Copyright (C) 2014-2016 Feel++ Consortium

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define BOOST_TEST_MODULE test_fit
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfit/fit.hpp>



FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testfit_suite )

BOOST_AUTO_TEST_CASE( test_main_fit )
{
    using namespace Feel;

    // tag::fit[]
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Xh = Pch<1>(mesh);
    auto T = Xh->element(); // the temperature (say)
    auto K = Xh->element(); // K(T) - the dependant of the temperature conductivity
    auto Kd= Xh->element(); // K'(T)
    auto T_K = Xh->element(); // T_ = true 
    auto T_Kd= Xh->element(); // 
    auto D_K = Xh->element(); // D_ = difference
    auto D_Kd= Xh->element(); // 
    T.on(_range=elements(mesh), _expr=Px());
    T_K.on(_range=elements(mesh),_expr=(5*Px()+sin(Px())));
    T_Kd.on(_range=elements(mesh),_expr=(5+cos(Px())));
    //double f(double t = 0.) { return 5. * t + sin(t); }
    auto f = [](double x = 0.) { return 5. * x + sin(x); };
#if 1
    auto e = exporter(_mesh = mesh );
#endif
    std::string datafilename = (fs::current_path()/"data.txt").string();
    if ( Environment::worldComm().isMasterRank() )
    {
        // Generates the datafile
        // we assume an unitsquare as mesh
        std::ofstream datafile( datafilename );
        for(double t = -1; t < 2; t+=0.32)
            datafile << t << "\t" << f(t) << "\n";
        datafile.close();
    }
    Environment::worldComm().barrier();

    std::vector<std::string> interpTypeRange = { "P0" , "P1", "Spline", "Akima" };
    for(int k = 0; k < interpTypeRange.size(); ++k )
    {
        std::string const& interpType = interpTypeRange[k];
        BOOST_TEST_MESSAGE( boost::format("test %1%")% interpType );
        // evaluate K(T) with the interpolation from the datafile
        K.on(_range=elements(mesh), _expr=fit( idv(T), datafilename, interpType ) );
        Kd.on(_range=elements(mesh), _expr=fitDiff( idv(T), datafilename, interpType ) );

        D_K.on(_range=elements(mesh),_expr=vf::abs(idv(K)-idv(T_K)));
        D_Kd.on(_range=elements(mesh),_expr=vf::abs(idv(Kd)-idv(T_Kd)));

        auto max_K = D_K.max();
        auto max_Kd= D_Kd.max();
#if 1
        e->step(k)->add("T",T);
        e->step(k)->add("K",K);
        e->step(k)->add("Kd",Kd);
        e->step(k)->add("T_K",T_K);
        e->step(k)->add("T_Kd",T_Kd);
        e->step(k)->add("D_K",D_K);
        e->step(k)->add("D_Kd",D_Kd);
        e->save();
#endif
        /// Note : the error has nothing to do with the mesh size but the step on the datafile
        switch( InterpolationTypeMap[interpType] )
        {
        case InterpolationType::P0: //P0 interpolation
        {
            BOOST_CHECK_SMALL(max_K, 0.95);
            BOOST_CHECK_SMALL(max_Kd, 6.0001); // the derivative is null
            break;
        }
        case InterpolationType::P1: // P1 interpolation
        {
            BOOST_CHECK_SMALL(max_K, 0.01);
            BOOST_CHECK_SMALL(max_Kd, 0.15);
            break;
        }
        case InterpolationType::Spline: // CSpline interpolation
        {
            BOOST_CHECK_SMALL(max_K, 0.01);
            BOOST_CHECK_SMALL(max_Kd, 0.15);
            break;
        }
        case InterpolationType::Akima: // Akima interpolation
        {
            BOOST_CHECK_SMALL(max_K, 0.016);
            BOOST_CHECK_SMALL(max_Kd, 0.03);
            break;
        }
        }

        BOOST_TEST_MESSAGE( boost::format("test %1% done")% interpType );
    }
   // end::fit[]
}

BOOST_AUTO_TEST_SUITE_END()


