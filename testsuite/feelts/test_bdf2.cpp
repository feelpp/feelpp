/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-28

  Copyright (C) 2011-2014 Feel++ Consortium

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
//#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE test_bdf2
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
//using Feel::project;

template<int Dim>
class Test:
    public Simget
{
public :

    void run() override
    {
        auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim,1>>,
                                    _desc=domain( _name=( boost::format( "%1%-%2%" ) % soption(_name="gmsh.domain.shape") % Dim ).str() ,
                                                  _dim=Dim ) );

        auto Xh = Pch<2>( mesh );
        auto u = Xh->element();
        auto ue = Xh->element();
        auto v = Xh->element();
        auto solution = Xh->element();
        auto mybdf = bdf( _space=Xh, _name="mybdf" );

        double alpha = doption(_name="parameters.alpha");
        double beta = doption(_name="parameters.beta");

        std::string g = soption(_name="functions.g");
        auto ue_g= expr( g );
        auto fe = -laplacian( ue_g ) + diff( ue_g,"t" );

        auto e = exporter(_mesh=mesh);

        // initialize bdf with with known values prior to mybdf->initialTime()
        // (initialTime included)
        for( auto time : mybdf->priorTimes() )
        {
            if( Environment::worldComm().isMasterRank() )
                std::cout << "Initialize prior times (from timeInitial()) : " << time.second << "s index: " << time.first << "\n";

            ue_g.setParameterValues( {
                    {"t", time.second},
                    {"alpha", alpha},
                    {"beta", beta} } );
            ue = project( _space=Xh, _expr=ue_g );
            mybdf->setUnknown( time.first, ue );
        }

        fe.setParameterValues( {
                {"t", mybdf->timeInitial()},
                {"alpha", alpha},
                {"beta", beta} } );

        solution.on( _range=elements(mesh), _expr=ue_g );
        ue.on(_range=elements(mesh), _expr=ue_g );
        auto error = vf::project( _space=Xh, _expr=idv(ue)-idv(solution) );

        e->step(0)->add("exact",ue);
        e->step(0)->add("solution",solution);
        e->step(0)->add("error",error);

        e->save();
        double maxerror = error.linftyNorm();
        if( Environment::worldComm().isMasterRank() )
        {
            std::cout << "max error at time " << mybdf->timeInitial() << "s :" << std::setprecision(16) << maxerror << "\n";
        }

        for ( mybdf->start();  mybdf->isFinished() == false; mybdf->next(solution) )
        {
            // update time value in expression
            ue_g.setParameterValues( {{"t", mybdf->time()}} );
            fe.setParameterValues( {{"t", mybdf->time()}} );

            auto bdf_poly = mybdf->polyDeriv();

            auto ft = form1(_test=Xh);
            ft = integrate( _range=elements(mesh), _expr=(fe+idv(bdf_poly))*id(u) );
            auto at = form2( _test=Xh, _trial=Xh );
            at = integrate( _range= elements( mesh ),
                            _expr= gradt( u )*trans( grad( v ) ) + mybdf->polyDerivCoefficient(0)*idt(u)*id(u));

            at += on(_range=boundaryfaces(mesh),_element=solution,_rhs=ft,_expr=ue_g);
            at.solve( _solution=solution, _rhs=ft );

            // project onto space for error estimation
            ue = project( _space=Xh, _expr=ue_g );

            // compute max error which should be 0
            auto error = vf::project( _space=Xh, _expr=idv(ue)-idv(solution) );
            double maxerror = error.linftyNorm();
            if( Environment::worldComm().isMasterRank() )
                std::cout << "max error at time " << mybdf->time() << "s :" << std::setprecision(16) << maxerror << "\n";
            BOOST_CHECK_SMALL( maxerror,1e-9 );

            e->step(mybdf->time())->add("exact",ue);
            e->step(mybdf->time())->add("solution",solution);
            e->step(mybdf->time())->add("error",error);
            e->save();
        }
    }
};




FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( bdf2 )

BOOST_AUTO_TEST_CASE( test_1 )
{
    Test<2> test;
    test.run();
}

BOOST_AUTO_TEST_SUITE_END()
