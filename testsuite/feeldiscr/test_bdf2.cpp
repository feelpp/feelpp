/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-12-28

  Copyright (C) 2011 - 2014 Feel++ Consortium

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
/**
   \file test_bdf2.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-12-28
*/

//#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_bdf2
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>

/** use Feel namespace */
using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_bdf2" ,
                     "test_bdf2" ,
                     "0.2",
                     "nD(n=2,3) test bdf2",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}

template<int Dim>
class Test:
    public Simget
{
public :

    void run()
    {
        auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim,1>>,
                                    _desc=domain( _name=( boost::format( "%1%-%2%" ) % option(_name="gmsh.domain.shape").template as<std::string>() % Dim ).str() ,
                                                  _dim=Dim ) );


        auto Xh = Pch<1>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();
        auto solution = Xh->element();
        auto mybdf = bdf( _space=Xh, _name="mybdf" );
        double mu0=0.5;
        //stifness matrix
        auto a = form2( _test=Xh, _trial=Xh );
        a = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) + mybdf->polyDerivCoefficient(0)*idt(u)*id(u));
        //we have a robin condition
        a += integrate( _range= markedfaces( mesh, "Dirichlet" ), _expr= mu0 * idt( u )*id( v ) );
        //mass matrix
        auto m = form2( _test=Xh, _trial=Xh);
        m = integrate( _range= elements( mesh ), _expr= idt( u )* id( v ) );
        //Rhs
        auto f = form1( Xh );
        f = integrate( _range=markedfaces( mesh,"Dirichlet" ), _expr= mu0 * id( v ) );

        auto ft = form1(_test=Xh);

        double bdf_coeff;
        auto vec_bdf_poly = backend()->newVector( Xh );
        mybdf->start(solution);

        for ( ;  mybdf->isFinished() == false; mybdf->next(solution) )
        {
            auto bdf_poly = mybdf->polyDeriv();
            ft.zero();
            ft = integrate( _range=elements(mesh), _expr=idv(bdf_poly)*id(u) );
            ft += f;
            a.solve( _solution=solution, _rhs=ft );
        }

        //check that we obtain the same result in sequential or in parallel
        double l2_inner_prod = m( solution, solution );
        if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            std::cout<<"l2_inner_prod : "<<std::setprecision(16) <<l2_inner_prod<<std::endl;
    }
};




#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( bdf2 )

BOOST_AUTO_TEST_CASE( test_1 )
{
    Test<2> test;
    test.run();
}

BOOST_AUTO_TEST_SUITE_END()
#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=feel_options() );
    Test<2>  test;
    test.run();
}
#endif





