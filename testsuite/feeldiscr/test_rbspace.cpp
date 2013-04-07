/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-04-07

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
   \file test_rbspace.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-07
*/

#define BOOST_TEST_MODULE test_rbspace
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feeldiscr/reducedbasisspace.hpp>

/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testrbspace( "RBSpace test options" );
    testrbspace.add_options()
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testrbspace.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_rbspace" ,
                     "test_rbspace" ,
                     "0.2",
                     "nD(n=1,2,3) test context of functionspace",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}



template<int Dim, int Order>
class Model:
    public Simget,
    public boost::enable_shared_from_this< Model<Dim,Order> >
{

public :

    Model()
        :
        Simget()
    {
    }

    void run()
    {
        auto mesh=unitHypercube<Dim>();
        auto Xh = Pch<Order>( mesh );
        LOG(INFO)<<"nDof : "<<Xh->nDof();

        auto RbSpace = RbSpacePch<Order>( this , mesh );
    }


};


/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( rbspace )

BOOST_AUTO_TEST_CASE( test_1 )
{
    Application app;
    app.add( new Model<2,1>() );
    app.run();
}

BOOST_AUTO_TEST_SUITE_END()


