/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 Aug 2015
 
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
#define BOOST_TEST_MODULE test_convolve
#include <feel/feelcore/testsuite.hpp>


#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelvf/convolve.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Options" );
    options.add_options()
        ( "N",po::value<int>()->default_value( 10 ),"N" )
        ;
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_convolve" ,
                     "test_convolve" ,
                     "0.2",
                     "test_convolve",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
    
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
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % soption(_name="gmsh.domain.shape") % Dim ).str() ,
                                                      _dim=Dim ) );

            auto Xh = Pch<1>( mesh );
            auto u = Xh->element();
            auto f=expr(soption("functions.f"));
            u.on(_range=elements(mesh),_expr=f);
            tic();
            auto v = convolve(_range=elements(mesh), _expr=idv(u)*inner(_e1v-P()),_space=Xh);
            toc("convolve", FLAGS_v>0);
            auto e = exporter( _mesh=mesh );
            e->add("v",v);
            e->save();
        }
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( test_convolve )


BOOST_AUTO_TEST_CASE( test )
{
    Test<3> test;
    test.run();
}


BOOST_AUTO_TEST_SUITE_END()


