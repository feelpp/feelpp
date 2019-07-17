/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 juin 2015
 
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
#define BOOST_TEST_MODULE test_integrals
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/norm2.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/cross.hpp>
#include <feel/feelvf/print.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Laplacian Options" );
    options.add_options()
        ( "N",po::value<int>()->default_value( 10 ),"number of integrals to compute" )
        ( "f1",po::value<std::string>()->default_value( "" ),"test function 1D" )
        ( "f2",po::value<std::string>()->default_value( "" ),"test function 2D" )
        ( "f3",po::value<std::string>()->default_value( "" ),"test function 3D" )
        ;
    options.add( feel_options() );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_integrals" ,
                     "test_integrals" ,
                     "0.2",
                     "nD(n=2,3) test integrals",
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

            auto Xh = Pchv<2>( mesh );
            auto u = Xh->element();
            u.on(_range=elements(mesh),_expr=vec(cst(1.),cst(1.),cst(1.)));
            
            Eigen::Matrix<double,3,1> n;
            n << 1,2,3;
            
            std::vector<Eigen::Matrix<double,3,1>> x( ioption("N"),  n);


            tic();
            auto v = integrate(_range=elements(mesh),_expr=_e1v,_quad=_Q<1>()).evaluate( x );
            if ( Environment::numberOfProcessors() == 1 )
                std::for_each( v.begin(), v.end(), 
                               [&n]( Eigen::Matrix<double,3,1> const& nn ) 
                               { 
                                   BOOST_CHECK_SMALL( (n-nn).norm(), 1e-11); 
                               });
            toc("lambda integral _e1=vec(x,y,z)",FLAGS_v>0);
            
            int i=0;
            for( auto e : x )
            {
                tic();
                auto vv = integrate(_range=elements(mesh),_expr=vec(cst(e(0)),cst(e(1)),cst(e(2))),_quad=_Q<1>()).evaluate();
                toc("integral vec(x,y,z)",FLAGS_v>0);
                BOOST_CHECK_SMALL( (vv-v[i]).norm(), 1e-11);
                ++i;
            }

            tic();
            v = integrate(_range=elements(mesh),_expr=cross(idv(u),_e1v-P()),_quad=_Q<3>()).evaluate( x );
            toc("lambda integral u x (_1-X)", FLAGS_v>0);
            i = 0;
            for( auto e : x )
            {
                tic();
                auto vv = integrate(_range=elements(mesh),_expr=cross(idv(u),vec(cst(e(0)),cst(e(1)),cst(e(2)))-P()),_quad=_Q<3>()).evaluate();
                toc("integral u x ((x,y,z)-X) ",FLAGS_v>0);
                BOOST_CHECK_SMALL( (v[i]-vv).norm(), 1e-10 );
                ++i;
            }
            tic();
            auto mid = mat<3,3>( cst(1.), cst(0.), cst(0.),
                                 cst(0.), cst(1.), cst(0.),
                                 cst(0.), cst(0.), cst(1.) );
            v = integrate(_range=elements(mesh),_expr=cross(trans(trans(idv(u))*mid),_e1v-P()),_quad=_Q<3>()).evaluate( x );
            toc("lambda integral u x (_1-X)", FLAGS_v>0);
            i = 0;
            for( auto e : x )
            {
                tic();
                auto vv = integrate(_range=elements(mesh),_expr=cross(idv(u),vec(cst(e(0)),cst(e(1)),cst(e(2)))-P()),_quad=_Q<3>()).evaluate();
                toc("integral u x ((x,y,z)-X) ",FLAGS_v>0);
                BOOST_CHECK_SMALL( (v[i]-vv).norm(), 1e-10 );
                ++i;
            }
        }
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( inner_suite )


BOOST_AUTO_TEST_CASE( test_3 )
{
    Test<3> test;
    test.run();
}


BOOST_AUTO_TEST_SUITE_END()

