/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 août 2015
 
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
#define BOOST_TEST_MODULE test_interpolation_grad
#include <testsuite/testsuite.hpp>

#include <feel/feel.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Options" );
    options.add_options()
        ( "f",po::value<std::vector<std::string>>()->default_value( {"1"} ),
          "list of functions to test the Grad interpolation operator (default : {\"1\"}" )
        ;
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_interpolation_grad" ,
                     "test_interpolation_grad" ,
                     "0.2",
                     "test_interpolation_grad",
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
            auto mesh = loadMesh( _mesh=new Mesh<Simplex<Dim>>() );
            auto Xh = Pch<1>(mesh);
            auto Gh = Pdhv<0>(mesh);
            auto Igrad = Grad( _domainSpace = Xh, _imageSpace=Gh );
            auto e = exporter(_mesh=mesh);
            
            int i = 0;
            for( auto f : vsoption( _name="f" ) )
            {
                
                std::string n { str(format("u%1%")%i) };
                auto u = Xh->element( expr(f), n, f );
                std::string Ign { str(format("Igrad_u%1%")%i) };
                e->add(n,u);
                e->add(Ign,Igrad(u));

                std::string gn = str(boost::format("grad_u%1%")%i);
                auto v = Gh->element(trans(grad<Dim>(expr(f))), gn, str(grad<Dim>(expr(f))));
                
                e->add(gn,v);
                                     
            }
            e->save();
        }
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( test_interpolation_grad )


BOOST_AUTO_TEST_CASE( test )
{
    Test<2> test;
    test.run();
}


BOOST_AUTO_TEST_SUITE_END()


