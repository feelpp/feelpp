/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-11

  Copyright (C) 2010-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/mean.hpp>
#include <feel/feelvf/trans.hpp>
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description thermalfinoptions( "Thermalfin options" );
    thermalfinoptions.add_options()
    // physical coeff
    ( "k0", Feel::po::value<double>()->default_value( 1 ), "k0 diffusion parameter" )
    ( "k1", Feel::po::value<double>()->default_value( 1 ), "k1 diffusion parameter" )
    ( "k2", Feel::po::value<double>()->default_value( 1 ), "k2 diffusion parameter" )
    ( "k3", Feel::po::value<double>()->default_value( 1 ), "k3 diffusion parameter" )
    ( "Bimin", Feel::po::value<double>()->default_value( 0.01 ), "minimum value of Biot number" )
    ( "Bimax", Feel::po::value<double>()->default_value( 1 ), "maximum value of Biot number" )

    ( "N", Feel::po::value<int>()->default_value( 1 ), "number of samples withing parameter space" )

    ;
    return thermalfinoptions;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "thermalfin" ,
                           "thermalfin" ,
                           "0.3",
                           "Heat transfer in a 2D thermal fin",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

int main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );


    Environment::changeRepository( boost::format( "%1%/%2%/" )
                                   % Environment::about().appName()
                                   % option( _name="gmsh.hsize" ).as<double>() );

    /*
     * First we create the mesh
     */
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = Pch<2>( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    // flux a bottom of fin
    auto f = form1( Xh );
    f = integrate( markedfaces( mesh, "Tflux" ), id( v ) );


    double Bimin = doption(_name="Bimin");
    double Bimax = doption(_name="Bimax");
    int N = option(_name="N").as<int>();

    auto e = exporter( _mesh=mesh );
    for ( int i = 0; i < N; ++i )
    {
        int Nb = 2;
        Eigen::MatrixXd k( Nb, Nb );
        k( 0 , 0 ) = doption(_name="k0");
        k( 1 , 0 ) = doption(_name="k1");
        k( 0 , 1 ) = doption(_name="k2");
        k( 1 , 1 ) = doption(_name="k3");
        std::map<std::string,double> material{ {"Mat0", k(0,0)}, {"Mat1",k(1,0)}, {"Mat2",k(0,1)},  {"Mat3",k(1,1)} };

        auto a = form2( Xh, Xh );
        for ( auto marker : material )
        {
            LOG(INFO) << "Material " << marker.first << " with conductivity " << marker.second;
            a += integrate( markedelements( mesh, marker.first ),
                            marker.second*gradt( u )*trans( grad( v ) ) );
        }

        double Bi;

        if ( N == 1 )
            Bi = 0.1;
        else
            Bi = math::exp( math::log( Bimin )+double( i )*( math::log( Bimax )-math::log( Bimin ) )/double( N-1 ) );

        LOG(INFO) << "Bi = " << Bi << "\n";

        a += integrate( markedfaces( mesh, "Tfourier" ), Bi*idt( u )*id( v ) );


        a.solve( _solution=u, _rhs=f );

        double moy_u = mean( markedfaces( mesh, "Tflux" ), idv( u ) )( 0,0 );
        if (Environment::worldComm().isMasterRank())
        {
            std::cout.precision( 5 );
            std::cout << std::setw( 5 ) << k( 0,0 ) << " "
                      << std::setw( 5 ) << k( 1,0 ) << " "
                      << std::setw( 5 ) << k( 0,1 ) << " "
                      << std::setw( 5 ) << k( 1,1 ) << " "
                      << std::setw( 5 ) << Bi << " "
                      << std::setw( 10 ) << moy_u << "\n";
        }
        e->step( i )->setMesh( mesh );
        e->step( i )->add( "Temperature", u );
        e->save();




    }
}
