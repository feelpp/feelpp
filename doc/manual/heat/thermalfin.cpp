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
#include <feel/feel.hpp>

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
    return thermalfinoptions.add( Feel::feel_options() );
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


namespace Feel
{
/**
 * Thermal fin application
 *
 */
class ThermalFin : public Application
{
public:

    /**
     * run the convergence test
     */
    void run();
}; // ThermalFin


void
ThermalFin::run()
{
    Environment::changeRepository( boost::format( "%1%/%2%/" )
                                   % Environment::about().appName()
                                   % option( _name="gmsh.hsize" ).as<double>() );

    /*
     * First we create the mesh
     */
    auto mesh = loadMesh( _mesh=new mesh_type );

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = Pch<2>( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    // flux a bottom of fin
    auto f = form1( Xh );
    f = integrate( markedfaces( mesh,1 ), id( v ) );

    /*
     * Construction of the left hand side
     */
    Dcst = M_backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, Dcst ) = integrate( markedelements( mesh,1 ), ( gradt( u )*trans( grad( v ) ) ) );

    Dcst->close();

    double Bimin = this->vm()["Bimin"].as<double>();
    double Bimax = this->vm()["Bimax"].as<double>();
    int N = this->vm()["N"].as<int>();

    for ( int i = 0; i < N; ++i )
    {

        int Nb = 2;
        Eigen::MatrixXd k( Nb, Nb );
        k( 0 , 0 ) = this->vm()["k0"].as<double>();
        k( 1 , 0 ) = this->vm()["k1"].as<double>();
        k( 0 , 1 ) = this->vm()["k2"].as<double>();
        k( 1 , 1 ) = this->vm()["k3"].as<double>();

        auto a = form2( Xh, Xh );
        for ( int r = 0; r < Nb; ++r )
            for ( int c = 0; c < Nb; ++c )
            {

                a += integrate( markedelements( mesh,Nb*c+r+1 ),
                                k( r,c )*gradt( u )*trans( grad( v ) ) );
            }

        double Bi;

        if ( N == 1 )
            Bi = 0.1;

        else
            Bi = math::exp( math::log( Bimin )+value_type( i )*( math::log( Bimax )-math::log( Bimin ) )/value_type( N-1 ) );

        LOG(INFO) << "Bi = " << Bi << "\n";

        a += integrate( markedfaces( mesh,2 ), Bi*idt( u )*id( v ) );


        a.solve( _solution=u, _rhs=Fcst );

        double moy_u = ( integrate( markedfaces( mesh,1 ), idv( u ) ).evaluate()( 0,0 ) /
                         integrate( markedfaces( mesh,1 ), constant( 1.0 ) ).evaluate()( 0,0 ) );
        std::cout.precision( 5 );
        std::cout << std::setw( 5 ) << k( 0,0 ) << " "
                  << std::setw( 5 ) << k( 1,0 ) << " "
                  << std::setw( 5 ) << k( 0,1 ) << " "
                  << std::setw( 5 ) << k( 1,1 ) << " "
                  << std::setw( 5 ) << Bi << " "
                  << std::setw( 10 ) << moy_u << "\n";

        e->step( i )->setMesh( U.functionSpace()->mesh() );
        e->step( i )->add( "Temperature", U );
        e->save();




    }
} // ThermalFin::run

} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    /* define and run application */
    ThermalFin thermalfin;

    thermalfin.run();
}





