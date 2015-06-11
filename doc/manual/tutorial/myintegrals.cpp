/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
/// [all]
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelfilters/exporter.hpp>



using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description options( "DAR options" );
    options.add_options()
        ( "nbcores", Feel::po::value<int>()->default_value( 1 ), "Number of cores" )
        ;
    return options;
}

int
main( int argc, char** argv )
{
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about( _name="myintegrals" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    /// [mesh]
    // create the mesh (specify the dimension of geometric entity)
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>> );
    /// [mesh]

    /// [expression]
    // our function to integrate
#if 0
    auto g = expr( soption(_name="functions.g") );
#endif
    /// [expression]

#if defined(FEELPP_HAS_OPENMP)

    omp_set_num_threads(ioption(_name="nbcores"));
    std::cout << "Using " << ioption(_name="nbcores") << " threads" << std::endl;
    std::cout << boption( _name="parallel.cpu.enable" ) << " " << soption( _name="parallel.cpu.impl") << " " << ioption( _name="parallel.cpu.restrict") << std::endl;

    struct timespec ts1;
    struct timespec ts2;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

#pragma omp parallel
{
    int id = omp_get_thread_num();

    struct timespec tts1;
    struct timespec tts2;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tts1);

    /// [integrals]
    // compute integral of g (global contribution): \(\int_{\partial \Omega} g\)
    auto intf_1 = integrate( _range = elements( mesh ),
                                 _expr = cst(1.0) ).evaluate();

    clock_gettime(CLOCK_MONOTONIC_RAW, &tts2);
    double tt1 = (double)(tts1.tv_sec) + (double)(tts1.tv_nsec) / (1000000000.0); 
    double tt2 = (double)(tts2.tv_sec) + (double)(tts2.tv_nsec) / (1000000000.0); 
    std::cout << Environment::worldComm().globalRank() << "|" << id << " elapsed (MONOTONIC_RAW)=" <<  (tt2 - tt1) << std::endl;

    std::cout << "intf_1" << " = " << intf_1  << std::endl;
}

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    double t1 = (double)(ts1.tv_sec) + (double)(ts1.tv_nsec) / (1000000000.0); 
    double t2 = (double)(ts2.tv_sec) + (double)(ts2.tv_nsec) / (1000000000.0); 
    std::cout << Environment::worldComm().globalRank() << " elapsed (MONOTONIC_RAW)=" <<  (t2 - t1) << std::endl;

#else

    std::cout << "Harts is not configured" << std::endl;

#endif


#if 0
    // compute integral g on boundary: \( \int_{\partial \Omega} g \)
    auto intf_2 = integrate( _range = boundaryfaces( mesh ),
                             _expr = g ).evaluate();

    // compute integral of grad f (global contribution): \( \int_{\Omega} \nabla g \)
    auto grad_g = grad<2>(g);
    auto intgrad_f = integrate( _range = elements( mesh ),
                                _expr = grad_g ).evaluate();
#endif

    // only the process with rank 0 prints to the screen to avoid clutter
    /*
    if ( Environment::isMasterRank() )
    {
        std::cout << "intf_1" << " = " << intf_1  << std::endl;
    }
    */
#if 0
    if ( Environment::isMasterRank() )
        std::cout << "int_Omega " << g << " = " << intf_1  << std::endl
                  << "int_{boundary of Omega} " << g << " = " << intf_2 << std::endl
                  << "int_Omega grad " << g << " = "
                  << "int_Omega  " << grad_g << " = "
                  << intgrad_f  << std::endl;
#endif
    /// [integrals]
}
/// [all]
