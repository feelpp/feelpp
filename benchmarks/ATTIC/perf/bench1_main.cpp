/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2010 Universite Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2010-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file bench1_main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-06-14
 */
#include <bench1.hpp>

namespace Feel
{
AboutData
makeAbout()
{
    AboutData about( "bench1" ,
                     "bench1" ,
                     "0.2",
                     "assembly performance",
                     AboutData::License_LGPL,
                     "Copyright (c) 2005,2006 EPFL"
                     "Copyright (c) 2006-2010 UniversitÃÂ© Joseph Fourier (Grenoble 1)"
                   );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

po::options_description
makeOptions()
{
    po::options_description desc( "Specific options" );
    desc.add_options()
    ( "dim", po::value<int>()->default_value( 1 ), "dimension (1,2,3)" )
    ( "hsize", po::value<double>()->default_value( 0.1 ), "element size" )
    ( "shape", po::value<std::string>()->default_value( "simplex" ), "type of domain shape: simplex, hypercube ellipsoid" )
    ;
    return desc.add( feel_options() );
}


}
int
main( int argc, char** argv )
{
    Feel::Bench1 bench1( argc, argv, Feel::makeAbout(), Feel::makeOptions() );

#ifdef FEELPP_HAS_TBB
    int n = tbb::task_scheduler_init::default_num_threads();
#else
    int n = 1 ;
#endif

    for ( int p=1; p<=n; ++p )
    {
        Feel::fs::path cur = Feel::fs::current_path();
        std::ostringstream os;
        os << "thread_" << p;
        Feel::fs::path pth = cur / os.str();
        Feel::fs::create_directory( pth );
        ::chdir( pth.string().c_str() );
        std::cout << "benchmark starts with nthreads = " << p << "\n";
#ifdef FEELPP_HAS_TBB
        tbb::task_scheduler_init init( p );

        tbb::tick_count t0 = tbb::tick_count::now();
#endif
        bench1.run();
#ifdef FEELPP_HAS_TBB
        tbb::tick_count t1 = tbb::tick_count::now();
        double t = ( t1-t0 ).seconds();
        ::chdir( cur.string().c_str() );
        std::cout << "benchmark stops with nthreads = " << p << " in " << t << " seconds\n";
#endif
    }
}

