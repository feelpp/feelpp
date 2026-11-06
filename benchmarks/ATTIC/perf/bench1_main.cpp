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

    bench1.run();
}
