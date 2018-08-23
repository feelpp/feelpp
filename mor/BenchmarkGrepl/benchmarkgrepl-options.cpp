/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include "benchmarkgrepl-options.hpp"

namespace Feel
{
po::options_description
makeBenchmarkGreplLinearEllipticOptions()
{
    po::options_description bgoptions( "BenchmarkGreplLinearElliptic options" );
    bgoptions.add_options()
        ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
        ( "hsize", Feel::po::value<double>()->default_value( 1e-1 ), "hsize")
        ( "trainset-eim-size", Feel::po::value<int>()->default_value( 40 ), "EIM trainset is built using a equidistributed grid 40 * 40 by default")
        ( "gamma", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
        ;
    return bgoptions;
}

AboutData
makeBenchmarkGreplLinearEllipticAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Benchmark Grepl",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2011-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}

po::options_description
makeBenchmarkGreplNonlinearEllipticOptions()
{
    po::options_description bgoptions( "BenchmarkGreplNonlinearElliptic options" );
    bgoptions.add_options()
        ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
        ( "hsize", Feel::po::value<double>()->default_value( 1e-1 ), "hsize")
        ( "trainset-eim-size", Feel::po::value<int>()->default_value( 15 ), "EIM trainset is built using a equidistributed grid 15*15 by default")
        ( "gamma", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
        ( "use-deim", Feel::po::value<bool>()->default_value( false ), "use deim or eim if false (default)" )
        ;
    return bgoptions;
}
AboutData
makeBenchmarkGreplNonlinearEllipticAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Benchmark Grepl",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2011-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    return about;
}

}
