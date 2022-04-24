/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_BENCHMARKGREPL_OPTIONS_HPP
#define FEELPP_BENCHMARKGREPL_OPTIONS_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>

namespace Feel
{

po::options_description
makeBenchmarkGreplLinearEllipticOptions();
AboutData
makeBenchmarkGreplLinearEllipticAbout( std::string const& str = "benchmarkGrepl" );

po::options_description
makeBenchmarkGreplNonlinearEllipticOptions();
AboutData
makeBenchmarkGreplNonlinearEllipticAbout( std::string const& str = "benchmarkGrepl" );

}

#endif
