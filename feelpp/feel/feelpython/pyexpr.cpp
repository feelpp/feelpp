//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 18 Sep 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <iostream>
#include <fmt/core.h>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelio.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <pybind11/stl.h>

namespace Feel
{

std::vector<std::string> lookups_ =
    {
        "$cfgdir/../../python/",
        "$cfgdir/../python/",
        "$datadir/testcases/python/",
        "$datadir/python/",
        "$top_srcdir/feelpp/quickstart/python/",
        "$top_srcdir/feelpp/feel/feelpython/" };

void
// pyexprFromFile( std::string const& pyfilename, std::vector<std::string> const& vars, std::map<std::string,std::string> const& _locals )
pyexprFromFile( std::string const& pyfilename, std::map<std::string, std::map<std::string, std::string>>& _locals )
{
    using Feel::cout;

    std::map<std::string, std::string> r;

    if ( Environment::isMasterRank() )
    {
        try
        {
            std::string sympytoginac_f = Environment::findFile( "sympy2ginac.py", lookups_ );
            // std::cout << "sympytoginac_f = " << sympytoginac_f << std::endl;
            py::module::import( "sys" ).attr( "path" ).cast<py::list>().append( fs::path( sympytoginac_f ).parent_path().string() );
            // py::print(py::module::import("sys").attr("path"));
            py::dict locals = py::cast( _locals );
            // py::print(locals);
            // std::cout << "eval_f = " << Environment::findFile( pyfilename.c_str(), lookups_ ) << std::endl;
            py::eval_file( Environment::findFile( pyfilename.c_str(), lookups_ ), py::globals(), locals );

#if 1
            for ( auto l : _locals )
            {
                std::cout << "l: " << l.first << std::endl;
                for ( auto n : l.second )
                {
                    std::string v = l.first + "['" + n.first + "']";
                    std::string cmd = v + "= sympytoginac(" + v + " );";
                    py::exec( cmd, py::globals(), locals );
                    _locals[l.first][n.first] = locals[l.first.c_str()][n.first.c_str()].cast<std::string>();
                }
            }
#endif
        }
        catch ( const pybind11::error_already_set& ex )
        {
            std::cerr << fmt::format( "[feelpp.pybind11.error_already_set] python interpreter failed : {}", ex.what() ) << std::endl;
            throw;
        }
    }
    mpi::broadcast( Environment::worldComm(), _locals, 0 );

    // return r;
}

FEELPP_NO_EXPORT std::map<std::string, std::string> clean_locals( std::map<std::string, std::string>& locals )
{
    using namespace boost::algorithm;
    std::map<std::string, std::string> ls;
    for ( auto const& [key, value] : locals )
    {
        std::vector<std::string> v;
        split( v, value, is_any_of( ":" ) );
        ls[key] = v[0];
    }
    return ls;
}
// std::map<std::string,std::string>
void
// pyexprFromFile( std::string const& pyfilename, std::vector<std::string> const& vars, std::map<std::string,std::string> const& _locals )
pyexprFromFile( std::string const& pyfilename, std::map<std::string, std::string>& _locals )
{
    using Feel::cout;

    std::map<std::string, std::string> r;

    if ( Environment::isMasterRank() )
    {
        try
        {
            std::string sympytoginac_f = Environment::findFile( "sympy2ginac.py", lookups_ );
            // std::cout << "sympytoginac_f = " << sympytoginac_f << std::endl;
            py::module::import( "sys" ).attr( "path" ).cast<py::list>().append( fs::path( sympytoginac_f ).parent_path().string() );
            // py::print(py::module::import("sys").attr("path"));
            py::dict locals = py::cast( clean_locals( _locals ) );
            // py::print(locals);
            // std::cout << "eval_f = " << Environment::findFile( pyfilename.c_str(), lookups_ ) << std::endl;
            py::eval_file( Environment::findFile( pyfilename.c_str(), lookups_ ), py::globals(), locals );

#if 0
        for( auto const& [key,value] : _locals )
        {
            std::cout << key << "," << value << std::endl;
        }
#endif
            for ( auto const& [key, value] : _locals )
            {
                if ( key.empty() ) continue;
                std::string cmd = key + "= sympytoginac(" + key + " );";
                // std::cout << "cmd:" << cmd << std::endl;
                py::exec( cmd, py::globals(), locals );
                _locals[key] = locals[key.c_str()].cast<std::string>();
            }
        }
        catch ( const pybind11::error_already_set& ex )
        {
            std::cerr << fmt::format( "[feelpp.pybind11.error_already_set] python interpreter failed : {}", ex.what() ) << std::endl;
            throw;
        }
    }
    mpi::broadcast( Environment::worldComm(), _locals, 0 );

    // return r;
}

std::map<std::string, std::string>
pyexpr( std::string const& pycode, std::vector<std::string> const& vars, std::map<std::string, std::string> const& locals )
{
    using Feel::cout;

    std::map<std::string, std::string> r;

    if ( Environment::isMasterRank() )
    {
        try
        {
            py::module::import( "sys" ).attr( "path" ).cast<py::list>().append( Environment::expand( "$top_srcdir/feelpp/feel/feelpython/" ) );
            py::dict _locals = py::cast( locals );
            py::exec( pycode.c_str(), py::globals(), _locals );

            for ( std::string const& l : vars )
            {
                // std::string cmd = l + "= toginac(sympify(" + l + "), [x] if len(" + l + ".free_symbols)==0 else " + l + ".free_symbols );";
                std::string cmd = l + "=sympytoginac( " + l + " );";
                py::exec( cmd, py::globals(), _locals );
                r[l] = _locals[l.c_str()].cast<std::string>();
            }
        }
        catch ( const pybind11::error_already_set& ex )
        {
            std::cerr << fmt::format( "[feelpp.pybind11.error_already_set] python interpreter failed : {}", ex.what() ) << std::endl;
            throw;
        }
    }
    mpi::broadcast( Environment::worldComm(), r, 0 );

    return r;
}

} // namespace Feel
