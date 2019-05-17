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



#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <pybind11/stl.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelio.hpp>

namespace Feel {

void
//pyexprFromFile( std::string const& pyfilename, std::vector<std::string> const& vars, std::map<std::string,std::string> const& _locals )
pyexprFromFile( std::string const& pyfilename, std::map<std::string,std::map<std::string,std::string>> & _locals )
{
    using Feel::cout;


    std::map<std::string,std::string> r;
    
    if ( Environment::isMasterRank() )
    {
        py::scoped_interpreter guard{}; // start the interpreter and keep it alive
        py::module::import("sys").attr("path").cast<py::list>().append(Environment::expand("$top_srcdir/feelpp/feel/feelpython/"));
        py::print(py::module::import("sys").attr("path"));
        py::dict locals = py::cast(_locals);
        py::print(locals);
        py::eval_file( pyfilename.c_str(), py::globals(), locals );
        
#if 1
        for( auto l : _locals )
        {
            std::cout << "l: " << l.first << std::endl;
            for( auto n : l.second )
            {
                std::string v = l.first+"['"+n.first+"']";
                std::string cmd =  v + "= sympytoginac(" + v +" );";
                py::exec(cmd, py::globals(), locals );
                _locals[l.first][n.first] = locals[l.first.c_str()][n.first.c_str()].cast<std::string>();
            }
        }
#endif
    }
    mpi::broadcast( Environment::worldComm(), _locals, 0 );
    
    //return r;
}

//std::map<std::string,std::string> 
void
//pyexprFromFile( std::string const& pyfilename, std::vector<std::string> const& vars, std::map<std::string,std::string> const& _locals )
pyexprFromFile( std::string const& pyfilename, std::map<std::string,std::string> & _locals )
{
    using Feel::cout;


    std::map<std::string,std::string> r;
    
    if ( Environment::isMasterRank() )
    {
        py::scoped_interpreter guard{}; // start the interpreter and keep it alive
        py::module::import("sys").attr("path").cast<py::list>().append(Environment::expand("$top_srcdir/feelpp/feel/feelpython/"));
        py::print(py::module::import("sys").attr("path"));
        py::dict locals = py::cast(_locals);
        py::print(locals);
        py::eval_file( pyfilename.c_str(), py::globals(), locals );

#if 0
        for( auto const& [key,value] : _locals )
        {
            std::cout << key << "," << value << std::endl;
        }
#endif
        for( auto const& [key,value] : _locals )
        {
            if ( key.empty() ) continue;
            std::string cmd = key + "= sympytoginac(" + key +" );";
            //std::cout << "cmd:" << cmd << std::endl;
            py::exec(cmd, py::globals(), locals );
            _locals[key] = locals[key.c_str()].cast<std::string>();
        }
    }
    mpi::broadcast( Environment::worldComm(), _locals, 0 );
    
    //return r;
}

std::map<std::string,std::string> 
pyexpr( std::string const& pycode, std::vector<std::string> const& vars, std::map<std::string,std::string> const& locals )
{
    using Feel::cout;


    std::map<std::string,std::string> r;
    
    if ( Environment::isMasterRank() )
    {
        py::scoped_interpreter guard{}; // start the interpreter and keep it alive
        py::module::import("sys").attr("path").cast<py::list>().append(Environment::expand("$top_srcdir/feelpp/feel/feelpython/"));
        py::dict _locals=py::cast(locals);
        py::exec(pycode.c_str(), py::globals(),_locals);

        for( auto l : vars )
        {
            std::string cmd = l + "= toginac(sympify(" + l + "), [x] if len(" + l + ".free_symbols)==0 else " + l + ".free_symbols );";
            py::exec(cmd, py::globals(), _locals );
            r[l] = _locals[l.c_str()].cast<std::string>();
        }
    }
    mpi::broadcast( Environment::worldComm(), r, 0 );
    
    return r;
}




}
