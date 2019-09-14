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
//! @author Vincent Chabannes <vincent.chabannes@cemosis.fr>
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 02 Sep 2019
//! @copyright 2019 Feel++ Consortium
//!

#include <feel/feelcrb/options.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>

#include <feel/feelcrb/crbonlinerun.hpp>


int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description crbonlinerunoptions( "crb online run options" );
	crbonlinerunoptions.add_options()
        ( "plugin.dir", po::value<std::string>()->default_value(Info::libdir()) , "plugin directory" )
        ( "plugin.name", po::value<std::string>(), "CRB online code name" )
        ( "plugin.libname", po::value<std::string>(), "CRB online libname" )
        ( "plugin.dbid", po::value<std::string>(), "CRB online code id" )
        ( "plugin.last", po::value<int>()->default_value( 2 ), "use last created(=1) or modified(=2) or not (=0)" )
        ( "plugin.db", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )
        ( "parameter", po::value<std::vector<double> >()->multitoken(), "database filename" )
        ( "sampling.size", po::value<int>()->default_value( 10 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
        ( "rb-dim", po::value<int>()->default_value( -1 ), "reduced basis dimension used (-1 use the max dim)" )
        ( "query", po::value<std::string>(), "query string for mongodb DB feelpp.crbdb" )
        ( "compare", po::value<std::string>(), "compare results from query in mongodb DB feelpp.crbdb" )
        ( "list", "list registered DB in mongoDB  in feelpp.crbdb" )
	 	;
	po::options_description crbonlinerunliboptions( "crb online run lib options" );
#if 1
    crbonlinerunliboptions.add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"))
        ;
#endif
	Environment env( _argc=argc, _argv=argv,
                     _desc=crbonlinerunoptions,
                     _desc_lib=crbonlinerunliboptions.add( feel_options() ),
                     _about=about(_name="crbonlinerun",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string dirname = Environment::expand( soption(_name="plugin.dir") );
    std::string pluginlibname = "";
    if ( Environment::vm().count("plugin.libname") )
        pluginlibname = soption(_name="plugin.libname");
    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    int loadLast = ioption(_name="plugin.last");
                         
    if ( Environment::vm().count( "list" )  )
    {
        runCrbOnlineList();
        return 0;
    }
    if ( Environment::vm().count( "compare" )  )
    {
        runCrbOnlineCompare( dirname, pluginlibname, loadFiniteElementDatabase, loadLast  );
        return 0;
    }
    if ( Environment::vm().count( "query" ) == 0 )
        runCrbOnline( { loadPlugin( "thermalbuilding" ) } );
    else
    {
        runCrbOnlineQuery(  dirname, pluginlibname, loadFiniteElementDatabase, loadLast );
    }

    return 0;
}
