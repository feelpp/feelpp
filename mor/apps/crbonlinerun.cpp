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
//! @date 11 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <boost/dll.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/enumerate.hpp>
#include <feel/feelcore/table.hpp>
#include <feel/feelfilters/loadcsv.hpp>
#include <feel/feelmor/mormodels.hpp>
#include <feel/feelmor/options.hpp>

#include <fmt/ranges.h>
#include <iostream>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/builder/basic/array.hpp>
#endif

bool runCrbOnline( std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugin );
/*
 ( "crbmodel.root", po::value<std::string>(), "CRB online code root repository" )
            ( "crbmodel.name", po::value<std::string>(), "CRB online code name" )
            ( "crbmodel.attribute", po::value<std::string>()->default_value( "last_modified" ), "last_created, last_modified, id, name" )
            ( "crbmodel.db.id", po::value<std::string>(), "CRB online code id" )
            ( "crbmodel.db.name", po::value<std::string>(), "CRB online code name" )
            ( "crbmodel.db.last", po::value<std::string>()->default_value( "modified" ), "use created or modified" )
            ( "crbmodel.db.load", po::value<std::string>()->default_value( "rb" ), "load rb, fe or all (fe and rb)" )
            ( "crbmodel.db.root_directory", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )

*/

std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin( std::string const& name, std::string const& id )
{
    using namespace Feel;
    namespace dll=boost::dll;
    std::string dirname = Environment::expand( soption(_name="plugin.dir") );
    std::string pluginname = name;

    std::string pluginlibname = "";
    if ( Environment::vm().count("plugin.libname") )
        pluginlibname = soption(_name="plugin.libname");

    auto plugin = factoryCRBPlugin( pluginname, pluginlibname, dirname );
    std::cout << "Loaded the plugin " << plugin->name() << std::endl;
    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    if(boption(_name="export-solution"))
        loadFiniteElementDatabase=true;
    std::string jsonfilename = (fs::path(Environment::expand( soption(_name="plugin.db") )) / fs::path(pluginname) / fs::path(id) / (pluginname+".crb.json")).string() ;

    std::cout << " . using db " << jsonfilename << std::endl;

    plugin->loadDB( jsonfilename, (loadFiniteElementDatabase)? crb::load::all : crb::load::rb );

    return plugin;
}

void runCrbOnlineList()
{
#if defined ( FEELPP_HAS_MONGOCXX )
    mongocxx::client conn{mongocxx::uri{}};

    using namespace bsoncxx::builder::basic;
    document document{};

    auto collection = conn["feelpp"]["crbdb"];

    std::cout << "document " << bsoncxx::to_json(document.view()) << std::endl;
    auto cursor = collection.find(document.extract());


    for (auto&& doc : cursor) {
        std::string n = doc["crbmodel"]["name"].get_utf8().value.to_string();
        std::string d = doc["crb"]["dimension"].get_utf8().value.to_string();
        std::string o = doc["crb"]["output-index"].get_utf8().value.to_string();
        std::cout << " . " << n << " " << d << " " << o << std::endl;
    }
#else
    std::cout << "Feel++ was not compiled with Mongo C++ support" << std::endl;
#endif
}

void runCrbOnlineQuery()
{
#if defined ( FEELPP_HAS_MONGOCXX )
    mongocxx::client conn{mongocxx::uri{}};

    using namespace bsoncxx::builder::basic;
    document document{};


    auto collection = conn["feelpp"]["crbdb"];

    typedef std::vector< std::string > split_vector_type;

    split_vector_type SplitVec; // #2: Search for tokens
    std::string q = Feel::soption("query");
    boost::split( SplitVec, q, boost::is_any_of(","), boost::token_compress_on );
    array arr{};
    for( auto const& p: SplitVec )
    {
        split_vector_type s; // #2: Search for tokens
        boost::split( s, p, boost::is_any_of(":"), boost::token_compress_on );
        arr.append( [&s](sub_document subdoc) { subdoc.append(kvp(s[0], s[1] ) ); });
    }
    document.append( kvp( "$and", arr ) );
    std::cout << "document " << bsoncxx::to_json(document.view()) << std::endl;
    auto cursor = collection.find(document.extract());


    for (auto&& doc : cursor) {
        //std::cout << bsoncxx::to_json(doc) << std::endl;
        LOG(INFO) << "crbmodel.name: " << doc["crbmodel"]["name"].get_utf8().value.to_string() << std::endl;
        LOG(INFO) << "uuid: " << doc["uuid"].get_utf8().value.to_string() << std::endl;
        runCrbOnline( { loadPlugin( doc["crbmodel"]["name"].get_utf8().value.to_string(),doc["uuid"].get_utf8().value.to_string() ) } );
    }
#else
    std::cout << "Feel++ was not compiled with Mongo C++ support" << std::endl;
#endif
}

void runCrbOnlineCompare()
{
#if defined ( FEELPP_HAS_MONGOCXX )
    mongocxx::client conn{mongocxx::uri{}};

    using namespace bsoncxx::builder::basic;
    document document{};


    auto collection = conn["feelpp"]["crbdb"];

    typedef std::vector< std::string > split_vector_type;

    split_vector_type SplitVec; // #2: Search for tokens
    std::string q = Feel::soption("compare");
    boost::split( SplitVec, q, boost::is_any_of(","), boost::token_compress_on );
    array arr{};
    for( auto const& p: SplitVec )
    {
        split_vector_type s; // #2: Search for tokens
        boost::split( s, p, boost::is_any_of(":"), boost::token_compress_on );
        arr.append( [&s](sub_document subdoc) { subdoc.append(kvp(s[0], s[1] ) ); });
    }
    document.append( kvp( "$and", arr ) );
    std::cout << "document " << bsoncxx::to_json(document.view()) << std::endl;
    auto cursor = collection.find(document.extract());

    std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugins;

    for (auto&& doc : cursor) {
        //std::cout << bsoncxx::to_json(doc) << std::endl;
        LOG(INFO) << "crbmodel.name: " << doc["crbmodel"]["name"].get_utf8().value.to_string() << std::endl;
        LOG(INFO) << "uuid: " << doc["uuid"].get_utf8().value.to_string() << std::endl;
        plugins.push_back( loadPlugin( doc["crbmodel"]["name"].get_utf8().value.to_string(),doc["uuid"].get_utf8().value.to_string() ) );
    }
    runCrbOnline( plugins );
#else
    std::cout << "Feel++ was not compiled with Mongo C++ support" << std::endl;
#endif
}


std::string
loadModelName( std::string const& filename )
{
    using namespace Feel;
    if ( !fs::exists( filename ) )
    {
        LOG(INFO) << "Could not find " << filename << std::endl;
        return std::string("");
    }

    auto json_str_wo_comments = removeComments(readFromFile(filename));
    //LOG(INFO) << "json file without comment:" << json_str_wo_comments;

    boost::property_tree::ptree ptree;
    std::istringstream istr( json_str_wo_comments );
    boost::property_tree::read_json( istr, ptree );

    auto const& ptreeCrbModel = ptree.get_child( "crbmodel" );
    std::string modelName = ptreeCrbModel.template get<std::string>( "name" );
    return modelName;
}

namespace Feel {

auto
setupSampling( MORModels const& m )
{
    auto muspace = m.parameterSpace();

    std::ostringstream ostrmumin, ostrmumax;
    auto mumin = muspace->min();
    auto mumax = muspace->max();
    for ( uint16_type d = 0; d < muspace->dimension(); ++d )
    {
        ostrmumin << mumin( d ) << " ";
        ostrmumax << mumax( d ) << " ";
    }
    std::cout << "dimension of parameter space : " << muspace->dimension() << "\n";
    std::cout << "min element in parameter space : " << ostrmumin.str() << "\n";
    std::cout << "max element in parameter space : " << ostrmumax.str() << "\n";

    auto mysampling = muspace->sampling();

    std::vector<double> inputParameter;
    if ( Environment::vm().count( "parameter" ) )
    {
        auto inputParameterParsed = Environment::vm()["parameter"].as<std::vector<std::string>>();

        if ( inputParameterParsed.size() == 1 )
        {
            std::vector<std::string> stringParsedSplitted;
            boost::split( stringParsedSplitted, inputParameterParsed.front(), boost::is_any_of( " " ), boost::token_compress_on );
            inputParameterParsed = stringParsedSplitted;
        }

        for ( std::string const& paramParsed : inputParameterParsed )
            inputParameter.push_back( std::stod( paramParsed ) );
    }
    else if ( Environment::vm().count( "parameter.filename" ) )
    {
        std::string fname = Environment::vm()["parameter.filename"].as<std::string>();
        auto r = loadXYFromCSV( fname, muspace->parameterNames() );
        auto mu = muspace->element();
        for ( auto const& p : r )
        {
            mu = p;
            mysampling->push_back( mu );
        }
    }
    // inputParameter = Environment::vm()["parameter"].as<std::vector<double> >();
    if ( !inputParameter.empty() )
    {
        CHECK( inputParameter.size() == muspace->dimension() ) << "parameter has a wrong size : " << inputParameter.size() << " but must be " << muspace->dimension() << ":" << inputParameter;
        auto mu = muspace->element();
        for ( uint16_type d = 0; d < muspace->dimension(); ++d )
            mu( d ) = inputParameter[d];
        mysampling->push_back( mu );
    }
    else if ( mysampling->empty() )
    {
        int nSample = ioption( _name = "sampling.size" );
        std::string sampler = soption( "sampling.type" );
        mysampling->sample( nSample, sampler );
    }
    return mysampling;
}
}
int main(int argc, char**argv )
{
    using namespace Feel;
    try {
        po::options_description crbonlinerunoptions( "crb online run options" );
        crbonlinerunoptions.add_options()
            ( "plugin.dir", po::value<std::string>()->default_value( Info::libdir().string() ) , "plugin directory" )
            ( "morjson", po::value<std::string>(), "filename describing the mor models associated to the case study" )
# if 0
            ( "crbmodel.root", po::value<std::string>(), "CRB online code root repository" )
            ( "crbmodel.name", po::value<std::string>(), "CRB online code name" )
            ( "crbmodel.attribute", po::value<std::string>()->default_value( "last_modified" ), "last_created, last_modified, id, name" )
            ( "crbmodel.db.id", po::value<std::string>(), "CRB online code id" )
            ( "crbmodel.db.name", po::value<std::string>(), "CRB online code name" )
            ( "crbmodel.db.last", po::value<std::string>()->default_value( "modified" ), "use created or modified" )
            ( "crbmodel.db.load", po::value<std::string>()->default_value( "rb" ), "load rb, fe or all (fe and rb)" )
            ( "crbmodel.db.root_directory", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )
#endif
            //( "plugin.name", po::value<std::string>(), "CRB online code name" )
            //( "plugin.libname", po::value<std::string>(), "CRB online libname" )
            //( "plugin.dbid", po::value<std::string>(), "CRB online code id" )
            //( "plugin.last", po::value<int>()->default_value( 2 ), "use last created(=1) or modified(=2) or not (=0)" )
            //( "plugin.db", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )
            //( "parameter", po::value<std::vector<double> >()->multitoken(), "database filename" )
            ( "parameter", po::value<std::vector<std::string> >()->multitoken(), "database filename" )
            ( "parameter.filename", po::value<std::string>(), "parameters from csv file" )
            ( "sampling.size", po::value<int>()->default_value( 10 ), "size of sampling" )
            ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
            ( "rb-dim", po::value<int>()->default_value( -1 ), "reduced basis dimension used (-1 use the max dim)" )
            ( "output_results.save.path", po::value<std::string>(), "output_results.save.path" )
            ( "output_results.precision", po::value<int>()->default_value( 6 ), "float precision for output results")
            ( "output_results.print", po::value<bool>()->default_value( true ), "print results in shell")

            ( "query", po::value<std::string>(), "query string for mongodb DB feelpp.crbdb" )
            ( "compare", po::value<std::string>(), "compare results from query in mongodb DB feelpp.crbdb" )
            ( "list", "list registered DB in mongoDB  in feelpp.crbdb" )
            ( "export-solution", po::value<bool>()->default_value(false), "export the solutions for visualization")
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
        Environment env(_argc=argc, _argv=argv,
                        _desc=crbonlinerunoptions,
                        _desc_lib=crbonlinerunliboptions.add( feel_options() ),
                        _about=about(_name="crbonlinerun",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        std::ifstream is( Environment::expand( soption(_name="morjson" ) ) );
        nl::json js = nl::json::parse(is);
        LOG(INFO) << "js = " << js.dump(4) << std::endl;

        MORModels ms{ js };
        ms.load();
        MORTable table( ms );
        ms.addObserver( std::make_shared<MORTable>( table ) );
        ms.run(setupSampling( ms ), nl::json{ { "N", ioption(_name="rb-dim") }, { "print_rb_matrix", true }, { "tolerance", 1e-2 } } );
#if 0
        if ( Environment::vm().count( "list" )  )
        {
            runCrbOnlineList();
            return 0;
        }
        if ( Environment::vm().count( "compare" )  )
        {
            runCrbOnlineCompare();
            return 0;
        }
        if ( Environment::vm().count( "query" ) == 0 )
        {
            std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugins;
            for( auto const& j : jmodels )
            {
                plugins.push_back( loadPlugin( j ) );
            }

            runCrbOnline( plugins );
        }
        else
        {
            runCrbOnlineQuery();
        }
#endif        
    }
    catch( ... )
    {
        handleExceptions();
    }
    return 0;
}
