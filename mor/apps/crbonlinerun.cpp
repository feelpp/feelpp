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

#include <feel/feelcore/table.hpp>
#include <feel/feelfilters/loadcsv.hpp>
#include <feel/feelmor/crbmodeldb.hpp>
#include <feel/feelmor/crbplugin_interface.hpp>
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
bool
runCrbOnline( std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugin )
{
    using namespace Feel;

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    if(boption(_name="export-solution"))
        loadFiniteElementDatabase=true;

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;//Feel::doption(Feel::_name="crb.online-tolerance");
    bool print_rb_matrix = false;//boption(_name="crb.print-rb-matrix");
    auto muspace = plugin[0]->parameterSpace();

    std::ostringstream ostrmumin,ostrmumax;
    auto mumin=muspace->min();
    auto mumax=muspace->max();
    for ( uint16_type d = 0; d < muspace->dimension(); ++d)
    {
        ostrmumin << mumin(d) << " ";
        ostrmumax << mumax(d) << " ";
    }
    std::cout << "dimension of parameter space : " << muspace->dimension() << "\n";
    std::cout << "min element in parameter space : "<< ostrmumin.str() << "\n";
    std::cout << "max element in parameter space : "<< ostrmumax.str() << "\n";


    auto mysampling = muspace->sampling();

    std::vector<double> inputParameter;
    if ( Environment::vm().count("parameter"))
    {
        auto inputParameterParsed = Environment::vm()["parameter"].as<std::vector<std::string> >();

        if ( inputParameterParsed.size() == 1 )
        {
            std::vector<std::string > stringParsedSplitted;
            boost::split( stringParsedSplitted, inputParameterParsed.front(), boost::is_any_of(" "), boost::token_compress_on );
            inputParameterParsed = stringParsedSplitted;
        }

        for ( std::string const& paramParsed : inputParameterParsed )
            inputParameter.push_back( std::stod(paramParsed) );
    }
    else if( Environment::vm().count( "parameter.filename" ) )
    {
        std::string fname = Environment::vm()["parameter.filename"].as<std::string>();
        auto r = loadXYFromCSV( fname, muspace->parameterNames() );
        auto mu = muspace->element();
        for(auto const& p : r)
        {
            mu = p;
            mysampling->push_back( mu );
        }
    }
    //inputParameter = Environment::vm()["parameter"].as<std::vector<double> >();
    if ( !inputParameter.empty() )
    {
        CHECK( inputParameter.size() == muspace->dimension() ) << "parameter has a wrong size : "<< inputParameter.size() << " but must be " << muspace->dimension() << ":"<<inputParameter;
        auto mu = muspace->element();
        for ( uint16_type d = 0; d < muspace->dimension(); ++d)
            mu(d) = inputParameter[d];
        mysampling->push_back( mu );
    }
    else if ( mysampling->empty() )
    {
        int nSample = ioption(_name="sampling.size");
        std::string sampler = soption("sampling.type");
        mysampling->sample( nSample, sampler );
    }

    if ( loadFiniteElementDatabase )
        plugin[0]->initExporter();

    int rbDim = ioption(_name="rb-dim");
    int nSamples = mysampling->size();

    Feel::Table tableOutputResults;
    tableOutputResults.format().setFloatingPointPrecision( ioption(_name="output_results.precision") );
    std::vector<std::string> tableRowHeader = muspace->parameterNames();
    tableRowHeader.push_back( "output" );
    tableRowHeader.push_back( "errorBound" );
    tableRowHeader.push_back( "time(s)" );
    tableOutputResults.add_row( tableRowHeader );
    tableOutputResults.format().setFirstRowIsHeader( true );

    std::vector<double> tableRowValues(tableRowHeader.size());

    for ( int k = 0; k < nSamples; ++k )
    {
        
        auto const& mu = (*mysampling)[k];
        std::ostringstream ostrmu;
        for ( uint16_type d = 0; d < muspace->dimension(); ++d)
            ostrmu << mu(d) << " ";
        // std::cout << "--------------------------------------\n";
        // std::cout << "mu["<<k<<"] : " << ostrmu.str() << "\n";
        //auto mu = crb->Dmu()->element();
        //std::cout << "input mu\n" << mu << "\n";
        for( auto const& p : plugin )
        {
            tic();
            auto crbResult = p->run( mu, time_crb, online_tol, rbDim, print_rb_matrix);
            double t = toc("rb_online", FLAGS_v>0);

            double resOuptut = crbResult.output();
            double resErrorBound = crbResult.errorbound();
            int curRowValIndex = 0;
            // add values of parameter mu to table
            for ( uint16_type d = 0; d < muspace->dimension(); ++d)
                tableRowValues[curRowValIndex++] = mu(d);

            // add values of output and error bound to table
            tableRowValues[curRowValIndex++] = resOuptut;
            tableRowValues[curRowValIndex++] = resErrorBound;
            tableRowValues[curRowValIndex++] = t;
            tableOutputResults.add_row( tableRowValues );

            if ( loadFiniteElementDatabase )
            {
                p->exportField( (boost::format("sol-%1%")%k).str(), crbResult );
            }
        }
    }

    bool printResults = true;
    if ( printResults )
        std::cout << tableOutputResults << std::endl;
    bool saveResults = true;

    std::string outputResultPath = "output.csv";
    if ( Environment::vm().count("output_results.save.path") )
        outputResultPath = soption(_name="output_results.save.path");
    outputResultPath = Environment::expand( outputResultPath );
    if ( !fs::exists( fs::path(outputResultPath).parent_path() ) && !fs::path(outputResultPath).parent_path().empty() )
        fs::create_directories( fs::path(outputResultPath).parent_path() );

    if ( saveResults )
    {
        std::ofstream ofs( outputResultPath );
        tableOutputResults.exportCSV( ofs );
    }
    if ( loadFiniteElementDatabase )
        plugin[0]->saveExporter();

    return true;

}

std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin()
{
    using namespace Feel;

    std::string crbmodelName = Environment::expand( soption(_name = "crbmodel.name") );
    CRBModelDB crbmodelDB{ crbmodelName, uuids::nil_uuid() };

    std::string attribute = soption(_name = "crbmodel.attribute" );
    std::string attribute_data;
    if ( attribute == "id"  || attribute == "name")
    {
        attribute_data = Environment::expand( soption(_name = fmt::format("crbmodel.db.{}",attribute) ) );
    }
    else if ( attribute == "last_created" || attribute == "last_modified" )
    {
        std::vector<std::string> split_;
        boost::split(split_, attribute, boost::is_any_of("_"));
        attribute_data = split_[1];
    }
    else
    {
        throw std::runtime_error( "no crbmodel selection, crbmodel.db.id or crbmodel.db.last should be defined" );
    }
    auto meta = crbmodelDB.loadDBMetaData( attribute, attribute_data );
    std::string pluginlibdir = Environment::expand( soption(_name = "plugin.dir") );
    std::cout << "-- crbmodelDB::dbRepository() = " << crbmodelDB.dbRepository() << std::endl;
    std::cout << "-- plugin libdir = " << pluginlibdir << std::endl;

    if(boption(_name="export-solution"))
        return crbmodelDB.loadDBPlugin( meta, "all", pluginlibdir );
    else
        return crbmodelDB.loadDBPlugin( meta, soption(_name="crbmodel.db.load" ), pluginlibdir );
}

int main(int argc, char**argv )
{
    using namespace Feel;
    try {
        po::options_description crbonlinerunoptions( "crb online run options" );
        crbonlinerunoptions.add_options()
            ( "plugin.dir", po::value<std::string>()->default_value( Info::libdir().string() ) , "plugin directory" )

            ( "crbmodel.name", po::value<std::string>(), "CRB online code name" )
            ( "crbmodel.attribute", po::value<std::string>()->default_value( "last_modified" ), "last_created, last_modified, id, name" )
            ( "crbmodel.db.id", po::value<std::string>(), "CRB online code id" )
            ( "crbmodel.db.name", po::value<std::string>(), "CRB online code name" )
            ( "crbmodel.db.last", po::value<std::string>()->default_value( "modified" ), "use created or modified" )
            ( "crbmodel.db.load", po::value<std::string>()->default_value( "rb" ), "load rb, fe or all (fe and rb)" )
            ( "crbmodel.db.root_directory", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )

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
        Environment env( _argc=argc, _argv=argv,
                        _desc=crbonlinerunoptions,
                        _desc_lib=crbonlinerunliboptions.add( feel_options() ),
                        _about=about(_name="crbonlinerun",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

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
            runCrbOnline( { loadPlugin() } );
        }
        else
        {
            runCrbOnlineQuery();
        }
    }
    catch( ... )
    {
        handleExceptions();
    }
    return 0;
}
