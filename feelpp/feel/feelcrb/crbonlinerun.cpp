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

#include <feel/feelcrb/options.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>

#include <iostream>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/builder/basic/array.hpp>
#endif

namespace Feel {

bool runCrbOnline( std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugin );

std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin( std::string const& name, std::string const& id, std::string const& dirname,
            std::string const& libname,
            bool loadFiniteElementDatabase,
            int loadLast )
{
    using namespace Feel;
    namespace dll=boost::dll;
    std::string pluginname = name;


    auto plugin = factoryCRBPlugin( pluginname, libname, dirname );
    std::cout << "Loaded the plugin " << plugin->name() << std::endl;

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

void runCrbOnlineQuery(  std::string const& dirname, std::string const& libname, bool loadFiniteElementDatabase, int loadLast  )
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
        runCrbOnline( { loadPlugin( doc["crbmodel"]["name"].get_utf8().value.to_string(),doc["uuid"].get_utf8().value.to_string(),
                                    dirname, libname, loadFiniteElementDatabase, loadLast ) } );
    }
#else
    std::cout << "Feel++ was not compiled with Mongo C++ support" << std::endl;
#endif
}

void runCrbOnlineCompare( std::string const& dirname, std::string const& libname, bool loadFiniteElementDatabase, int loadLast )
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
        plugins.push_back( loadPlugin( doc["crbmodel"]["name"].get_utf8().value.to_string(),
                                       doc["uuid"].get_utf8().value.to_string(),
                                       dirname, libname, loadFiniteElementDatabase, loadLast ) );

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

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;//Feel::doption(Feel::_name="crb.online-tolerance");
    bool print_rb_matrix = false;//boption(_name="crb.print-rb-matrix");
    auto muspace = plugin[0]->parameterSpace();

    std::ostringstream ostrmumin,ostrmumax;
    auto mumin=muspace->min();
    auto mumax=muspace->max();
    for ( uint16_type d=0;d<muspace->dimension();++d)
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
        inputParameter = Environment::vm()["parameter"].as<std::vector<double> >();
    if ( !inputParameter.empty() )
    {
        CHECK( inputParameter.size() == muspace->dimension() ) << "parameter has a wrong size : "<< inputParameter.size() << " but must be " << muspace->dimension();
        auto mu = muspace->element();
        for ( uint16_type d=0;d<muspace->dimension();++d)
            mu(d)=inputParameter[d];
        mysampling->push_back( mu );
    }
    else
    {
        int nSample = ioption(_name="sampling.size");
        std::string sampler = soption("sampling.type");
        mysampling->sample( nSample, sampler );
    }

    if ( loadFiniteElementDatabase )
        plugin[0]->initExporter();

    int rbDim = ioption(_name="rb-dim");
    int nSamples = mysampling->size();
    for ( int k=0;k<nSamples;++k )
    {
        auto const& mu = (*mysampling)[k];
        std::ostringstream ostrmu;
        for ( uint16_type d=0;d<muspace->dimension();++d)
            ostrmu << mu(d) << " ";
        std::cout << "--------------------------------------\n";
        std::cout << "mu["<<k<<"] : " << ostrmu.str() << "\n";
        //auto mu = crb->Dmu()->element();
        //std::cout << "input mu\n" << mu << "\n";
        for( auto const& p : plugin )
        {
            auto crbResult = p->run( mu, time_crb, online_tol, rbDim, print_rb_matrix);
            auto resOuptut = boost::get<0>( crbResult );
            auto resError = boost::get<0>( boost::get<6>( crbResult ) );
            std::cout << "output " << resOuptut.back() << " " << resError.back() << "\n";

            if ( loadFiniteElementDatabase )
            {
                p->exportField( (boost::format("sol-%1%")%k).str(), crbResult );
            }
        }
    }
    if ( loadFiniteElementDatabase )
        plugin[0]->saveExporter();
    
    return true;

}

std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin( std::string const& pluginname,
            std::string const& dirname,
            std::string const& libname,
            bool loadFiniteElementDatabase,
            int loadLast )
{
    using namespace Feel;
    namespace dll=boost::dll;

    std::string pluginlibname = "";
#if 0
    if ( Environment::vm().count("plugin.libname") )
        pluginlibname = soption(_name="plugin.libname");
#endif
    std::cout << "try loading the plugin " << pluginname << std::endl;
    auto plugin = factoryCRBPlugin( pluginname, pluginlibname, dirname );
    std::cout << "Loaded the plugin " << plugin->name() << std::endl;

    std::string jsonfilename;
    if ( loadLast )
    {
        plugin->loadDBLast( static_cast<crb::last>(loadLast), (loadFiniteElementDatabase)? crb::load::all : crb::load::rb );
    }
    else
    {
        std::string plugindbid = Environment::expand( soption(_name="plugin.dbid") );
        std::string jsonfilename = (fs::path(Environment::expand( soption(_name="plugin.db") )) / fs::path(pluginname) / fs::path(plugindbid) / (pluginname+".crb.json")).string() ;
        std::cout << " . using db " << jsonfilename << std::endl;

        plugin->loadDB( jsonfilename, (loadFiniteElementDatabase)? crb::load::all : crb::load::rb );
    }
    return plugin;
}

} // Feel
