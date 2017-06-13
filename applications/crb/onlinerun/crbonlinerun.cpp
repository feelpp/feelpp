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

#include <feel/feelcrb/options.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>

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
runCrbOnline()
{
    using namespace Feel;
    namespace dll=boost::dll;
    std::string dirname = Environment::expand( soption(_name="plugin.dir") );
    std::string pluginname = Environment::expand( soption(_name="plugin.name") );
    std::string plugindbid = Environment::expand( soption(_name="plugin.dbid") );
    std::string jsonfilename = (fs::path(Environment::expand( soption(_name="plugin.db") )) / fs::path(pluginname) / fs::path(plugindbid) / (pluginname+".crb.json")).string() ;
    
    boost::function<crbpluginapi_create_t> creator;
    fs::path pname = fs::path(dirname) / ("libfeelpp_crb_" + pluginname + ".so");
    
    creator = boost::dll::import_alias<crbpluginapi_create_t>(pname,
                                                              "create_crbplugin",
                                                              dll::load_mode::append_decorations );
    auto plugin = creator();
    std::cout << "Loaded the plugin " << plugin->name() << std::endl
              << " . from " << pname.string() << std::endl
              << " . using db " << jsonfilename << std::endl;
    plugin->loadDB( jsonfilename );

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;//Feel::doption(Feel::_name="crb.online-tolerance");
    bool print_rb_matrix = false;//boption(_name="crb.print-rb-matrix");
    auto muspace = plugin->parameterSpace();

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

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    if ( loadFiniteElementDatabase )
        plugin->initExporter();

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
        auto crbResult = plugin->run( mu, time_crb, online_tol, -1, print_rb_matrix);
        auto resOuptut = boost::get<0>( crbResult );
        auto resError = boost::get<0>( boost::get<6>( crbResult ) );
        std::cout << "output " << resOuptut.back() << "\n";
        std::cout << "err " << resError.back() << "\n";


        if ( loadFiniteElementDatabase )
        {
            plugin->exportField( (boost::format("sol-%1%")%k).str(), crbResult );
        }
   }
    if ( loadFiniteElementDatabase )
        plugin->saveExporter();

    return true;

}


int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description crbonlinerunoptions( "crb online run options" );
	crbonlinerunoptions.add_options()
        ( "plugin.dir", po::value<std::string>()->default_value(Info::libdir()) , "plugin directory" )
        ( "plugin.name", po::value<std::string>(), "CRB online code name" )
        ( "plugin.dbid", po::value<std::string>(), "CRB online code id" )
        ( "plugin.db", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )
        ( "parameter", po::value<std::vector<double> >()->multitoken(), "database filename" )
        ( "sampling.size", po::value<int>()->default_value( 10 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
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

    runCrbOnline();

    return 0;
}
