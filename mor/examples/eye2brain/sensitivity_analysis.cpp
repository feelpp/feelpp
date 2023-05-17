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
//! @author Thomas Saigre <saigre@math.unistra.fr>
//! @date 01 April 2023 🐟
//! @copyright 2023 Feel++ Consortium
//!
#include <boost/dll.hpp>
#include <boost/algorithm/string/split.hpp>
#include <openturns/OT.hxx>

#include <feel/feelcore/table.hpp>
#include <feel/feelmor/options.hpp>
#include <feel/feelmor/crbplugin_interface.hpp>
#include <feel/feelmor/crbmodeldb.hpp>

#include <iostream>
#include <ctime>
#include <execution>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/builder/basic/array.hpp>
#endif

// #include <omp.h>
#include "tqdm/tqdm.h"
#include "results.hpp"
#include "FunctionalChaos.hpp"
// #include "MartinezSensitivity.hpp"

using element_t = Feel::ParameterSpaceX::element_type;
using plugin_ptr_t = std::shared_ptr<Feel::CRBPluginAPI>;
using parameter_space_ptr_t = std::shared_ptr<Feel::ParameterSpaceX>;


const std::vector<std::string> NAMES = {"h_bl", "h_amb", "T_bl", "T_amb", "E", "k_lens"};
const std::vector<double> MINS       = {50    , 8      , 308.3 , 283.15 , 20 , 0.21    };
const std::vector<double> MAXS       = {110   , 100    , 312   , 303.15 , 320, 0.544   };
const size_t SIZE = 6;


std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin()
{
    using namespace Feel;

    std::string crbmodelName = Environment::expand( soption(_name="crbmodel.name") );
    CRBModelDB crbmodelDB{ crbmodelName, uuids::nil_uuid() };

    std::string attribute = soption(_name="crbmodel.attribute" );
    std::string attribute_data;
    if ( attribute == "id"  || attribute == "name")
    {
        attribute_data = Environment::expand( soption(_name=fmt::format("crbmodel.db.{}",attribute) ) );
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
    Feel::cout << "-- crbmodelDB::dbRepository()=" << crbmodelDB.dbRepository() << std::endl;

    return crbmodelDB.loadDBPlugin( meta, soption(_name="crbmodel.db.load" ) );
}



inline Feel::AboutData makeAbout()
{
    Feel::AboutData about( "sensitivity_analysis",
                     "SA" ,
                     "0.1",
                     "Sensitivity analysis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2023 Feel++ Consortium" );

    about.addAuthor( "Thomas Saigre", "developer", "saigre@math.unistra.fr", "" );
    return about;
}

/**
 * @brief Generate of composed distribution for the eye model
 *
 * @return OT::ComposedDistribution
 */
OT::ComposedDistribution composedFromModel()
{
    using namespace Feel;

    OT::Collection<OT::Distribution> marginals( SIZE );

    for (size_t d = 0; d < SIZE; ++d)
    {
        OT::Distribution dist;
        std::string name = NAMES[d];

        if (name == "h_bl")
        {
            double s = 0.15; double mu = log(65) - 0.5*s*s;
            dist = OT::TruncatedDistribution(OT::LogNormal(mu, s, 0), OT::Interval(50, 120));
        }
        else if (name == "h_amb")
        {
            double s = 1; double mu = log(10) - 0.5*s*s;
            dist = OT::TruncatedDistribution(OT::LogNormal(mu, s, 8), OT::Interval(8, 100));
        }
        else if (name == "E")
        {
            // double s = 0.7; double mu = log(40.) - 0.5*s*s;
            // dist = OT::TruncatedDistribution(OT::LogNormal(mu, s, 20), OT::Interval(20, 130));
            double Emin = 20, Emax = 320;
            dist = OT::Uniform( MINS[d], MAXS[d] );
        }
        else
        {
            dist = OT::Uniform( MINS[d], MAXS[d] );
        }

        dist.setDescription( {name} );
        marginals[d] = dist;
        Feel::cout << tc::blue << "Distribution " << d << " (" << name << ") = " << dist << tc::reset << std::endl;
    }

    return OT::ComposedDistribution( marginals );
}

/**
 * @brief Generate the output sample from a given input sample
 *
 * @param input Sample of input parameters
 * @param plugin loaded plugin
 * @param time_crb collection of timers
 * @param online_tol online tolerance
 * @param rbDim size of the reduced basis
 * @return OT::Sample
 */
OT::Sample output(OT::Sample const& input, plugin_ptr_t const& plugin, Eigen::VectorXd &time_crb, double online_tol, int rbDim)
{
    size_t n = input.getSize();
    OT::Sample output(n, 1);
    {
        parameter_space_ptr_t Dmu = plugin->parameterSpace();
        std::vector<std::string> names = Dmu->parameterNames();
        Feel::cout << "Start to compute outputs, sampling of size " << n << std::endl;

        double k_lens, h_amb, h_bl, h_r=6, T_amb, T_bl, E;

        for (size_t i: tqdm::range(n))          // std::for_each
        {
            element_t mu = Dmu->element();
            OT::Point X = input[i];
            // X = [ h_bl, h_amb, T_bl, T_amb, E, k_lens ]
            k_lens = X[5];
            h_amb = X[1];
            h_bl = X[0];
            T_amb = X[3];
            T_bl = X[2];
            E = X[4];

            // "mu0": "k_lens:k_lens",
            mu.setParameter(0, k_lens);
            // "mu1": "h_amb:h_amb",
            mu.setParameter(1, h_amb);
            // "mu2": "h_bl:h_bl",
            mu.setParameter(2, h_bl);
            // "mu3": "h_r:h_r"
            mu.setParameter(3, h_r);
            // "mu4": "1",
//          mu.setParameter(4, 1);
            // "mu4": "h_amb * T_amb + h_r * T_amb - E:h_amb:T_amb:E",
            mu.setParameter(5, h_amb*T_amb + h_r*T_amb - E);
            // "mu5": "h_bl * T_bl:h_bl:T_bl"
            mu.setParameter(6, h_bl*T_bl);
            // Feel::cout << "mu = " << mu << std::endl;

            Feel::CRBResults crbResult = plugin->run( mu, time_crb, online_tol, rbDim, false );
            output[i] = OT::Point( crbResult.output() );
        }
        Feel::cout << "output computed" << std::endl;
    }
    return output;
}


void computeSobolIndicesBootstrap(std::vector<plugin_ptr_t> plugin, OT::ComposedDistribution composed_distribution, size_t sampling_size,
    OT::Sample input_sample, OT::Sample output_sample, std::vector<std::string> tableRowHeader, Eigen::VectorXd &time_crb, double online_tol, int rbDim)
{
    using namespace Feel;
    size_t dim = composed_distribution.getDimension();
    size_t bootstrap_size = ioption(_name="algo.bootstrap-size");
    Feel::cout << tc::bold << tc::red << "Run polynomial chaos and bootstrap : sampling of size " << sampling_size
        << " (bootstrap size " << bootstrap_size << ")" << tc::reset << std::endl;

    OT::Collection<OT::Distribution> marginals(dim);
    for ( size_t d = 0; d < dim; ++d )
        marginals[d] = composed_distribution.getMarginal(d);
    auto basis = OT::OrthogonalProductPolynomialFactory( marginals );
    OT::UnsignedInteger total_degree = 3;

    Results res( dim, tableRowHeader, "polynomial-chaos-bootstrap", sampling_size );


    // Check the meta-model
    if ( boption(_name="algo.check-meta-model"))
    {
        Feel::cout << "Compute Sparse Least Squares Chaos" << std::endl;
        tic();
        OT::FunctionalChaosResult polynomialChaosResult =
            computeSparseLeastSquaresChaos(input_sample, output_sample, basis, total_degree, composed_distribution);
        toc("computeSparseLeastSquaresChaos");
        tic();
        OT::Function metaModel = polynomialChaosResult.getMetaModel();
        OT::UnsignedInteger n_valid = 1000;
        OT::Sample X_test = composed_distribution.getSample(n_valid);
        OT::Sample Y_test = output(X_test, plugin[0], time_crb, online_tol, rbDim);
        checkMetaModel( X_test, Y_test, metaModel );
        toc("checkMetaModel");
    }

    OT::Graph graph = computeAndDrawSobolIndices( res, input_sample, output_sample, basis, total_degree, composed_distribution, bootstrap_size=bootstrap_size);
    res.print();
    res.exportValues( soption( _name="save.path" ) + "-bootstrap.json" );
    //graph.draw("sobol-indices.png");
    Feel::cout << "Files are saved in " << fs::current_path() << std::endl;
}


/**
 * @brief Compute sobol indices
 *
 * @param plugin std::vector containing the plugin from load_plugin
 * @param sampling_size size of the input sample used for computation of sobol indices
 * @param rbDim size of the reduced basis
 * @param computeSecondOrder boolean to compute second order sobol indices
 */
void runSensitivityAnalysis( std::vector<plugin_ptr_t> plugin, size_t sampling_size, int rbDim, bool computeSecondOrder=true )
{
    using namespace Feel;

    Feel::cout << "Running sensisivity analysis with a sample of size " << sampling_size << std::endl;

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;               //Feel::doption(Feel::_name="crb.online-tolerance");
    bool print_rb_matrix = false;           //boption(_name="crb.print-rb-matrix");
    parameter_space_ptr_t muspace = plugin[0]->parameterSpace();

    OT::ComposedDistribution composed_distribution = composedFromModel();


    double adapt_tol = doption(_name="adapt.tol");

    OT::Sample input_sample = composed_distribution.getSample(sampling_size);
    tic();
    OT::Sample output_sample = output(input_sample, plugin[0], time_crb, online_tol, rbDim);
    toc("output sample");

    Feel::cout << "input_sample =\n" << input_sample << std::endl;
    Feel::cout << "output_sample =\n" << output_sample << std::endl;


    computeSobolIndicesBootstrap(plugin, composed_distribution, sampling_size, input_sample, output_sample, NAMES, time_crb, online_tol, rbDim);

}

int main( int argc, char** argv )
{
    using namespace Feel;
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

        ( "parameter", po::value<std::vector<std::string> >()->multitoken(), "database filename" )
        ( "sampling.size", po::value<int>()->default_value( 2000 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
        ( "rb-dim", po::value<int>()->default_value( -1 ), "reduced basis dimension used (-1 use the max dim)" )
        ( "output_results.save.path", po::value<std::string>(), "output_results.save.path" )

        ( "algo.poly", po::value<bool>()->default_value(true), "use polynomial chaos" )
        ( "algo.bootstrap", po::value<bool>()->default_value(true), "use polynomial chaos and bootstrap" )
        ( "algo.nrun", po::value<int>()->default_value(5), "number to run algorithm" )
        ( "adapt.tol", po::value<double>()->default_value(0.01), "tolerance for adaptative algoritmh" )
        ( "algo.bootstrap-size", po::value<int>()->default_value(100), "bootstrap size for sensitivity analysis" )
        ( "algo.check-meta-model", po::value<bool>()->default_value(false), "Check the metamodel" )
        ( "save.path", po::value<std::string>()->default_value( "sensitivity" ))

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

    Environment env( _argc = argc, _argv = argv,
                     _desc = crbonlinerunoptions,
                     _desc_lib = crbonlinerunliboptions.add( feel_options() ),
                     _about = makeAbout() );

    OT::RandomGenerator::SetSeed( ::time(NULL) );
    plugin_ptr_t plugin = loadPlugin();
    int rbDim = ioption(_name="rb-dim");
    // runCrbOnline( { plugin } );
    runSensitivityAnalysis( { plugin }, ioption(_name="sampling.size"), rbDim, false );

    Feel::cout << tc::green << "Done ✓" << tc::reset << std::endl;
    return 0;
}
