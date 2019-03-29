//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
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
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/feel.hpp>
#include <feel/feelopt/nlopt.hpp>

#include "biotsavartbase.hpp"
#include "biotsavart-electric-alpha.hpp"
#include <electric-alpha.hpp>
#include <feel/feelmodels/electric/electric.hpp>

int iter=0;

int main(int argc, char**argv )
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using Feel::cout;

    using biotsavart_type = BiotSavartAlphaElectricCRB<AlphaElectric>;
    using electric_tb_type = Electric<Simplex<3,1>, Lagrange<1, Scalar,Continuous,PointSetFekete> >;

    po::options_description nloptoptions( "NLOpt options" );
    nloptoptions.add_options()
        ( "nlopt.algo", po::value<std::string>()->default_value( "LN_NEWUOA" ), "NLOPT algorithm [LN_NEWUOA,LD_LBFGS]" )
        ( "biotsavart.do-opt", po::value<bool>()->default_value(false),
          "do or not the optimization" )
        ( "biotsavart.use-bg-field", po::value<bool>()->default_value(false),
          "use a background field for the optimization" )
        ;

    nloptoptions.add(biotsavart_type::makeOptions("biotsavart"));
    nloptoptions.add(electricity_options("toolbox"));

    Environment env( _argc=argc, _argv=argv,
                     _desc=nloptoptions,
                     _about=about(_name="biotsavart_alpha_electro",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    // List of NLOPT algorithm
    const boost::unordered_map< const std::string, ::nlopt::algorithm >& authAlgo = boost::assign::map_list_of
        ("LN_NEWUOA", ::nlopt::LN_NEWUOA )
        ("LN_COBYLA", ::nlopt::LN_COBYLA )
        ("LN_BOBYQA", ::nlopt::LN_BOBYQA )
        ("LD_LBFGS", ::nlopt::LD_LBFGS )
        ("LD_MMA", ::nlopt::LD_MMA )
        ("LD_SLSQP", ::nlopt::LD_SLSQP );

    auto BS = BiotSavartAlphaElectricCRB<AlphaElectric>::New(crb::stage::offline, "biotsavart");
    BS->initModel();
    if( !boption("biotsavart.do-opt") )
    {
        auto mu = BS->param0();
        BS->computeFE(mu);
        auto V0 = BS->potentialFE();
        auto B0 = BS->magneticFluxFE();
        mu = BS->paramFromProperties();
        BS->computeFE(mu);
        auto VFE = BS->potentialFE();
        auto BFE = BS->magneticFluxFE();
        auto rangeB = BFE.functionSpace()->dof()->meshSupport()->rangeElements();
        double normBFE = normL2( rangeB, idv(BFE) );

        boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
        fs::ofstream file( "cvg.dat" );
        if( file && Environment::isMasterRank() )
            file << fmter % "M" % "errB" % "relErrB";

        int size = BS->dimension();
        std::vector<std::vector<double> > errs(size, std::vector<double>(2));
        for( int m = 1; m <= size; ++m)
        {
            Feel::cout << "M = " << m << std::endl;
            BS->online(mu, m);
            BS->expand();
            auto B = BS->magneticFlux();
            errs[m-1][0] = normL2( rangeB, idv(B)-idv(BFE) );
            errs[m-1][1] = errs[m-1][0]/normBFE;
            if( Environment::isMasterRank() )
                file << fmter % m % errs[m-1][0] % errs[m-1][1];
        }

        file.close();
        Feel::cout << fmter % "M" % "errB" % "relErrB";
        for( int m = 0; m < size; ++m)
            Feel::cout << fmter % (m+1) % errs[m][0] % errs[m][1];

        auto alpha = BS->alpha(mu);
        auto V = BS->potential();
        auto B = BS->magneticFlux();

        auto eCM = exporter(_mesh=V.functionSpace()->mesh(), _name="biotsavart");

        eCM->add("alpha", alpha);
        eCM->add( "V", V);
        eCM->add( "B", B);
        eCM->add( "VFE", VFE);
        eCM->add( "BFE", BFE );
        eCM->add( "V0", V0);
        eCM->add( "B0", B0);
        eCM->save();
    }
    else
    {
        auto bBg = BS->spaceMgn()->element();
        if( boption("biotsavart.use-bg-field") )
        {
            auto tb = std::make_shared<electric_tb_type>("toolbox", false);
            tb->setMesh(BS->mesh());
            tb->init();
            tb->printAndSaveInfo();
            tb->solve();
            auto jBg = tb->fieldCurrentDensity();
            auto bsbg = BiotSavartBase(tb->spaceElectricField(), BS->spaceMgn());
            bBg = bsbg.computeMagneticField(jBg);
        }

        int N = BS->nbParameters();

        auto salgo = soption("nlopt.algo");
        cout << "NLOP algorithm: " << salgo << "\n";
        auto algo = authAlgo.at(salgo);
        opt::OptimizationNonLinear opt( algo, N );

        // lower and upper bounds
        auto muMin = BS->parameterSpace()->min();
        auto muMax = BS->parameterSpace()->max();

        std::vector<double> lb(muMin.data(), muMin.data()+muMin.size());
        std::vector<double> ub(muMax.data(), muMax.data()+muMax.size());

        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);

        // stopping criteria
        opt.set_maxeval( ioption("nlopt.maxeval") );
        opt.set_xtol_rel( doption("nlopt.xtol_rel") );
        opt.set_ftol_rel( doption("nlopt.ftol_rel") );
        opt.set_xtol_abs( doption("nlopt.xtol_abs") );
        opt.set_ftol_abs( doption("nlopt.ftol_abs") );

        // Objective function.
        auto myfunc = [&]( const std::vector<double> &x, std::vector<double> &grad, void *my_func_data )->double
            {
                iter++;
                auto mu = BS->newParameter();
                mu = Eigen::Map<Eigen::VectorXd const>( x.data(), mu.size());
                BS->online(mu);
                auto B = BS->magneticFlux();
                // BS->computeFE(mu);
                // auto B = BS->magneticFluxFE();
                if( boption("biotsavart.use-bg-field") )
                    B += bBg;
                return BS->homogeneity(B);
            };

        opt.set_min_objective( myfunc, nullptr );

        // optimization
        double minf;
        tic();
        std::vector<double> x(N);
        Eigen::VectorXd::Map(x.data(), N) = muMin;

        ::nlopt::result result = opt.optimize(x, minf);

        double optiTime = toc("optimization", false);

        auto mu = BS->newParameter();
        mu = Eigen::Map<Eigen::VectorXd>( x.data(), mu.size());

        cout << iter << " iterations in " << optiTime << " (" << optiTime/iter << "/iter)" << std::endl;

        auto alpha = BS->alpha(mu);
        auto V = BS->potential();
        auto B = BS->magneticFlux();
        BS->computeFE(mu);
        auto VFE = BS->potentialFE();
        auto BFE = BS->magneticFluxFE();
        auto mu0 = BS->param0();
        BS->computeFE(mu0);
        auto V0 = BS->potentialFE();
        auto B0 = BS->magneticFluxFE();

        // export

        auto eCM = exporter(_mesh=V.functionSpace()->mesh(), _name="biotsavart");

        eCM->add("alpha", alpha);
        eCM->add( "V", V);
        eCM->add( "B", B);
        eCM->add( "VFE", VFE);
        eCM->add( "BFE", BFE );
        eCM->add( "V0", V0);
        eCM->add( "B0", B0);
        if( boption("biotsavart.use-bg-field") )
            eCM->add("Bbg", bBg);
        eCM->save();

        auto homoFE = BS->homogeneity(BFE);
        auto homo0 = BS->homogeneity(B0);
        std::stringstream ss;
        ss << "homogeneity = " << homoFE << " for parameter\n" << mu << "\n"
           << "default homogeneity = " << homo0 << "\n"
           << "Evaluation number: " << iter << "\n";

        switch( result )
        {
        case ::nlopt::FAILURE:
            Feel::cerr << "NLOPT Generic Failure!" << "\n";
            break;
        case ::nlopt::INVALID_ARGS:
            Feel::cerr << "NLOPT Invalid arguments!" << "\n";
            break;
        case ::nlopt::OUT_OF_MEMORY:
            Feel::cerr << "NLOPT Out of memory!" << "\n";
            break;
        case ::nlopt::ROUNDOFF_LIMITED:
            Feel::cerr << "NLOPT Roundoff limited!" << "\n";
            break;
        case ::nlopt::FORCED_STOP:
            Feel::cerr << "NLOPT Forced stop!" << "\n";
            break;
        case::nlopt::SUCCESS:
            cout << "NLOPT coefficient found! (status " << result << ") \n"
                 << ss.str();
            break;
        case ::nlopt::STOPVAL_REACHED:
            cout << "NLOPT Stop value reached!" << "\n";
            cout << "NLOPT coefficient found! (status " << result << ") \n"
                 << ss.str();
            break;
        case ::nlopt::FTOL_REACHED:
            cout << "NLOPT ftol reached!" << "\n";
            cout << "NLOPT coefficient found! (status " << result << ") \n"
                 << ss.str();
            break;
        case ::nlopt::XTOL_REACHED:
            cout << "NLOPT xtol reached!" << "\n";
            cout << "NLOPT coefficient found! (status " << result << ") \n"
                 << ss.str();
            break;
        case ::nlopt::MAXEVAL_REACHED:
            cout << "NLOPT Maximum number of evaluation reached!" << "\n";
            cout << "NLOPT coefficient found! (status " << result << ") \n"
                 << ss.str();
            break;
        case ::nlopt::MAXTIME_REACHED:
            cout << "NLOPT Maximum time reached" << "\n";
            cout << "NLOPT coefficient found! (status " << result << ") \n"
                 << ss.str();
            break;
        }
    }
    return 0;
}
