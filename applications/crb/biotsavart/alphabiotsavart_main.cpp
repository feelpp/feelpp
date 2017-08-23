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

#include "alphabiotsavart.hpp"
#include <alphathermoelectric.hpp>

int iter=0;

int main(int argc, char**argv )
{
    using namespace Feel;
    using Feel::cout;

    po::options_description nloptoptions( "NLOpt options" );
    nloptoptions.add_options()
        // ( "Tcritic", po::value<double>()->default_value( 400 ), "critical value for the temperature")
        ( "nlopt.algo", po::value<std::string>()->default_value( "LN_NEWUOA" ), "NLOPT algorithm [LN_NEWUOA,LD_LBFGS]" )
        ;

    nloptoptions.add(biotsavartOptions());
    nloptoptions.add(crbOptions());
    nloptoptions.add(crbSEROptions());
    nloptoptions.add(eimOptions());
    nloptoptions.add(eimOptions());
    nloptoptions.add(podOptions());
    nloptoptions.add(backend_options("backend-primal"));
    nloptoptions.add(backend_options("backend-dual"));
    nloptoptions.add(backend_options("backend-l2"));
    nloptoptions.add(bdf_options("ThermoElectricCRB"));
    nloptoptions.add(makeOptions());

    Environment env( _argc=argc, _argv=argv,
                     _desc=nloptoptions,
                     _about=about(_name="biotsavart_alpha",
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

    auto BS = AlphaBiotSavartCRB<AlphaThermoelectric>::New();
    BS->init();
    if( !boption("biotsavart.do-opt") )
        BS->runBS();
    else
    {
        // BS->offline();
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
                mu=Eigen::Map<Eigen::VectorXd const>( x.data(), mu.size());
                // BS->online(mu);
                // auto B = BS->magneticFlux();
                BS->computeFE(mu);
                auto B = BS->magneticFluxFE();
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

        // export

        auto eC = Exporter<Mesh<Simplex<3> > >::New( "conductor");
        eC->setMesh(BS->meshCond());
        auto eM = Exporter<Mesh<Simplex<3> > >::New( "magneto");
        eM->setMesh(BS->meshMgn());

        auto VTFe = BS->potentialTemperatureFE();
        auto VFe = VTFe.template element<0>();
        auto TFe = VTFe.template element<1>();
        auto BFe = BS->magneticFluxFE();
        auto alphaStr = BS->alpha(mu);
        Feel::cout << "alpha: " <<alphaStr << std::endl;
        auto alphaExpr = expr(BS->alpha(mu));
        auto Xh = BS->spaceCond();
        auto Vh = Xh->template functionSpace<0>();
        auto alpha = Vh->element();
        alpha.on( elements( BS->meshMgn()), alphaExpr);
        eC->add("alpha", alpha);
        eC->add("VFe", VFe);
        eC->add("TFe", TFe);
        eM->add("BFe", BFe);

        mu = BS->param0();
        BS->computeFE(mu);
        auto VT0 = BS->potentialTemperatureFE();
        auto V0 = VTFe.template element<0>();
        auto T0 = VTFe.template element<1>();
        auto B0 = BS->magneticFluxFE();
        eC->add("V0", V0);
        eC->add("T0", T0);
        eM->add("B0", B0);

        eC->save();
        eM->save();

        auto homoFE = BS->homogeneity(BFe);
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
