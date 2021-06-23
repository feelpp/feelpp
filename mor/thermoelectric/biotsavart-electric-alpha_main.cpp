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


#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelopt/nlopt.hpp>

#include "biotsavartbase.hpp"
#include "biotsavart-electric-alpha.hpp"
#include <electric-alpha.hpp>
#include <feel/feelmodels/electric/electric.hpp>

int iter=0;

// std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> >
// computeStats(std::vector<std::vector<double> > const& err)
// {
//     int M = err.size();
//     std::vector<double> min(M), max(M), mean(M), stdev(M);
//     if( M > 0 )
//     {
//         int size = err[0].size();
//         for(int m = 0; m < M; ++m)
//         {
//             min[m] = *std::min_element(err[m].begin(), err[m].end());
//             max[m] = *std::max_element(err[m].begin(), err[m].end());
//             double s = std::accumulate(err[m].begin(), err[m].end(), 0.0);
//             mean[m] = s/size;
//             double accum = std::accumulate(err[m].begin(), err[m].end(), 0.0,
//                                            [s,size](double a, double b) {
//                                                return a + (b-s/size)*(b-s/size);
//                                        });
//             stdev[m] = accum/size;
//         }
//     }
//     return std::make_tuple(min,max,mean,stdev);
// }

// void writeErrors(fs::ofstream& out, std::vector<std::vector<double> > const& err)
// {
//     if( out && Environment::isMasterRank() )
//     {
//         int M = err.size();
//         int size = err[0].size();
//         out << std::setw(5) << "M";
//         for(int i = 0; i < size; ++i)
//             out << std::setw(24) << "mu_" << i;
//         out << std::endl;
//         for(int m = 0; m < M; ++m)
//         {
//             out << std::setw(5) << m+1;
//             for(int i = 0; i < size; ++i)
//                 out << std::setw(25) << err[m][i];
//             out << std::endl;
//         }
//         out.close();
//     }
// }

// void writeStats(std::ostream& out, std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > const& stat, std::string const& base = "M")
// {
//     auto min = std::get<0>(stat);
//     auto max = std::get<1>(stat);
//     auto mean = std::get<2>(stat);
//     auto stdev = std::get<3>(stat);
//     if( Environment::isMasterRank() )
//     {
//         int M = min.size();
//         out << std::setw(5) << base << std::setw(25) << "min" << std::setw(25) << "max"
//             << std::setw(25) << "mean" << std::setw(25) << "stdev" << std::endl;
//         for(int m = 0; m < M; ++m)
//             out << std::setw(5) << m+1 << std::setw(25) << min[m] << std::setw(25) << max[m]
//                 << std::setw(25) << mean[m] << std::setw(25) << stdev[m] << std::endl;
//     }
// }

int main(int argc, char**argv )
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using Feel::cout;

    using biotsavart_type = BiotSavartAlphaElectricCRB<AlphaElectric>;
    using electric_tb_type = Electric<Simplex<3,1>, Lagrange<1, Scalar,Continuous,PointSetFekete> >;
    using sampling_type = typename biotsavart_type::crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;

    po::options_description nloptoptions( "NLOpt options" );
    nloptoptions.add_options()
        ( "nlopt.algo", po::value<std::string>()->default_value( "LN_NEWUOA" ), "NLOPT algorithm [LN_NEWUOA,LD_LBFGS]" )
        ( "biotsavart.do-opt", po::value<bool>()->default_value(false),
          "do or not the optimization" )
        ( "biotsavart.use-bg-field", po::value<bool>()->default_value(false),
          "use a background field for the optimization" )
        ( "biotsavart.size", po::value<int>()->default_value(1), "size of convergence sampling" )
        ( "biotsavart.use-rb-cvg", po::value<bool>()->default_value(true), "" )
        ;

    nloptoptions.add(biotsavart_type::makeOptions("biotsavart"));
    nloptoptions.add(electricity_options("toolbox"));

    Environment env( _argc=argc, _argv=argv,
                     _desc=nloptoptions,
                     _about=about(_name="biotsavart_alpha_electro",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto BS = BiotSavartAlphaElectricCRB<AlphaElectric>::New(crb::stage::offline, "biotsavart");
    BS->initModel();
    if( !boption("biotsavart.do-opt") )
    {
        auto sampling = sampling_type( BS->parameterSpace());
        int size = ioption("biotsavart.size");
        if( size == -1 )
        {
            auto mus = BS->deim()->mus();
            sampling.setElements(mus);
            size = mus.size();
        }
        else if( size == 0 )
        {
            auto mu = BS->paramFromProperties();
            sampling.setElements({mu});
            size = 1;
        }
        else
        {
            sampling.clear();
            sampling.randomize( size, true );
        }

        auto Xh = BS->spaceMgn();
        auto rangeB = Xh->dof()->meshSupport()->rangeElements();
        auto BFE = Xh->element();
        auto BRB = Xh->element();
        // auto BRBn = Xh->element();
        auto Vh = BS->spaceCond();
        auto rangeV = Vh->dof()->meshSupport()->rangeElements();
        auto VFE = Vh->element();
        auto VRB = Vh->element();

        int M = BS->dimension();
        std::vector<double> errsV(size);
        std::vector<double> errsVR(size);
        std::vector<double> errsB(size);
        std::vector<double> errsBR(size);
        std::vector<double> errsH(size);
        std::vector<double> errsHR(size);
        // std::vector<std::vector<double> > errs(M, std::vector<double>(size));
        // std::vector<std::vector<double> > errsRel(M, std::vector<double>(size));
        // std::vector<std::vector<double> > errsH(M, std::vector<double>(size));
        // std::vector<std::vector<double> > errsRelH(M, std::vector<double>(size));
        // int N = BS->crbDimension();
        // std::vector<std::vector<double> > errsV(N, std::vector<double>(size));
        // std::vector<std::vector<double> > errsRelV(N, std::vector<double>(size));

        Feel::cout << "start convergence study with " << size << " parameters" << std::endl;
        int i = 0;
        for( auto const& mu : sampling )
        {
            if( boption("biotsavart.use-rb-cvg") )
            {
                BFE = *BS->assembleForDEIM(mu,0);
            }
            else
            {
                BS->computeFE(mu);
                VFE = BS->potentialFE();
                BFE = BS->magneticFluxFE();
            }
            double normV = normL2( rangeV, idv(VFE) );
            double normB = normL2( rangeB, idv(BFE) );
            double homoFE = BS->homogeneity(BFE);
            BS->online(mu);
            VRB = BS->potential();
            BRB = BS->magneticFlux();
            double homoRB = BS->homogeneity(BRB);
            errsV[i] = normL2(rangeV, idv(VRB)-idv(VFE) );
            errsVR[i] = errsV[i]/normV;
            errsB[i] = normL2(rangeB, idv(BRB)-idv(BFE) );
            errsBR[i] = errsB[i]/normB;
            errsH[i] = std::abs(homoRB-homoFE);
            errsHR[i] = errsH[i]/homoFE;
            // for( int m = 0; m < M; ++m)
            // {
            //     BS->online(mu, m+1);
            //     BRB = BS->magneticFlux();
            //     double homoRB = BS->homogeneity(BRB);
            //     errs[m][i] = normL2( rangeB, idv(BRB)-idv(BFE) );
            //     errsRel[m][i] = errs[m][i]/normB;
            //     errsH[m][i] = std::abs(homoRB-homoFE);
            //     errsRelH[m][i] = errs[m][i]/homoFE;
            // }
            // for( int n = 0; n < N; ++n)
            // {
            //     BS->computeRB(mu, n+1);
            //     BRBn = BS->magneticFlux();
            //     errsV[n][i] = normL2( rangeB, idv(BRBn)-idv(BFE) );
            //     errsRelV[n][i] = errsV[n][i]/normB;
            // }
            ++i;
        }

        double s;
        double minV = *std::min_element(errsVR.begin(), errsVR.end());
        double maxV = *std::max_element(errsVR.begin(), errsVR.end());
        s = std::accumulate(errsVR.begin(), errsVR.end(), 0.0);
        double meanV = s/size;
        double minB = *std::min_element(errsBR.begin(), errsBR.end());
        double maxB = *std::max_element(errsBR.begin(), errsBR.end());
        s = std::accumulate(errsBR.begin(), errsBR.end(), 0.0);
        double meanB = s/size;
        double minH = *std::min_element(errsHR.begin(), errsHR.end());
        double maxH = *std::max_element(errsHR.begin(), errsHR.end());
        s = std::accumulate(errsHR.begin(), errsHR.end(), 0.0);
        double meanH = s/size;

        Feel::cout << std::setw(5) << "V" << std::setw(24) << "min"
                   << std::setw(24) << "max" << std::setw(24) << "mean" << std::endl;
        Feel::cout << std::setw(5) << " " << std::setw(24) << minV
                   << std::setw(24) << maxV << std::setw(24) << meanV << std::endl;
        Feel::cout << std::setw(5) << "B" << std::setw(24) << "min"
                   << std::setw(24) << "max" << std::setw(24) << "mean" << std::endl;
        Feel::cout << std::setw(5) << " " << std::setw(24) << minB
                   << std::setw(24) << maxB << std::setw(24) << meanB << std::endl;
        Feel::cout << std::setw(5) << "H" << std::setw(24) << "min"
                   << std::setw(24) << "max" << std::setw(24) << "mean" << std::endl;
        Feel::cout << std::setw(5) << " " << std::setw(24) << minH
                   << std::setw(24) << maxH << std::setw(24) << meanH << std::endl;

        // auto stats = computeStats(errsRel);

        // fs::ofstream cvgErr( "err.dat" );
        // fs::ofstream cvgErrR( "errR.dat" );
        // fs::ofstream cvgStat( "stat.dat" );
        // writeErrors(cvgErr, errs);
        // writeErrors(cvgErrR, errsRel);
        // if( cvgStat )
        // {
        //     writeStats(cvgStat, stats);
        //     cvgStat.close();
        // }
        // if( Environment::isMasterRank() )
        //     writeStats(std::cout, stats);

        // auto statsH = computeStats(errsRelH);

        // fs::ofstream cvgErrH( "errH.dat" );
        // fs::ofstream cvgErrRH( "errRH.dat" );
        // fs::ofstream cvgStatH( "statH.dat" );
        // writeErrors(cvgErrH, errsH);
        // writeErrors(cvgErrRH, errsRelH);
        // if( cvgStatH )
        // {
        //     writeStats(cvgStatH, statsH);
        //     cvgStatH.close();
        // }
        // if( Environment::isMasterRank() )
        //     writeStats(std::cout, statsH);

        // auto statsV = computeStats(errsRelV);

        // fs::ofstream cvgErrV( "errV.dat" );
        // fs::ofstream cvgErrRV( "errRV.dat" );
        // fs::ofstream cvgStatV( "statV.dat" );
        // writeErrors(cvgErrV, errsV);
        // writeErrors(cvgErrRV, errsRelV);
        // if( cvgStatV )
        // {
        //     writeStats(cvgStatV, statsV);
        //     cvgStatV.close();
        // }
        // if( Environment::isMasterRank() )
        //     writeStats(std::cout, statsV);

        auto mu = sampling.back();
        auto alpha = BS->alpha(mu);
        // BS->computeFE(mu);
        // BS->expandV();
        // auto VRB = BS->potential();
        // auto VFE = BS->potentialFE();
        mu = BS->param0();
        BS->computeFE(mu);
        auto V0 = BS->potentialFE();
        auto B0 = BS->magneticFluxFE();

        auto eCM = exporter(_mesh=BRB.functionSpace()->mesh(), _name="biotsavart");

        eCM->add("alpha", alpha);
        eCM->add( "V_RB", VRB);
        eCM->add( "B_RB", BRB);
        eCM->add( "V_FE", VFE);
        eCM->add( "B_FE", BFE );
        eCM->add( "V_0", V0);
        eCM->add( "B_0", B0);
        eCM->save();
    }
    else
    {
        auto tb = std::make_shared<electric_tb_type>("toolbox");
        auto bBg = BS->spaceMgn()->element();
        if( boption("biotsavart.use-bg-field") )
        {
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
        auto algo = nloptAlgoMap.at(salgo);
        opt::OptimizationNonLinear opt( algo, N );

        // lower and upper bounds
        auto muMin = BS->parameterSpace()->min();
        auto muMax = BS->parameterSpace()->max();
        auto mu0 = BS->paramFromProperties();

        std::vector<double> lb(muMin.data(), muMin.data()+muMin.size());
        std::vector<double> ub(muMax.data(), muMax.data()+muMax.size());

        Feel::cout << "lb: [";
        for(int j = 0; j < lb.size()-1; ++j)
            Feel::cout << lb[j] << ", ";
        Feel::cout << lb.back() << "]" << std::endl;
        Feel::cout << "ub: [";
        for(int j = 0; j < ub.size()-1; ++j)
            Feel::cout << ub[j] << ", ";
        Feel::cout << ub.back() << "]" << std::endl;
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
                auto h = BS->homogeneity(B);
                Feel::cout << "iter " << iter << " homogeneity for parameter mu=[";
                for(int j = 0; j < x.size()-1; ++j)
                    Feel::cout << x[j] << ", ";
                Feel::cout << x.back() << "] = " << h << std::endl;
                return h;
            };

        opt.set_min_objective( myfunc, nullptr );

        // optimization
        double minf;
        tic();
        std::vector<double> x(N);
        Eigen::VectorXd::Map(x.data(), N) = mu0;

        Feel::cout << "start optimization with mu = " << x << std::endl;

        ::nlopt::result result = opt.optimize(x, minf);

        double optiTime = toc("optimization", false);

        auto mu = BS->newParameter();
        mu = Eigen::Map<Eigen::VectorXd>( x.data(), mu.size());

        Feel::cout << nloptResultMap.at(result) << std::endl
                   << iter << " iterations in " << optiTime
                   << " (" << optiTime/iter << "/iter)" << std::endl;

        auto alpha = BS->alpha(mu);
        auto V = BS->potential();
        auto B = BS->magneticFlux();
        if( boption("biotsavart.use-bg-field") )
            B += bBg;
        auto div = integrate(_range=elements(V.functionSpace()->mesh()),
                             _expr= -58e6*laplacianv(V)).evaluate()(0,0);
        BS->computeFE(mu);
        auto VFE = BS->potentialFE();
        auto BFE = BS->magneticFluxFE();
        if( boption("biotsavart.use-bg-field") )
            BFE += bBg;
        auto divFE = integrate(_range=elements(V.functionSpace()->mesh()),
                               _expr= -58e6*laplacianv(V)).evaluate()(0,0);
        BS->computeFE(mu0);
        auto V0FE = BS->potentialFE();
        auto B0FE = BS->magneticFluxFE();
        if( boption("biotsavart.use-bg-field") )
            B0FE += bBg;
        BS->online(mu0);
        auto B0 = BS->magneticFlux();

        if( result > 0 )
        {
            auto homo = BS->homogeneity(B);
            auto homoFE = BS->homogeneity(BFE);
            auto homo0 = BS->homogeneity(B0);
            auto homo0FE = BS->homogeneity(B0FE);
            Feel::cout << "homogeneity = " << homo << " (" << homoFE
                       << ") for parameter\n" << mu << "\n"
                       << "default homogeneity = " << homo0 << " (" << homo0FE << ")\n"
                       << "Evaluation number: " << iter << std::endl;
            Feel::cout << "divergence = " << div << " (" << divFE << ")"<< std::endl;
        }

        // export

        auto eCM = exporter(_mesh=V.functionSpace()->mesh(), _name="biotsavart");

        eCM->add("alpha", alpha);
        eCM->add( "V", V);
        eCM->add( "B", B);
        eCM->add( "VFE", VFE);
        eCM->add( "BFE", BFE );
        eCM->add( "V0", V0FE);
        eCM->add( "B0FE", B0FE);
        eCM->add( "B0", B0);
        if( boption("biotsavart.use-bg-field") )
        {
            eCM->add("Bbg", bBg);
            eCM->add("Vbg", tb->fieldElectricPotential());
            eCM->add("jbg", tb->fieldCurrentDensity());
        }
        eCM->save();

    }
    return 0;
}
