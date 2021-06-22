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

#include "thermoelectric-nl.hpp"
#include <feel/feel.hpp>
#include <feel/feelopt/nlopt.hpp>
#include <feel/feelcrb/ser.hpp>


using namespace Feel;

int iter=0;

int main( int argc, char** argv)
{
    using rb_model_type = ThermoElectricNL;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using ser_type = SER<crb_type>;
    using ser_ptrtype = std::shared_ptr<ser_type>;
    using vectorN_type = Eigen::VectorXd;

    auto paramopt = rb_model_type::makeOptions();
    paramopt.add_options()
        ( "use-ser", po::value<bool>()->default_value(false), "force the use of SER")
        ( "nlopt.algo", po::value<std::string>()->default_value( "LN_COBYLA" ), "NLOPT algorithm [LN_NEWUOA,LD_LBFGS]" )
        ( "online.dbid", po::value<std::string>()->default_value(""), "id of the db to load" )
        ( "online.sampling-size", po::value<int>()->default_value(1), "number of parameters" )
        ( "online.do-cvg", po::value<bool>()->default_value(true), "" )
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=paramopt
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("PoissonCRB")) );

    auto crb = crb_type::New(soption("thermoelectric.basename"), crb::stage::offline);
    auto crbmodel = crb->model();
    auto model = crb->model()->model();
    auto ser = std::make_shared<ser_type>(crb, crbmodel);
    if( boption("use-ser") )
        ser->run();
    else
        crb->offline();
    // crb->loadDBFromId( soption("online.dbid"), crb::load::all );
    model->initOutputsPoints();

    int N = crb->dimension();
    std::vector<vectorN_type> uNs(1, vectorN_type(N)), uNolds(1, vectorN_type(N));
    std::vector<double> outputs(1, 0);

    int P = crbmodel->parameterSpace()->dimension();
    auto salgo = soption("nlopt.algo");
    Feel::cout << "NLOP algorithm: " << salgo << "\n";
    auto algo = nloptAlgoMap.at(salgo);
    opt::OptimizationNonLinear opt( algo, P );

    // lower and upper bounds
    auto muMin = crb->Dmu()->min();
    auto muMax = crb->Dmu()->max();
    auto mu0 = model->parameterProperties();

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

    // std::vector<double> vexp = {6.7, 8.1, 22.6, 13.2, 16.0, 17.7, 20.3, 22.9, 25.1, 27.1, 28.6, 26.1, 24.2};
    // // std::vector<int> vi = {1,2,4,5,6,7,8,9,10,11,12,13,14};
    // std::vector<int> vi = {9};
    // Objective function.
    auto myfunc = [&]( const std::vector<double> &x, std::vector<double> &grad, void *my_func_data )->double
                      {
                          iter++;
                          auto mu = crb->Dmu()->element();
                          mu = Eigen::Map<Eigen::VectorXd const>( x.data(), mu.size());
                          crb->fixedPointPrimal(N, mu, uNs, uNolds, outputs);
                          auto Vs = model->computeOutputsPointsElectro(uNs[0]);
                          double V10 = std::abs(Vs[1]-Vs[0]);
                          double r = std::abs(V10-25.1);
                          Feel::cout << iter << " mu = " << mu.toString() << " U9 = " << Vs[0]
                                     << " U10 = " << Vs[1] << " U10 = " << V10
                                     << " erreur = " << r << std::endl;
                          // r += (Vs[0]+156.38395364759077-22.9)*(Vs[0]+156.38395364759077+22.9)/(22.9*22.9);
                          // for( int i = 0; i < vi.size(); ++i )
                          //     r += (Vs[vi[i]]-Vs[vi[i]-1]-vexp[i])*(Vs[vi[i]]-Vs[vi[i]-1]-vexp[i])/(vexp[i]*vexp[i]);
                          return r;
                      };

    opt.set_min_objective( myfunc, nullptr );

    // optimization
    double minf;
    tic();
    std::vector<double> x(P);
    Eigen::VectorXd::Map(x.data(), P) = mu0;

    ::nlopt::result result = opt.optimize(x, minf);

    double optiTime = toc("optimization", false);

    auto mu = crb->Dmu()->element();
    mu = Eigen::Map<Eigen::VectorXd>( x.data(), mu.size());

    Feel::cout << nloptResultMap.at(result) << std::endl
               << iter << " iterations in " << optiTime
               << " (" << optiTime/iter << "/iter)" << std::endl;
    Feel::cout << "mu = " << mu.toString() << std::endl;

    crb->fixedPointPrimal(N, mu, uNs, uNolds, outputs);
    auto Vs = model->computeOutputsPointsElectro(uNs[0]);
    Feel::cout << Vs << std::endl;

    return 0;
}

