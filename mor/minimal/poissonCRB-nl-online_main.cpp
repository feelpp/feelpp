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

#include "poissonCRB-nl.hpp"

using namespace Feel;

void writeErrors(fs::ofstream& out, std::vector<std::vector<double> > const& err)
{
    if( out && Environment::isMasterRank() )
    {
        int N = err.size();
        int size = err[0].size();
        out << std::setw(5) << "N";
        for(int i = 0; i < size; ++i)
            out << std::setw(24) << "mu_" << i;
        out << std::endl;
        for(int n = 0; n < N; ++n)
        {
            out << std::setw(5) << n+1;
            for(int i = 0; i < size; ++i)
                out << std::setw(25) << err[n][i];
            out << std::endl;
        }
    }
}

void printStats(fs::ofstream& out, std::vector<double > const& min, std::vector<double > const& max, std::vector<double > const& mean)
{
    int N = min.size();
    if( out && Environment::isMasterRank() )
    {
        out << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
                << std::setw(25) << "mean" << std::endl;
        for(int n = 0; n < N; ++n)
            out << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                    << std::setw(25) << mean[n] << std::endl;
    }
    Feel::cout << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
               << std::setw(25) << "mean" << std::endl;
    for(int n = 0; n < N; ++n)
        Feel::cout << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                   << std::setw(25) << mean[n] << std::endl;
}

int main( int argc, char** argv)
{
    using rb_model_type = PoissonNL;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;
    using vectorN_type = Eigen::VectorXd;

    auto opt = rb_model_type::makeOptions();
    opt.add_options()
        ( "online.dbid", po::value<std::string>()->default_value(""), "id of the db to load" )
        ( "online.sampling-size", po::value<int>()->default_value(1), "number of parameters" )
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("PoissonCRB")) );

    auto crb = crb_type::New("poissonmodel-nl_crb", crb::stage::offline);
    crb->loadDBFromId( soption("online.dbid"), crb::load::all );
    auto model = crb->model()->model();

    sampling_ptrtype sampling( new sampling_type( model->parameterSpace() ) );
    int size = ioption("online.sampling-size");
    sampling->clear();
    sampling->randomize( size, true );

    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);

    auto Xh = model->functionSpace();
    auto rangeU = elements(Xh->mesh());
    auto VTFE = Xh->element();
    auto VFE = VTFE.template element<0>();
    auto TFE = VTFE.template element<1>();
    auto VTRB = Xh->element();
    auto VRB = VTRB.template element<0>();
    auto TRB = VTRB.template element<1>();

    std::vector<std::vector<std::vector<double> > > errs(2, std::vector<std::vector<double> >(N, std::vector<double>(size)));
    std::vector<std::vector<std::vector<double> > > errsRel(2, std::vector<std::vector<double> >(N, std::vector<double>(size)));
    Feel::cout << "start convergence study with " << size << " parameters" << std::endl;
    int i = 0;
    for( auto const& mu : *sampling )
    {
        VTFE = model->solve(mu);
        auto normV = normL2( rangeU, idv(VFE) );
        auto normT = normL2( rangeU, idv(TFE) );
        for(int n = 0; n < N; ++n)
        {
            crb->fixedPointPrimal(n+1, mu, uNs, uNolds, outputs);
            vectorN_type uN = uNs[0];
            VTRB = crb->expansion( uN, n+1 );
            errs[0][n][i] = normL2( rangeU, idv(VRB)-idv(VFE) );
            errs[1][n][i] = normL2( rangeU, idv(TRB)-idv(TFE) );
            errsRel[0][n][i] = errs[0][n][i]/normV;
            errsRel[1][n][i] = errs[1][n][i]/normT;
        }
        ++i;
    }

    std::vector<std::vector<double> > min(2, std::vector<double>(N)), max(2, std::vector<double>(N)), mean(2, std::vector<double>(N));
    for(int n = 0; n < N; ++n)
    {
        for( int i = 0; i < 2; ++i )
        {
            min[i][n] = *std::min_element(errsRel[i][n].begin(), errsRel[i][n].end());
            max[i][n] = *std::max_element(errsRel[i][n].begin(), errsRel[i][n].end());
            double s = std::accumulate(errsRel[i][n].begin(), errsRel[i][n].end(), 0.0);
            mean[i][n] = s/size;
        }
    }

    fs::ofstream cvgErrV( "errV.dat" ), cvgErrT( "errT.dat" );
    fs::ofstream cvgErrVR( "errVR.dat" ), cvgErrTR( "errTR.dat" );
    fs::ofstream cvgStatV( "statV.dat" ), cvgStatT( "statT.dat" );
    writeErrors(cvgErrV, errs[0]);
    writeErrors(cvgErrTR, errsRel[0]);
    writeErrors(cvgErrT, errs[1]);
    writeErrors(cvgErrTR, errsRel[1]);
    printStats(cvgStatV, min[0], max[0], mean[0]);
    printStats(cvgStatT, min[1], max[1], mean[1]);
    cvgErrV.close();
    cvgErrVR.close();
    cvgErrT.close();
    cvgErrTR.close();
    cvgStatV.close();
    cvgStatT.close();

    auto mesh = model->mesh();
    auto e = exporter(mesh);
    e->add("VFE", VFE);
    e->add("VRB", VRB);
    e->add("TFE", TFE);
    e->add("TRB", TRB);
    e->save();

    return 0;
}
