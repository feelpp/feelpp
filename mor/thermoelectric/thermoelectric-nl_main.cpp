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
#include <feel/feelcrb/ser.hpp>

#define DO_CVG_TEST 0

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
    using rb_model_type = ThermoElectricNL;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using ser_type = SER<crb_type>;
    using ser_ptrtype = std::shared_ptr<ser_type>;
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;
    using vectorN_type = Eigen::VectorXd;

    auto opt = rb_model_type::makeOptions();
    opt.add_options()
        ( "online.dbid", po::value<std::string>()->default_value(""), "id of the db to load" )
        ( "online.sampling-size", po::value<int>()->default_value(1), "number of parameters" )
        ( "online.do-cvg", po::value<bool>()->default_value(true), "" )
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

    auto crb = crb_type::New(soption("thermoelectric.basename"), crb::stage::offline);
    auto crbmodel = crb->model();
    auto model = crb->model()->model();
    auto ser = std::make_shared<ser_type>(crb, crbmodel);
    if( boption("crb.rebuild-database") || boption("eim.rebuild-database") )
        ser->run();
    else
        crb->offline();
    // crb->loadDBFromId( soption("online.dbid"), crb::load::all );
    // model->initOutputsPoints();

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
    double outFE=1, outRB=1;

#if 0
    if( !boption("online.do-cvg") )
    {
        auto eim_k = model->scalarContinuousEim()[0];
        auto eim_sigma = model->scalarContinuousEim()[1];
        auto eim_grad = model->scalarDiscontinuousEim()[0];
        auto sigmaMax = model->parameterSpace()->min().parameterNamed("sigma");

        auto mu = model->parameterSpace()->element();
        mu << 2.65644e-08, 275.225, 0.0034501, 0.0494472, 0.0920968, 50501.4;

        VTFE = model->solve(mu);
        auto normV = normL2( rangeU, idv(VFE) );
        auto normT = normL2( rangeU, idv(TFE) );
        auto sigma0 = mu.parameterNamed("sigma");
        auto alpha = mu.parameterNamed("alpha");
        auto L = mu.parameterNamed("L");
        auto sigma = sigma0/(cst(1.)+alpha*(idv(TFE)-cst(293.)) );
        auto nSigma = normL2(rangeU,sigma/sigmaMax);
        auto qSigma = eim_sigma->q();
        auto betaSigmaFE = eim_sigma->beta(mu, VTFE);
        auto sigmaFE = Feel::expansion(qSigma,betaSigmaFE);
        auto k = sigma*L*idv(TFE);
        auto nK = normL2(rangeU,k);
        auto qK = eim_k->q();
        auto betaKFE = eim_k->beta(mu, VTFE);
        auto kFE = Feel::expansion(qK,betaKFE);
        auto gradgrad = inner(gradv(VFE));
        auto nGrad = normL2(rangeU,gradgrad);
        auto qGrad = eim_grad->q();
        auto betaGradFE = eim_grad->beta(mu, VTFE);
        auto gradFE = Feel::expansion(qGrad,betaGradFE);

        auto mesh = model->mesh();
        auto e = exporter(mesh);
        for(int n = 0; n < 2; ++n)
        {
            crb->computeProjectionInitialGuess(mu, n+1, uNs[0]);
            auto Init = crb->expansion( uNs[0], n+1 );
            auto VInit = Init.template element<0>();
            auto TInit = Init.template element<1>();
            Feel::cout << "solving for N=" << n+1 << "...";
            crb->fixedPointPrimal(n+1, mu, uNs, uNolds, outputs);
            Feel::cout << "solved" << std::endl;
            vectorN_type uN = uNs[0];
            VTRB = crb->expansion( uN, n+1 );
            auto errV = normL2(rangeU, idv(VRB)-idv(VFE));
            auto errT = normL2(rangeU, idv(TRB)-idv(TFE));
            Feel::cout << "["<<n+1<<"] error V= " << errV << "\t rel=" << errV/normV << std::endl;
            Feel::cout << "["<<n+1<<"] error T= " << errT << "\t rel=" << errT/normT << std::endl;
            auto betaSigma = eim_sigma->beta(mu, uN);
            auto sigmaRB = Feel::expansion(qSigma,betaSigma);
            auto errSigma = normL2(rangeU,idv(sigmaRB)-sigma/sigmaMax);
            Feel::cout << "["<<n+1<<"] error sigma=" << errSigma << "\t rel=" << errSigma/nSigma << std::endl;
            auto betaK = eim_k->beta(mu, uN);
            auto kRB = Feel::expansion(qK,betaK);
            auto errK = normL2(rangeU,idv(kRB)-k);
            Feel::cout << "["<<n+1<<"] error k=" << errK << "\t rel=" << errK/nK << std::endl;
            auto betaGrad = eim_grad->beta(mu, uN);
            auto gradRB = Feel::expansion(qGrad,betaGrad);
            auto errGrad = normL2(rangeU,idv(gradRB)-gradgrad);
            Feel::cout << "["<<n+1<<"] error grad=" << errGrad << "\t rel=" << errGrad/nGrad << std::endl;
            e->step(n)->add("VFE", VFE);
            e->step(n)->add("VRB", VRB);
            e->step(n)->add("TFE", TFE);
            e->step(n)->add("TRB", TRB);
            e->step(n)->add("sigma", sigma);
            e->step(n)->add("sigmaRB", sigmaRB);
            e->step(n)->add("k", k);
            e->step(n)->add("kRB", kRB);
            e->step(n)->add("grad", gradgrad);
            e->step(n)->add("gradRB", gradRB);
            e->step(n)->add("sigmaFE", sigmaFE);
            e->step(n)->add("kFE", kFE);
            e->step(n)->add("gradFE", gradFE);
            e->step(n)->add("Vinit", VInit);
            e->step(n)->add("Tinit", TInit);
            e->save();
        }
    }
    else
#endif
    {
        sampling_ptrtype sampling( new sampling_type( model->parameterSpace() ) );
        int size = ioption("online.sampling-size");
        sampling->clear();
        sampling->randomize( size, true );

        auto modelOutputs = model->modelProperties()->outputs().outputsOfType("point");
        auto ctxFeV = Xh->template functionSpace<0>()->context();
        auto ctxFeT = Xh->template functionSpace<1>()->context();
        for( auto const& [name,output] : modelOutputs )
        {
            node_type t(3);
            auto coord = expr<3,1>(output.getString("coord")).evaluate();
            t(0) = coord(0); t(1) = coord(1); t(2) = coord(2);
            if( output.getString("field") == "electric-potential" )
                ctxFeV.add( t );
            else if( output.getString("field") == "temperature" )
                ctxFeT.add( t );
        }

        int nbErr = 2 + modelOutputs.size() + (ioption("crb.output-index") != 0 ? 1 : 0);
        std::vector<std::vector<std::vector<double> > > errs(nbErr, std::vector<std::vector<double> >(N, std::vector<double>(size)));
        std::vector<std::vector<std::vector<double> > > errsRel(nbErr, std::vector<std::vector<double> >(N, std::vector<double>(size)));
        Feel::cout << "start convergence study with " << size << " parameters" << std::endl;
        int i = 0;
        for( auto const& mu : *sampling )
        {
            Feel::cout << "solving for mu = " << mu.toString() << std::endl;
            VTFE = model->solve(mu);
            // auto evV = evaluateFromContext( _context=ctxFeV, _expr=idv(VFE));
            // auto evT = evaluateFromContext( _context=ctxFeT, _expr=idv(TFE));
            if( ioption("crb.output-index") != 0 )
                outFE = model->output(ioption("crb.output-index"), mu, VTFE);
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
                // auto evRbV = model->computeOutputsPointsElectro(uN);
                // auto evRbT = model->computeOutputsPointsThermo(uN);
                // for( int j = 0; j < evRbV.size(); ++j )
                // {
                //     errs[j+2][n][i] = std::abs(evV(j)-evRbV(j));
                //     errsRel[j+2][n][i] = errs[j+2][n][i]/std::abs(evV(j));
                // }
                // for( int j = 0; j < evRbT.size(); ++j )
                // {
                //     errs[j+2+evRbV.size()][n][i] = std::abs(evT(j)-evRbT(j));
                //     errsRel[j+2+evRbV.size()][n][i] = errs[j+2+evRbV.size()][n][i]/std::abs(evT(j));
                // }
                if( ioption("crb.output-index") != 0 )
                {
                    outRB = outputs[0];
                    errs[2+modelOutputs.size()][n][i] = std::abs(outFE-outRB);
                    errsRel[2+modelOutputs.size()][n][i] = errs[2+modelOutputs.size()][n][i]/std::abs(outFE);
                }
            }
            ++i;
        }

        std::vector<std::vector<double> > min(nbErr, std::vector<double>(N)), max(nbErr, std::vector<double>(N)), mean(nbErr, std::vector<double>(N));
        for(int n = 0; n < N; ++n)
        {
            for( int i = 0; i < nbErr; ++i )
            {
                min[i][n] = *std::min_element(errsRel[i][n].begin(), errsRel[i][n].end());
                max[i][n] = *std::max_element(errsRel[i][n].begin(), errsRel[i][n].end());
                double s = std::accumulate(errsRel[i][n].begin(), errsRel[i][n].end(), 0.0);
                mean[i][n] = s/size;
            }
        }

        std::vector<std::string> values({"V","T"});
        for( auto const& [name,o] : modelOutputs )
            values.push_back("P_"+name);
        if( ioption("crb.output-index") != 0 )
            values.push_back("O");
        for( int i = 0; i < values.size(); ++i )
        {
            Feel::cout << values[i] << std::endl;
            fs::ofstream cvgErr( "err"+values[i]+".dat" );
            fs::ofstream cvgErrR( "err"+values[i]+"R.dat" );
            fs::ofstream stat( "stat"+values[i]+".dat" );
            writeErrors(cvgErr, errs[i]);
            writeErrors(cvgErrR, errsRel[i]);
            cvgErr.close();
            cvgErrR.close();
            printStats(stat, min[i], max[i], mean[i]);
            stat.close();
        }

        auto mesh = model->mesh();
        auto e = exporter(mesh);
        e->add("VFE", VFE);
        e->add("VRB", VRB);
        e->add("TFE", TFE);
        e->add("TRB", TRB);
        e->save();
    }

    return 0;
}
