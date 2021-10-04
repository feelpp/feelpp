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
#include "biotsavart-maxwell.hpp"
#include <feel/feelcrb/ser.hpp>

#define DO_CVG_TEST 0
#define USE_BS 0

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

template<int Dim, int Order, int G_Order>
int runApplicationThermoElectricNL()
{
    using rb_model_type = ThermoElectricNL<Dim, Order, G_Order>;
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

    auto crb = crb_type::New(soption("thermoelectric.basename"), crb::stage::offline);
    auto crbmodel = crb->model();
    auto model = crb->model()->model();
    auto ser = std::make_shared<ser_type>(crb, crbmodel);
    if( boption("crb.rebuild-database") || boption("eim.rebuild-database") )
        ser->run();
    else
        crb->offline();
    // crb->loadDBFromId( soption("online.dbid"), crb::load::all );
    model->initOutputsPoints();

    int N = crb->dimension();
    auto wn = crb->wn();
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

#if USE_BS
    auto sigmaEim = model->scalarContinuousEim()[1];
    int M = sigmaEim->mMax();
    auto center = vdoption("online.biotsavart.center");
    GeoTool::Node c(center[0],center[1],center[2]);
    auto radius = doption("online.biotsavart.radius");
    auto BS = BiotSavart<Dim>(model->mesh(), c, radius);
    BS.init();
    std::vector<decltype(BS)::magneticfield_element_type> bn;
    Feel::cout << "computing " << N << " times " << M << " reduced basis for BS" << std::endl;
    for(int i = 0; i< N; ++i )
    {
        for(int m = 0; m < M; ++m )
        {
            BS.compute(idv(sigmaEim->q(m))*trans(gradv(wn[i]->template element<0>())), true, false );
            bn.push_back(BS.magneticField());
        }
    }

    auto Bh = BS.mgnSpace();
    auto rangeB = elements(Bh->mesh());
    auto BFE = Bh->element();
    auto BRB = Bh->element();

#endif

#if DO_CVG_TEST
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

        auto modelOutputs = model->modelProperties()->outputs();
        auto modelOutputsPoints = modelOutputs.outputsOfType("point");
        auto ctxFeV = Xh->template functionSpace<0>()->context();
        auto ctxFeT = Xh->template functionSpace<1>()->context();
        for( auto const& [name,output] : modelOutputsPoints )
        {
            node_type t(Dim);
            auto coord = expr<Dim,1>(output.getString("coord")).evaluate();
            for(int i = 0; i < Dim; ++i )
                t(i) = coord(i);
            if( output.getString("field") == "electric-potential" )
                ctxFeV.add( t );
            else if( output.getString("field") == "temperature" )
                ctxFeT.add( t );
        }

        int nbErr = 2 + modelOutputsPoints.size() + (ioption("crb.output-index") != 0 ? 1 : 0);
#if USE_BS
        nbErr += 1;
#endif
        std::vector<std::vector<std::vector<double> > > errs(nbErr, std::vector<std::vector<double> >(N, std::vector<double>(size)));
        std::vector<std::vector<std::vector<double> > > errsRel(nbErr, std::vector<std::vector<double> >(N, std::vector<double>(size)));
        Feel::cout << "start convergence study with " << size << " parameters" << std::endl;
        int i = 0;
        for( auto const& mu : *sampling )
        {
            Feel::cout << i << ": solving for mu = " << mu.toString() << std::endl;
            VTFE = model->solve(mu);
            auto evV = evaluateFromContext( _context=ctxFeV, _expr=idv(VFE));
            auto evT = evaluateFromContext( _context=ctxFeT, _expr=idv(TFE));
            if( ioption("crb.output-index") != 0 )
            {
                int j = 1;
                for( auto const& [key, output] : modelOutputs )
                {
                    if( output.type() == "averageTemp" )
                    {
                        if( j == ioption("crb.output-index") )
                        {
                            if( output.dim() == Dim )
                                outFE = mean(_range=markedelements(model->mesh(), output.markers()),
                                             _expr=idv(TFE) )(0,0);
                            else if( output.dim() == Dim-1)
                                outFE = mean(_range=markedfaces(model->mesh(), output.markers()),
                                             _expr=idv(TFE) )(0,0);
                            break;
                        }
                        ++j;
                    }
                    else if( output.type() == "intensity" )
                    {
                        if( j == ioption("crb.output-index") )
                        {
                            auto mat = model->elecMaterials().at(output.getString("material"));
                            auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
                            auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
                            auto sigma = sigma0/(cst(1.)+alpha*(idv(TFE)-cst(293.)));
                            outFE = integrate( _range=markedfaces(model->mesh(), output.markers()),
                                               _expr=-sigma*gradv(VFE)*vf::N()).evaluate()(0,0);
                            break;
                        }
                        ++j;
                    }
                }
            }
            auto normV = normL2( rangeU, idv(VFE) );
            auto normT = normL2( rangeU, idv(TFE) );

#if USE_BS
            auto sigma0 = mu.parameterNamed("sigma");
            auto alpha = mu.parameterNamed("alpha");
            auto sigma = sigma0/(cst(1.)+alpha*(idv(TFE)-cst(293.)));
            BS.compute(sigma*trans(gradv(VFE)),true, false);
            BFE = BS.magneticField();
            auto normB = normL2( _range=elements(BS.mesh()), _expr=idv(BFE));
#endif
            for(int n = 0; n < N; ++n)
            {
                int j = 0;
                crb->fixedPointPrimal(n+1, mu, uNs, uNolds, outputs);
                vectorN_type uN = uNs[0];
                VTRB = crb->expansion( uN, n+1 );
                errs[j][n][i] = normL2( rangeU, idv(VRB)-idv(VFE) );
                errsRel[j][n][i] = errs[j][n][i]/normV;
                j++;
                errs[j][n][i] = normL2( rangeU, idv(TRB)-idv(TFE) );
                errsRel[j][n][i] = errs[j][n][i]/normT;
                j++;
#if USE_BS
                auto sigmaE = sigmaEim->beta(mu, uN);
                std::vector<double> betaB;
                for(int ii = 0; ii < n+1; ++ii )
                    for(int m = 0; m < M; ++m )
                        betaB.push_back(50000*sigmaE[m]*uN[ii]);
                BRB = Feel::expansion(bn, betaB, betaB.size());
                errs[j][n][i] = normL2( _range=elements(BS.mesh()), _expr=idv(BRB)-idv(BFE) );
                errsRel[j][n][i] = errs[j][n][i]/normB;
                j++;
#endif

                auto evRbV = model->computeOutputsPointsElectro(uN);
                auto evRbT = model->computeOutputsPointsThermo(uN);
                for( int k = 0; k < evRbV.size(); ++k,++j )
                {
                    errs[j][n][i] = std::abs(evV(k)-evRbV(k));
                    errsRel[j][n][i] = std::abs(evV(k)) < 1e-7 ? errs[j][n][i] : errs[j][n][i]/std::abs(evV(k));
                }
                for( int k = 0; k < evRbT.size(); ++k,++j )
                {
                    errs[j][n][i] = std::abs(evT(k)-evRbT(k));
                    errsRel[j][n][i] = std::abs(evT(k)) < 1e-7 ? errs[j][n][i] : errs[j][n][i]/std::abs(evT(k));
                }
                if( ioption("crb.output-index") != 0 )
                {
                    outRB = outputs[0];
                    errs[j][n][i] = std::abs(outFE-outRB);
                    errsRel[j][n][i] = std::abs(outFE) < 1e-7 ? errs[j][n][i] : errs[j][n][i]/std::abs(outFE);
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
#if USE_BS
        values.push_back("B");
#endif
        for( auto const& [name,o] : modelOutputsPoints )
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

#if USE_BS
        auto e2 = exporter(_mesh=BS.mesh(),_name="bs");
        e2->add("BFE",BFE);
        e2->add("BRB",BRB);
        e2->save();
#endif
        double maxErr = 0;
        for( int i = 0; i < nbErr; ++i )
            if( mean[i][N-1] > maxErr )
                maxErr = mean[i][N-1];
        return maxErr > 1e-4;
    }
    return 0;
}


int main( int argc, char** argv)
{
    auto opt = ThermoElectricNLBase::makeOptions();
    opt.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1 ")
        ( "online.dbid", po::value<std::string>()->default_value(""), "id of the db to load" )
        ( "online.sampling-size", po::value<int>()->default_value(1), "number of parameters" )
        ( "online.do-cvg", po::value<bool>()->default_value(true), "" )
        ( "online.biotsavart.center", po::value<std::vector<double>>()->default_value({{0,0,0}}),"")
        ( "online.biotsavart.radius", po::value<double>()->default_value(1), "")
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(biotsavart_options())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("PoissonCRB")) );

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1>, hana::int_c<1> ));
    int status = 1;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)),
                    [&discretization,&dimension,&status]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            constexpr int _gorder = std::decay_t<decltype(hana::at_c<2>( hana::at_c<1>(d) ))>::value;
                            if ( dimension == _dim && discretization == _discretization )
                                status = runApplicationThermoElectricNL<_dim,_torder,_gorder>();
                        } );

    return status;
}
