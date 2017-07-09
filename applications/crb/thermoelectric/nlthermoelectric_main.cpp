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

#include "nlthermoelectric.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="nl-thermoelectric"),
                     _desc=makeOptions()
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("PoissonCRB")) );

    using rb_model_type = NLThermoelectric;
    using rb_model_ptrtype = boost::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = boost::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = boost::shared_ptr<crb_type>;
    using wn_type = typename crb_type::wn_type;
    using vectorN_type = Eigen::VectorXd;
    using export_vector_wn_type = typename crb_type::export_vector_wn_type;

    // init
    rb_model_ptrtype model = boost::make_shared<rb_model_type>();
    crb_model_ptrtype crbModel = boost::make_shared<crb_model_type>(model);
    crb_ptrtype crb = boost::make_shared<crb_type>("nlthermoelectric", crbModel);

    // offline
    crb->offline();

    // parameter from option
    double sigma = doption("thermoelectric.sigma");
    double alpha = doption("thermoelectric.alpha");
    double L = doption("thermoelectric.L");
    double current = doption("thermoelectric.current");
    double h = doption("thermoelectric.h");
    double Tw = doption("thermoelectric.Tw");
    auto paramSpace = crbModel->parameterSpace();
    auto mu = paramSpace->element();
    mu.setParameterNamed("sigma", sigma);
    mu.setParameterNamed("current", current);
    mu.setParameterNamed("alpha", alpha);
    mu.setParameterNamed("k", L);
    mu.setParameterNamed("h", h);
    mu.setParameterNamed("Tw", Tw);
    Feel::cout << "using parameter:" << std::endl << mu << std::endl;

    // online
    int N = ioption("thermoelectric.N");
    if( N > crb->dimension() || N < 1 )
        N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);
    crb->fixedPointPrimal(N, mu, uNs, uNolds, outputs);
    vectorN_type uN = uNs[0];
    Feel::cout << "uN:" << std::endl << uN << std::endl;

    // export basis
    wn_type WN = crb->wn();
    std::vector<wn_type> WN_vec = std::vector<wn_type>();
    WN_vec.push_back(WN);
    std::vector<std::string> name_vec = std::vector<std::string>(1, "primal");
    export_vector_wn_type exportWn = boost::make_tuple(WN_vec, name_vec);
    crb->exportBasisFunctions(exportWn);

    auto e = Exporter<rb_model_type::mesh_type>::New( "thermoelectric" );
    e->setMesh( model->mesh() );

    // export
    auto u = Feel::expansion( WN, uN, uN.size() );
    auto V = u.element<0>();
    auto T = u.element<1>();
    e->add( "V", V );
    e->add( "T", T );

    // export FE
    auto VTFE = model->solve(mu);
    auto VFE = VTFE.element<0>();
    auto TFE = VTFE.element<1>();
    e->add( "VFE", VFE );
    e->add( "TFE", TFE );
#if 0
    // export EIM basis
    auto eimSigma = model->scalarContinuousEim()[0];
    auto spaceSigma = eimSigma->functionSpace();
    auto qSigma = eimSigma->q();
    auto betaSigma = eimSigma->beta(mu, VTFE);
    auto sigmaEim = spaceSigma->element();
    for( int m = 0; m < qSigma.size(); ++m )
    {
        e->add( (boost::format("qSigma_%1%" ) % m).str(), qSigma[m] );
        sigmaEim += vf::project( _range=elements(spaceSigma->mesh()),
                                 _space=spaceSigma,
                                 _expr=betaSigma[m]*idv(qSigma[m]) );
    }
    e->add( "sigma", sigmaEim );

    // auto eimJoule = model->scalarDiscontinuousEim()[0];
    // auto spaceJoule = eimJoule->functionSpace();
    // auto qJoule = eimJoule->q();
    // auto betaJoule = eimJoule->beta(mu, VTFE);
    // auto joule = spaceJoule->element();
    // for( int m = 0; m < qJoule.size(); ++m )
    // {
    //     e->add( (boost::format("qJoule_%1%" ) % m).str(), qJoule[m] );
    //     joule += vf::project( elements(spaceJoule->mesh()), spaceJoule, betaJoule[m]*idv(qJoule[m]));
    // }
    // e->add( "joule", joule );

    auto eimJoule = model->scalarDiscontinuousEim()[0];
    auto spaceJoule = eimJoule->functionSpace();
    auto qJoule = eimJoule->q();
    auto betaJoule = eimJoule->beta(mu, VTFE);
    auto joule = spaceJoule->element();
    for( int m = 0; m < qJoule.size(); ++m )
    {
        e->add( (boost::format("qJoule_%1%" ) % m).str(), qJoule[m] );
        joule += vf::project( _range=elements(spaceJoule->mesh()),
                              _space=spaceJoule,
                              _expr=betaJoule[m]*idv(qJoule[m]) );
    }
    e->add( "joule", joule );

    // auto j = vf::project( _range=elements(model->mesh()), _space=eimG->functionSpace(),
    //                       _expr=mu.parameterNamed("sigma")*inner(gradv(V)) );
    // e->add( "j2", j );
#endif
    e->save();

    return 0;
}
