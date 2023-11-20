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

#include "thermic.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    po::options_description opt("options");
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt.add(makeThermoElectricOptions())
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("ThermoElectricCRB")) );

    using rb_model_type = Thermic;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using wn_type = typename crb_type::wn_type;
    using vectorN_type = Eigen::VectorXd;
    using export_vector_wn_type = typename crb_type::export_vector_wn_type;
    using mesh_type = rb_model_type::mesh_type;
    using mesh_ptrtype = rb_model_type::mesh_ptrtype;

    auto mesh = loadMesh( _mesh=new mesh_type,
                          _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES);

    // init
    rb_model_ptrtype model = std::make_shared<rb_model_type>(mesh);
    crb_model_ptrtype crbModel = std::make_shared<crb_model_type>(model);
    crb_ptrtype crb = crb_type::New("thermic", crbModel, crb::stage::offline);

    // offline
    crb->offline();

    // online
    auto mu = model->paramFromProperties();
    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);

    auto e = exporter(mesh, _name="thermic");

    // auto T = crb->expansion(uN,N);
    // e->add("T",T);
    // for(auto i =0; i < N; i++)
    // {
    //     auto rb=crbModel->rBFunctionSpace()->primalBasisElement(i);
    //     auto norm = normL2( elements(model->mesh()), idv(rb));
    //     Feel::cout << "||rb_" << i << "|| = " << norm << std::endl;
    //     e->add((boost::format("rb_%1%") % i ).str(), rb);
    // }
    // e->save();
    auto TFE = model->solve(mu);
    auto normT = normL2( elements(model->mesh()), idv(TFE) );
    boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
    std::ofstream file( "cvg.dat" );
    if( file && Environment::isMasterRank() )
    {
        file << fmter % "N" % "errT" % "relErrT";
        Feel::cout << fmter % "N" % "errT" % "relErrT";
    }
    for(int n = 1; n <= N; ++n)
    {
        crb->fixedPointPrimal(n, mu, uNs, uNolds, outputs);
        vectorN_type uN = uNs[0];
        auto TRB = crb->expansion( uN, n );
        auto errTRB = normL2( elements(model->mesh()), idv(TRB)-idv(TFE) );
        auto errRel = errTRB/normT;
        if( Environment::isMasterRank() )
            file << fmter % n % errTRB % errRel;
        Feel::cout << fmter % n % errTRB % errRel;
        e->step(n)->add("TFE", TFE);
        e->step(n)->add("TRB", TRB);
        e->save();
    }
    file.close();

    return 0;
}
