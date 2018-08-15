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

#include <feel/feelcrb/crbsaddlepoint.hpp>
#include <feel/feelcrb/crbmodelsaddlepoint.hpp>
#include "thermoelectric-simple.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    po::options_description opt("options");
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt.add(makeThermoElectricOptions())
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(crbSaddlePointOptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("ThermoElectricCRB")) );

    using rb_model_type = ThermoElectric;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModelSaddlePoint<rb_model_type>;
    // using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRBSaddlePoint<crb_model_type>;
    // using crb_type = CRB<crb_model_type>;
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
    crb_model_ptrtype crbModel = std::make_shared<crb_model_type>(model, crb::stage::offline);
    crb_ptrtype crb = crb_type::New("thermoelectric", crbModel, crb::stage::offline);

    // offline
    crb->offline();

    // online
    auto mu = model->paramFromProperties();
    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);

    auto e = exporter(mesh, _name="thermoelectric");

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
    auto VTFE = model->solve(mu);
    auto VFE = VTFE.template element<0>();
    auto TFE = VTFE.template element<1>();
    auto normV = normL2( elements(model->mesh()), idv(VFE) );
    auto normT = normL2( elements(model->mesh()), idv(TFE) );
    boost::format fmter("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5%\n");
    fs::ofstream file( "cvg.dat" );
    if( file && Environment::isMasterRank() )
    {
        file << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT";
        Feel::cout << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT";
    }
    for(int n = 1; n <= N; ++n)
    {
        crb->fixedPointPrimal(n, mu, uNs, uNolds, outputs);
        vectorN_type uN = uNs[0];
        auto VTRB = crb->expansion( uN, n );
        auto VRB = VTRB.template element<0>();
        auto TRB = VTRB.template element<1>();
        auto errVRB = normL2( elements(model->mesh()), idv(VRB)-idv(VFE) );
        auto errRelV = errVRB/normV;
        auto errTRB = normL2( elements(model->mesh()), idv(TRB)-idv(TFE) );
        auto errRelT = errTRB/normT;
        if( Environment::isMasterRank() )
            file << fmter % n % errVRB % errRelV % errTRB % errRelT;
        Feel::cout << fmter % n % errVRB % errRelV % errTRB % errRelT;
        e->step(n)->add("TFE", TFE);
        e->step(n)->add("TRB", TRB);
        e->step(n)->add("VFE", VFE);
        e->step(n)->add("VRB", VRB);
        e->save();
    }
    file.close();

    return 0;
}
