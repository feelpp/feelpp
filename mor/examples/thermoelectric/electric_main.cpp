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

#include "electric.hpp"

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

    using rb_model_type = Electric;
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
    crb_ptrtype crb = crb_type::New("electric", crbModel, crb::stage::offline);

    // offline
    crb->offline();

    // online
    auto mu = model->paramFromProperties();
    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);

    auto e = exporter(mesh, _name="electric");

    auto VFE = model->solve(mu);
    auto normV = normL2( elements(model->mesh()), idv(VFE) );
    boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
    std::ofstream file( "cvg.dat" );
    if( file && Environment::isMasterRank() )
    {
        file << fmter % "N" % "errV" % "relErrV";
        Feel::cout << fmter % "N" % "errV" % "relErrV";
    }
    for(int n = 1; n <= N; ++n)
    {
        crb->fixedPointPrimal(n, mu, uNs, uNolds, outputs);
        vectorN_type uN = uNs[0];
        auto VRB = crb->expansion( uN, n );
        auto errVRB = normL2( elements(model->mesh()), idv(VRB)-idv(VFE) );
        auto errRel = errVRB/normV;
        if( Environment::isMasterRank() )
            file << fmter % n % errVRB % errRel;
        Feel::cout << fmter % n % errVRB % errRel;
        e->step(n)->add("VFE", VFE);
        e->step(n)->add("VRB", VRB);
        e->save();
    }
    file.close();

    return 0;
}
