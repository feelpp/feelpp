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

#if defined THERMOELECTRIC_SADDLEPOINT
#include <feel/feelmor/crbsaddlepoint.hpp>
#include <feel/feelmor/crbmodelsaddlepoint.hpp>
#endif

#include "thermoelectric-nonlinear.hpp"
#include <feel/feelmor/ser.hpp>

using namespace Feel;

int main( int argc, char** argv)
{
    po::options_description opt("options");
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="thermoelectric-nonlinear"),
                     _desc=opt.add(makeThermoElectricOptions())
                     .add(crbOptions())
                     .add(crbSEROptions())
#if defined THERMOELECTRIC_SADDLEPOINT
                     .add(crbSaddlePointOptions())
#endif
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("ThermoElectricCRB")) );

    using rb_model_type = ThermoElectric;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;

#if defined THERMOELECTRIC_SADDLEPOINT
    using crb_model_type = CRBModelSaddlePoint<rb_model_type>;
#else
    using crb_model_type = CRBModel<rb_model_type>;
#endif
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;

#if defined THERMOELECTRIC_SADDLEPOINT
    using crb_type = CRBSaddlePoint<crb_model_type>;
#else
    using crb_type = CRB<crb_model_type>;
#endif
    using crb_ptrtype = std::shared_ptr<crb_type>;

    using ser_type = SER<crb_type>;
    using ser_ptrtype = std::shared_ptr<ser_type>;

    using wn_type = typename crb_type::wn_type;
    using vectorN_type = Eigen::VectorXd;
    using export_vector_wn_type = typename crb_type::export_vector_wn_type;
    using mesh_type = rb_model_type::mesh_type;
    using mesh_ptrtype = rb_model_type::mesh_ptrtype;

    auto mesh = loadMesh( _mesh=new mesh_type,
                          _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES);

    // init
    crb::stage stage = crb::stage::offline;
    rb_model_ptrtype model = std::make_shared<rb_model_type>(mesh);
    crb_model_ptrtype crbModel = std::make_shared<crb_model_type>(model, stage);
    crb_ptrtype crb = crb_type::New(model->modelName(), crbModel, stage );
    ser_ptrtype ser = std::make_shared<ser_type>( crb, crbModel );

    // offline
    ser->run();

    // online
    auto mu = model->paramFromProperties();
    int N = crb->dimension();
    int timeSteps = 1;
    int output_index = ioption("crb.output-index");

    auto e = exporter(mesh, _name="thermoelectric");

    auto VTFE = model->solve(mu);
    double truthOutput = 0;
    if( output_index != 0 )
        truthOutput = model->output( output_index, mu, VTFE, false);

    auto VFE = VTFE.template element<0>();
    auto TFE = VTFE.template element<1>();
    auto normV = normL2( VFE.functionSpace()->template rangeElements<0>(), idv(VFE) );
    auto normT = normL2( TFE.functionSpace()->template rangeElements<0>(), idv(TFE) );

    boost::format fmter;
    if( output_index == 0 )
        fmter = boost::format("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5%\n");
    else
        fmter = boost::format("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5% %|70t|%6%\n");
    std::string fileName = "cvg.dat";
    fs::ofstream file( fileName );
    if( file && Environment::isMasterRank() )
    {
        if( output_index == 0 )
        {
            file << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT";
            Feel::cout << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT";
        }
        else
        {
            file << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT" % "errOutput";
            Feel::cout << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT" % "errOutput";
        }
    }

    std::vector<vectorN_type> uNs(timeSteps), uNolds(timeSteps), uNdus(timeSteps), uNduolds(timeSteps);
    std::vector<double> outputs(timeSteps, 0);

    for(int n = 1; n <= N; ++n)
    {
        uNs[0].setZero(n);
        uNolds[0].setZero(n);
        uNdus[0].setZero(n);
        uNduolds[0].setZero(n);
        crb->fixedPoint(n, mu, uNs, uNdus, uNolds, uNduolds, outputs, 0, boption("crb.print-rb-matrix"));
        vectorN_type uN = uNs[0];
        double output = 0;
        if( output_index > 0)
            output = outputs[0];

        auto VTRB = crb->expansion( uN, n );
        auto VRB = VTRB.template element<0>();
        auto TRB = VTRB.template element<1>();

        auto errVRB = normL2( VFE.functionSpace()->template rangeElements<0>(), idv(VRB)-idv(VFE) );
        auto errRelV = errVRB/normV;
        auto errTRB = normL2( TFE.functionSpace()->template rangeElements<0>(), idv(TRB)-idv(TFE) );
        auto errRelT = errTRB/normT;
        double errOutput = 0;
        if( output_index > 0 )
            errOutput = std::abs(output-truthOutput);

        if( Environment::isMasterRank() )
        {
            if( output_index == 0 )
                file << fmter % n % errVRB % errRelV % errTRB % errRelT;
            else
                file << fmter % n % errVRB % errRelV % errTRB % errRelT % errOutput;
        }
        if( output_index == 0 )
            Feel::cout << fmter % n % errVRB % errRelV % errTRB % errRelT;
        else
            Feel::cout << fmter % n % errVRB % errRelV % errTRB % errRelT % errOutput;

        e->step(n)->add("TFE", TFE);
        e->step(n)->add("TRB", TRB);
        e->step(n)->add("VFE", VFE);
        e->step(n)->add("VRB", VRB);
        e->save();
    }
    file.close();

    return 0;
}
