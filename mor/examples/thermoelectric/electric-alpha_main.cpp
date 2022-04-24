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

#include "electric-alpha.hpp"

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
        out.close();
    }
}

int main( int argc, char** argv)
{
    po::options_description opt("options");
    opt.add_options()
        ( "alphaelectric.size", po::value<int>()->default_value(1), "size of convergence sampling" )
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt.add(Feel::AlphaElectric::makeOptions("electric")));

    using rb_model_type = AlphaElectric;
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
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;

    auto mesh = loadMesh( _mesh=new mesh_type,
                          _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES);

    // init
    rb_model_ptrtype model = std::make_shared<rb_model_type>(mesh, "electric");
    crb_model_ptrtype crbModel = std::make_shared<crb_model_type>(model);
    crb_ptrtype crb = crb_type::New(model->dbBasename(), crbModel, crb::stage::offline);

    // offline
    crb->offline();

    // online
    sampling_ptrtype sampling( new sampling_type( crbModel->parameterSpace() ) );
    int size = ioption("alphaelectric.size");
    sampling->clear();
    sampling->randomize( size, true );

    auto mu = model->paramFromProperties();
    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);

    auto Xh = model->functionSpace();
    auto rangeV = Xh->dof()->meshSupport()->rangeElements();
    auto VFE = Xh->element();
    auto VRB = Xh->element();

    std::vector<std::vector<double> > errs(N, std::vector<double>(size));
    std::vector<std::vector<double> > errsRel(N, std::vector<double>(size));
    Feel::cout << "start convergence study with " << size << " parameters" << std::endl;
    int i = 0;
    for( auto const& mu : *sampling )
    {
        VFE = model->solve(mu);
        auto normV = normL2( rangeV, idv(VFE) );
        for(int n = 0; n < N; ++n)
        {
            crb->fixedPointPrimal(n+1, mu, uNs, uNolds, outputs);
            vectorN_type uN = uNs[0];
            VRB = crb->expansion( uN, n+1 );
            errs[n][i] = normL2( rangeV, idv(VRB)-idv(VFE) );
            errsRel[n][i] = errs[n][i]/normV;
        }
        ++i;
    }

    std::vector<double> min(N), max(N), mean(N), stdev(N);
    for(int n = 0; n < N; ++n)
    {
        min[n] = *std::min_element(errsRel[n].begin(), errsRel[n].end());
        max[n] = *std::max_element(errsRel[n].begin(), errsRel[n].end());
        double s = std::accumulate(errsRel[n].begin(), errsRel[n].end(), 0.0);
        mean[n] = s/size;
        double accum = std::accumulate(errsRel[n].begin(), errsRel[n].end(), 0.0,
                                       [s,size](double a, double b) {
                                           return a + (b-s/size)*(b-s/size);
                                       });
        stdev[n] = accum/size;
    }

    fs::ofstream cvgErr( "err.dat" );
    fs::ofstream cvgErrR( "errR.dat" );
    fs::ofstream cvgStat( "stat.dat" );
    writeErrors(cvgErr, errs);
    writeErrors(cvgErrR, errsRel);
    if( cvgStat && Environment::isMasterRank() )
    {
        cvgStat << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
                << std::setw(25) << "mean" << std::setw(25) << "stdev" << std::endl;
        for(int n = 0; n < N; ++n)
            cvgStat << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                    << std::setw(25) << mean[n] << std::setw(25) << stdev[n] << std::endl;
        cvgStat.close();
    }
    Feel::cout << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
               << std::setw(25) << "mean" << std::setw(25) << "stdev" << std::endl;
    for(int n = 0; n < N; ++n)
        Feel::cout << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                   << std::setw(25) << mean[n] << std::setw(25) << stdev[n] << std::endl;

    auto e = exporter(mesh);
    e->add("VFE", VFE);
    e->add("VRB", VRB);
    e->save();

    return 0;
}
