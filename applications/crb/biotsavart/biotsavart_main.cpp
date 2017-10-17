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

#include "biotsavart.hpp"
#include "thermoelectric-linear.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=biotsavartOptions()
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("ThermoElectricCRB"))
                     .add(makeThermoElectricOptions()) );

    BiotSavartCRB<ThermoElectric> bs = BiotSavartCRB<ThermoElectric>();
    //bs.runBS();

    auto mu = bs.paramFromOption();

    bs.offline();

    auto eC = Exporter<Mesh<Simplex<3> > >::New("conductor");
    eC->setMesh(bs.meshCond());
    auto eM = Exporter<Mesh<Simplex<3> > >::New("magneto");
    eM->setMesh(bs.meshMgn());

    bs.computeFE(mu);
    auto VTFe = bs.potentialTemperatureFE();
    auto VFe = VTFe.template element<0>();
    auto normV = normL2( elements(bs.meshCond()), idv(VFe) );
    auto TFe = VTFe.template element<1>();
    auto normT = normL2( elements(bs.meshCond()), idv(TFe) );
    auto BFe = bs.magneticFluxFE();
    auto normB = normL2( elements(bs.meshMgn()), idv(BFe) );

    int N = bs.nMax();
    std::vector<double> errV(N-1), errT(N-1), errB(N-1), relErrV(N-1), relErrT(N-1), relErrB(N-1);
    for( int n = 1; n < N; ++n )
    {
        bs.computeUn(mu, n);
        auto uN = bs.uN();
        bs.computeB(uN);
        auto VT = bs.potentialTemperature();
        auto V = VT.template element<0>();
        auto T = VT.template element<1>();
        auto B = bs.magneticFlux();
        errV[n-1] = normL2( elements(bs.meshCond()), idv(VFe)-idv(V) );
        relErrV[n-1] = errV[n-1]/normV;
        errT[n-1] = normL2( elements(bs.meshCond()), idv(TFe)-idv(T) );
        relErrT[n-1] = errT[n-1]/normT;
        errB[n-1] = normL2( elements(bs.meshMgn()), idv(BFe)-idv(B) );
        relErrB[n-1] = errB[n-1]/normB;

        eC->step(n)->add("VFE", VFe);
        eC->step(n)->add("TFE", TFe);
        eC->step(n)->add("V", V);
        eC->step(n)->add("T", T);
        eM->step(n)->add("BFe", BFe);
        eM->step(n)->add("B", B);
    }

    eC->save();
    eM->save();

    boost::format fmter("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5% %|70t|%6% %|84t|%7%\n");
    fs::ofstream file( "cvg.dat" );
    if( file && Environment::isMasterRank() )
    {
        file << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT" % "errB" % "relErrB";
        Feel::cout << fmter % "N" % "errV" % "relErrV" % "errT" % "relErrT" % "errB" % "relErrB";
        for( int n = 0; n < N-1; ++n )
        {
            file << fmter % (n+1) % errV[n] % relErrV[n] % errT[n] % relErrT[n] % errB[n] % relErrB[n];
            Feel::cout << fmter % (n+1) % errV[n] % relErrV[n] % errT[n] % relErrT[n] % errB[n] % relErrB[n];
        }
    }


    return 0;
}
