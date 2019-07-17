/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Christophe Trophime

   Date 2016-02-18

   Copyright (C) 2016 Feel++ Consortium

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   */

#include <vector>
#include <Eigen/StdVector>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/function.hpp>
#include <feel/feelvf/minmax.hpp>
#include <feel/feelfilters/unitsquare.hpp>

namespace Feel
{
struct b_ana
{
    static const size_type context = vm::JACOBIAN|vm::POINT|vm::NORMAL;
    typedef double value_type;
    typedef value_type evaluate_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 1;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    using B_t = Eigen::Matrix<value_type,Eigen::Dynamic,1>;

    b_ana()
        {}
    template<typename Geo_t>
    void init( Geo_t g )
        {
            B.resize(g->nPoints());
            std::for_each( B.begin(), B.end(), [&g]( auto& m ) { m = B_t::Zero( g->N() );  });
        }
    template<typename Geo_t>
    void update( Geo_t g )
        {
            for(int q = 0; q < g->nPoints();++q )
            {
                auto const& x=g->xReal( q );
                double r = math::sqrt(x[0]*x[0] + x[1]*x[1]);
                double z = x[2];
                double theta = math::atan2(x[1], x[0]);

                auto Cos_Theta = std::cos(theta);
                auto Sin_Theta = std::sin(theta);

                B[q](0) = -Sin_Theta/r;
                B[q](1) = Cos_Theta/r;
                B[q](2) = 0;
            }

        }
    template<typename Geo_t>
    void update( Geo_t g, uint16_type f )
        {}

    double evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return B[q](c1);
        }

    B_t const& evalq( uint16_type q ) const
        {
            return B[q]; //.row(q);
        }

    std::vector<B_t> B;
};

} // Feel


int main(int argc, char**argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="divergence",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new  Mesh<Simplex<3>>); //unitHypercube()

    auto L2 = Pdh<0>(mesh);
    auto H1= Pchv<1>(mesh);
    auto Hcurl = Ned1h<0>(mesh);
    auto Hdiv = Dh<0>(mesh);

    auto b_H1 = H1->element();
    auto b_Hdiv = Hdiv->element();
    auto b_Hcurl = Hcurl->element();

    b_ana bfield;
    b_H1.on( _range=elements(mesh), _expr=idf2(bfield));
    b_Hcurl.on( _range=elements(mesh), _expr=idf2(bfield));
    b_Hdiv.on( _range=elements(mesh), _expr=idf2(bfield));

    auto div_H1 = L2->element();
    auto div_Hdiv = L2->element();
    auto div_Hcurl = L2->element();
    div_H1.on( _range=elements(mesh), _expr=divv(b_H1));
    div_Hcurl.on( _range=elements(mesh), _expr=divv(b_Hcurl));
    div_Hdiv.on( _range=elements(mesh), _expr=divv(b_Hdiv));

    auto m = minmax( _range=elements(mesh), _pset=_Q<2>(), _expr=idv(div_H1));
    auto mmin = m.min();
    auto mmax = m.max();
    auto p_mmin = m.argmin();
    auto p_mmax = m.argmax();

    auto n = minmax( _range=elements(mesh), _pset=_Q<2>(), _expr=idv(div_Hcurl));
    auto nmin = n.min();
    auto nmax = n.max();
    auto p_nmin = n.argmin();
    auto p_nmax = n.argmax();

    auto p = minmax( _range=elements(mesh), _pset=_Q<2>(), _expr=idv(div_Hdiv));
    auto pmin = p.min();
    auto pmax = p.max();
    auto p_pmin = p.argmin();
    auto p_pmax = p.argmax();

    Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "{", "}");
    Feel::cout << "******************************** \n";
    Feel::cout << "Space   |\tMin |\tMax :\n";
    Feel::cout << "******************************** \n";
    Feel::cout << "H1      |\t";
    Feel::cout << mmin << p_mmin.format(CleanFmt) << "|\t";
    Feel::cout << mmax << p_mmax.format(CleanFmt) << "\n";
    Feel::cout << "H(curl) |\t";
    Feel::cout << mmin << p_nmin.format(CleanFmt) << "|\t";
    Feel::cout << mmax << p_nmax.format(CleanFmt) << "\n";
    Feel::cout << "H(div) |\t";
    Feel::cout << mmin << p_nmin.format(CleanFmt) << "|\t";
    Feel::cout << mmax << p_nmax.format(CleanFmt) << "\n";
    Feel::cout << "******************************** \n";

    auto postproc = exporter(_mesh=mesh);
    postproc->step(0)->add( "b_H1", b_H1 );
    postproc->step(0)->add( "div_H1", div_H1 );
    postproc->step(0)->add( "b_Hcurl", b_Hdiv );
    postproc->step(0)->add( "div_Hcurl", div_Hcurl );
    postproc->step(0)->add( "b_Hdiv", b_Hcurl );
    postproc->step(0)->add( "div_Hdiv", div_Hdiv );
    postproc->save();
}
