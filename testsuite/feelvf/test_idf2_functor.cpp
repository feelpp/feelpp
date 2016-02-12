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
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/function.hpp>



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

                B[q](0) = r*Cos_Theta;
                B[q](1) = r*Sin_Theta;
                B[q](2) = z;
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
                     _about=about(_name="myfunctor",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

    auto Vh= Pchv<3>(mesh);
    auto b = Vh->element();
    b_ana bfield;
    b.on( _range=elements(mesh), _expr=idf2(bfield));

}
