/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
   -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
/**
   \file test_projector.cpp
   \author Vincent Doyeux <vdoyeux at gmail.com>
   \date 2015-06-12
 */

#define BOOST_TEST_MODULE projector
#include <feel/feelcore/testsuite.hpp>

#if !defined(USE_BOOST_TEST)
#include <feel/feelfilters/exporter.hpp>
#endif

#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

/* ProjectorType are given in enum.hpp:
enum ProjectorType
{
    NODAL=-1,
    L2=0,
    H1=1,
    DIFF=2,
    HDIV=3,
    HCURL=4,
    LIFT=5,
    CIP=6
};
 */

/// Test the projectors
///     \tparam DIM         Topological dimension.
///     \tparam H_ORDER     Space polynomial order..
template< int DIM, int H_ORDER >
class TestProjector
{

public:

    // Init the test, choose the projection, the mesh size and the diffusion coeff if needed
    ///     \param projType     The projection type to be tested.
    ///     \param projName     A name used for the exports.
    ///     \param epsilon      The smoothing parameter used in DIFF projection.
    TestProjector( ProjectorType projType=L2, std::string projName="L2", double epsilon=0.001 )
        {
            auto mesh = unitHypercube<DIM>();
            auto Xh = Pch<H_ORDER>(mesh);
            const double h = doption(_name="gmsh.hsize");

            auto f_expr = sin(2*M_PI*Px())*cos(2*M_PI*Py())*cos(2*M_PI*Pz());

            LOG(INFO) << "making projector "<< projName
                    << "\ntype =  " << projType
                    << "\ndim = "<< DIM
                    << "\norder = "<<H_ORDER<<"\n";

            auto proj = projector( Xh, Xh,
                                   Backend<double>::build( soption( _name="backend" ) ),
                                   projType,
                                   epsilon );

            LOG(INFO) << "projecting\n";
            auto f_proj = proj->project( _expr=f_expr );

            const double errL2 = normL2(_range=elements(mesh), _expr=idv(f_proj)-f_expr );
            LOG(INFO) << "projection error: "<< errL2 <<"\n";

#if defined (USE_BOOST_TEST)
            /*A more accurate convergence rate could be find, but this should be enough
             to detect bugs. For better convergence tests with projectors,
             see benchmarks/perf/curvature  */
            const double error_order = (projType==L2) ? std::pow(h, H_ORDER) : h;
            BOOST_CHECK_SMALL(errL2, error_order);
#else
            auto f_nod = vf::project(_space=Xh, _range=elements(mesh), _expr=f_expr);
            auto f_error = vf::project(_space=Xh, _range=elements(mesh), _expr=idv(f_nod)-idv(f_proj));

            const std::string suffix = projName + "_" + std::to_string(DIM) + "d_o" + std::to_string(H_ORDER);
            auto exp = exporter(_mesh=mesh, _name="projector_"+suffix);
            exp->step(0)->add("f_"+suffix, f_proj);
            exp->step(0)->add("fnod_"+suffix, f_nod);
            exp->step(0)->add("error_"+suffix, f_error);
            exp->save();
#endif
        }

};//class TestProjector


#if defined(USE_BOOST_TEST)

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( projector )
const std::vector< std::pair<ProjectorType, std::string> > projToTest = { {L2, "L2"},
                                                                    {DIFF, "DIFF"} };

BOOST_AUTO_TEST_CASE(test_proj_2d_o1)
{
    for (auto const& projtt : projToTest)
        TestProjector<2, 1> tp(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_CASE(test_proj_2d_o2)
{
    for (auto const& projtt : projToTest)
        TestProjector<2, 2> tp(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_CASE(test_proj_3d_o1)
{
    for (auto const& projtt : projToTest)
        TestProjector<3, 1> tp(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_CASE(test_proj_3d_o2)
{
    for (auto const& projtt : projToTest)
        TestProjector<3, 2> tp(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_SUITE_END()

#else

int main(int argc, char** argv)
{
    Feel::Environment env( _argc=argc, _argv=argv );

    const std::vector< std::pair<ProjectorType, std::string> > projToTest = { {L2, "L2"},
                                                                        {DIFF, "DIFF"} };

    for (auto const& projtt : projToTest )
    {
        // dim 2
        TestProjector<2, 1> tp2do1(projtt.first, projtt.second);
        TestProjector<2, 2> tp2do2(projtt.first, projtt.second);

        // dim  3
        TestProjector<3, 1> tp3do1(projtt.first, projtt.second);
        TestProjector<3, 2> tp3do2(projtt.first, projtt.second);
    }

    return 0;
}
#endif
