//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
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
//! @author Vincent Doyeux <vdoyeux at gmail.com>
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 06 Dec 2015
//! @copyright 2015-2017 Feel++ Consortium
//!
#define BOOST_TEST_MODULE projector
#include <testsuite/testsuite.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/projector.hpp>

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
// Init the test, choose the projection, the mesh size and the diffusion coeff if needed
///     \param projType     The projection type to be tested.
///     \param projName     A name used for the exports.
///     \param epsilon      The smoothing parameter used in DIFF projection.

template< int DIM, int H_ORDER >
void TestProjector( ProjectorType projType=L2, std::string projName="L2", double epsilon=0.001 )
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


template< int DIM, int H_ORDER >
void TestProjectorMatrix(double epsilon=0.001 )
{
    auto mesh = unitHypercube<DIM>();
    auto Mh = Pchms<H_ORDER>(mesh);
    const double h = doption(_name="gmsh.hsize");

    auto f_expr = sin(2*M_PI*Px())*cos(2*M_PI*Py())*cos(2*M_PI*Pz());
    auto fm_expr= f_expr*eye<DIM,DIM>( );
    LOG(INFO) << "making projector "<< "L2"
              << "\ntype =  "  << L2
              << "\ndim = "<< DIM
              << "\norder = "<<H_ORDER<<"\n";

    auto projm = projector( Mh, Mh,
                            Backend<double>::build( soption( _name="backend" ) ),
                            L2,
                            epsilon );
    return ;
    LOG(INFO) << "projecting\n";
    auto fm_proj = projm->project( _expr=fm_expr );
    const double errL2m = normL2(_range=elements(mesh), _expr=idv(fm_proj)-fm_expr );
    LOG(INFO) << "projection error matrix: "<< errL2m <<"\n";

#if defined (USE_BOOST_TEST)
    /*A more accurate convergence rate could be find, but this should be enough
     to detect bugs. For better convergence tests with projectors,
     see benchmarks/perf/curvature  */
    const double error_order = std::pow(h, H_ORDER);
    BOOST_CHECK_SMALL(errL2m, error_order);
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

inline
po::options_description
makeOptions()
{
    po::options_description o( "FESpace context options" );
    return o;
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_projector" ,
                     "test_projector" ,
                     "0.2",
                     "nD(n=1,2,3) test context of functionspace",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( projector )
const std::vector< std::pair<ProjectorType, std::string> > projToTest = { {L2, "L2"},
                                                                    {DIFF, "DIFF"} };


BOOST_AUTO_TEST_CASE(test_proj_2d_o1)
{
    for (auto const& projtt : projToTest)
        TestProjector<2, 1>(projtt.first, projtt.second);
}
BOOST_AUTO_TEST_CASE(test_proj_2d_o1m)
{
    TestProjectorMatrix<2, 1>();
}


BOOST_AUTO_TEST_CASE(test_proj_2d_o2)
{
    for (auto const& projtt : projToTest)
        TestProjector<2, 2>(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_CASE(test_proj_3d_o1)
{
    for (auto const& projtt : projToTest)
        TestProjector<3, 1>(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_CASE(test_proj_3d_o2)
{
    for (auto const& projtt : projToTest)
        TestProjector<3, 2>(projtt.first, projtt.second);
}

BOOST_AUTO_TEST_SUITE_END()


