/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 28 Mar 2016

 Copyright (C) 2016 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#define BOOST_TEST_MODULE test_eigen3_tensor
#include <feel/feelcore/testsuite.hpp>

#include <Eigen/CXX11/Tensor>

using namespace Feel;


using Eigen::DefaultDevice;
using Eigen::Tensor;

typedef Tensor<double, 1>::DimensionPair DimPair;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( eigen3_tensor )

BOOST_AUTO_TEST_CASE( tensor1 )
{
    Tensor<double, 3> mat1(2, 2, 2);
    Tensor<double, 2> mat2(2, 2);

    mat1.setRandom();
    mat2.setRandom();

    Tensor<double, 3 > mat4(2,2,2);
    mat4.setZero();
    Eigen::array<DimPair, 1> dims3 = {{DimPair(2, 0)}};
    mat4 = mat1.contract(mat2, dims3);

    std::cout << "mat1=" << mat1 << "\n";
    std::cout << "mat2=" << mat2 << "\n";
    std::cout << "mat4=" << mat4 << "\n";
    BOOST_CHECK_CLOSE(mat4(0,0,0), mat1(0,0,0)*mat2(0,0) +mat1(0,0,1)*mat2(1,0), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(0,0,1), mat1(0,0,0)*mat2(0,1) +mat1(0,0,1)*mat2(1,1), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(0,1,0), mat1(0,1,0)*mat2(0,0) +mat1(0,1,1)*mat2(1,0), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(0,1,1), mat1(0,1,0)*mat2(0,1) +mat1(0,1,1)*mat2(1,1), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(1,0,0), mat1(1,0,0)*mat2(0,0) +mat1(1,0,1)*mat2(1,0), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(1,0,1), mat1(1,0,0)*mat2(0,1) +mat1(1,0,1)*mat2(1,1), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(1,1,0), mat1(1,1,0)*mat2(0,0) +mat1(1,1,1)*mat2(1,0), 1e-13 );
    BOOST_CHECK_CLOSE(mat4(1,1,1), mat1(1,1,0)*mat2(0,1) +mat1(1,1,1)*mat2(1,1), 1e-13 );
}


BOOST_AUTO_TEST_SUITE_END()
