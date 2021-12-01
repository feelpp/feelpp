/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 09 Feb 2020

 Copyright (C) 2020 Feel++ Consortium

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
#include <cmath>
#define BOOST_TEST_MODULE vector testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/vectorublas.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( tensor )

BOOST_AUTO_TEST_CASE( test_tensor_1 )
{
    using namespace Eigen;
    using Eigen::Tensor;

    typedef Tensor<float, 1>::DimensionPair DimPair;

    Tensor<float, 3> mat1(1, 2, 3);
    Tensor<float, 2> mat2(2, 2);

    mat1.setRandom();
    mat2.setRandom();

    Tensor<float, 3> mat3(1, 2, 3);
    mat3.setZero();
    Eigen::array<DimPair, 1> dims = {{DimPair(1, 1)}};
    typedef TensorEvaluator<decltype(mat1.contract(mat2, dims)), DefaultDevice> Evaluator;
    Evaluator eval(mat1.contract(mat2, dims), DefaultDevice());
    std::cout << "eval.numdims=" << Evaluator::NumDims << std::endl;
    std::cout << "eval=" << eval.dimensions()[0] << "," << eval.dimensions()[1] << "," << eval.dimensions()[2] << std::endl;
    mat3 = mat1.contract(mat2, dims);
    std::cout << "mat1=" << mat1 << std::endl;
    std::cout << "mat2=" << mat2 << std::endl;
    std::cout << "mat3=" << mat3 << std::endl;

    TensorFixedSize<float, Eigen::Sizes<2,2>> mat4;
    mat4 = mat2;
}
BOOST_AUTO_TEST_SUITE_END()

