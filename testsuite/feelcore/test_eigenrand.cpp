/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 Date: 6 Jul 2023

 Copyright (C) 2023 Feel++ Consortium

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
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <chrono>
#include <random>

#include <EigenRand/EigenRand>
#include <Eigen/Dense>
#define BOOST_TEST_MODULE test_eigenrand
#include <feel/feelcore/testsuite.hpp>


BOOST_AUTO_TEST_SUITE( eigenrand )

BOOST_AUTO_TEST_CASE( test_0 )
{
    using namespace Eigen;
    // Initialize random number generator with seed=42 for following codes.
    // Or you can use C++11 RNG such as std::mt19937 or std::ranlux48.
    Rand::P8_mt19937_64 urng{ 42 };
    
    // this will generate 4x4 real matrix with range [-1, 1]
    MatrixXf mat = Rand::balanced<MatrixXf>(4, 4, urng);
    std::cout << mat << std::endl;
    
    // this will generate 10x10 real 2d array on the normal distribution
    ArrayXXf arr = Rand::normal<ArrayXXf>(10, 10, urng);
    std::cout << arr << std::endl;
}
BOOST_AUTO_TEST_CASE( test_1 )
{
    using namespace Eigen;
    Rand::P8_mt19937_64 urng{ 42 };
    
    MatrixXf mat{ 10, 10 };
    // this will generate a random matrix in MatrixXf type with the shape (10, 10)
    // note: it doesn't change mat at all.
    Rand::balancedLike(mat, urng);
    
    // if you want to assign a random matrix into itself, use assignment operator.
    mat = Rand::balancedLike(mat, urng);
    std::cout << mat << std::endl;
}

/**
 * @brief Vectorization over Parameters
 * EigenRand's random number generators typically accept scalar parameters. However, certain generators can generate random numbers efficiently 
 * for an array of parameters in an element-wise manner. You can see the full list of distributions which support the vectorization over parameters 
 * at list_of_supported_distribution.
 */
BOOST_AUTO_TEST_CASE( test_2 )
{
    using namespace Eigen;
    Rand::P8_mt19937_64 urng{ 42 };
    
    ArrayXf a{ 10 }, b{ 10 }, c{ 10 };
        a << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
        b << 10, 12, 14, 16, 18, 20, 22, 24, 26, 28;
    
    // You can use two array parameters.
    // The shape of two parameters should be equal in this case.
    c = Rand::uniformReal(urng, a, b);
    std::cout << c << std::endl;
    // c[0] is generated in the range [a[0], b[0]), 
    // c[1] is generated in the range [a[1], b[1]) ...
    
    // Or you can provide one parameter as a scalar
    // In this case, a scalar parameter is broadcast to the shape of the array parameter.
    c = Rand::uniformReal(urng, -5, b);
    std::cout << c << std::endl;
    // c[0] is generated in the range [-5, b[0]), 
    // c[1] is generated in the range [-5, b[1]) ...
    
    c = Rand::uniformReal(urng, a, 11);
    std::cout << c << std::endl;
    // c[0] is generated in the range [a[0], 11), 
    // c[1] is generated in the range [a[1], 11) ...
}

/**
 * @brief Efficient Reusable Generator
 * In the example above, functions, such as Eigen::Rand::balancedLike, Eigen::Rand::normal and so on, 
 * creates a generator internally each time to be called. If you want to generate random matrices from the same distribution, 
 * consider using Generator classes as following:
 */
BOOST_AUTO_TEST_CASE( test_3 )
{
    using namespace Eigen;
    Rand::P8_mt19937_64 urng{ 42 };
    // constructs generator for normal distribution with mean=1.0, stdev=2.0
    Rand::NormalGen<float> norm_gen{ 1.0, 2.0 };
    
    // Generator classes have a template function `generate`.
    // 10 by 10 random matrix will be assigned to `mat`.
    MatrixXf mat = norm_gen.template generate<MatrixXf>(10, 10, urng);
    std::cout << mat << std::endl;
    
    // Generator classes also have `generateLike`. 
    mat = norm_gen.generateLike(mat, urng);
    std::cout << mat << std::endl;
}

/**
 * @brief Drawing samples from Multivariate Distribution
 * EigenRand provides generators for some multivariate distributions.
 * 
 */
BOOST_AUTO_TEST_CASE( test_4 )
{
    using namespace Eigen;
    Rand::P8_mt19937_64 urng{ 42 };

    Vector4f mean{ 0, 1, 2, 3 };
    Matrix4f cov;
    cov << 1, 1, 0, 0,
            1, 2, 0, 0,
            0, 0, 3, 1,
            0, 0, 1, 2;
    {
    // constructs MvNormalGen with Scalar=float, Dim=4
    Rand::MvNormalGen<float, 4> gen1{ mean, cov };

    // or you can use `make-` helper function. It can deduce the type of generator to be created.
    auto gen2 = Rand::makeMvNormalGen(mean, cov);

    // generates one sample ( shape (4, 1) )
    Vector4f sample = gen1.generate(urng);

    // generates 10 samples ( shape (4, 10) )
    Matrix<float, 4, -1> samples = gen1.generate(urng, 10);
    // or you can just use `MatrixXf` type
    }

    {
    // construct MvWishartGen with Scalar=float, Dim=4, df=4
    auto gen3 = Rand::makeWishartGen(4, cov);

    // generates one sample ( shape (4, 4) )
    Matrix4f sample = gen3.generate(urng);

    // generates 10 samples ( shape (4, 40) )
    Matrix<float, 4, -1> samples = gen3.generate(urng, 10);
    // or you can just use `MatrixXf` type
    }
}

BOOST_AUTO_TEST_SUITE_END()