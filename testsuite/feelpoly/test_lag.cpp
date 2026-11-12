/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE test_lag
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/mp11/algorithm.hpp>

#include "testlagrange.hpp"

#include <array>
#include <type_traits>

using namespace Feel;
namespace bdata = boost::unit_test::data;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Reference dimensions included in the ordinary nodal-identity dataset. */
constexpr std::array<unsigned, 3> dimensions{ 1, 2, 3 };
/** Static degrees exercised for scalar and vector-valued elements. */
constexpr std::array<unsigned, 3> orders{ 1, 2, 3 };
/** Boolean choices reused for cell-family and scalar/vector datasets. */
constexpr std::array<bool, 2> alternatives{ false, true };

/** @brief Dispatch dataset values to one shared typed Lagrange assertion helper.
 * @tparam MaxDimension Largest dimension instantiated by this bounded dataset.
 * @tparam Orders MP11 list of the selected compile-time degrees.
 * @tparam IncludeVectorial Whether vector-valued element types are instantiated.
 * @param dimension Reference dimension, in `[1, MaxDimension]`.
 * @param orderIndex Index into Orders, rather than a polynomial degree.
 * @param tensorProduct Select a hypercube instead of a simplex.
 * @param vectorial Select a vector-valued instead of scalar Lagrange basis.
 */
template<unsigned MaxDimension, typename Orders, bool IncludeVectorial>
void checkConfiguration( unsigned dimension, unsigned orderIndex, bool tensorProduct, bool vectorial )
{
    BOOST_REQUIRE_GE( dimension, 1 );
    BOOST_REQUIRE_LE( dimension, MaxDimension );
    BOOST_REQUIRE_LT( orderIndex, boost::mp11::mp_size<Orders>::value );
    boost::mp11::mp_with_index<MaxDimension>( dimension - 1, [&]( auto dimensionTag )
    {
        constexpr unsigned Dim = decltype( dimensionTag )::value + 1;
        boost::mp11::mp_with_index<boost::mp11::mp_size<Orders>::value>( orderIndex, [&]( auto orderTag )
        {
            constexpr unsigned Order = boost::mp11::mp_at_c<Orders, decltype( orderTag )::value>::value;
            unsigned const index = IncludeVectorial ? 2 * tensorProduct + vectorial : tensorProduct;
            boost::mp11::mp_with_index<IncludeVectorial ? 4 : 2>( index, [&]( auto configurationTag )
            {
                constexpr unsigned configuration = decltype( configurationTag )::value *
                                                   ( IncludeVectorial ? 1 : 2 );
                using cell_type = std::conditional_t<( configuration & 2 ) != 0,
                                                      Hypercube<Dim>, Simplex<Dim>>;
                using family_type = std::conditional_t<( configuration & 1 ) != 0,
                    Lagrange<Order, Vectorial>, Lagrange<Order, Scalar>>;
                using fe_type = typename family_type::template apply<Dim, Dim, double, cell_type>::type;
                // Match the existing high-order static/runtime parity tolerances.
                TestLagrange<fe_type>{ Order >= 5 ? 1e-8 : 1e-11 }();
            } );
        } );
    } );
}
} // namespace

BOOST_AUTO_TEST_SUITE( lagrange_testsuite )

/** @test Low-order scalar/vector Lagrange bases interpolate and differentiate. */
BOOST_DATA_TEST_CASE( nodal_identity_and_differentiation,
    bdata::make( dimensions ) * bdata::make( orders ) *
    bdata::make( alternatives ) * bdata::make( alternatives ),
    dimension, order, tensorProduct, vectorial )
{
    using order_types = boost::mp11::mp_list_c<unsigned, 1, 2, 3>;
    checkConfiguration<3, order_types, true>( dimension, order - 1, tensorProduct, vectorial );
}

/** @test Degree-five/ten scalar bases retain nodal identity in one/two dimensions. */
BOOST_DATA_TEST_CASE( high_order_scalar_identity_and_differentiation,
    bdata::make( std::array<unsigned, 2>{ 1, 2 } ) *
    bdata::make( std::array<unsigned, 2>{ 5, 10 } ) * bdata::make( alternatives ),
    dimension, order, tensorProduct )
{
    using order_types = boost::mp11::mp_list_c<unsigned, 5, 10>;
    checkConfiguration<2, order_types, false>( dimension, order == 5 ? 0 : 1, tensorProduct, false );
}

/** @test Degree-five tetrahedron and hexahedron bases share the same assertions. */
BOOST_DATA_TEST_CASE( degree_five_3d_scalar_identity_and_differentiation,
    bdata::make( alternatives ), tensorProduct )
{
    using order_types = boost::mp11::mp_list_c<unsigned, 5>;
    checkConfiguration<3, order_types, false>( 3, 0, tensorProduct, false );
}

BOOST_AUTO_TEST_SUITE_END()
