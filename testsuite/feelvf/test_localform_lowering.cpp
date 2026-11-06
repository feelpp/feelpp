/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_localform_lowering
#include <feel/feelcore/testsuite.hpp>

#include <algorithm>
#include <cmath>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>

#include <feel/feelvf/vf_eval.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/detail/localform.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{

/** @brief Require equality of complete assembled matrices through zero norms. */
template<typename OptimizedForm, typename ReferenceForm>
void checkCompleteMatrixParity( OptimizedForm const& optimized,
                                ReferenceForm const& reference,
                                double tolerance = 1e-12 )
{
    optimized.matrixPtr()->close();
    reference.matrixPtr()->close();
    BOOST_REQUIRE_EQUAL( optimized.matrixPtr()->size1(),
                         reference.matrixPtr()->size1() );
    BOOST_REQUIRE_EQUAL( optimized.matrixPtr()->size2(),
                         reference.matrixPtr()->size2() );
    auto const scale = std::max(
        { 1., std::abs( optimized.matrixPtr()->l1Norm() ),
          std::abs( reference.matrixPtr()->l1Norm() ),
          std::abs( optimized.matrixPtr()->linftyNorm() ),
          std::abs( reference.matrixPtr()->linftyNorm() ) } );
    reference.matrixPtr()->addMatrix( -1., optimized.matrixPtr() );
    reference.matrixPtr()->close();
    BOOST_CHECK_SMALL( reference.matrixPtr()->l1Norm()/scale, tolerance );
    BOOST_CHECK_SMALL( reference.matrixPtr()->linftyNorm()/scale, tolerance );
}

/** @brief Check canonical scalar mass and stiffness matrices at one order. */
template<int Order, typename MeshType>
void checkScalarLocalMatrices( std::shared_ptr<MeshType> const& mesh )
{
    auto space = Pch<Order>( mesh );
    auto trial = space->element();
    auto test = space->element();
    auto coefficient = cst( 1.75 ) + Px();

    auto massOptimized = form2( _trial = space, _test = space );
    auto massReference = form2( _trial = space, _test = space );
    auto mass = coefficient*idt( trial )*id( test );
    massOptimized = integrate(
        _range = elements( mesh ), _expr = mass );
    massReference = integrate(
        _range = elements( mesh ),
        _expr = mass, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( massOptimized, massReference );

    auto stiffnessOptimized = form2( _trial = space, _test = space );
    auto stiffnessReference = form2( _trial = space, _test = space );
    auto stiffness = coefficient*inner( gradt( trial ), grad( test ) );
    stiffnessOptimized = integrate(
        _range = elements( mesh ), _expr = stiffness,
        _local_form_lowering = local_form_lowering::enabled );
    stiffnessReference = integrate(
        _range = elements( mesh ),
        _expr = stiffness, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( stiffnessOptimized, stiffnessReference );
}

} // namespace

BOOST_AUTO_TEST_SUITE( localform_lowering_suite )

BOOST_AUTO_TEST_CASE( scalar_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );

    auto expr_mass = cst( 2.5 ) * idt( u ) * id( v );
    auto expr_source = cst( 3.0 ) * id( v );
    auto expr_unsupported = id( v ) + cst( 1.0 );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_linear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass,
                           _local_form_lowering=local_form_lowering::disabled );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source,
                           _local_form_lowering=local_form_lowering::disabled );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto measure = integrate( _range=elements( mesh ), _expr=cst( 1.0 ) ).evaluate()( 0, 0 );
    auto expected_mass = 2.5 * measure;
    auto expected_source = 3.0 * measure;

    BOOST_CHECK_CLOSE( a_generic( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_lowered( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_chi_weighted_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );

    auto indicator = chi( Px() > cst( 0.5 ) );
    auto expr_mass = indicator * cst( 4.0 ) * idt( u ) * id( v );
    auto expr_source = indicator * cst( 1.25 ) * id( v );
    auto expr_nonlinear = chi( id( v ) > cst( 0.5 ) ) * id( v );

    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( indicator )>::supported );
    static_assert( vf::detail::scalar_localform_factors<decltype( indicator * cst( 4.0 ) )>::supported );
    static_assert( vf::detail::scalar_localform_factors<decltype( indicator * cst( 4.0 ) * idt( u ) )>::supported );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerNonlinear =
        vf::detail::can_lower_scalar_linear_localform_v<decltype( expr_nonlinear )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerNonlinear );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass,
                           _local_form_lowering=local_form_lowering::disabled );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source,
                           _local_form_lowering=local_form_lowering::disabled );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto indicator_measure = integrate( _range=elements( mesh ), _expr=indicator ).evaluate()( 0, 0 );
    auto expected_source = 1.25 * indicator_measure;

    BOOST_TEST( indicator_measure > 0.45 );
    BOOST_TEST( indicator_measure < 0.55 );
    BOOST_TEST( a_generic( one, one ) > 0.0 );
    BOOST_TEST( a_lowered( one, one ) > 0.0 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_field_coefficient_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto Mh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    auto beta = Mh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );
    beta.on( _range=elements( mesh ), _expr=Px() + cst( 2. ) );

    auto expr_mass = idv( beta ) * idt( u ) * id( v );
    auto expr_source = ( idv( beta ) + cst( 0.5 ) ) * id( v );
    auto expr_unsupported = idv( beta ) * id( v ) * id( v );

    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( idv( beta ) )>::supported );
    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( idv( beta ) )>::kind ==
                   vf::detail::localform_scalar_expr_kind::scalar_value_leaf );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_bilinear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass,
                           _local_form_lowering=local_form_lowering::disabled );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source,
                           _local_form_lowering=local_form_lowering::disabled );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto expected_mass = integrate( _range=elements( mesh ), _expr=idv( beta ) ).evaluate()( 0, 0 );
    auto expected_source = integrate( _range=elements( mesh ), _expr=idv( beta ) + cst( 0.5 ) ).evaluate()( 0, 0 );

    BOOST_CHECK_CLOSE( a_generic( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_lowered( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_geometry_coefficient_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );

    auto expr_mass = Px() * idt( u ) * id( v );
    auto expr_source = ( Px() + cst( 0.5 ) ) * id( v );
    auto expr_unsupported = Px() * id( v ) * id( v );

    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( Px() )>::supported );
    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( Px() )>::kind ==
                   vf::detail::localform_scalar_expr_kind::geometry_scalar_leaf );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_bilinear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass,
                           _local_form_lowering=local_form_lowering::disabled );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source,
                           _local_form_lowering=local_form_lowering::disabled );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto expected_mass = integrate( _range=elements( mesh ), _expr=Px() ).evaluate()( 0, 0 );
    auto expected_source = integrate( _range=elements( mesh ), _expr=Px() + cst( 0.5 ) ).evaluate()( 0, 0 );

    BOOST_CHECK_CLOSE( a_generic( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_lowered( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_mixed_coefficient_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto Mh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    auto beta = Mh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );
    beta.on( _range=elements( mesh ), _expr=Px() + cst( 2. ) );

    auto indicator = chi( Px() > cst( 0.5 ) );
    auto mixed_coefficient = indicator * ( idv( beta ) + cst( 0.5 ) );
    auto expr_mass = mixed_coefficient * idt( u ) * id( v );
    auto expr_source = mixed_coefficient * id( v );
    auto expr_unsupported = indicator * ( id( v ) + cst( 0.5 ) ) * id( v );

    using coefficient_support = vf::detail::scalar_localform_coefficient_expr<decltype( mixed_coefficient )>;
    static_assert( coefficient_support::supported );
    static_assert( vf::detail::is_localform_coeff_ir_v<typename coefficient_support::coefficient_expr_type> );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_linear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass,
                           _local_form_lowering=local_form_lowering::disabled );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source,
                           _local_form_lowering=local_form_lowering::disabled );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    BOOST_TEST( a_generic( one, one ) > 0.0 );
    BOOST_TEST( a_lowered( one, one ) > 0.0 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_TEST( l_generic( one ) > 0.0 );
    BOOST_TEST( l_lowered( one ) > 0.0 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

/** Verify complete scalar local-tensor assembly for representative orders. */
BOOST_AUTO_TEST_CASE( scalar_mass_stiffness_complete_matrix_orders )
{
    auto mesh = unitSquare();
    checkScalarLocalMatrices<1>( mesh );
    checkScalarLocalMatrices<2>( mesh );
    checkScalarLocalMatrices<3>( mesh );
}

/** Verify fixed and runtime P2 spaces use the same canonical contractions. */
BOOST_AUTO_TEST_CASE( scalar_runtime_order_complete_matrix )
{
    auto mesh = unitSquare();
    auto space = Pch<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto trial = space->element();
    auto test = space->element();
    auto coefficient = cst( 2. ) + Py();

    auto optimized = form2( _trial = space, _test = space );
    auto reference = form2( _trial = space, _test = space );
    auto stiffness = coefficient*inner( gradt( trial ), grad( test ) );
    auto loweredStiffness =
        vf::detail::lower_scalar_h1_stiffness_localform( stiffness );
    optimized = integrate(
        _range = elements( mesh ), _expr = loweredStiffness );
    reference = integrate(
        _range = elements( mesh ),
        _expr = stiffness, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( optimized, reference );
}

/** Verify tensor diffusion on curved geometry for static and runtime P2. */
BOOST_AUTO_TEST_CASE( anisotropic_scalar_h1_static_runtime_curved )
{
    auto affineMesh = unitSquare();
    using curved_mesh_type = Mesh<Simplex<2, 2>>;
    auto mesh = createGMSHMesh(
        _mesh = new curved_mesh_type,
        _desc = domain( _name = "localform-anisotropic-g2",
                        _shape = "ellipsoid", _dim = 2, _order = 2,
                        _xmin = -1., _xmax = 1.,
                        _ymin = -1., _ymax = 1.,
                        _h = 0.35, _straighten = false ),
        _h = 0.35, _straighten = false, _verbose = false );

    auto coefficient = mat<2, 2>(
        cst( 2. ) + Px()*Px(), cst( 0.25 ) + cst( 0.1 )*Py(),
        cst( 0.25 ) + cst( 0.1 )*Py(), cst( 1.5 ) + Py()*Py() );

    auto checkSpace = [&]( auto const& currentMesh, auto const& space )
    {
        auto trial = space->element();
        auto test = space->element();
        auto anisotropic =
            gradt( trial )*coefficient*trans( grad( test ) );
        static_assert(
            vf::detail::can_lower_anisotropic_scalar_h1_stiffness_localform_v<
                decltype( anisotropic )> );
        auto lowered =
            vf::detail::lower_anisotropic_scalar_h1_stiffness_localform(
                anisotropic );
        using lowered_expression_type =
            typename std::decay_t<decltype( lowered )>::expression_type;
        using space_type = typename std::decay_t<decltype( space )>::element_type;
        using basis_type = typename space_type::reference_element_type;
        static_assert(
            lowered_expression_type::test_workspace_rows ==
            ( basis_type::is_order_dynamic ? Eigen::Dynamic :
              static_cast<int>( basis_type::nLocalDof ) ) );
        auto optimized = form2( _trial = space, _test = space );
        auto reference = form2( _trial = space, _test = space );
        optimized = integrate(
            _range = elements( currentMesh ), _expr = anisotropic,
            _local_form_lowering = local_form_lowering::enabled,
            _geomap = GeomapStrategyType::GEOMAP_HO );
        reference = integrate(
            _range = elements( currentMesh ),
            _expr = anisotropic,
            _local_form_lowering = local_form_lowering::disabled,
            _geomap = GeomapStrategyType::GEOMAP_HO );
        checkCompleteMatrixParity( optimized, reference, 2e-12 );
    };

    checkSpace( affineMesh, Pch<1>( affineMesh ) );
    checkSpace( mesh, Pch<2>( mesh ) );
    checkSpace( mesh, Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ) );
}

/** Verify vector, H(curl), H(div), and mixed canonical matrices. */
BOOST_AUTO_TEST_CASE( compatible_space_complete_matrices )
{
    auto mesh = unitSquare();

    auto vectorSpace = Pchv<1>( mesh );
    auto vectorTrial = vectorSpace->element();
    auto vectorTest = vectorSpace->element();
    auto vectorMass = inner( idt( vectorTrial ), id( vectorTest ) );
    auto vectorOptimized = form2( _trial = vectorSpace, _test = vectorSpace );
    auto vectorReference = form2( _trial = vectorSpace, _test = vectorSpace );
    vectorOptimized = integrate(
        _range = elements( mesh ), _expr = vectorMass,
        _local_form_lowering = local_form_lowering::enabled );
    vectorReference = integrate(
        _range = elements( mesh ),
        _expr = vectorMass, _local_form_lowering = local_form_lowering::disabled );
    auto vectorAutomatic = form2( _trial = vectorSpace, _test = vectorSpace );
    vectorAutomatic = integrate(
        _range = elements( mesh ), _expr = vectorMass );
    checkCompleteMatrixParity( vectorOptimized, vectorAutomatic );
    checkCompleteMatrixParity( vectorOptimized, vectorReference );

    auto hcurlSpace = Ned1h<0>( mesh );
    auto hcurlTrial = hcurlSpace->element();
    auto hcurlTest = hcurlSpace->element();
    auto hcurlMass = inner( idt( hcurlTrial ), id( hcurlTest ) );
    auto curlCurl = inner( curlt( hcurlTrial ), curl( hcurlTest ) );
    auto hcurlMassOptimized = form2( _trial = hcurlSpace, _test = hcurlSpace );
    auto hcurlMassReference = form2( _trial = hcurlSpace, _test = hcurlSpace );
    hcurlMassOptimized = integrate(
        _range = elements( mesh ), _expr = hcurlMass,
        _local_form_lowering = local_form_lowering::enabled );
    hcurlMassReference = integrate(
        _range = elements( mesh ),
        _expr = hcurlMass, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( hcurlMassOptimized, hcurlMassReference );

    auto curlOptimized = form2( _trial = hcurlSpace, _test = hcurlSpace );
    auto curlReference = form2( _trial = hcurlSpace, _test = hcurlSpace );
    curlOptimized = integrate(
        _range = elements( mesh ), _expr = curlCurl,
        _local_form_lowering = local_form_lowering::enabled );
    curlReference = integrate(
        _range = elements( mesh ),
        _expr = curlCurl, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( curlOptimized, curlReference );

    auto hdivSpace = RTh<0>( mesh );
    auto pressureSpace = Pdh<0>( mesh );
    auto hdivTrial = hdivSpace->element();
    auto hdivTest = hdivSpace->element();
    auto pressureTrial = pressureSpace->element();
    auto pressureTest = pressureSpace->element();
    auto hdivMass = inner( idt( hdivTrial ), id( hdivTest ) );
    auto divDiv = divt( hdivTrial )*div( hdivTest );

    auto hdivMassOptimized = form2( _trial = hdivSpace, _test = hdivSpace );
    auto hdivMassReference = form2( _trial = hdivSpace, _test = hdivSpace );
    hdivMassOptimized = integrate(
        _range = elements( mesh ), _expr = hdivMass,
        _local_form_lowering = local_form_lowering::enabled );
    hdivMassReference = integrate(
        _range = elements( mesh ),
        _expr = hdivMass, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( hdivMassOptimized, hdivMassReference );

    auto divOptimized = form2( _trial = hdivSpace, _test = hdivSpace );
    auto divReference = form2( _trial = hdivSpace, _test = hdivSpace );
    divOptimized = integrate(
        _range = elements( mesh ), _expr = divDiv,
        _local_form_lowering = local_form_lowering::enabled );
    divReference = integrate(
        _range = elements( mesh ),
        _expr = divDiv, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( divOptimized, divReference );

    auto mixed = divt( hdivTrial )*id( pressureTest );
    auto mixedOptimized = form2( _trial = hdivSpace, _test = pressureSpace );
    auto mixedReference = form2( _trial = hdivSpace, _test = pressureSpace );
    mixedOptimized = integrate(
        _range = elements( mesh ), _expr = mixed,
        _local_form_lowering = local_form_lowering::enabled );
    mixedReference = integrate(
        _range = elements( mesh ),
        _expr = mixed, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( mixedOptimized, mixedReference );

    auto adjointMixed = idt( pressureTrial )*div( hdivTest );
    auto adjointOptimized = form2(
        _trial = pressureSpace, _test = hdivSpace );
    auto adjointReference = form2(
        _trial = pressureSpace, _test = hdivSpace );
    adjointOptimized = integrate(
        _range = elements( mesh ), _expr = adjointMixed,
        _local_form_lowering = local_form_lowering::enabled );
    adjointReference = integrate(
        _range = elements( mesh ),
        _expr = adjointMixed, _local_form_lowering = local_form_lowering::disabled );
    checkCompleteMatrixParity( adjointOptimized, adjointReference );
}

BOOST_AUTO_TEST_SUITE_END()
