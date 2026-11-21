/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_dynamic_order_math_parity.cpp
 * @brief Dynamic-order mathematical parity checks (static vs dynamic)
 */
#define BOOST_TEST_MODULE test_dynamic_order_math_parity
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>
#include <feel/feelvf/vf.hpp>

#include <string>
#include <type_traits>
#include <utility>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using mesh_type = Mesh<Simplex<2,1>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;

mesh_ptrtype
makeOneElementSimplexMesh( std::string const& name, bool referenceDomain )
{
    auto const domainType = referenceDomain ? Gmsh::GMSH_REFERENCE_DOMAIN : Gmsh::GMSH_REAL_DOMAIN;
    auto desc = std::make_shared<GmshSimplexDomain>( 2, 1, domainType );
    desc->setPrefix( name );
    desc->usePhysicalNames( true );
    if ( referenceDomain )
    {
        desc->setCharacteristicLength( 2.0 );
    }
    else
    {
        desc->setX( { 0.2, 1.4 } );
        desc->setY( { -0.1, 1.1 } );
        desc->setCharacteristicLength( 1.2 );
    }
    return createGMSHMesh( _mesh = new mesh_type, _desc = desc, _h = desc->h() );
}

template<int Order>
auto manufacturedScalarExpr()
{
    if constexpr ( Order == 0 )
        return cst( 1.0 );
    else if constexpr ( Order == 1 )
        return cst( 1.0 ) + 2.0 * Px() - 0.5 * Py();
    else if constexpr ( Order == 2 )
        return cst( 1.0 ) + Px() + Py() + Px() * Px() + Px() * Py() + Py() * Py();
    else if constexpr ( Order == 3 )
        return cst( 1.0 ) + Px() * Px() * Px() + Px() * Px() * Py() + Px() * Py() * Py() + Py() * Py() * Py();
    else
        return cst( 1.0 ) + pow( Px(), 4 ) + pow( Py(), 4 ) + Px() * Px() * Py() * Py();
}

template<int Order>
auto manufacturedVectorExpr()
{
    if constexpr ( Order == 0 )
        return vec( cst( 0.75 ), cst( -1.25 ) );
    else if constexpr ( Order == 1 )
        return vec( cst( 0.75 ) + 2.0 * Px() - Py(),
                    cst( -1.25 ) + Px() + 3.0 * Py() );
    else if constexpr ( Order == 2 )
        return vec( cst( 0.75 ) + Px() * Px() + Px() * Py(),
                    cst( -1.25 ) + Py() * Py() + 2.0 * Px() * Py() );
    else if constexpr ( Order == 3 )
        return vec( cst( 0.75 ) + pow( Px(), 3 ) + Px() * Py() * Py(),
                    cst( -1.25 ) + pow( Py(), 3 ) + Px() * Px() * Py() );
    else
        return vec( cst( 0.75 ) + pow( Px(), 4 ) + Px() * Py() * Py() * Py(),
                    cst( -1.25 ) + pow( Py(), 4 ) + Px() * Px() * Px() * Py() );
}

template<int Order>
auto manufacturedTensorExpr()
{
    if constexpr ( Order == 0 )
        return mat<2,2>( cst( 1.0 ), cst( -2.0 ), cst( 0.5 ), cst( 3.0 ) );
    else if constexpr ( Order == 1 )
        return mat<2,2>( cst( 1.0 ) + Px(),
                         cst( -2.0 ) + Py(),
                         cst( 0.5 ) + Px() - Py(),
                         cst( 3.0 ) + Px() + Py() );
    else if constexpr ( Order == 2 )
        return mat<2,2>( cst( 1.0 ) + Px() * Px() + Px() * Py(),
                         cst( -2.0 ) + Py() * Py(),
                         cst( 0.5 ) + Px() * Py(),
                         cst( 3.0 ) + Px() * Px() + Py() * Py() );
    else if constexpr ( Order == 3 )
        return mat<2,2>( cst( 1.0 ) + pow( Px(), 3 ),
                         cst( -2.0 ) + Px() * Px() * Py(),
                         cst( 0.5 ) + Px() * Py() * Py(),
                         cst( 3.0 ) + pow( Py(), 3 ) );
    else
        return mat<2,2>( cst( 1.0 ) + pow( Px(), 4 ),
                         cst( -2.0 ) + Px() * Px() * Py() * Py(),
                         cst( 0.5 ) + Px() * pow( Py(), 3 ),
                         cst( 3.0 ) + pow( Py(), 4 ) );
}

template<int Order>
void
checkScalarParity( mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto expr = manufacturedScalarExpr<Order>();

    BOOST_REQUIRE_EQUAL( VhDynamic->order(), Order );
    BOOST_REQUIRE_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto uStatic = project( _space = VhStatic, _expr = expr );
    auto uDynamic = project( _space = VhDynamic, _expr = expr );

    auto const errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - expr, _quad = _Q<20>() );
    auto const errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - expr, _quad = _Q<20>() );
    auto const errParity = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<20>() );

    BOOST_TEST_CONTEXT( "Pch static/dynamic parity P" << Order )
    {
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errParity, tol );
    }
}

template<int Order>
void
checkVectorParity( mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pdhv<Order>( mesh );
    auto VhDynamic = Pdhv<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto expr = manufacturedVectorExpr<Order>();

    BOOST_REQUIRE_EQUAL( VhDynamic->order(), Order );
    BOOST_REQUIRE_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto uStatic = VhStatic->element();
    auto uDynamic = VhDynamic->element();
    uStatic.on( _range = elements( mesh ), _expr = expr );
    uDynamic.on( _range = elements( mesh ), _expr = expr );

    auto const errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - expr, _quad = _Q<20>() );
    auto const errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - expr, _quad = _Q<20>() );
    auto const errParity = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<20>() );

    BOOST_TEST_CONTEXT( "Pdhv static/dynamic parity P" << Order )
    {
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errParity, tol );
    }
}

template<int Order>
void
checkTensorParity( mesh_ptrtype const& mesh, double tol )
{
    auto MhStatic = Pdhm<Order>( mesh );
    auto MhDynamic = Pdhm<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto expr = manufacturedTensorExpr<Order>();

    BOOST_REQUIRE_EQUAL( MhDynamic->order(), Order );
    BOOST_REQUIRE_EQUAL( MhStatic->nDof(), MhDynamic->nDof() );

    auto mStatic = MhStatic->element();
    auto mDynamic = MhDynamic->element();
    mStatic.on( _range = elements( mesh ), _expr = expr );
    mDynamic.on( _range = elements( mesh ), _expr = expr );

    auto const errStatic = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - expr, _quad = _Q<20>() );
    auto const errDynamic = normL2( _range = elements( mesh ), _expr = idv( mDynamic ) - expr, _quad = _Q<20>() );
    auto const errParity = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - idv( mDynamic ), _quad = _Q<20>() );

    BOOST_TEST_CONTEXT( "Pdhm static/dynamic parity P" << Order )
    {
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errParity, tol );
    }
}

template<int... Orders, typename Fn>
void
forEachOrder( std::integer_sequence<int, Orders...>, Fn&& fn )
{
    ( fn( std::integral_constant<int, Orders>{} ), ... );
}

template<int... Orders>
void
runDynamicOrderChecks( mesh_ptrtype const& mesh )
{
    constexpr double tolScalar = 1e-10;
    constexpr double tolVector = 1e-10;
    constexpr double tolTensor = 1e-10;

    forEachOrder( std::integer_sequence<int, Orders...>{},
                  [&]( auto orderC )
                  {
                      constexpr int order = decltype( orderC )::value;
                      checkScalarParity<order>( mesh, tolScalar );
                      checkVectorParity<order>( mesh, tolVector );
                      checkTensorParity<order>( mesh, tolTensor );
                  } );
}
/** @brief Check real modal space construction, DG cardinality and a constant field. */
template<int Degree,class Cell>
void checkModalConstantField(std::shared_ptr<Mesh<Cell>> const& mesh)
{
        using static_space = FunctionSpace<Mesh<Cell>,bases<OrthonormalPolynomialSet<Degree>>>;
        using dynamic_space = FunctionSpace<Mesh<Cell>,bases<OrthonormalPolynomialSet<Dynamic>>>;
        auto fixed = static_space::New(_mesh=mesh);
        auto runtime = dynamic_space::New(_mesh=mesh,_runtime_order=RuntimeOrder{Degree});
        BOOST_REQUIRE_EQUAL(fixed->nDof(),runtime->nDof());
        BOOST_CHECK_EQUAL(runtime->order(),Degree);
        constexpr auto count=Cell::is_hypercube ? (Degree+1)*(Degree+1) : (Degree+1)*(Degree+2)/2;
        BOOST_CHECK_EQUAL(runtime->fe()->localDof(),count);
        BOOST_CHECK_EQUAL(runtime->nDof(),mesh->numGlobalElements()*count);
        auto makeConstant = [&mesh](auto const& space)
        {
            auto field=space->element();
            field.zero();
            // A normalized constant mode has value 1/sqrt(reference measure).
            for(auto const& entry : elements(mesh))
            {
                auto const& element=unwrap_ref(entry);
                auto index=space->dof()->localToGlobal(element.id(),0,0).index();
                field(index)=std::sqrt(Cell::is_hypercube ? 4.0 : 2.0);
            }
            field.close();
            return field;
        };
        auto fixedConstant=makeConstant(fixed);
        auto runtimeConstant=makeConstant(runtime);
        BOOST_CHECK_SMALL(normL2(_range=elements(mesh),_expr=idv(fixedConstant)-cst(1.0)),2e-11);
        BOOST_CHECK_SMALL(normL2(_range=elements(mesh),_expr=idv(runtimeConstant)-cst(1.0)),2e-11);
        BOOST_CHECK_SMALL(normL2(_range=elements(mesh),_expr=idv(fixedConstant)-idv(runtimeConstant)),2e-11);
}
} // namespace

BOOST_AUTO_TEST_SUITE( dynamic_order_math_parity_suite )

BOOST_AUTO_TEST_CASE( dynamicOrderReferenceSimplexMathParity )
{
    auto mesh = makeOneElementSimplexMesh( "dynamic_order_reference", true );
    runDynamicOrderChecks<0,1,2,3,4>( mesh );
}

BOOST_AUTO_TEST_CASE( dynamicOrderAffineSimplexMathParity )
{
    auto mesh = makeOneElementSimplexMesh( "dynamic_order_affine", false );
    runDynamicOrderChecks<0,1,2,3,4>( mesh );
}

/** @test FunctionSpace validates modal order through the public descriptor policy. */
BOOST_AUTO_TEST_CASE( modalFunctionSpaceValidatesDescriptorOrder )
{
    using modal_mesh_type = Mesh<Simplex<2, 1>>;
    using static_space_type = FunctionSpace<modal_mesh_type, bases<OrthonormalPolynomialSet<2>>>;
    using dynamic_space_type = FunctionSpace<modal_mesh_type, bases<OrthonormalPolynomialSet<Dynamic>>>;
    static_assert( static_space_type::defaultBasisOrder().value == 2 );
    static_assert( dynamic_space_type::defaultBasisOrder().value == 0 );
    BOOST_CHECK_EQUAL( static_space_type::validatedBasisOrder( RuntimeOrder{ 2 } ).value, 2 );
    BOOST_CHECK_THROW( (void)static_space_type::validatedBasisOrder( RuntimeOrder{ 3 } ),
                       std::invalid_argument );
    for ( int degree : { 0, 1, 2, 3 } )
        BOOST_TEST_CONTEXT( "modal degree=" << degree )
            BOOST_CHECK_EQUAL( dynamic_space_type::validatedBasisOrder( RuntimeOrder{ degree } ).value, degree );
}

/** @test Modal FunctionSpaces own cell DOFs on simplices and hypercubes. */
BOOST_AUTO_TEST_CASE( modalFunctionSpaceConstructionAndConstantField )
{
    auto triangle=makeOneElementSimplexMesh("modal_constant_field",false);
    auto description=std::make_shared<GmshHypercubeDomain>(2,1,2,true);
    description->setPrefix("modal_tensor_product_constant_field");
    description->setCharacteristicLength(1.0);
    auto quadrilateral=createGMSHMesh(_mesh=new Mesh<Hypercube<2>>,_desc=description,_h=1.0);
    forEachOrder(std::integer_sequence<int,0,1,2,3>{},[&](auto order)
    {
        checkModalConstantField<decltype(order)::value>(triangle);
        checkModalConstantField<decltype(order)::value>(quadrilateral);
    });
}

BOOST_AUTO_TEST_SUITE_END()
