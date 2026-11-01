/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_hdiv_runtime_order_parity.cpp
 * @brief Runtime-order constructor parity checks for H(div)/H(curl) wrappers.
 */
#define BOOST_TEST_MODULE test_hdiv_runtime_order_parity
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/neh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelvf/vf.hpp>

#include <stdexcept>
#include <string>

#include "fe_parity_checks.hpp"

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<int Dim>
using simplex_mesh_type = Mesh<Simplex<Dim,1>>;

template<int Dim>
using simplex_mesh_ptrtype = std::shared_ptr<simplex_mesh_type<Dim>>;

template<int Dim>
simplex_mesh_ptrtype<Dim>
makeOneElementSimplexMesh( std::string const& name )
{
    auto desc = std::make_shared<GmshSimplexDomain>( Dim, 1, Gmsh::GMSH_REFERENCE_DOMAIN );
    desc->setPrefix( name );
    desc->usePhysicalNames( true );
    desc->setCharacteristicLength( 2.0 );
    return createGMSHMesh( _mesh = new simplex_mesh_type<Dim>, _desc = desc, _h = desc->h() );
}

template<typename StaticSpacePtrType, typename RuntimeSpacePtrType, typename MeshPtrType, typename ExprType>
void
checkOnParity( StaticSpacePtrType const& staticSpace,
               RuntimeSpacePtrType const& runtimeSpace,
               MeshPtrType const& mesh,
               ExprType const& expr,
               double tol,
               bool compareOrder = true )
{
    BOOST_REQUIRE_EQUAL( staticSpace->nDof(), runtimeSpace->nDof() );
    if ( compareOrder )
        BOOST_CHECK_EQUAL( staticSpace->order(), runtimeSpace->order() );

    auto uStatic = staticSpace->element();
    auto uRuntime = runtimeSpace->element();
    uStatic.on( _range = elements( mesh ), _expr = expr );
    uRuntime.on( _range = elements( mesh ), _expr = expr );

    auto const parityError = normL2( _range = elements( mesh ),
                                     _expr = idv( uStatic ) - idv( uRuntime ),
                                     _quad = _Q<12>() );
    BOOST_CHECK_SMALL( parityError, tol );
}

template<typename StaticMatrix, typename RuntimeMatrix>
void
checkMatrixParity( StaticMatrix const& staticMatrix,
                   RuntimeMatrix const& runtimeMatrix,
                   double tolerance )
{
    BOOST_REQUIRE_EQUAL( staticMatrix.size1(), runtimeMatrix.size1() );
    BOOST_REQUIRE_EQUAL( staticMatrix.size2(), runtimeMatrix.size2() );
    for ( std::size_t i = 0; i < staticMatrix.size1(); ++i )
        for ( std::size_t j = 0; j < staticMatrix.size2(); ++j )
            BOOST_CHECK_SMALL(
                staticMatrix( i, j )-runtimeMatrix( i, j ), tolerance );
}

template<typename StaticSpacePtrType, typename RuntimeSpacePtrType,
         typename MeshPtrType>
void
checkReferenceFiniteElementParity( StaticSpacePtrType const& staticSpace,
                                   RuntimeSpacePtrType const& runtimeSpace,
                                   MeshPtrType const& mesh )
{
    auto const& staticFe = *staticSpace->basis();
    auto const& runtimeFe = *runtimeSpace->basis();
    BOOST_REQUIRE_EQUAL( staticFe.localDofCount(), runtimeFe.localDofCount() );
    checkMatrixParity( staticFe.primal().coeff(),
                       runtimeFe.primal().coeff(), 1e-12 );
    checkMatrixParity( staticFe.dual()( staticFe.primal() ),
                       runtimeFe.dual()( runtimeFe.primal() ), 1e-11 );
    checkMatrixParity( staticFe.coeff(), runtimeFe.coeff(), 1e-10 );
    checkMatrixParity( staticFe.points(), runtimeFe.points(), 1e-12 );
    for ( uint16_type localDof = 0;
          localDof < staticFe.localDofCount(); ++localDof )
    {
        auto const staticAttachment = staticFe.dofAttachment( localDof );
        auto const runtimeAttachment = runtimeFe.dofAttachment( localDof );
        BOOST_CHECK_EQUAL( staticAttachment.entityDim,
                           runtimeAttachment.entityDim );
        BOOST_CHECK_EQUAL( staticAttachment.entityId,
                           runtimeAttachment.entityId );
        BOOST_CHECK_EQUAL( staticAttachment.ordinal,
                           runtimeAttachment.ordinal );
        BOOST_CHECK_EQUAL( staticFe.dofFunctionalKind( localDof ),
                           runtimeFe.dofFunctionalKind( localDof ) );
    }
    BOOST_CHECK( !staticFe.dofAttachment(
                     staticFe.localDofCount() ).isValid() );
    BOOST_CHECK( !runtimeFe.dofAttachment(
                     runtimeFe.localDofCount() ).isValid() );

    Feel::test_detail::checkFiniteElementTransformParity(
        staticSpace, runtimeSpace, mesh );
}

template<typename SpacePtrType, typename MeshPtrType, typename ExprType>
void
checkInterpolationExactness( SpacePtrType const& space,
                             MeshPtrType const& mesh,
                             ExprType const& expr,
                             double tol )
{
    auto u = space->element();
    u.on( _range = elements( mesh ), _expr = expr );

    auto const error = normL2( _range = elements( mesh ),
                               _expr = idv( u ) - expr,
                               _quad = _Q<14>() );
    BOOST_CHECK_SMALL( error, tol );
}

template<typename StaticSpacePtrType, typename RuntimeSpacePtrType, typename MeshPtrType, typename ExprType>
void
checkBoundaryTangentialInterpolation(
    StaticSpacePtrType const& staticSpace,
    RuntimeSpacePtrType const& runtimeSpace,
    MeshPtrType const& mesh,
    ExprType const& expr,
    double tol )
{
    auto uStatic = staticSpace->element();
    auto uRuntime = runtimeSpace->element();
    uStatic.on( _range = boundaryfaces( mesh ), _expr = expr );
    uRuntime.on( _range = boundaryfaces( mesh ), _expr = expr );

    const auto parityError = normL2(
        _range = boundaryfaces( mesh ),
        _expr = cross( idv( uStatic )-idv( uRuntime ), N() ),
        _quad = _Q<14>() );
    const auto staticError = normL2(
        _range = boundaryfaces( mesh ),
        _expr = cross( idv( uStatic )-expr, N() ),
        _quad = _Q<14>() );
    const auto runtimeError = normL2(
        _range = boundaryfaces( mesh ),
        _expr = cross( idv( uRuntime )-expr, N() ),
        _quad = _Q<14>() );
    BOOST_CHECK_SMALL( parityError, tol );
    BOOST_CHECK_SMALL( staticError, tol );
    BOOST_CHECK_SMALL( runtimeError, tol );
}

template<int Dim, int Order>
auto
rtExactExpr()
{
    if constexpr ( Dim == 2 && Order == 0 )
        return vec( cst( 1.25 ), cst( -0.75 ) );
    else if constexpr ( Dim == 3 && Order == 0 )
        return vec( cst( 1.25 ), cst( -0.75 ), cst( 0.5 ) );
    else if constexpr ( Dim == 2 && Order == 1 )
        return vec( cst( 1.0 ) + 2.0*Px() - 0.25*Py(),
                    cst( -0.5 ) + Px() + 3.0*Py() );
    else if constexpr ( Dim == 3 && Order == 1 )
        return vec( cst( 1.0 ) + 2.0*Px() - 0.25*Py(),
                    cst( -0.5 ) + Px() + 3.0*Py(),
                    cst( 0.25 ) - Py() + 0.5*Pz() );
    else if constexpr ( Dim == 2 && Order == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py() );
    else if constexpr ( Dim == 3 && Order == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py(),
                    cst( 0.25 ) + Pz() + Px()*Pz() );
    else if constexpr ( Dim == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py() + Px()*Py()*Py() );
    else
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py() + Px()*Py()*Py(),
                    cst( 0.25 ) + Pz() + Px()*Pz() + Px()*Py()*Pz() );
}

template<int Dim, int Order>
auto
bdmExactExpr()
{
    if constexpr ( Dim == 2 && Order == 0 )
        return vec( cst( 1.0 ) + 0.5*Px() - Py(),
                    cst( -0.25 ) + 2.0*Py() + 0.125*Px() );
    else if constexpr ( Dim == 3 && Order == 0 )
        return vec( cst( 1.0 ) + 0.5*Px() - Py(),
                    cst( -0.25 ) + 2.0*Py() + 0.125*Px(),
                    cst( 0.5 ) - Px() + 1.5*Pz() );
    else if constexpr ( Dim == 2 && Order == 1 )
        return vec( cst( 0.75 ) + Px() + Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() );
    else if constexpr ( Dim == 3 && Order == 1 )
        return vec( cst( 0.75 ) + Px() + Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px(),
                    cst( 0.5 ) + Pz() + Px()*Pz() );
    else if constexpr ( Dim == 2 && Order == 2 )
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py() );
    else if constexpr ( Dim == 3 && Order == 2 )
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py(),
                    cst( 0.5 ) + Pz() + Px()*Pz() + Py()*Pz() + Px()*Py()*Pz() );
    else if constexpr ( Dim == 2 )
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py() + Px()*Px()*Py()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py() + Px()*Px()*Px()*Py() );
    else
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py() + Px()*Px()*Py()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py() + Px()*Px()*Px()*Py(),
                    cst( 0.5 ) + Pz() + Px()*Pz() + Py()*Pz() + Px()*Py()*Pz() + Px()*Px()*Py()*Pz() );
}

template<int Dim, int Order>
void
checkRtStaticDynamicParityAndExactness( simplex_mesh_ptrtype<Dim> const& mesh )
{
    auto expr = rtExactExpr<Dim, Order>();
    auto rtStatic = RTh<Order>( mesh );
    auto rtDynamic = RTh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "RT dim=" << Dim << " order=" << Order )
    {
        BOOST_CHECK_EQUAL( rtDynamic->basis()->order(), Order );
        BOOST_CHECK_EQUAL( rtDynamic->basis()->localDofPerComponent(), rtStatic->basis()->localDofPerComponent() );
        checkReferenceFiniteElementParity( rtStatic, rtDynamic, mesh );
        checkOnParity( rtStatic, rtDynamic, mesh, expr, 1e-11, false );
        checkInterpolationExactness( rtStatic, mesh, expr, 1e-10 );
        checkInterpolationExactness( rtDynamic, mesh, expr, 1e-10 );
    }
}

template<int Dim, int Order>
void
checkBdmStaticDynamicParityAndExactness( simplex_mesh_ptrtype<Dim> const& mesh )
{
    auto expr = bdmExactExpr<Dim, Order>();
    auto bdmStatic = BDMh<Order>( mesh );
    auto bdmDynamic = BDMh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "BDM dim=" << Dim << " order=" << Order )
    {
        BOOST_CHECK_EQUAL( bdmDynamic->basis()->order(), Order );
        BOOST_CHECK_EQUAL( bdmDynamic->basis()->localDofPerComponent(), bdmStatic->basis()->localDofPerComponent() );
        checkReferenceFiniteElementParity( bdmStatic, bdmDynamic, mesh );
        checkOnParity( bdmStatic, bdmDynamic, mesh, expr, 1e-11, false );
        checkInterpolationExactness( bdmStatic, mesh, expr, 1e-10 );
        checkInterpolationExactness( bdmDynamic, mesh, expr, 1e-10 );
    }
}

template<int Dim, int Order>
void
checkNedelecStaticDynamicParityAndExactness( simplex_mesh_ptrtype<Dim> const& mesh )
{
    // Vector P_Order is a strict polynomial subspace of first-kind NED1_Order
    // and exercises every admitted scalar degree without relying on a
    // lower-order field for the higher-order cases.
    auto expr = rtExactExpr<Dim, Order>();
    auto nedStatic = Neh<Order>( mesh );
    auto nedDynamic = Neh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "NED1 dim=" << Dim << " order=" << Order )
    {
        BOOST_CHECK_EQUAL( nedDynamic->basis()->order(), Order );
        BOOST_CHECK_EQUAL( nedDynamic->basis()->localDofPerComponent(),
                           nedStatic->basis()->localDofPerComponent() );
        checkReferenceFiniteElementParity( nedStatic, nedDynamic, mesh );
        checkOnParity( nedStatic, nedDynamic, mesh, expr, 2e-10, false );
        checkInterpolationExactness( nedStatic, mesh, expr, 2e-9 );
        checkInterpolationExactness( nedDynamic, mesh, expr, 2e-9 );
    }
}

template<typename SpacePtrType>
void
checkNedelecOrderOneConstantLocalDofs( SpacePtrType const& space,
                                       simplex_mesh_ptrtype<3> const& mesh )
{
    std::array<double,3> const field{ 1.25, -0.75, 0.5 };
    auto u = space->element();
    u.on( _range = elements( mesh ),
          _expr = vec( cst( field[0] ), cst( field[1] ), cst( field[2] ) ) );

    auto const& fe = *space->basis();
    auto referenceFaceSurfaceJacobian = [&fe]( uint16_type face )
    {
        auto const vertices = fe.referenceConvex().faceVertices( face );
        Eigen::Vector3d tangent0, tangent1;
        for ( int component = 0; component < 3; ++component )
        {
            tangent0( component ) = vertices( component, 1 )-vertices( component, 0 );
            tangent1( component ) = vertices( component, 2 )-vertices( component, 0 );
        }
        return tangent0.cross( tangent1 ).norm();
    };
    IMGeneral<3, double, Simplex> faceQuadrature(
        Feel::detail::hcurlTriangularFaceMomentQuadratureOrder( 0 ) );
    auto const& element = unwrap_ref( *elements( mesh ).begin() );
    for ( uint16_type localDof = 0; localDof < fe.localDofPerComponent(); ++localDof )
    {
        auto const attachment = fe.dofAttachment( localDof );
        double expected = 0;
        if ( attachment.entityDim == 1 )
        {
            auto const p0 = element.point( element.eToP( attachment.entityId, 0 ) ).node();
            auto const p1 = element.point( element.eToP( attachment.entityId, 1 ) ).node();
            for ( int component = 0; component < 3; ++component )
                expected += field[component]*( p1( component )-p0( component ) );
        }
        else if ( attachment.entityDim == 2 )
        {
            auto const p0 = element.point( element.fToP( attachment.entityId, 0 ) ).node();
            auto const p1 = element.point( element.fToP(
                attachment.entityId, static_cast<uint16_type>( attachment.ordinal+1 ) ) ).node();
            double tangentValue = 0;
            for ( int component = 0; component < 3; ++component )
                tangentValue += field[component]*( p1( component )-p0( component ) );
            expected = ublas::sum( faceQuadrature.weights( attachment.entityId ) )*
                       tangentValue/referenceFaceSurfaceJacobian( attachment.entityId );
        }
        else
            continue;

        BOOST_TEST_CONTEXT( "NED1 order-one constant local dof " << localDof )
            BOOST_CHECK_SMALL( u.orientedLocalCoefficient( element.id(), localDof )-expected, 2e-10 );
    }
}

template<typename SpacePtrType>
void
checkNedelecOrderOneAffineLocalDofs( SpacePtrType const& space,
                                     simplex_mesh_ptrtype<3> const& mesh )
{
    auto field = []( std::array<double,3> const& x )
    {
        return std::array<double,3>{
            1.0+2.0*x[0]-0.25*x[1],
            -0.5+x[0]+3.0*x[1],
            0.25-x[1]+0.5*x[2] };
    };
    auto u = space->element();
    u.on( _range = elements( mesh ), _expr = rtExactExpr<3, 1>() );

    auto const& fe = *space->basis();
    auto referenceFaceSurfaceJacobian = [&fe]( uint16_type face )
    {
        auto const vertices = fe.referenceConvex().faceVertices( face );
        Eigen::Vector3d tangent0, tangent1;
        for ( int component = 0; component < 3; ++component )
        {
            tangent0( component ) = vertices( component, 1 )-vertices( component, 0 );
            tangent1( component ) = vertices( component, 2 )-vertices( component, 0 );
        }
        return tangent0.cross( tangent1 ).norm();
    };
    auto const& element = unwrap_ref( *elements( mesh ).begin() );
    auto mapToPhysical = [&element]( auto const& referencePoint )
    {
        std::array<double,3> physical{};
        auto const x0 = element.point( 0 ).node();
        for ( int component = 0; component < 3; ++component )
            physical[component] = x0( component );
        for ( int direction = 0; direction < 3; ++direction )
        {
            auto const xi = element.point( direction+1 ).node();
            const double barycentric = 0.5*( referencePoint( direction )+1.0 );
            for ( int component = 0; component < 3; ++component )
                physical[component] += barycentric*( xi( component )-x0( component ) );
        }
        return physical;
    };

    IMGeneral<3, double, Simplex> faceQuadrature(
        Feel::detail::hcurlTriangularFaceMomentQuadratureOrder( 0 ) );
    for ( uint16_type localDof = 0; localDof < fe.localDofPerComponent(); ++localDof )
    {
        auto const attachment = fe.dofAttachment( localDof );
        double expected = 0;
        if ( attachment.entityDim == 1 )
        {
            auto const physicalPoint = mapToPhysical( ublas::column( fe.points(), localDof ) );
            auto const value = field( physicalPoint );
            auto const p0 = element.point( element.eToP( attachment.entityId, 0 ) ).node();
            auto const p1 = element.point( element.eToP( attachment.entityId, 1 ) ).node();
            for ( int component = 0; component < 3; ++component )
                expected += value[component]*( p1( component )-p0( component ) );
        }
        else if ( attachment.entityDim == 2 )
        {
            auto const p0 = element.point( element.fToP( attachment.entityId, 0 ) ).node();
            auto const p1 = element.point( element.fToP(
                attachment.entityId, static_cast<uint16_type>( attachment.ordinal+1 ) ) ).node();
            const double surfaceJacobian =
                referenceFaceSurfaceJacobian( attachment.entityId );
            auto const& points = faceQuadrature.fpoints( attachment.entityId );
            for ( uint16_type q = 0; q < faceQuadrature.nPointsOnFace( attachment.entityId ); ++q )
            {
                auto const value = field( mapToPhysical( ublas::column( points, q ) ) );
                double tangentValue = 0;
                for ( int component = 0; component < 3; ++component )
                    tangentValue += value[component]*( p1( component )-p0( component ) );
                expected += faceQuadrature.weight( attachment.entityId, q )*
                            tangentValue/surfaceJacobian;
            }
        }
        else
            continue;

        BOOST_TEST_CONTEXT( "NED1 order-one affine local dof " << localDof )
            BOOST_CHECK_SMALL( u.orientedLocalCoefficient( element.id(), localDof )-expected, 2e-10 );
    }

    IMGeneral<3, double, Simplex> referenceQuadrature( 8 );
    auto const nodalValues = fe.evaluate( referenceQuadrature.points() );
    Eigen::MatrixXd sampleMatrix(
        3*referenceQuadrature.nPoints(), fe.localDofPerComponent() );
    Eigen::VectorXd pullbackValues( 3*referenceQuadrature.nPoints() );
    Eigen::Matrix3d jacobian;
    auto const x0 = element.point( 0 ).node();
    for ( int direction = 0; direction < 3; ++direction )
    {
        auto const xi = element.point( direction+1 ).node();
        for ( int component = 0; component < 3; ++component )
            jacobian( component, direction ) = 0.5*( xi( component )-x0( component ) );
    }
    for ( int q = 0; q < referenceQuadrature.nPoints(); ++q )
    {
        auto const physicalValue = field(
            mapToPhysical( ublas::column( referenceQuadrature.points(), q ) ) );
        Eigen::Vector3d physical( physicalValue[0], physicalValue[1], physicalValue[2] );
        pullbackValues.segment<3>( 3*q ) = jacobian.transpose()*physical;
        for ( uint16_type basis = 0; basis < fe.localDofPerComponent(); ++basis )
            for ( int component = 0; component < 3; ++component )
                sampleMatrix( 3*q+component, basis ) = nodalValues( 3*basis+component, q );
    }
    Eigen::VectorXd const expectedLocal =
        sampleMatrix.colPivHouseholderQr().solve( pullbackValues );
    BOOST_CHECK_SMALL( ( sampleMatrix*expectedLocal-pullbackValues ).norm(), 2e-10 );
    std::vector<bool> visitedLocalDof( fe.localDofPerComponent(), false );
    for ( auto const& localEntry : space->dof()->localDof( element.id() ) )
    {
        const uint16_type localDof = localEntry.first.localDof();
        BOOST_REQUIRE_LT( localDof, visitedLocalDof.size() );
        BOOST_CHECK( !visitedLocalDof[localDof] );
        visitedLocalDof[localDof] = true;
    }
    BOOST_CHECK( std::ranges::all_of( visitedLocalDof, std::identity{} ) );
    for ( uint16_type localDof = 0; localDof < fe.localDofPerComponent(); ++localDof )
        BOOST_TEST_CONTEXT( "NED1 order-one pullback local dof " << localDof )
            BOOST_CHECK_SMALL(
                u.orientedLocalCoefficient( element.id(), localDof )-expectedLocal( localDof ),
                2e-10 );

    typename std::remove_cvref_t<decltype( u )>::matrix_node_type physicalPoints( 3, 1 );
    auto const physicalPoint = mapToPhysical(
        ublas::column( referenceQuadrature.points(), 0 ) );
    for ( int component = 0; component < 3; ++component )
        physicalPoints( component, 0 ) = physicalPoint[component];
    typename std::remove_cvref_t<decltype( u )>::id_array_type evaluated( boost::extents[1] );
    evaluated[0].setZero();
    u.idInterpolate(
        physicalPoints, evaluated, false,
        typename std::remove_cvref_t<decltype( u )>::matrix_node_type{} );
    auto const expectedValue = field( physicalPoint );
    for ( int component = 0; component < 3; ++component )
        BOOST_CHECK_SMALL( evaluated[0]( component, 0 )-expectedValue[component], 2e-10 );
}

} // namespace

BOOST_AUTO_TEST_SUITE( hdiv_runtime_order_parity_suite )

BOOST_AUTO_TEST_CASE( rt_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh<2>( "rt-runtime-parity" );
    auto expr = vec( cst( 1.0 ) + Px(), cst( -0.5 ) + 2.0 * Py() );

    auto rt0Static = RTh<0>( mesh );
    auto rt0Runtime = RTh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( rt0Static, rt0Runtime, mesh, expr, 1e-12 );
}

BOOST_AUTO_TEST_CASE( rt_dynamic_low_order_space_parity )
{
    auto mesh2d = makeOneElementSimplexMesh<2>( "rt-dynamic-low-order-parity-2d" );
    auto mesh3d = makeOneElementSimplexMesh<3>( "rt-dynamic-low-order-parity-3d" );
    checkRtStaticDynamicParityAndExactness<2, 0>( mesh2d );
    checkRtStaticDynamicParityAndExactness<2, 1>( mesh2d );
    checkRtStaticDynamicParityAndExactness<2, 2>( mesh2d );
    checkRtStaticDynamicParityAndExactness<2, 3>( mesh2d );
    checkRtStaticDynamicParityAndExactness<3, 0>( mesh3d );
    checkRtStaticDynamicParityAndExactness<3, 1>( mesh3d );
    checkRtStaticDynamicParityAndExactness<3, 2>( mesh3d );
    checkRtStaticDynamicParityAndExactness<3, 3>( mesh3d );
}

BOOST_AUTO_TEST_CASE( rt_static_higher_order_simplex_space_layout )
{
    auto mesh = makeOneElementSimplexMesh<2>( "rt-static-simplex-higher-order-layout" );

    auto rt1Static = RTh<1>( mesh );
    auto rt2Static = RTh<2>( mesh );

    BOOST_CHECK_EQUAL( rt1Static->basis()->localDofPerComponent(), 8 );
    BOOST_CHECK_EQUAL( rt2Static->basis()->localDofPerComponent(), 15 );
    BOOST_CHECK_EQUAL( rt1Static->nLocalDof(), 8 );
    BOOST_CHECK_EQUAL( rt2Static->nLocalDof(), 15 );
    BOOST_CHECK_EQUAL( rt1Static->basis()->nDofPerEdge, 2 );
    BOOST_CHECK_EQUAL( rt2Static->basis()->nDofPerEdge, 3 );
    BOOST_CHECK_EQUAL( rt1Static->basis()->nDofPerFace, 2 );
    BOOST_CHECK_EQUAL( rt2Static->basis()->nDofPerFace, 6 );

    auto rt1Dynamic = RTh<Dynamic>( mesh, RuntimeOrder{ 1 } );
    auto rt2Dynamic = RTh<Dynamic>( mesh, RuntimeOrder{ 2 } );

    BOOST_CHECK_EQUAL( rt1Dynamic->basis()->localDofPerComponent(), 8 );
    BOOST_CHECK_EQUAL( rt2Dynamic->basis()->localDofPerComponent(), 15 );
    BOOST_CHECK_EQUAL( rt1Dynamic->nLocalDof(), 8 );
    BOOST_CHECK_EQUAL( rt2Dynamic->nLocalDof(), 15 );
    BOOST_CHECK_EQUAL( rt1Dynamic->basis()->dofPerEdge(), 2 );
    BOOST_CHECK_EQUAL( rt2Dynamic->basis()->dofPerEdge(), 3 );
    BOOST_CHECK_EQUAL( rt1Dynamic->basis()->dofPerFace(), 2 );
    BOOST_CHECK_EQUAL( rt2Dynamic->basis()->dofPerFace(), 6 );
}

BOOST_AUTO_TEST_CASE( bdm_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh<2>( "bdm-runtime-parity" );
    auto mesh3d = makeOneElementSimplexMesh<3>( "bdm-runtime-parity-3d" );
    auto expr = bdmExactExpr<2, 0>();

    auto bdm0Static = BDMh<0>( mesh );
    auto bdm0Runtime = BDMh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( bdm0Static, bdm0Runtime, mesh, expr, 1e-12 );
    checkBdmStaticDynamicParityAndExactness<2, 0>( mesh );
    checkBdmStaticDynamicParityAndExactness<2, 1>( mesh );
    checkBdmStaticDynamicParityAndExactness<2, 2>( mesh );
    checkBdmStaticDynamicParityAndExactness<2, 3>( mesh );
    checkBdmStaticDynamicParityAndExactness<3, 0>( mesh3d );
    checkBdmStaticDynamicParityAndExactness<3, 1>( mesh3d );
    checkBdmStaticDynamicParityAndExactness<3, 2>( mesh3d );
    checkBdmStaticDynamicParityAndExactness<3, 3>( mesh3d );
}

BOOST_AUTO_TEST_CASE( neh_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh<2>( "neh-runtime-parity" );
    auto mesh3d = makeOneElementSimplexMesh<3>( "neh-runtime-parity-3d" );
    auto expr = vec( -Py(), Px() );

    auto neh0Static = Neh<0>( mesh );
    auto neh0Runtime = Neh<0>( mesh, RuntimeOrder{ 0 } );
    auto neh0Dynamic = Neh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto neh1Static = Neh<1>( mesh );
    auto neh1Dynamic = Neh<Dynamic>( mesh, RuntimeOrder{ 1 } );
    auto neh2Static = Neh<2>( mesh );
    auto neh2Dynamic = Neh<Dynamic>( mesh, RuntimeOrder{ 2 } );
    checkOnParity( neh0Static, neh0Runtime, mesh, expr, 1e-12 );
    checkOnParity( neh0Static, neh0Dynamic, mesh, expr, 1e-12 );
    checkOnParity( neh1Static, neh1Dynamic, mesh, expr, 1e-11 );
    checkOnParity( neh2Static, neh2Dynamic, mesh, expr, 1e-11 );
    BOOST_CHECK_EQUAL( neh0Dynamic->basis()->order(), 0 );
    BOOST_CHECK_EQUAL( neh0Dynamic->basis()->localDofCount(), neh0Static->basis()->localDofCount() );
    BOOST_CHECK_EQUAL( neh1Dynamic->basis()->order(), 1 );
    BOOST_CHECK_EQUAL( neh1Dynamic->basis()->localDofCount(), neh1Static->basis()->localDofCount() );
    BOOST_CHECK_EQUAL( neh2Dynamic->basis()->order(), 2 );
    BOOST_CHECK_EQUAL( neh2Dynamic->basis()->localDofCount(), neh2Static->basis()->localDofCount() );

    checkNedelecStaticDynamicParityAndExactness<2, 0>( mesh );
    checkNedelecStaticDynamicParityAndExactness<2, 1>( mesh );
    checkNedelecStaticDynamicParityAndExactness<2, 2>( mesh );
    checkNedelecStaticDynamicParityAndExactness<2, 3>( mesh );
    checkNedelecStaticDynamicParityAndExactness<3, 0>( mesh3d );
    checkNedelecStaticDynamicParityAndExactness<3, 1>( mesh3d );
    checkNedelecStaticDynamicParityAndExactness<3, 2>( mesh3d );
    checkNedelecStaticDynamicParityAndExactness<3, 3>( mesh3d );

    checkNedelecOrderOneConstantLocalDofs( Neh<1>( mesh3d ), mesh3d );
    checkNedelecOrderOneConstantLocalDofs(
        Neh<Dynamic>( mesh3d, RuntimeOrder{ 1 } ), mesh3d );
    checkNedelecOrderOneAffineLocalDofs( Neh<1>( mesh3d ), mesh3d );
    checkNedelecOrderOneAffineLocalDofs(
        Neh<Dynamic>( mesh3d, RuntimeOrder{ 1 } ), mesh3d );

    checkBoundaryTangentialInterpolation(
        Neh<1>( mesh3d ), Neh<Dynamic>( mesh3d, RuntimeOrder{ 1 } ),
        mesh3d, rtExactExpr<3, 1>(), 2e-9 );
    checkBoundaryTangentialInterpolation(
        Neh<2>( mesh3d ), Neh<Dynamic>( mesh3d, RuntimeOrder{ 2 } ),
        mesh3d, rtExactExpr<3, 2>(), 2e-9 );
    checkBoundaryTangentialInterpolation(
        Neh<3>( mesh3d ), Neh<Dynamic>( mesh3d, RuntimeOrder{ 3 } ),
        mesh3d, rtExactExpr<3, 3>(), 2e-9 );

    auto mesh3dMulti = unitCube( 0.8 );
    checkNedelecStaticDynamicParityAndExactness<3, 1>( mesh3dMulti );
}

BOOST_AUTO_TEST_SUITE_END()
