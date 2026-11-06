// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

/** @file test_polynomial_order_operations.cpp
 * @brief Mathematical identities across polynomial order and result policies.
 */
#define BOOST_TEST_MODULE polynomial_order_operations
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/operations.hpp>
#include <boost/mp11.hpp>
#include <boost/test/data/test_case.hpp>
#include <optional>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Require equal matrix extents and uniformly small entrywise error. */
template<class A, class B>
void checkMatrix( A const& a, B const& b )
{
    BOOST_REQUIRE_EQUAL( a.size1(), b.size1() );
    BOOST_REQUIRE_EQUAL( a.size2(), b.size2() );
    for ( std::size_t i = 0; i < a.size1(); ++i )
        for ( std::size_t j = 0; j < a.size2(); ++j )
            BOOST_CHECK_SMALL( a(i,j) - b(i,j), 2e-10 );
}

/** @brief An analytic polynomial and its derivatives, independent of basis kernels. */
double analytic( int degree, double x, double y, unsigned dx = 0, unsigned dy = 0 )
{
    if ( dx + dy == 0 )
        return 1 + (degree >= 1 ? 0.5*x - 0.75*y : 0)
                 + (degree >= 2 ? x*x + x*y + 2*y*y : 0)
                 + (degree >= 3 ? x*x*x + x*y*y : 0);
    if ( dx == 1 && dy == 0 )
        return (degree >= 1 ? 0.5 : 0) + (degree >= 2 ? 2*x+y : 0)
               + (degree >= 3 ? 3*x*x+y*y : 0);
    if ( dx == 0 && dy == 1 )
        return (degree >= 1 ? -0.75 : 0) + (degree >= 2 ? x+4*y : 0)
               + (degree >= 3 ? 2*x*y : 0);
    if ( dx == 2 ) return (degree >= 2 ? 2 : 0) + (degree >= 3 ? 6*x : 0);
    if ( dy == 2 ) return (degree >= 2 ? 4 : 0) + (degree >= 3 ? 2*x : 0);
    return (degree >= 2 ? 1 : 0) + (degree >= 3 ? 2*y : 0);
}

/** @brief Exact integral from monomial moments on the reference cell. */
template<class Cell>
double analyticIntegral( int degree )
{
    constexpr double d=Cell::nDim;
    double mean=1;
    double volume=std::pow(2.0,d);
    if constexpr(Cell::is_hypercube)
    {
        if(degree>=2) mean+=(Cell::nDim>1 ? 3.0 : 1.0)/3.0;
    }
    else
    {
        volume/=std::tgamma(d+1);
        double x=2/(d+1)-1;
        double xx=1-4/(d+1)+8/((d+1)*(d+2));
        double xy=1-4/(d+1)+4/((d+1)*(d+2));
        double xxx=-1+6/(d+1)-24/((d+1)*(d+2))+48/((d+1)*(d+2)*(d+3));
        double xyy=-1+6/(d+1)-16/((d+1)*(d+2))+16/((d+1)*(d+2)*(d+3));
        if(degree>=1) mean+=(Cell::nDim>1 ? -0.25 : 0.5)*x;
        if(degree>=2) mean+=xx+(Cell::nDim>1 ? xy+2*xx : 0);
        if(degree>=3) mean+=xxx+(Cell::nDim>1 ? xyy : 0);
    }
    return volume*mean;
}

/** @brief A static/runtime order and reference-cell test configuration. */
template<class Order, class Cell>
struct OperationCase
{
    /** Authoritative order policy under test. */
    static constexpr int order = Order::value;
    /** Reference-cell geometry. */
    using cell_type = Cell;
};
/** Configurations exercise the same client code with both order policies. */
using operation_cases = boost::mp11::mp_product<OperationCase,
    boost::mp11::mp_list_c<int, 0, 1, 2, 3, Dynamic>,
    boost::mp11::mp_list<Simplex<1>, Hypercube<1>, Simplex<2>, Hypercube<2>, Simplex<3>, Hypercube<3>>>;

/** @brief Check extraction and differential identities against an analytic polynomial. */
template<class Case>
void checkOperations( int degree )
{
    using cell_type = typename Case::cell_type;
    constexpr unsigned dim = cell_type::nDim;
    using fe_type = typename Lagrange<Case::order, Scalar, Continuous, PointSetEquiSpaced>::
        template Apply<dim, dim, double, cell_type>::type;
    auto fe = makeOrderAware<fe_type>( RuntimeOrder{degree} );
    using matrix_type = typename fe_type::matrix_type;
    using polynomial_type = typename fe_type::polynomial_type;
    static_assert( polynomial_type::order_spec == Case::order );
    static_assert( fe_type::primal_space_type::order_spec == Case::order );
    static_assert( fe_type::dual_space_type::order_spec == Case::order );
    static_assert( orderIsDynamic<typename fe_type::pointset_type> == fe_type::is_order_dynamic );
    BOOST_CHECK_EQUAL(fe.dual().pointSet().nPoints(),fe.localDof());
    if constexpr(fe_type::is_order_dynamic)
        BOOST_CHECK_EQUAL(fe.dual().pointSet().order(),degree);
    matrix_type nodalWeights( 1, fe.points().size2() );
    for ( std::size_t q = 0; q < nodalWeights.size2(); ++q )
        nodalWeights(0,q) = analytic( degree, fe.points()(0,q), dim > 1 ? fe.points()(1,q) : 0.0 );
    polynomial_type polynomial( fe, nodalWeights );
    matrix_type points(dim,3);
    for(unsigned d=0;d<dim;++d)
        for(unsigned q=0;q<3;++q) points(d,q)=-0.81;
    points(0,0)=-0.73; points(0,1)=-0.21; points(0,2)=-0.44;
    if constexpr(dim>1)
    {
        points(1,0)=-0.61; points(1,1)=-0.56; points(1,2)=-0.18;
    }
    auto y=[&points](std::size_t q) { return dim>1 ? points(1,q) : 0.0; };
    auto values = polynomial( points );
    auto set = polynomial.toSet();
    BOOST_CHECK_EQUAL( polynomial.order(), degree );
    BOOST_CHECK_EQUAL( set.order(), degree );
    checkMatrix( values, set.evaluate( points ) );
    auto basisValues = fe.evaluate( points );
    auto extracted = fe.polynomial(0);
    checkMatrix( extracted( points ), matrix_type( ublas::project( basisValues,
        ublas::range(0,1), ublas::range(0,points.size2()) ) ) );
    auto subset = fe.polynomials( {0} );
    checkMatrix( subset.evaluate( points ), extracted( points ) );
    checkMatrix( fe.polynomialsRange(0,1).evaluate(points), subset.evaluate(points) );
    checkMatrix( fe.polynomialsUpToDimension(1).evaluate(points), subset.evaluate(points) );
    BOOST_CHECK_THROW( fe.polynomial( fe.polynomialDimension() ), std::out_of_range );
    BOOST_CHECK_THROW( fe.polynomialsRange(1,0), std::out_of_range );

    auto gradient = set.gradient();
    BOOST_CHECK_EQUAL( gradient.order(), degree );
    BOOST_CHECK_EQUAL( gradient.polynomialDimension(), 1 );
    auto gradientValues = gradient.evaluate(points);
    for ( unsigned i = 0; i < dim; ++i )
    {
        auto const& op = polynomial.d(i);
        BOOST_REQUIRE_EQUAL( op.size1(), polynomial.coeff().size2() );
        BOOST_REQUIRE_EQUAL( op.size2(), polynomial.coeff().size2() );
        BOOST_CHECK( &op == &polynomial.d(i) );
        auto first = polynomial.derivate(i,points);
        checkMatrix( first, polynomial.derivative(i)(points) );
        checkMatrix( first, set.derivate(i).evaluate(points) );
        for ( std::size_t q = 0; q < points.size2(); ++q )
        {
            BOOST_CHECK_SMALL( first(0,q) - (i<2 ? analytic(degree,points(0,q),y(q),i==0,i==1) : 0.0), 2e-10 );
            BOOST_CHECK_SMALL( gradientValues(i,q) - first(0,q), 2e-10 );
        }
        for ( unsigned j = 0; j < dim; ++j )
        {
            auto second = polynomial.derivate(i,j,points);
            checkMatrix( second, polynomial.derivative(i).derivative(j)(points) );
            checkMatrix( second, set.derivate(i,j,points) );
            checkMatrix( polynomial.d(i,j), set.d(i,j) );
            for ( std::size_t q = 0; q < points.size2(); ++q )
                BOOST_CHECK_SMALL( second(0,q)-(i<2 && j<2 ? analytic(degree,points(0,q),y(q),(i==0)+(j==0),(i==1)+(j==1)) : 0.0), 2e-10 );
        }
    }
    for ( std::size_t q = 0; q < points.size2(); ++q )
    {
        BOOST_CHECK_SMALL( values(0,q)-analytic(degree,points(0,q),y(q)), 2e-10 );
        typename polynomial_type::node_type point(dim);
        for(unsigned d=0;d<dim;++d) point(d)=points(d,q);
        BOOST_CHECK_SMALL( polynomial(point)(0,0)-values(0,q), 2e-10 );
        BOOST_CHECK_SMALL( set.evaluate(point)(0,0)-values(0,q), 2e-10 );
        BOOST_CHECK_SMALL( fe(0,point)-extracted(point)(0,0), 2e-10 );
    }
    auto copied = polynomial;
    auto moved = std::move(copied);
    checkMatrix( moved(points), values );
    polynomial_type assigned;
    assigned = moved;
    checkMatrix( assigned(points), values );
    checkMatrix( (polynomial-assigned)(points), matrix_type(ublas::zero_matrix<double>(1,points.size2())) );
    BOOST_CHECK_SMALL(set.firstMoments()(0,0)-analyticIntegral<cell_type>(degree),2e-10);
    BOOST_CHECK_SMALL(Feel::integrate(polynomial)-analyticIntegral<cell_type>(degree),2e-10);
    checkMatrix(Feel::dx(polynomial)(points),polynomial.derivative(0)(points));
    if constexpr(dim>1) checkMatrix(Feel::dy(set).evaluate(points),set.derivative(1).evaluate(points));
    if constexpr(dim>2) checkMatrix(Feel::dz(polynomial)(points),polynomial.derivative(2)(points));
    checkMatrix(Feel::gradient(polynomial)(points),gradientValues);
    checkMatrix(Feel::hessian(polynomial)(points),Feel::hessian(set).evaluate(points));
    auto span=Feel::unite(set,set);
    BOOST_CHECK_EQUAL(span.order(),degree);
    BOOST_CHECK_EQUAL(span.polynomialDimension(),1);
    double normSquared=0;
    for(std::size_t k=0;k<span.coeff().size2();++k) normSquared+=span.coeff()(0,k)*span.coeff()(0,k);
    BOOST_CHECK_SMALL(normSquared-1.0,2e-10);
    if constexpr(dim==3)
    {
        matrix_type nodalVector(3,fe.points().size2());
        for(std::size_t q=0;q<nodalVector.size2();++q)
        {
            auto const x=fe.points()(0,q), y=fe.points()(1,q), z=fe.points()(2,q);
            nodalVector(0,q)=degree>0 ? x+y : 1;
            nodalVector(1,q)=degree>0 ? 2*y+z : 2;
            nodalVector(2,q)=degree>0 ? 3*z+x : 3;
        }
        matrix_type vectorCoefficients=ublas::prod(nodalVector,fe.coeff());
        Polynomial<typename fe_type::primal_space_type,Vectorial,matrix_type,Case::order> vector(fe,vectorCoefficients,true);
        auto divergence=Feel::divergence(vector)(points);
        auto curl=Feel::curl(vector)(points);
        auto transposedCurl=Feel::curlt(vector.toSet()).evaluate(points);
        for(std::size_t q=0;q<points.size2();++q)
        {
            BOOST_CHECK_SMALL(divergence(0,q)-(degree>0 ? 6.0 : 0.0),2e-10);
            for(unsigned c=0;c<3;++c)
            {
                BOOST_CHECK_SMALL(curl(c,q)+(degree>0 ? 1.0 : 0.0),2e-10);
                BOOST_CHECK_SMALL(transposedCurl(c,q)-(degree>0 ? 1.0 : 0.0),2e-10);
            }
        }
    }

    using set_type = decltype(set);
    auto owner = std::make_shared<set_type>(set);
    typename set_type::PreCompute precompute(owner,points);
    for ( std::size_t q=0; q<points.size2(); ++q )
        BOOST_CHECK_SMALL( precompute.phi(0,0,0,q)-values(0,q), 2e-10 );
}
}

/** @test Both order policies support identical mathematical client code. */
BOOST_AUTO_TEST_CASE_TEMPLATE( polynomialResultsPreserveOrderAndIdentities, Case, operation_cases )
{
    if constexpr ( Case::order == Dynamic )
        for ( int degree : {0,1,2,3} )
        {
            BOOST_TEST_CONTEXT("runtime degree=" << degree) { checkOperations<Case>(degree); }
        }
    else checkOperations<Case>(Case::order);
}

/** @test Static conversions retain their degree even with a different evaluator seed. */
BOOST_AUTO_TEST_CASE( polynomialConversionRetainsStaticDegree )
{
    using modal_type = Feel::detail::OrthonormalPolynomialSet<2,Dynamic,2>;
    modal_type modal{RuntimeOrder{3}};
    PolynomialSet<typename modal_type::basis_type,Scalar,Dynamic> dynamic(modal);
    auto fixed = dynamic.toStatic<3>();
    BOOST_REQUIRE(fixed);
    static_assert(decltype(fixed)::value_type::order_spec == 3);
    static_assert(std::is_empty_v<Feel::detail::PolynomialBasisState<typename modal_type::basis_type,1>>);
    BOOST_CHECK_EQUAL(fixed->order(),3);
    BOOST_CHECK_EQUAL(fixed->d(0).size1(),10);
    BOOST_CHECK(!dynamic.toStatic<2>());
    auto points = PointSetEquiSpaced<Simplex<2>,3,double>().points();
    checkMatrix(fixed->evaluate(points),dynamic.evaluate(points));
    checkMatrix(fixed->polynomial(0)(points),dynamic.polynomial(0)(points));
    auto roundtrip=fixed->toDynamic();
    BOOST_CHECK(decltype(roundtrip)::is_order_dynamic);
    checkMatrix(roundtrip.evaluate(points),dynamic.evaluate(points));
    using fixed_type=PolynomialSet<typename modal_type::basis_type,Scalar,2>;
    BOOST_CHECK_THROW((void)fixed_type{dynamic},std::invalid_argument);
    auto retained=dynamic;
    auto const* retainedOperator=&retained.d(0);
    dynamic.setOrder(1);
    BOOST_CHECK_EQUAL(dynamic.d(0).size1(),3);
    BOOST_CHECK_EQUAL(retainedOperator->size1(),10);
    BOOST_CHECK_EQUAL(retained.order(),3);
}

namespace
{
/** @brief Verify component extraction, modal layout and Jacobian row conventions. */
template<template<uint16_type> class Shape, template<int,int,int> class Cell>
void checkComponents()
{
    using modal_type=Feel::detail::OrthonormalPolynomialSet<2,Dynamic,2,Shape,double,0,Cell>;
    for(int degree : {0,1,3})
    {
        modal_type modal{RuntimeOrder{degree}};
        auto const count=modal.localDof();
        BOOST_CHECK_EQUAL(modal.localDofCount(),modal.nComponents*count);
        BOOST_CHECK_EQUAL(modal.localDofCountOnFacet(),0);
        BOOST_CHECK_EQUAL(modal.localDofCountOnEntity(2,0,true),count);
        for(uint16_type c=0;c<modal.nComponents;++c)
            for(uint16_type parent=0;parent<count;++parent)
            {
                auto id=modal.localDofId(parent,c);
                auto layout=modal.localDofLayout(id);
                BOOST_CHECK_EQUAL(layout.component,c);
                BOOST_CHECK_EQUAL(layout.parentLocalDofId,parent);
                BOOST_CHECK_EQUAL(layout.attachment.entityDim,2);
                BOOST_CHECK(!modal.dofHasRepresentativePoint(id));
            }
        BOOST_CHECK_THROW(modal.localDofId(count),std::out_of_range);
        auto one=modal.polynomial(0);
        typename modal_type::matrix_type coefficients(modal.nComponents,count);
        for(std::size_t c=0;c<coefficients.size1();++c)
            for(std::size_t k=0;k<count;++k) coefficients(c,k)=(c+1)*0.17+(k+1)*0.03;
        one.setCoefficient(coefficients,true);
        auto points=PointSetEquiSpaced<Cell<2,1,2>,2,double>().points();
        auto values=one(points);
        auto set=one.toSet();
        for(uint16_type c=0;c<modal.nComponents;++c)
        {
            auto component=one[c];
            BOOST_CHECK_EQUAL(component.order(),degree);
            checkMatrix(component(points),typename modal_type::matrix_type(ublas::project(values,
                ublas::range(c,c+1),ublas::range(0,values.size2()))));
            if constexpr(modal_type::is_vectorial)
                checkMatrix(set[c].evaluate(points),component(points));
            else
                checkMatrix(set(c/modal.nComponents2,c%modal.nComponents2).evaluate(points),component(points));
        }
        if constexpr(modal_type::is_vectorial)
        {
            auto jacobian=set.gradient();
            BOOST_CHECK_EQUAL(jacobian.polynomialDimension(),1);
            auto entries=jacobian.evaluate(points);
            for(uint16_type d=0;d<2;++d)
            {
                auto derivative=one.derivate(d,points);
                for(uint16_type c=0;c<modal.nComponents;++c)
                    for(std::size_t q=0;q<points.size2();++q)
                        BOOST_CHECK_SMALL(entries(c*2+d,q)-derivative(c,q),2e-10);
            }
        }
    }
}
}

/** @test Vector and tensor results retain components and cell-owned modal DOFs. */
BOOST_AUTO_TEST_CASE(modalComponentsAndLocalLayout)
{
    checkComponents<Vectorial,Simplex>();
    checkComponents<Vectorial,Hypercube>();
    checkComponents<Tensor2,Simplex>();
    checkComponents<Tensor2,Hypercube>();
}
