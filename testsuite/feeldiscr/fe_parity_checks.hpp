/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_TESTSUITE_FE_PARITY_CHECKS_HPP
#define FEELPP_TESTSUITE_FE_PARITY_CHECKS_HPP 1

#include <algorithm>
#include <cmath>

namespace Feel::test_detail
{

/** Compare every field of two local DoF-transform descriptors. */
template<typename StaticTransform, typename RuntimeTransform>
void
checkDofTransformParity( StaticTransform const& staticTransform,
                         RuntimeTransform const& runtimeTransform )
{
    BOOST_CHECK_EQUAL( static_cast<int>( staticTransform.kind ),
                       static_cast<int>( runtimeTransform.kind ) );
    BOOST_CHECK_EQUAL( staticTransform.sign, runtimeTransform.sign );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        staticTransform.permutation.begin(), staticTransform.permutation.end(),
        runtimeTransform.permutation.begin(), runtimeTransform.permutation.end() );
    auto const staticMatrixRow = staticTransform.matrixRowView();
    auto const runtimeMatrixRow = runtimeTransform.matrixRowView();
    BOOST_REQUIRE_EQUAL( staticMatrixRow.size(), runtimeMatrixRow.size() );
    for ( std::size_t column = 0; column < staticMatrixRow.size(); ++column )
        BOOST_CHECK_SMALL(
            staticMatrixRow[column]-runtimeMatrixRow[column], 1e-12 );
}

/** Compare all local FE transforms on the first locally owned mesh element. */
template<typename StaticSpacePtrType, typename RuntimeSpacePtrType,
         typename MeshPtrType>
void
checkFiniteElementTransformParity( StaticSpacePtrType const& staticSpace,
                                   RuntimeSpacePtrType const& runtimeSpace,
                                   MeshPtrType const& mesh )
{
    auto const& staticFe = *staticSpace->basis();
    auto const& runtimeFe = *runtimeSpace->basis();
    for ( auto const& elementWrapper : elements( mesh ) )
    {
        auto const& element = unwrap_ref( elementWrapper );
        for ( uint16_type localDof = 0;
              localDof < staticFe.localDofCount(); ++localDof )
            checkDofTransformParity(
                staticFe.dofTransform( element, localDof ),
                runtimeFe.dofTransform( element, localDof ) );
        return;
    }
}

/** Compare complete assembled operators without modifying either input form. */
template<typename StaticForm, typename RuntimeForm>
void
checkAssembledMatrixParity( StaticForm const& staticForm,
                            RuntimeForm const& runtimeForm,
                            double tolerance = 1e-11 )
{
    staticForm.matrixPtr()->close();
    runtimeForm.matrixPtr()->close();
    BOOST_REQUIRE_EQUAL( staticForm.matrixPtr()->size1(),
                         runtimeForm.matrixPtr()->size1() );
    BOOST_REQUIRE_EQUAL( staticForm.matrixPtr()->size2(),
                         runtimeForm.matrixPtr()->size2() );
    auto const scale = std::max(
        { 1., std::abs( staticForm.matrixPtr()->l1Norm() ),
          std::abs( runtimeForm.matrixPtr()->l1Norm() ),
          std::abs( staticForm.matrixPtr()->linftyNorm() ),
          std::abs( runtimeForm.matrixPtr()->linftyNorm() ) } );
    auto difference = runtimeForm.matrixPtr()->clone();
    difference->addMatrix( -1., staticForm.matrixPtr() );
    difference->close();
    BOOST_CHECK_SMALL( difference->l1Norm()/scale, tolerance );
    BOOST_CHECK_SMALL( difference->linftyNorm()/scale, tolerance );
}

/** Compare complete assembled right-hand sides without modifying either form. */
template<typename StaticForm, typename RuntimeForm>
void
checkAssembledVectorParity( StaticForm const& staticForm,
                            RuntimeForm const& runtimeForm,
                            double tolerance = 1e-11 )
{
    staticForm.vectorPtr()->close();
    runtimeForm.vectorPtr()->close();
    auto staticVector = [&]
    {
        if constexpr ( requires { staticForm.vectorPtr()->getVector(); } )
            return staticForm.vectorPtr()->getVector();
        else
            return staticForm.vectorPtr();
    }();
    auto runtimeVector = [&]
    {
        if constexpr ( requires { runtimeForm.vectorPtr()->getVector(); } )
            return runtimeForm.vectorPtr()->getVector();
        else
            return runtimeForm.vectorPtr();
    }();
    BOOST_REQUIRE_EQUAL( staticVector->size(), runtimeVector->size() );
    auto const scale = std::max(
        { 1., std::abs( staticVector->l1Norm() ),
          std::abs( runtimeVector->l1Norm() ),
          std::abs( staticVector->linftyNorm() ),
          std::abs( runtimeVector->linftyNorm() ) } );
    auto difference = runtimeVector->clone();
    difference->add( 1., runtimeVector );
    difference->add( -1., staticVector );
    difference->close();
    BOOST_CHECK_SMALL( difference->l1Norm()/scale, tolerance );
    BOOST_CHECK_SMALL( difference->linftyNorm()/scale, tolerance );
}

/** Compare block operators that use equivalent but distinct product maps. */
template<typename StaticForm, typename RuntimeForm>
void
checkAssembledMatrixNormParity( StaticForm const& staticForm,
                                RuntimeForm const& runtimeForm,
                                double tolerance = 1e-11 )
{
    staticForm.matrixPtr()->close();
    runtimeForm.matrixPtr()->close();
    BOOST_REQUIRE_EQUAL( staticForm.matrixPtr()->size1(),
                         runtimeForm.matrixPtr()->size1() );
    BOOST_REQUIRE_EQUAL( staticForm.matrixPtr()->size2(),
                         runtimeForm.matrixPtr()->size2() );
    auto const scale = std::max(
        { 1., std::abs( staticForm.matrixPtr()->l1Norm() ),
          std::abs( runtimeForm.matrixPtr()->l1Norm() ),
          std::abs( staticForm.matrixPtr()->linftyNorm() ),
          std::abs( runtimeForm.matrixPtr()->linftyNorm() ) } );
    BOOST_CHECK_SMALL(
        ( staticForm.matrixPtr()->l1Norm()-
          runtimeForm.matrixPtr()->l1Norm() )/scale,
        tolerance );
    BOOST_CHECK_SMALL(
        ( staticForm.matrixPtr()->linftyNorm()-
          runtimeForm.matrixPtr()->linftyNorm() )/scale,
        tolerance );
}

/** Require equivalent convergence behavior from static/runtime solves. */
template<typename StaticSolveData, typename RuntimeSolveData>
void
checkSolveParity( StaticSolveData const& staticSolve,
                  RuntimeSolveData const& runtimeSolve )
{
    BOOST_REQUIRE( staticSolve.isConverged() );
    BOOST_REQUIRE( runtimeSolve.isConverged() );
    auto const staticIterations = staticSolve.nIterations();
    auto const runtimeIterations = runtimeSolve.nIterations();
    BOOST_CHECK_LE( staticIterations > runtimeIterations
                        ? staticIterations-runtimeIterations
                        : runtimeIterations-staticIterations,
                    1 );
    auto const residualScale = std::max(
        { 1., std::abs( staticSolve.residual() ),
          std::abs( runtimeSolve.residual() ) } );
    BOOST_CHECK_SMALL(
        ( staticSolve.residual()-runtimeSolve.residual() )/residualScale,
        1e-8 );
}

} // namespace Feel::test_detail

#endif // FEELPP_TESTSUITE_FE_PARITY_CHECKS_HPP
