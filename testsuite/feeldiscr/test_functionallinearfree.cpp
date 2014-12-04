/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-04-29

  Copyright (C) 2013 Feel++ Consortium

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file test_functionallinearfree.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-29
*/

#define BOOST_TEST_MODULE test_functionallinearfree
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <feel/feeldiscr/fsfunctionallinearfree.hpp>
#include <feel/feeldiscr/fsfunctionallinearcomposite.hpp>


/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testfunctionallinearfree( "functionallinearfree options" );
    testfunctionallinearfree.add_options()
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testfunctionallinearfree.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_functionallinearfree" ,
                     "test_functionallinearfree" ,
                     "0.2",
                     "nD(n=1,2,3) test FsFunctionalLinearFree and FsFunctionalLinearComposite",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}



template<int Dim, int Order>
void
testFunctionalLinearFree()
{

    typedef Backend<double> backend_type;
#if 0
    auto backend = backend_type::build( BACKEND_PETSC );
#else
    auto backend = backend_type::build( soption( _name="backend" ) );
#endif

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto v = Xh->element();
    auto element = project( Xh, elements(mesh) , cos( Px() ) );

    auto expr = integrate( _range=elements(mesh), _expr=grad(v)*trans(grad(v)) );
    auto functionalfree = functionalLinearFree( _space=Xh, _expr=expr, _backend=backend );

    auto functional = functionalLinear( _space=Xh, _backend=backend );
    *functional = expr;

    double result = functional->operator()( element );
    double resultfree = functionalfree->operator()( element );

    auto vector_functional = backend->newVector( Xh );
    functional->containerPtr(vector_functional);
    auto vector_functionalfree = backend->newVector( Xh );
    functionalfree->containerPtr(vector_functionalfree);

    double norm_vector_functional = vector_functional->l2Norm();
    double norm_vector_functionalfree = vector_functionalfree->l2Norm();

    double epsilon=1e-13;
    BOOST_CHECK_SMALL( math::abs(result-resultfree), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vector_functional-norm_vector_functionalfree), epsilon );
}


template<int Dim, int Order>
void
testFunctionalLinearComposite()
{

    typedef Backend<double> backend_type;
#if 0
    auto backend = backend_type::build( BACKEND_PETSC );
#else
    auto backend = backend_type::build( soption( _name="backend" ) );
#endif

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto composite = functionalLinearComposite( _space=Xh,  _backend=backend );
    auto compositefree = functionalLinearComposite( _space=Xh, _backend=backend );

    auto v = Xh->element();
    auto element = project( Xh, elements(mesh) , cos( Px() ) );

    //operators
    auto expr1 = integrate( _range=elements(mesh), _expr=grad(v)*trans(grad(v)) );
    auto expr2 = integrate( _range=elements(mesh), _expr=id(v)*trans(id(v)) );

    auto functionalfree1 = functionalLinearFree(  _space=Xh, _expr=expr1, _backend=backend );
    auto functionalfree2 = functionalLinearFree(  _space=Xh, _expr=expr2, _backend=backend );

    auto functional1 = functionalLinear( _space=Xh, _backend=backend );
    *functional1 = expr1;
    auto functional2 = functionalLinear( _space=Xh, _backend=backend );
    *functional2 = expr2;

    auto expr = expr1 + expr2 ;

    auto functionalfree = functionalLinearFree( _space=Xh, _expr=expr, _backend=backend );
    auto functional = functionalLinear( _space=Xh, _backend=backend );
    *functional = expr;

    //fill operator composite
    composite->addElement( functional1 );
    composite->addElement( functional2 );
    compositefree->addElement( functionalfree1 );
    compositefree->addElement( functionalfree2 );


    //test apply function
    double result_composite = composite->operator()( element );
    double result_compositefree = compositefree->operator()( element );
    double result = functional->operator()( element );
    double result_free = functionalfree->operator()( element );

    double epsilon=1e-13;
    BOOST_CHECK_SMALL( math::abs(result_composite-result_compositefree), epsilon );
    BOOST_CHECK_SMALL( math::abs(result-result_compositefree), epsilon );
    BOOST_CHECK_SMALL( math::abs(result_free-result_compositefree), epsilon );


    //test access functions
    auto vector_composite1 = backend->newVector( Xh );
    auto vector_composite2 = backend->newVector( Xh );
    auto vector_compositefree1 = backend->newVector( Xh );
    auto vector_compositefree2 = backend->newVector( Xh );

    composite->vecPtr( 0 , vector_composite1 );
    composite->vecPtr( 1 , vector_composite2 );
    compositefree->vecPtr( 0 , vector_compositefree1 );
    compositefree->vecPtr( 1 , vector_compositefree2 );

    auto vector_free1 = backend->newVector( Xh );
    auto vector_free2 = backend->newVector( Xh );

    auto vector1 = functional1->containerPtr();
    functionalfree1->containerPtr( vector_free1 );
    auto vector2 = functional2->containerPtr();
    functionalfree2->containerPtr( vector_free2 );

    double norm_vector_composite1 = vector_composite1->l2Norm();
    double norm_vector_compositefree1 = vector_compositefree1->l2Norm();
    double norm_vector_composite2 = vector_composite2->l2Norm();
    double norm_vector_compositefree2 = vector_compositefree2->l2Norm();
    double norm_vector1 = vector1->l2Norm();
    double norm_vector2 = vector2->l2Norm();
    double norm_vectorfree1 = vector_free1->l2Norm();
    double norm_vectorfree2 = vector_free2->l2Norm();

    BOOST_CHECK_SMALL( math::abs(norm_vector_composite1 - norm_vector_compositefree1 ), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vector_composite2 - norm_vector_compositefree2 ), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vector1 - norm_vector_composite1 ), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vector2 - norm_vector_composite2 ), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vector2 - norm_vector_compositefree2 ), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vectorfree2 - norm_vector_compositefree2 ), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_vectorfree1 - norm_vector_compositefree1 ), epsilon );

    //test sum of all vectors
    bool scalars_are_one=true;
    auto vec_sum = backend->newVector( Xh );
    composite->sumAllVectors(vec_sum, scalars_are_one );
    auto vec_sum_free = backend->newVector( Xh );
    compositefree->sumAllVectors(vec_sum_free, scalars_are_one );
    auto vec_free = backend->newVector( Xh );
    functionalfree->containerPtr(vec_free);
    auto vec = functional->containerPtr();

    double norm_sum_composite = vec_sum->l2Norm();
    double norm_sum_compositefree = vec_sum_free->l2Norm();
    double norm_sum_functional = vec->l2Norm();
    double norm_sum_functionalfree = vec_free->l2Norm();

    BOOST_CHECK_SMALL( math::abs(norm_sum_composite - norm_sum_functional), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_sum_compositefree - norm_sum_functionalfree), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_sum_compositefree - norm_sum_functional), epsilon );

}
/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( functionallinearfree )

BOOST_AUTO_TEST_CASE( test_1 )
{
    testFunctionalLinearFree<2,1>();
    //testFunctionalLinearFree<3,1>();
}

BOOST_AUTO_TEST_CASE( test_2 )
{
    testFunctionalLinearComposite<2,1>();
    //testFunctioalLinearComposite<3,1>();
}

BOOST_AUTO_TEST_SUITE_END()


