/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-03-29

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
   \file test_operatorlinearfree.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-29
*/

#define BOOST_TEST_MODULE test_operatorlinearfree
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <feel/feeldiscr/operatorlinearfree.hpp>
#include <feel/feeldiscr/operatorlinearcomposite.hpp>
#include <feel/feelcrb/eim.hpp>


/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testoperatorlinearfree( "operatorlinearfree options" );
    testoperatorlinearfree.add_options()
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testoperatorlinearfree.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_operatorlinearfree" ,
                     "test_operatorlinearfree" ,
                     "0.2",
                     "nD(n=1,2,3) test OperatorLinearFree and OperatorLinearComposite",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}



template<int Dim, int Order>
void
testOperatorLinearFree()
{

    typedef Backend<double> backend_type;
#if 0
    auto backend = backend_type::build( BACKEND_PETSC );
#else
    auto backend = backend_type::build( soption( _name="backend" ) );
#endif

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto vector = project( Xh, elements(mesh) , cos( Px() ) );
    auto result = Xh->element();

    auto expr = integrate( _range=elements(mesh), _expr=gradt(u)*trans(grad(v)) );
    auto opfree = opLinearFree( _domainSpace=Xh, _imageSpace=Xh, _expr=expr, _backend=backend );

    auto op = opLinear( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );
    *op = expr;

    op->apply( vector , result );
    double norm_op = result.l2Norm();

    result.zero();
    opfree->apply( vector , result );
    double norm_opfree = result.l2Norm();

    auto matrix_op = backend->newMatrix( Xh , Xh );
    *matrix_op=op->mat();

    auto matrix_opfree = backend->newMatrix( Xh , Xh );
    opfree->matPtr(matrix_opfree);

    double mat = matrix_op->l1Norm();
    double matfree = matrix_opfree->l1Norm();

    double energy_op = op->energy( vector , vector );
    double energy_opfree = opfree->energy( vector , vector );

    double epsilon=1e-13;
    BOOST_CHECK_SMALL( math::abs(norm_op-norm_opfree), epsilon );
    BOOST_CHECK_SMALL( math::abs(mat-matfree), epsilon );
    BOOST_CHECK_SMALL( math::abs(energy_op-energy_opfree), epsilon );
}


template<int Dim, int Order>
void
testOperatorLinearComposite()
{

    typedef Backend<double> backend_type;
#if 0
    auto backend = backend_type::build( BACKEND_PETSC );
#else
    auto backend = backend_type::build( soption( _name="backend" ) );
#endif

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto composite = opLinearComposite( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );
    auto compositefree = opLinearComposite( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );

    auto u = Xh->element();
    auto v = Xh->element();
    auto vector = project( Xh, elements(mesh) , cos( Px() ) );
    auto result = Xh->element();

    //operators
    auto expr1 = integrate( _range=elements(mesh), _expr=gradt(u)*trans(grad(v)) );
    auto opfree1 = opLinearFree( _domainSpace=Xh, _imageSpace=Xh, _expr=expr1, _backend=backend );
    auto expr2 = integrate( _range=elements(mesh), _expr=idt(u)*trans(id(v)) );
    auto opfree2 = opLinearFree( _domainSpace=Xh, _imageSpace=Xh, _expr=expr2, _backend=backend );

    auto expr = expr1 + expr2 ;

    auto op1 = opLinear( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );
    *op1 = expr1;
    auto op2 = opLinear( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );
    *op2 = expr2;

    //fill operator composite
    composite->addElement( op1 );
    composite->addElement( op2 );
    compositefree->addElement( opfree1 );
    compositefree->addElement( opfree2 );

    //test apply function
    result.zero(); composite->apply( vector , result );
    double norm_op_composite = result.l2Norm();
    result.zero();
    auto op = opLinear( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );
    *op = expr;
    op->apply( vector , result );
    double norm_op = result.l2Norm();

    result.zero();compositefree->apply( vector , result );
    double norm_opfree_composite = result.l2Norm();
    result.zero();
    auto opfree = opLinearFree( _domainSpace=Xh, _imageSpace=Xh, _expr=expr, _backend=backend );
    opfree->apply( vector , result );
    double norm_opfree = result.l2Norm();

    double epsilon=1e-13;

    BOOST_CHECK_SMALL( math::abs(norm_op_composite-norm_op), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_opfree_composite-norm_opfree), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_op_composite-norm_opfree_composite), epsilon );


    //test access functions

    auto mat_operator1 = backend->newMatrix( _test=Xh , _trial=Xh );
    auto mat_operator2 = backend->newMatrix( _test=Xh , _trial=Xh );
    composite->matrixPtr(0 , mat_operator1 );
    composite->matrixPtr(1 , mat_operator2 );
    double norm_mat_operator1_comp = mat_operator1->l1Norm();
    double norm_mat_operator1 = op1->matPtr()->l1Norm();
    double norm_mat_operator2_comp = mat_operator2->l1Norm();
    double norm_mat_operator2 = op2->matPtr()->l1Norm();

    BOOST_CHECK_SMALL( math::abs(norm_mat_operator1_comp-norm_mat_operator1), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_mat_operator2_comp-norm_mat_operator2), epsilon );

    auto mat_operatorfree1 = backend->newMatrix( _test=Xh , _trial=Xh );
    auto mat_operatorfree2 = backend->newMatrix( _test=Xh , _trial=Xh );
    compositefree->matrixPtr( 0 , mat_operatorfree1 );
    compositefree->matrixPtr( 1 , mat_operatorfree2 );
    double norm_mat_operatorfree1_comp = mat_operatorfree1->l1Norm();
    double norm_mat_operatorfree2_comp = mat_operatorfree2->l1Norm();

    auto mat_free1 = backend->newMatrix( _test=Xh, _trial=Xh);
    auto mat_free2 = backend->newMatrix( _test=Xh, _trial=Xh);
    opfree1->matPtr( mat_free1 );
    opfree2->matPtr( mat_free2 );
    double norm_mat_operatorfree1 = mat_free1->l1Norm();
    double norm_mat_operatorfree2 = mat_free2->l1Norm();

    BOOST_CHECK_SMALL( math::abs(norm_mat_operator1_comp - norm_mat_operator1), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_mat_operator2_comp - norm_mat_operator2), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_mat_operatorfree2_comp - norm_mat_operatorfree2), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_mat_operatorfree1_comp - norm_mat_operatorfree1), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_mat_operator1_comp - norm_mat_operatorfree1_comp), epsilon );

    //test sum of all matrices

    bool scalars_are_one=true;
    auto mat_sum = backend->newMatrix( _test=Xh , _trial=Xh );
    composite->sumAllMatrices( mat_sum, scalars_are_one );
    auto mat_sum_free = backend->newMatrix( _test=Xh , _trial=Xh );
    compositefree->sumAllMatrices(mat_sum_free, scalars_are_one );
    opfree->matPtr(mat_free1);
    double norm_sum_composite = mat_sum->l1Norm();
    double norm_sum_compositefree = mat_sum_free->l1Norm();
    double norm_sum_op = op->matPtr()->l1Norm();
    double norm_sum_opfree = mat_free1->l1Norm();

    BOOST_CHECK_SMALL( math::abs(norm_sum_composite - norm_sum_op), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_sum_compositefree - norm_sum_opfree), epsilon );
    BOOST_CHECK_SMALL( math::abs(norm_sum_compositefree - norm_sum_op), epsilon );

}



template<int Dim, int Order>
void
testExpression()
{

    typedef Backend<double> backend_type;
#if 0
    auto backend = backend_type::build( BACKEND_PETSC );
#else
    auto backend = backend_type::build( soption( _name="backend" ) );
#endif

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();

    auto operatorlinear = opLinear( _domainSpace=Xh, _imageSpace=Xh, _backend=backend );

    auto matrix = backend->newMatrix( _test=Xh , _trial=Xh );
    double norm=0,normfree=0;

    typedef decltype( operatorlinear ) operator_type;
    typedef decltype( u ) element_type;

    std::vector< operator_type > operators_vector;
    std::vector< operator_type > operators_free_vector;
    std::vector< element_type > elements_vector;

    auto x = project( Xh, elements(mesh), Px() );         elements_vector.push_back( x );
    auto xx = project( Xh, elements(mesh), Px()*Px() );   elements_vector.push_back( xx );
    auto cosy = project( Xh, elements(mesh), cos(Py()) ); elements_vector.push_back( cosy );
    auto xy = project( Xh, elements(mesh), Px()*Py() );   elements_vector.push_back( xy );

    double last_value=0;

    double epsilon=1e-13;

    for(int i=0; i<4; i++)
    {
        auto expr = integrate( _range=elements(mesh),
                               _expr=-idv(elements_vector[i])*gradt(u)*trans(grad(v)) );

        //build operators
        auto opfree = opLinearFree( _domainSpace=Xh, _imageSpace=Xh, _expr=expr, _backend=backend );
        auto op = opLinear( _domainSpace=Xh, _imageSpace=Xh,  _backend=backend );
        *op = expr;

        //store operators
        operators_vector.push_back( op );
        operators_free_vector.push_back( opfree );

        //compute norm of associated matrix
        opfree->matPtr( matrix ); normfree=matrix->l1Norm();
        *matrix=op->mat();            norm=matrix->l1Norm();

        LOG(INFO)<<"during the construction loop, for i = "<<i<<" - norm : "<<norm<<" and normfree : "<<normfree;

        //this test is ok
        BOOST_CHECK_SMALL( norm-normfree , epsilon );

        if( i == 3 ) last_value = norm;
    }
    for(int i=0; i<operators_vector.size(); i++)
    {
        //compute norm of the associated matrix
        *matrix=operators_vector[i]->mat();            norm=matrix->l1Norm();
        operators_free_vector[i]->matPtr(matrix);  normfree=matrix->l1Norm();
        LOG(INFO)<<"outside the construction loop for i = "<<i<<" - norm : "<<norm<<" and normfree : "<<normfree;

        //this test is not ok, except for the last value of i
        BOOST_CHECK_SMALL( norm-normfree , epsilon );
    }


}

/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( operatorlinearfree )

BOOST_AUTO_TEST_CASE( test_1 )
{
    testOperatorLinearFree<2,1>();
    //testOperatorLinearFree<3,1>();
}

BOOST_AUTO_TEST_CASE( test_2 )
{
    testOperatorLinearComposite<2,1>();
    //testOperatorLinearComposite<3,1>();
}

BOOST_AUTO_TEST_CASE( test_3 )
{
    testExpression<2,1>();
}

BOOST_AUTO_TEST_SUITE_END()
