/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-03-14

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
   \file test_preassemble.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-23
*/

#define BOOST_TEST_MODULE test_preassemble
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <feel/feeldiscr/preassembleobject.hpp>

/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testpreassemble( "preassemble options" );
    testpreassemble.add_options()
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testpreassemble.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_preassemble" ,
                     "test_preassemble" ,
                     "0.2",
                     "nD(n=1,2,3) test pre-assemble matrix and vector",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}



template<int Dim, int Order>
void
testPreassemble()
{
    //solve the laplacian - without store PreAssemble objects

    typedef Backend<double> backend_type;
    auto backend = backend_type::build( BACKEND_PETSC );

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();

    auto solution1 = Xh->element();
    auto solution2 = Xh->element();

    double penaldir = 50; // penalization coefficient

    //first, without using preassemble objects
    auto A1 = backend->newMatrix( _test=Xh , _trial=Xh );
    form2( _trial=Xh, _test=Xh ,_matrix=A1 )  = integrate(_range=elements(mesh), _expr=gradt(u)*trans(grad(v)) );
    form2( _trial=Xh, _test=Xh ,_matrix=A1 ) += integrate( _range=boundaryfaces(mesh),
                                                           _expr=-gradt(u)*vf::N()*id(v) - grad(v)*vf::N()*idt(u) + penaldir*id(v)*idt(u)/hFace() );

    auto F1 = backend->newVector( Xh );
    form1( _test=Xh , _vector=F1) = integrate ( _range=elements(mesh), _expr=id(v) );

    backend->solve( _matrix=A1 , _solution=solution1, _rhs=F1 );


    //now, solve the same equation but using preassemble
    //fill preassemble objects :
    auto range = elements(mesh);
    auto expr = gradt(u)*trans(grad(v));
    auto preassemble_matrix1 = PreAssembleMatrix( range , expr , Xh, Xh );
    auto range2 = boundaryfaces(mesh);
    auto expr2  =  -gradt(u)*vf::N()*id(v) - grad(v)*vf::N()*idt(u) + penaldir*id(v)*idt(u)/hFace() ;
    auto preassemble_matrix2 = PreAssembleMatrix( range2 , expr2 , Xh, Xh );

    //use preassemble objects
    auto A2 = backend->newMatrix( _test=preassemble_matrix1->test() , _trial=preassemble_matrix1->trial() );
    form2( _test=preassemble_matrix1->test() , _trial=preassemble_matrix1->trial() , _matrix=A2 )  = integrate( preassemble_matrix1->range() , preassemble_matrix1->expr() );
    form2( _test=preassemble_matrix2->test() , _trial=preassemble_matrix2->trial() , _matrix=A2 ) += integrate( preassemble_matrix2->range() , preassemble_matrix2->expr() );

    auto rangeV = elements(mesh);
    auto exprV = id(v);
    auto preassemble_vector = PreAssembleVector(rangeV , exprV , Xh);
    auto F2 = backend->newVector( preassemble_vector->test() );
    form1( _test=preassemble_vector->test() , _vector=F2) = integrate ( _range=rangeV, _expr=exprV );

    backend->solve( _matrix=A2 , _solution=solution2, _rhs=F2 );

    // check

    double normA1 = A1->l1Norm();
    double normA2 = A2->l1Norm();

    double normSolution1 = solution1.l2Norm();
    double normSolution2 = solution2.l2Norm();

    double normF1 = F1->l2Norm();
    double normF2 = F2->l2Norm();

    BOOST_CHECK_SMALL( (normA1-normA2) , 1e-14 );
    BOOST_CHECK_SMALL( (normF1-normF2) , 1e-14 );
    BOOST_CHECK_SMALL( (normSolution1-normSolution2) , 1e-14 );

}

template<int Dim, int Order>
void
testPreassemble2()
{
    //test to store PreAssemble objects

    typedef Backend<double> backend_type;
    auto backend = backend_type::build( BACKEND_PETSC );

    typedef PreAssembleObjectBase object_base_type;
    typedef boost::shared_ptr<object_base_type> object_base_ptrtype;

    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();

    auto solution1 = Xh->element();
    auto solution2 = Xh->element();

    double penaldir = 50; // penalization coefficient

    //first, without using preassemble objects
    auto A1 = backend->newMatrix( _test=Xh , _trial=Xh );
    form2( _trial=Xh, _test=Xh ,_matrix=A1 )  = integrate(_range=elements(mesh), _expr=gradt(u)*trans(grad(v)) );
    form2( _trial=Xh, _test=Xh ,_matrix=A1 ) += integrate( _range=boundaryfaces(mesh),
                                                           _expr=-gradt(u)*vf::N()*id(v) - grad(v)*vf::N()*idt(u) + penaldir*id(v)*idt(u)/hFace() );

    auto F1 = backend->newVector( Xh );
    form1( _test=Xh , _vector=F1) = integrate ( _range=elements(mesh), _expr=id(v) );

    backend->solve( _matrix=A1 , _solution=solution1, _rhs=F1 );

    //now we store assemble objects
    std::vector< object_base_ptrtype > matrices;
    std::vector< object_base_ptrtype > vectors;

    for(int i=0; i<2; i++ )
    {
        if( i==0)
        {
            auto range = elements(mesh);
            auto expr = gradt(u)*trans(grad(v));
            auto preassemble = PreAssembleMatrix( range , expr , Xh, Xh );
            matrices.push_back( preassemble );
        }
        else
        {
            auto range = boundaryfaces(mesh);
            auto expr  =  -gradt(u)*vf::N()*id(v) - grad(v)*vf::N()*idt(u) + penaldir*id(v)*idt(u)/hFace() ;
            auto preassemble = PreAssembleMatrix( range , expr , Xh, Xh );
            matrices.push_back( preassemble );
        }
    }//end if loop over matrix contributions

    auto range = elements(mesh);
    auto expr = id(v);
    auto preassemble = PreAssembleVector(range , expr , Xh);
    vectors.push_back( preassemble );

    //then we use them
    auto A2 = backend->newMatrix( _test=Xh , _trial=Xh );
    A2->zero();
    for(int i=0; i<matrices.size(); i++)
    {
        //assemble matrix from matrices[i] and add it to A2
        A2->addMatrix( 1.0 , *matrices[i]->assembleMatrix() );
    }
    auto F2 = backend->newVector( Xh );
    F2->zero();
    for(int i=0; i<vectors.size(); i++)
    {
        //assemble vector from vectors[i] and add it to F2
        F2->add( *vectors[i]->assembleVector() );
    }

    backend->solve( _matrix=A2 , _solution=solution2, _rhs=F2 );

    // check
    double normA1 = A1->l1Norm();
    double normA2 = A2->l1Norm();

    double normSolution1 = solution1.l2Norm();
    double normSolution2 = solution2.l2Norm();

    double normF1 = F1->l2Norm();
    double normF2 = F2->l2Norm();

    BOOST_CHECK_SMALL( (normA1-normA2) , 1e-14 );
    BOOST_CHECK_SMALL( (normF1-normF2) , 1e-14 );
    BOOST_CHECK_SMALL( (normSolution1-normSolution2) , 1e-14 );

}


/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( preassemble )

BOOST_AUTO_TEST_CASE( test_1 )
{
    testPreassemble<2,1>();
}
BOOST_AUTO_TEST_CASE( test_2 )
{
    testPreassemble2<2,1>();
}

BOOST_AUTO_TEST_SUITE_END()


