/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-09-10

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
   \file test_fspace_context.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2012-09-10
*/

#define BOOST_TEST_MODULE test_matrix_destructor
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */
#include <feel/feelcore/pslogger.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description testmatrixdestructor( "test memory after matrix (shared ptr) destructor" );
    testmatrixdestructor.add_options()
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testmatrixdestructor.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_matrix_destructor" ,
                     "test_matrix_destructor" ,
                     "0.2",
                     "nD(n=1,2,3) test matrix_destructor",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}


void testMatrixDestructorWriteInfo(std::string file_name)
{
    std::string real_file_name = file_name + (boost::format("-%1%_%2%") %Environment::worldComm().globalSize() %Environment::worldComm().globalRank()).str();
    std::ifstream file ( real_file_name );
    std::string str;
    double memory0, memory1, memory2;
    for(int i=0; i<11; i++)
        file>>str;
    file >> memory0;
    for(int i=0; i<4; i++)
        file>>str;
    file >> memory1;
    for(int i=0; i<4; i++)
        file>>str;
    file >> memory2;


    double memory_taken = memory1 - memory0;
    double memory_free = memory1 - memory2;
    double lost_memory = memory2 - memory0;
    LOG(INFO)<<" ======= MEMORY INFO ========== "<<real_file_name;
    LOG( INFO ) << "memory taken by the matrix : "<<memory_taken;
    LOG( INFO ) << "memory free : "<<memory_free;
    LOG( INFO ) << "lost memory : "<<lost_memory;
}
std::string
format(std::string logm, int Dim, int Order )
{
    return (boost::format("\"%1% (%2%D,Order %3%)\"") % logm % Dim %Order).str();
}
template<int Dim, int Order>
void
testMatrixDestructor()
{
    auto mesh=unitHypercube<Dim>();
    //auto mesh=unitSquare();
    auto Xh = Pch<Order>( mesh );
    LOG(INFO)<<"nDof : "<<Xh->nDof();

    std::string str = ( boost::format("pslog-%1%D-P%2%") %Dim %Order ).str();
    PsLogger ps (str);
    ps.log(format("before matrix creation", Dim, Order) );
    auto matrix = backend()->newMatrix( Xh , Xh);
    ps.log(format("matrix created",Dim,Order));

    //call destructor
    matrix.reset();
    ps.log(format("matrix destroyed",Dim,Order));

    testMatrixDestructorWriteInfo(str);
}




/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( matrix_destructor )

BOOST_AUTO_TEST_CASE( test_1d_Order9 )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker checker("checker for petsc matrix leaks ");
#endif /* FEELPP_HAS_GPERFTOOLS */

    testMatrixDestructor<1,9>();

#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(checker.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
}

BOOST_AUTO_TEST_CASE( test_2d_Order6 )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker checker("checker for petsc matrix leaks");
#endif /* FEELPP_HAS_GPERFTOOLS */

    testMatrixDestructor<2,6>();

#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(checker.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
}
BOOST_AUTO_TEST_CASE( test_3d_Order3 )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker checker("checker for petsc matrix leaks");
#endif /* FEELPP_HAS_GPERFTOOLS */

    testMatrixDestructor<3,3>();

#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(checker.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
}

BOOST_AUTO_TEST_SUITE_END()
