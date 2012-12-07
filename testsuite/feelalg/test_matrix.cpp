/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_matrix.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>

using namespace Feel;


Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_matrix" ,
                           "test_matrix" ,
                           "0.1",
                           "matrix test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2005,2006 EPFL" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
#if defined(FEELPP_HAS_PETSC_H)
class test_matrixpetsc
    :
public Application
{
public:

    test_matrixpetsc( int argc, char** argv )
        :
        Application( argc, argv, makeAbout() )
    {
        std::cout << "id: " << this->comm().rank()  << "\n"
                  << "nprocs: " << this->comm().size() << "\n";
    }

    void operator()() const
    {
        int nprocs = this->comm().size();
        int m = 8;//10*Application::nProcess();
        int n = 8;//10*Application::nProcess();

        VectorPetsc<double> vec;
        MatrixPetsc<double> mat;
        std::cout << "is initialized ? " << mat.isInitialized() << "\n";
        mat.init( m*n, m*n, m*n/nprocs, m*n );
        std::cout << "is initialized ? " << mat.isInitialized() << "\n";

        /*
           Currently, all PETSc parallel matrix formats are partitioned by
           contiguous chunks of rows across the processors.  Determine which
           rows of the matrix are locally owned.
        */
        int Istart;
        int Iend;
        MatGetOwnershipRange( mat.mat(),&Istart,&Iend );
        VLOG(1) << "Istart = "<< Istart << "\n";
        VLOG(1) << "Iend = "<< Iend << "\n";

        /*
           Set matrix elements for the 2-D, five-point stencil in parallel.
           - Each processor needs to insert only elements that it owns
           locally (but any non-local elements will be sent to the
           appropriate processor during matrix assembly).
           - Always specify global rows and columns of matrix entries.

           Note: this uses the less common natural ordering that orders first
           all the unknowns for x = h then for x = 2h etc; Hence you see J = I +- n
           instead of J = I +- m as you might expect. The more standard ordering
           would first do all variables for y = h, then y = 2h etc.

        */

        for ( int I=Istart; I<Iend; I++ )
        {
            int J;
            double v = -1.0;
            int i = I/n;
            int j = I - i*n;
            VLOG(1) << "I= " << I << "\n";

            if ( i>0 )
            {
                J = I - n+1;
                VLOG(1) << "1 J= " << J << "\n";
                mat.set( I,J,v );
            }

            if ( i<m-1 )
            {
                J = I + n-1;
                VLOG(1) << "2 J= " << J << "\n";
                mat.set( I,J,v );
            }

            if ( j>0 )
            {
                J = I - 1;
                VLOG(1) << "3 J= " << J << "\n";
                mat.set( I,J,v );
            }

            if ( j<n-1 )
            {
                J = I + 1;
                VLOG(1) << "4 J= " << J << "\n";
                mat.set( I,J,v );
            }

            v = 4.0;
            mat.set( I,I,v );
        }

        VLOG(1) << "closing petsc matrix\n";
        mat.close();
        VLOG(1) << "closing petsc matrix done\n";

        VLOG(1) << "saving petsc matrix in matlab\n";
        //mat.printMatlab("m");
        //mat.printMatlab(std::string("/tmp/mat.m") );
        mat.printMatlab( "mat.m" );
        VLOG(1) << "saving petsc matrix in matlab done\n";

    }
};
#endif

#if defined(FEELPP_HAS_PETSC_H)
int main( int argc, char** argv )
{
    Feel::Environment env( argc,argv );
    test_matrixpetsc t( argc, argv );
    t();
    return EXIT_SUCCESS;
}
#else
int main( int /*argc*/, char** /*argv*/ )
{
    return EXIT_SUCCESS;
}
#endif
