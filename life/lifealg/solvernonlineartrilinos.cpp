/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
       Date: 2009-05-25

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file solvernonlineartrilinos.cpp
   \author Florent Vielfaure <florent.vielfaure@gmail.com>
   \date 2009-05-25
 */

#include <life/lifealg/glas.hpp>
#include <life/lifealg/vectorpetsc.hpp>
#include <life/lifealg/matrixpetsc.hpp>
#include <life/lifealg/solvernonlineartrilinos.hpp>

namespace Life
{

bool computeF( const Epetra_Vector & x, Epetra_Vector & f, NOX::Epetra::Interface::Required::FillType F ) {
    return true;
}
bool computeJecobian( const Epetra_Vector & x, Epetra_Operator & Jac ) {
    return true;
}
bool computePrecMatrix( const Epetra_Vector & x, Epetra_RowMatrix & M ) {
    return true;
}
bool computePreconditioner( const Epetra_Vector & x, Epetra_Operator & O ) {
    return true;
}

// SolverNonLinearTrilinos<> methods
template <typename T>
void SolverNonLinearTrilinos<T>::clear ()
{

}

template <typename T>
void SolverNonLinearTrilinos<T>::init ()
{

}

template <typename T>
std::pair<unsigned int, typename SolverNonLinearPetsc<T>::real_type>
SolverNonLinearTrilinos<T>::solve ( sparse_matrix_ptrtype&  jac_in,  // System Jacobian Matrix
                                    vector_ptrtype& x_in,    // Solution vector
                                    vector_ptrtype& r_in,    // Residual vector
                                    const double,              // Stopping tolerance
                                    const unsigned int)
{

}

template <typename T>
std::pair<unsigned int, typename SolverNonLinearPetsc<T>::real_type>
SolverNonLinearTrilinos<T>::solve ( dense_matrix_type&  jac_in,  // System Jacobian Matrix
                                    dense_vector_type& x_in,    // Solution vector
                                    dense_vector_type& r_in,    // Residual vector
                                    const double,              // Stopping tolerance
                                    const unsigned int)
{

}



//------------------------------------------------------------------
// Explicit instantiations
template class SolverNonLinearTrilinos<double>;




} // namespace Life
