/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2007-05-30

   Copyright (C) 2007, 2009 Universit� Joseph Fourier (Grenoble I)

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
   \file solverlineartrilinos.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2007-08-21
*/

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solverlineartrilinos.hpp>

namespace Feel
{

template <typename T>
void
SolverLinearTrilinos<T>::init()
{
    if ( !this->initialized() )
    {
        this->setInitialized( true );

        M_Solver.SetParameters( M_List, true );
    }

    // AztecOO defined a certain number of output parameters, and store them
    // in a double vector called status.
    double status[AZ_STATUS_SIZE];
    solver.GetAllAztecStatus( status );
}

template <typename T>
void
SolverLinearTrilinos<T>::setOptions( list_type _list )
{
    M_List = _list;
}

template <typename T>
list_type&
SolverLinearTrilinos<T>::getOptions()
{
    return M_List;
}

template <typename T>
std::pair<unsigned int, real_type>
SolverLinearTrilinos<T>::solve ( MatrixSparse<T>  const& matrix,
                                 Vector<T> & solution,
                                 Vector<T> const& rhs,
                                 const double tol,
                                 const unsigned int m_its )
{
    DVLOG(2) << "Matrix solver...\n";

    setRHS( rhs );
    setLHS( solution );
    setUserOperator( matrix );

    M_Solver.SetParameters( M_List, true );
    M_Solver.Iterate( m_its, tol );

    return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
}

template <typename T>
std::pair<unsigned int, real_type>
SolverLinearTrilinos<T>::solve (  MatrixSparse<T> const& matrix,
                                  MatrixSparse<T> const& preconditioner,
                                  Vector<T>& solution,
                                  Vector<T> const& rhs,
                                  const double tol,
                                  const unsigned int m_its )
{
    std::cout << "Matrix solver with preconditioner...\n";

    setRHS( rhs );
    setLHS( solution );
    setUserOperator( matrix );

    M_Solver.SetParameters( M_List, true );
    M_Solver.Iterate( m_its, tol );

    return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
}


} // Feel
