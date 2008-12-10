/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-04

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file solvereigen.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-04
 */
#include <life/lifecore/life.hpp>

#include <life/lifealg/solvereigen.hpp>
#include <life/lifealg/solvereigenslepc.hpp>

namespace Life
{
template <typename T>
SolverEigen<T>::SolverEigen()
    :
    M_eigen_solver_type    (ARNOLDI),
    M_eigen_problem_type   (NHEP),
    M_position_of_spectrum (LARGEST_MAGNITUDE),
    M_is_initialized       (false)
{
}

template <typename T>
SolverEigen<T>::SolverEigen( SolverEigen const& eis )
    :
    M_eigen_solver_type    ( eis.M_eigen_solver_type ),
    M_eigen_problem_type   ( eis.M_eigen_problem_type ),
    M_position_of_spectrum ( eis.M_position_of_spectrum ),
    M_is_initialized       ( eis.M_is_initialized )
{
}

template <typename T>
SolverEigen<T>::~SolverEigen()
{
    this->clear ();
}
template <typename T>
boost::shared_ptr<SolverEigen<T> >
SolverEigen<T>::build(const SolverPackage solver_package)
{
    // Build the appropriate solver
    switch (solver_package)
        {

        case SOLVERS_SLEPC:
            {

#if defined( HAVE_SLEPC ) && defined( HAVE_PETSC_H )
                solvereigen_ptrtype ap(new SolverEigenSlepc<T>);
                return ap;
#else
                std::cerr << "Slepc is not available/installed" << std::endl;
                throw std::invalid_argument( "invalid solver slepc package" );
#endif
            }
            break;

        default:
            std::cerr << "ERROR:  Unrecognized eigen solver package: "
                      << solver_package
                      << std::endl;
            throw std::invalid_argument( "invalid solver package" );
        }

    return solvereigen_ptrtype();
}


/*
 * Explicit instantiations
 */
template class SolverEigen<double>;
}
