/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-09-25

  Copyright (C) 2006 Université Joseph Fourier

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
   \file solverlinearepetra.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-09-25
 */
#ifndef __epetra_linear_solver_h__
#define __epetra_linear_solver_h__


#include <feel/feelcore/feel.hpp>


#include <feel/feelalg/solverlinear.hpp>
#include <feel/feelalg/matrixepetra.hpp>
#include <feel/feelalg/vectorepetra.hpp>

// AztexOO
#include <AztecOO.h>

// Ifpack
#include <Ifpack_CrsIct.h>

// ml
#include <ml_config.h>
#include <ml_RowMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Teuchos_ParameterList.hpp>


namespace Feel
{
namespace trilinos
{
boost::shared_ptr<ML_Epetra::MultiLevelPreconditioner>
createMLPreconditioner( Epetra_RowMatrix const& mat, po::variables_map const& vm )
{
    // A is an Epetra_RowMatrix derived class object
    // solver is an AztecOO object
    Teuchos::ParameterList MList;
    // default values for smoothed aggregation
    ML_Epetra::SetDefaults( "SA",MList );
    MList.set( "max levels", vm["max-levels"].as<int>() );
    MList.set( "increasing or decreasing",vm["increasing-or-decreasing"].as<std::string>() );
    MList.set( "aggregation: type", vm["aggregation-type"].as<std::string>() );
    MList.set( "coarse: type",vm["coarse-type"].as<std::string>() );

    return boost::shared_ptr<ML_Epetra::MultiLevelPreconditioner>( new ML_Epetra::MultiLevelPreconditioner( mat,
            MList,
            true ) );
}

} // trilinos
} // Feel


#endif // #ifdef __epetra_linear_solver_h__

