/* -*- mode: c++ -*-

This file is part of the Life library

Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
Date: 2006-03-06

Copyright (C) 2006 EPFL

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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file preconditionerML.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 21-05-2007
*/

#ifndef __preconditionerML_H
#define __preconditionerML_H 1

#include <lifeconfig.h>

#include <life/lifealg/matrixepetra.hpp>

#if defined( HAVE_TRILINOS_ML )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#include <life/lifealg/matrixepetra.hpp>
#include <ml_config.h>
#include <ml_RowMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
namespace Life
{
class PreconditionerML
{
public:

    typedef ML_Epetra::MultiLevelPreconditioner prec_type;
    typedef boost::shared_ptr<prec_type> prec_ptrtype;

    typedef MatrixSparse<double> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef MatrixEpetra epetra_sparse_matrix_type;
    typedef boost::shared_ptr<epetra_sparse_matrix_type> epetra_sparse_matrix_ptrtype;

    typedef Teuchos::ParameterList list_type;

    PreconditionerML( list_type options );

    PreconditionerML( PreconditionerML const& tc );

    int buildPreconditioner( sparse_matrix_ptrtype& A );

    prec_ptrtype getPrec();

  private:

    prec_ptrtype M_Prec;

    Teuchos::ParameterList M_List;

    std::string M_precType;

  };

} // Life
#endif /* HAVE_TRILINOS_ML */
#endif /* __preconditionerML_H */
