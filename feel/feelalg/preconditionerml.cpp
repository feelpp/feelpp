/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
Date: 2006-03-06

Copyright (C) 2006, 2009 EPFL

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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file preconditionerML.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 14-02-2008
*/


#include <feel/feelalg/preconditionerml.hpp>

namespace Feel
{
#if defined( FEELPP_HAS_TRILINOS_ML )

PreconditionerML::PreconditionerML( list_type options )
    :
    M_Prec(),
    M_List(),
    M_precType( "MultiGrid" )
{
    ML_Epetra::SetDefaults( "SA", M_List );
    /*
    M_List.set("max levels", 5);
    M_List.set("increasing or decreasing", "increasing");
    M_List.set("aggregation: type", "MIS");
    M_List.set("smoother: type","Aztec");
    M_List.set("smoother: pre or post", "both");
    M_List.set("coarse: type","Amesos-KLU");
    */

    M_List.set( "output", 0 );
}

PreconditionerML::PreconditionerML( PreconditionerML const& tc )
    :
    M_Prec( tc.M_Prec ),
    M_List( tc.M_List ),
    M_precType( tc.M_precType )
{}

int
PreconditionerML::buildPreconditioner( sparse_matrix_ptrtype& A )
{
    epetra_sparse_matrix_type* A_ptr = dynamic_cast<epetra_sparse_matrix_type*>( A.get() );

    M_Prec = boost::shared_ptr<prec_type>( new prec_type( A_ptr->mat(),
                                           M_List,
                                           true ) );

    return 0;
}

PreconditionerML::prec_ptrtype
PreconditionerML::getPrec()
{
    return M_Prec;
}

#endif // FEELPP_HAS_TRILINOS_ML
} // Feel
