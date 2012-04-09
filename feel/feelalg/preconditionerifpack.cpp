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
   \file preconditionerifpack.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 14-02-2008
*/


#include <feel/feelalg/preconditionerifpack.hpp>

namespace Feel
{
#if defined( FEELPP_HAS_TRILINOS_IFPACK )
PreconditionerIfpack::PreconditionerIfpack( std::string str  )
    :
    M_Prec(),
    M_List(),
    M_precType( str )
{
    M_List.set( "fact: drop tolerance", 1e-5 );
    M_List.set( "fact: level-of-fill",  1 );
    M_List.set( "partitioner: local parts", 4 );
    M_List.set( "partitioner: overlap", 4 );

    M_overlap = 4;

    M_precType = str;

    M_List.set( "amesos: solver type", "Amesos_Klu" );

}

PreconditionerIfpack::PreconditionerIfpack( list_type options, std::string str  )
    :
    M_Prec(),
    M_List(),
    M_precType( str )
{
    M_List.set( "fact: drop tolerance", options.get( "fact: drop tolerance", 1e-5 ) );
    M_List.set( "fact: level-of-fill",  options.get( "fact: level-of-fill" , 1   ) );
    M_List.set( "partitioner: local parts", options.get( "partitioner: local parts", 4 ) );
    M_List.set( "partitioner: overlap", options.get( "partitioner: overlap", 4 ) );

    M_overlap = options.get( "partitioner: overlap", 4 );

    M_precType = str;

    M_List.set( "amesos: solver type", "Amesos_Klu" );

}

PreconditionerIfpack::PreconditionerIfpack( PreconditionerIfpack const& tc )
    :
    M_Prec( tc.M_Prec ),
    M_List( tc.M_List ),
    M_precType( tc.M_precType ),
    M_overlap( tc.M_overlap )
{
}


void
PreconditionerIfpack::setAmesosSolver( std::string str )
{
    M_List.set( "amesos: solver type", str );
}

void
PreconditionerIfpack::setOptions( list_type options )
{
    //M_List = options;

    M_List.set( "fact: drop tolerance", options.get( "fact: drop tolerance", 1e-5 ) );
    M_List.set( "fact: level-of-fill",  options.get( "fact: level-of-fill" , 3   ) );
    M_List.set( "partitioner: local parts", options.get( "partitioner: local parts", 4 ) );
    M_List.set( "partitioner: overlap", options.get( "partitioner: overlap", 4 ) );

    M_overlap = options.get( "partitioner: overlap", 4 );
}


int
PreconditionerIfpack::initializePreconditioner( sparse_matrix_ptrtype const& A )
{
    Ifpack factory;

    epetra_sparse_matrix_type* A_ptr = const_cast<epetra_sparse_matrix_type*>( dynamic_cast<epetra_sparse_matrix_type const*>( A.get() ) );

    M_Prec.reset( factory.Create( M_precType, &A_ptr->mat(), M_overlap ) );

    IFPACK_CHK_ERR( M_Prec->SetParameters( M_List ) );
    IFPACK_CHK_ERR( M_Prec->Initialize() );

    return EXIT_SUCCESS;
}

int
PreconditionerIfpack::computePreconditioner()
{
    IFPACK_CHK_ERR( M_Prec->Compute() );

    return EXIT_SUCCESS;
}

int
PreconditionerIfpack::buildPreconditioner( sparse_matrix_ptrtype const& A )
{
    initializePreconditioner( A );
    int result = computePreconditioner();

    return result;
}

PreconditionerIfpack::prec_ptrtype
PreconditionerIfpack::getPrec()
{
    return M_Prec;
}

#endif // FEELPP_HAS_TRILINOS_IFPACK
} // Feel
