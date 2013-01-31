/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   Date: 2006-03-06

   Copyright (C) 2006 EPFL

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
   \file operatortrilinos.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 31-05-2007
*/

#ifndef __operator_trilinos_matrix_H
#define __operator_trilinos_matrix_H 1


#include <feel/feelalg/matrixepetra.hpp>
#include <feel/feelalg/vectorepetra.hpp>
// #include <feel/feelalg/operatortrilinos.hpp>

#include <feel/feelalg/solverlineartrilinos.hpp>

#if defined( FEELPP_HAS_TRILINOS_EPETRA )
#include "Epetra_ConfigDefs.h"
#include "Epetra_Operator.h"


#if defined( FEELPP_HAS_MPI )
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"




namespace Feel
{
class OperatorMatrix : public Epetra_Operator
{
public:

    typedef MatrixSparse<double> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef MatrixEpetra epetra_sparse_matrix_type;
    typedef boost::shared_ptr<epetra_sparse_matrix_type> epetra_sparse_matrix_ptrtype;

    typedef Vector<double> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    typedef VectorEpetra<double> epetra_vector_type;
    typedef boost::shared_ptr<epetra_vector_type> epetra_vector_ptrtype;

    typedef Epetra_Operator prec_type;
    typedef boost::shared_ptr<prec_type> prec_ptrtype;

    typedef SolverLinearTrilinos<double> solver_type;
    typedef boost::shared_ptr<solver_type> solver_ptrtype;

    typedef Teuchos::ParameterList list_type;

    typedef solver_type::real_type real_type;


    OperatorMatrix( sparse_matrix_ptrtype const& F, std::string _label, bool transpose = 0 )
        :
        M_F( F ),
        M_Matrix( dynamic_cast<epetra_sparse_matrix_type const*>( F.get() ) ),
        M_Prec( ),
        M_Solver( new solver_type() ),
        M_hasInverse( 1 ),
        M_hasApply( 1 ),
        M_useTranspose( transpose ),
        M_domainMap( M_Matrix->getDomainMap() ),
        M_rangeMap( M_Matrix->getRangeMap() ),
        M_label( _label ),
        M_maxiter( 1000 ),
        M_tol( 1e-10 )
    {
        DVLOG(2) << "Create operator " << Label() << " ...\n";
    }

    OperatorMatrix( sparse_matrix_ptrtype const& F,
                    list_type const& options,
                    std::string _label,
                    prec_ptrtype  Prec, bool transpose = 0 )
        :
        M_F( F ),
        M_Matrix( dynamic_cast<epetra_sparse_matrix_type const*>( F.get() ) ),
        M_Prec( Prec ),
        M_Solver( new solver_type() ),
        M_hasInverse( 1 ),
        M_hasApply( 1 ),
        M_useTranspose( transpose ),
        M_domainMap( M_Matrix->getDomainMap() ),
        M_rangeMap( M_Matrix->getRangeMap() ),
        M_label( _label ),
        M_maxiter( 1000 ),
        M_tol( 1e-10 )
    {
        DVLOG(2) << "Create operator " << Label() << " ...\n";

        M_Solver->setOptions( options );

        M_tol = M_Solver->getOptions().get( "tol", 1e-10 );
        M_maxiter = M_Solver->getOptions().get( "max_iter", 100 );


    }

    OperatorMatrix( const OperatorMatrix& tc )
        :
        M_F( tc.M_F ),
        M_Matrix( tc.M_Matrix ),
        M_Prec( tc.M_Prec ),
        M_Solver( tc.M_Solver ),
        M_hasInverse( tc.M_hasInverse ),
        M_hasApply( tc.M_hasApply ),
        M_useTranspose( tc.M_useTranspose ),
        M_domainMap( tc.M_domainMap ),
        M_rangeMap( tc.M_rangeMap ),
        M_label( tc.M_label ),
        M_maxiter( tc.M_maxiter ),
        M_tol( tc.M_tol )
    {
        DVLOG(2) << "Copy operator " << Label() << " ...\n";
    }

    bool hasInverse() const
    {
        return M_hasInverse;
    }

    bool hasApply() const
    {
        return M_hasApply;
    }


    int Apply( const vector_ptrtype& X, vector_ptrtype& Y ) const
    {
        Apply( dynamic_cast<epetra_vector_type const&>( *X ).vec(),
               dynamic_cast<epetra_vector_type&>( *Y ).vec() );

        return !hasApply();
    }

    int Apply( const Epetra_MultiVector & X, Epetra_MultiVector & Y ) const
    {
        DVLOG(2) << "Apply Operator " << Label() << "\n";

        M_Matrix->multiply( false,X,Y );

        DVLOG(2) << "Apply Operator " << Label() << " successfully\n";
        return !hasApply();
    }

    int ApplyInverse ( const vector_ptrtype& X, vector_ptrtype& Y ) const
    {
        FEELPP_ASSERT( hasInverse() ).error( "This operator cannot be inverted." );

        int result = ApplyInverse( dynamic_cast<epetra_vector_type const&>( *X ).vec(),
                                   dynamic_cast<epetra_vector_type&>( *Y ).vec() );

        return result;
    }

    int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        FEELPP_ASSERT( hasInverse() ).error( "This operator cannot be inverted." );

        std::pair<unsigned int, real_type> result;

        DVLOG(2) << "Apply Inverse Operator " << Label() << "\n";

        if ( M_Prec.get() != 0 )
            result = M_Solver->solve( *M_Matrix, M_Prec, Y, X, M_tol, M_maxiter );

        else
            result = M_Solver->solve( *M_Matrix, Y, X, M_tol, M_maxiter );

        DVLOG(2) << "Apply Inverse Operator " << Label() << " successfully\n";

        int result_integer = result.first;

        return result_integer;
    }

    // other function
    int SetUseTranspose( bool UseTranspose = 0 )
    {
        M_useTranspose = UseTranspose;

        return UseTranspose;
    }

    double NormInf() const
    {
        return M_Matrix->linftyNorm();
    }

    const char * Label () const
    {
        return( M_label.c_str() );
    }

    bool UseTranspose() const
    {
        return M_useTranspose;
    }

    bool HasNormInf () const
    {
        return( true );
    }

    const Epetra_Comm& Comm() const
    {
        return( M_Matrix->getDomainMap().Comm() );
    }

    const Epetra_Map& OperatorDomainMap() const
    {
        return( M_domainMap );
    }

    const Epetra_Map& OperatorRangeMap() const
    {
        return( M_rangeMap );
    }

    ~OperatorMatrix()
    {
        DVLOG(2) << "Destroyed matrix operator: " << this->Label() << " ...\n";
    };

private:

    sparse_matrix_ptrtype M_F;
    epetra_sparse_matrix_type const* M_Matrix;

    prec_ptrtype M_Prec;

    solver_ptrtype M_Solver;

    bool M_hasInverse, M_hasApply, M_useTranspose;

    Epetra_Map M_domainMap, M_rangeMap;

    std::string M_label;

    int M_maxiter;

    double M_tol;
};









template< typename operator_type >
class OperatorInverse : public Epetra_Operator
{
public:

    typedef boost::shared_ptr<operator_type> operator_ptrtype;

    // This constructor implements the F^-1 operator
    OperatorInverse( operator_ptrtype& F )
        :
        M_F( F )
    {
        this->setName();
        DVLOG(2) << "Create inverse operator " << this->Label() << "...\n";
    }

    OperatorInverse( const OperatorInverse& tc )
        :
        M_F( tc.M_F ),
        M_Label( tc.M_Label )
    {
        DVLOG(2) << "Copy inverse operator " << this->Label() << "...\n";
    }

    bool hasInverse() const
    {
        return M_F->hasApply();
    }

    bool hasApply() const
    {
        return M_F->hasInverse();
    }

    int Apply( const Epetra_MultiVector & X, Epetra_MultiVector & Y ) const
    {
        DVLOG(2) << "Apply matrix " << Label() << "\n";

        FEELPP_ASSERT( hasApply() ).error( "This operator cannot be applied." );
        M_F->ApplyInverse( X,Y );

        return !hasApply();
    }

    int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        DVLOG(2) << "ApplyInverse matrix " << Label() << "\n";

        FEELPP_ASSERT( hasInverse() ).error( "This operator cannot be inverted." );

        M_F->Apply( X,Y );

        return !hasInverse();
    }

    // other function
    int SetUseTranspose( bool UseTranspose )
    {
        return( false );
    }

    double NormInf() const
    {
        return( false );
    }

    const char * Label () const
    {
        return M_Label.c_str();
    }

    bool UseTranspose() const
    {
        return( false );
    }

    bool HasNormInf () const
    {
        return( false );
    }


    const Epetra_Comm & Comm() const
    {
        return( M_F->Comm() );
    }

    const Epetra_Map & OperatorDomainMap() const
    {
        return( M_F->OperatorRangeMap() );
    }

    const Epetra_Map & OperatorRangeMap() const
    {
        return( M_F->OperatorDomainMap() );
    }

    ~OperatorInverse()
    {
        DVLOG(2) << "Destroyed inverse operator: " << this->Label() << " ...\n";
    };

private:

    operator_ptrtype M_F;

    std::string M_Label;

    void setName()
    {
        std::string L = M_F->Label();
        L.append( ")" );
        std::string temp( "inv(" );
        temp.append( L );

        M_Label = temp;
    }
};






template< typename op1_type, typename op2_type >
class OperatorCompose : public Epetra_Operator
{
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;
    typedef boost::shared_ptr<op2_type> op2_ptrtype;

    // This constructor implements the (F o G) operator

    OperatorCompose()
        :
        M_F(),
        M_G()
    {
        std::string t( M_F->Label() );
        std::string u( M_G->Label() );

        t.append( "*" );
        t.append( u );
        M_Label = t;
    }

    OperatorCompose( op1_ptrtype& F, op2_ptrtype& G )
        :
        M_F( F ),
        M_G( G )
    {
        std::string t( F->Label() );
        std::string u( G->Label() );

        t.append( "*" );
        t.append( u );
        M_Label = t;
        DVLOG(2) << "Create operator " << Label() << " ...\n";
    }

    OperatorCompose( const OperatorCompose& tc )
        :
        M_F( tc.M_F ),
        M_G( tc.M_G ),
        M_Label( tc.M_Label )
    {
        DVLOG(2) << "Copy operator " << Label() << " ...\n";
    }

    bool hasInverse() const
    {
        return M_F->hasInverse() * M_G->hasInverse();
    }

    bool hasApply() const
    {
        return M_F->hasApply() * M_G->hasApply();
    }

    int Apply( const Epetra_MultiVector & X, Epetra_MultiVector & Y ) const
    {
        FEELPP_ASSERT( hasApply() ).error( "This operator cannot be applied." );

        DVLOG(2) << "Apply operator " << Label() << " ...\n";

        Epetra_MultiVector Z( M_G->OperatorRangeMap(), 1 );

        M_G->Apply( X,Z );
        M_F->Apply( Z,Y );

        return !hasApply();
    }

    int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        FEELPP_ASSERT( hasInverse() ).error( "This operator cannot be inverted." );

        DVLOG(2) << "Apply Inverse operator " << Label() << " ...\n";

        Epetra_MultiVector Z( X );

        M_F->ApplyInverse( X,Z );
        M_G->ApplyInverse( Z,Y );

        return hasInverse();
    }

    // other function
    int SetUseTranspose( bool UseTranspose )
    {
        return( false );
    }

    double NormInf() const
    {
        return 0;
    }

    const char * Label () const
    {

        return M_Label.c_str();
    }

    bool UseTranspose() const
    {
        return( false );
    }

    bool HasNormInf () const
    {
        return( false );
    }


    const Epetra_Comm & Comm() const
    {
        return( M_F->Comm() );
    }

    const Epetra_Map & OperatorDomainMap() const
    {
        return( M_G->OperatorDomainMap() );
    }

    const Epetra_Map & OperatorRangeMap() const
    {
        return( M_F->OperatorRangeMap() );
    }

    ~OperatorCompose()
    {
        DVLOG(2) << "Destroyed compose operator: " << this->Label() << " ...\n";
    };

private:

    op1_ptrtype M_F;
    op2_ptrtype M_G;

    std::string M_Label;
};













template< typename op1_type>
class OperatorScale : public Epetra_Operator
{
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;

    // This constructor implements the (\alpha F) operator
    OperatorScale()
        :
        M_F(),
        M_alpha( 0 )
    {
    }

    OperatorScale( op1_ptrtype& F )
        :
        M_F( F ),
        M_alpha( 1 )
    {
        std::string temp( "alpha" );
        std::string t = M_F->Label();

        temp.append( "." );
        temp.append( t );

        M_Label = temp;

        DVLOG(2) << "Create scale operator " << Label() << " ...\n";
    }

    OperatorScale( op1_ptrtype& F, double alpha )
        :
        M_F( F ),
        M_alpha( alpha )
    {
        std::string temp( "alpha" );
        std::string t = M_F->Label();

        temp.append( "." );
        temp.append( t );

        M_Label = temp;

        DVLOG(2) << "Create scale operator " << Label() << " ...\n";
    }

    OperatorScale( const OperatorScale& tc )
        :
        M_F( tc.M_F ),
        M_alpha( tc.M_alpha ),
        M_Label( tc.M_Label )
    {
        DVLOG(2) << "Copy scale operator " << Label() << " ...\n";
    }

    bool hasInverse() const
    {
        return M_F->hasApply()*( M_alpha != 0 );
    }

    bool hasApply() const
    {
        return M_F->hasApply();
    }



    int Apply( const Epetra_MultiVector & X, Epetra_MultiVector & Y ) const
    {
        DVLOG(2) << "Apply scale operator " << Label() << "\n";

        M_F->Apply( X,Y );

        Y.Scale( M_alpha );

        return !hasApply();
    }

    int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        DVLOG(2) << "ApplyInverse scale operator " << Label() << "\n";

        FEELPP_ASSERT( hasInverse() && ( M_alpha != 0 ) ).error( "This operator cannot be inverted." );

        Epetra_MultiVector Z( X );
        Z.Scale( 1./M_alpha );

        M_F->ApplyInverse( Z,Y );

        return !hasInverse();
    }

    // other function
    int SetUseTranspose( bool UseTranspose )
    {
        return( false );
    }

    double NormInf() const
    {
        return 0;
    }

    const char * Label () const
    {
        return M_Label.c_str();
    }

    bool UseTranspose() const
    {
        return( false );
    }

    bool HasNormInf () const
    {
        return( false );
    }


    const Epetra_Comm & Comm() const
    {
        return( M_F->Comm() );
    }

    const Epetra_Map & OperatorDomainMap() const
    {
        return( M_F->OperatorDomainMap() );
    }

    const Epetra_Map & OperatorRangeMap() const
    {
        return( M_F->OperatorRangeMap() );
    }

    ~OperatorScale()
    {
        //DVLOG(2) << "Destroyed scale operator: " << Label() << " ...\n";
    };

private:

    op1_ptrtype M_F;

    double M_alpha;

    std::string M_Label;
};






template< typename operator_type >
class OperatorFree : public Epetra_Operator
{
public:

    typedef boost::shared_ptr<operator_type> operator_ptrtype;

    typedef Epetra_Operator prec_type;
    typedef boost::shared_ptr<prec_type> prec_ptrtype;

    typedef SolverLinearTrilinos<double> solver_type;
    typedef boost::shared_ptr<solver_type> solver_ptrtype;

    typedef solver_type::real_type real_type;

    typedef Teuchos::ParameterList list_type;

    OperatorFree( operator_ptrtype F )
        :
        M_op( F ),
        M_Solver( new solver_type() ),
        M_hasInverse( 0 ),
        M_hasApply( F->hasApply() ),
        M_useTranspose( F->UseTranspose() ),
        M_label( F->Label() )
    {
        DVLOG(2) << "Create operator " << Label() << " ...\n";
    }

    OperatorFree( operator_ptrtype F, prec_ptrtype Prec, list_type options )
        :
        M_op( F ),
        M_Solver( new solver_type() ),
        M_hasInverse( 1 ),
        M_hasApply( F->hasApply() ),
        M_useTranspose( F->UseTranspose() ),
        M_label( F->Label() ),
        M_Prec( Prec )
    {
        DVLOG(2) << "Create operator " << Label() << " ...\n";
        M_tol = options.get( "tol", 1e-7 );
        M_maxiter = options.get( "max_iter", 100 );

        M_Solver->setOptions( options );
    }


    OperatorFree( const OperatorFree& tc )
        :
        M_op( tc.M_op ),
        M_hasInverse( tc.M_hasInverse ),
        M_hasApply( tc.M_hasApply ),
        M_useTranspose( tc.M_useTranspose ),
        M_Solver( tc.M_Solver ),
        M_label( tc.M_label ),
        M_Prec( tc.M_Prec )
    {
        DVLOG(2) << "Copy operator " << Label() << " ...\n";
    }


    bool hasInverse() const
    {
        return M_hasInverse;
    }

    bool hasApply() const
    {
        return M_hasApply;
    }


    int Apply( const Epetra_MultiVector & X, Epetra_MultiVector & Y ) const
    {
        DVLOG(2) << "Apply operator " << Label() << "\n";
        M_op->Apply( X,Y );
        DVLOG(2) << "Finished Apply operator " << Label() << "\n";

        return !hasApply();
    }

    int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        DVLOG(2) << "ApplyInverse operator " << Label() << "\n";

        FEELPP_ASSERT( hasInverse() ).error( "This operator cannot be inverted." );

        std::pair<unsigned int, real_type> result = M_Solver->solve( M_op, M_Prec, Y, X, M_tol, M_maxiter );

        DVLOG(2) << "Finished ApplyInverse operator " << Label() << "\n";
        return !hasInverse();
    }

    // other function
    int SetUseTranspose( bool UseTranspose = 0 )
    {
        M_useTranspose = UseTranspose;

        return UseTranspose;
    }

    double NormInf() const
    {
        return 0;
    }

    const char * Label () const
    {
        return( M_label.c_str() );
    }

    bool UseTranspose() const
    {
        return M_useTranspose;
    }

    bool HasNormInf () const
    {
        return( true );
    }

    const Epetra_Comm & Comm() const
    {
        return( M_op->OperatorDomainMap().Comm() );
    }

    const Epetra_Map & OperatorDomainMap() const
    {
        return( M_op->OperatorDomainMap() );
    }

    const Epetra_Map & OperatorRangeMap() const
    {
        return( M_op->OperatorRangeMap() );
    }

    ~OperatorFree()
    {
        DVLOG(2) << "Destroyed operator: " << Label() << " ...\n";
    };

private:

    operator_ptrtype M_op;

    solver_ptrtype M_Solver;

    bool M_hasInverse, M_hasApply, M_useTranspose;

    int M_maxiter;
    double M_tol;

    std::string M_label;

    prec_ptrtype M_Prec;
};








} // Feel

#endif // FEELPP_HAS_TRILINOS_EPETRA
#endif /* __operator_trilinos_matrix_H */
