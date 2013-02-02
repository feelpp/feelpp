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
   \file solverlineartrilinos.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 21-05-2007
*/

#ifndef __trilinos_linear_solver_H
#define __trilinos_linear_solver_H 1



#include <feel/feelalg/solverlinear.hpp>

#include <feel/feelalg/backendtrilinos.hpp>

#if defined( FEELPP_HAS_TRILINOS_AZTECOO )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#include <AztecOO_config.h>
#include <AztecOO.h>
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION

namespace Feel
{

class BackendTrilinos;
class OperatorMatrix;

template< typename T = double >
class SolverLinearTrilinos : public SolverLinear<T>
{
    typedef SolverLinear<T> super;

public:

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;

    //typedef BackendTrilinos backend_type;

    typedef MatrixEpetra sparse_matrix_type;
    typedef VectorEpetra<value_type> vector_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    typedef Teuchos::ParameterList list_type;
    //typedef typename backend_type::list_type list_type;


    SolverLinearTrilinos()
    {
    }

    SolverLinearTrilinos( po::variables_map const& vm )
        :
        super( vm )
    {
    }

    ~SolverLinearTrilinos()
    {
        this->clear();
    }

    void clear()
    {
        if ( this->initialized() )
        {
            this->setInitialized( false );
        }
    }

    void init()
    {
        if ( !this->initialized() )
        {
            this->setInitialized( true );

            M_Solver.SetParameters( M_List, true );

            // AztecOO defined a certain number of output parameters, and store them
            // in a double vector called status.
            double status[AZ_STATUS_SIZE];
            M_Solver.GetAllAztecStatus( status );
        }
    }

    void setOptions( list_type _list )
    {
        M_List = _list;
    }

    list_type& getOptions()
    {
        return M_List;
    }

    //std::pair<unsigned int, real_type>
    boost::tuple<bool,unsigned int, real_type>
    solve ( MatrixSparse<T>  const& matrix,
            Vector<T> & solution,
            Vector<T> const& rhs,
            const double tol,
            const unsigned int m_its,
            bool transpose = false )
    {
        DVLOG(2) << "Matrix solver...\n";

        setRHS( rhs );
        setLHS( solution );
        setUserOperator( matrix );

        M_Solver.SetParameters( M_List, true );
        M_Solver.Iterate( m_its, tol );

        //return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
#warning todo!
        return boost::make_tuple( true, M_Solver.NumIters(), M_Solver.TrueResidual() );
    }

    //std::pair<unsigned int, real_type>
    boost::tuple<bool,unsigned int, real_type>
    solve (  MatrixSparse<T> const& matrix,
             MatrixSparse<T> const& preconditioner,
             Vector<T>& solution,
             Vector<T> const& rhs,
             const double tol,
             const unsigned int m_its,
             bool transpose = false )
    {
        DVLOG(2) << "Matrix solver with preconditioner...\n";

        setRHS( rhs );
        setLHS( solution );
        setUserOperator( matrix );

        M_Solver.SetParameters( M_List, true );
        M_Solver.Iterate( m_its, tol );

        //return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
#warning todo!
        return boost::make_tuple( true, M_Solver.NumIters(), M_Solver.TrueResidual() );

    }

    std::pair<unsigned int, real_type>
    solve ( boost::shared_ptr<OperatorMatrix>  const& op,
            Vector<T> & solution,
            Vector<T> const& rhs,
            const double tol,
            const unsigned int m_its,
            bool transpose = false )
    {
        DVLOG(2) << "Operator solver...\n";

        setRHS( rhs );
        setLHS( solution );
        setUserOperator( op );

        M_Solver.SetParameters( M_List, true );
        M_Solver.Iterate( m_its, tol );

        return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
    }

    template< typename operator_type >
    std::pair<unsigned int, real_type>
    solve ( operator_type  const& op,
            Epetra_MultiVector & solution,
            Epetra_MultiVector const& rhs,
            const double tol,
            const unsigned int m_its )
    {
        DVLOG(2) << "Operator solver...\n";

        setRHS( rhs );
        setLHS( solution );
        setUserOperator( op );

        M_Solver.SetParameters( M_List, true );
        M_Solver.Iterate( m_its, tol );

        return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
    }


    template< typename op1_type, typename op2_type >
    std::pair<unsigned int, real_type>
    solve ( op1_type const& op1,
            op2_type const& op2,
            Epetra_MultiVector & solution,
            Epetra_MultiVector const& rhs,
            const double tol,
            const unsigned int m_its )
    {
        DVLOG(2) << "Operator solver with preconditioner...\n";

        setRHS( rhs );
        setLHS( solution );
        setUserOperator( op1 );

        if ( op2.get() != 0 )
            M_Solver.SetPrecOperator( &( *op2 ) );

        M_Solver.SetParameters( M_List, true );
        M_Solver.Iterate( m_its, tol );

        return std::make_pair( M_Solver.NumIters(), M_Solver.TrueResidual() );
    }

private:

    AztecOO M_Solver;

    list_type M_List;

    // sets the user specified operator
    template< typename operator_type >
    void setUserOperator( operator_type const& D )
    {
        setUserOperator( mpl::bool_<
                         mpl::or_< boost::is_same< operator_type, MatrixSparse<T> >,
                         boost::is_same< operator_type, MatrixEpetra > >::value >(),
                         D );
    }

    void setUserOperator( mpl::bool_<true>, MatrixSparse<T> const& A )
    {
        DVLOG(2) << "Set matrix operator...\n";
        sparse_matrix_type* A_ptr = const_cast<sparse_matrix_type *>( dynamic_cast< sparse_matrix_type const*>( &A ) );

        M_Solver.SetUserMatrix( &A_ptr->mat() );
    }


    // case the operator is passed as an Epetra_Operator
    template< typename operator_type >
    void setUserOperator( mpl::bool_<false>, operator_type const& A )
    {
        DVLOG(2) << "Set epetra operator...\n";

        M_Solver.SetUserOperator( &( *A ) );
    }


    void setRHS( Vector<T> const& x )
    {
        VectorEpetra<T>* aux = const_cast< VectorEpetra<T>* >( dynamic_cast<VectorEpetra<T> const*>( &x ) );

        M_Solver.SetRHS( &aux->vec() );
    }

    void setLHS( Vector<T>& x )
    {
        VectorEpetra<T>* aux = dynamic_cast< VectorEpetra<T>* >( &x );

        M_Solver.SetLHS( &aux->vec() );
    }

    void setRHS( Epetra_MultiVector const& x )
    {
        Epetra_MultiVector* aux = const_cast< Epetra_MultiVector* >( dynamic_cast<Epetra_MultiVector const*>( &x ) );

        M_Solver.SetRHS( &( *aux ) );
    }

    void setLHS( Epetra_MultiVector& x )
    {
        M_Solver.SetLHS( &x );
    }
};


} // Feel

#endif /* FEELPP_HAS_TRILINOS_AZTECOO */
#endif /* __trilinos_linear_solver_H */


