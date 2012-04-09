/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
       Date: 2009-05-25

  Copyright (C) 2009-2011 Universite Joseph Fourier (Grenoble I)

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
   \file solvernonlineartrilinos.cpp
   \author Florent Vielfaure <florent.vielfaure@gmail.com>
   \date 2009-05-25
 */

#include <NOX.H>
#include <NOX_Epetra.H>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/solvernonlineartrilinos.hpp>

#if defined( FEELPP_HAS_TRILINOS )
namespace Feel
{

class SolverNonLinearTrilinosInterface
    :
public NOX::Epetra::Interface::Required,
public NOX::Epetra::Interface::Jacobian
{

private:
    SolverNonLinearTrilinos<double>* solver;

public:
    SolverNonLinearTrilinosInterface( SolverNonLinearTrilinos<double> * sol )
    {
        solver=sol;
    }

    bool computeF( const Epetra_Vector & x,
                   Epetra_Vector & f,
                   NOX::Epetra::Interface::Required::FillType ft )
    {
        //printf("Entering computeF...\n");
        boost::shared_ptr<Vector<double> > X( new VectorEpetra<double>( &x ) );
        boost::shared_ptr<Vector<double> > F( new VectorEpetra<double>( &x ) );

        if ( solver->residual != NULL ) solver->residual ( X, F );

        f=*( dynamic_cast<VectorEpetra<double>*>( F.get() )->epetraVector() );
        //std::cout << "X.use_count()=" << X.use_count() << std::endl;
        //std::cout << "F.use_count()=" << F.use_count() << std::endl;
        return true;
    }
    bool computeJacobian( const Epetra_Vector & x,
                          Epetra_Operator & Jac )
    {
        //printf("Entering computeJacobian...\n");
        boost::shared_ptr<Vector<double> > X( new VectorEpetra<double>( &x ) );
        //boost::shared_ptr<MatrixEpetra> M_Jac;
        boost::shared_ptr<MatrixSparse<double> > M_Jac;

        if ( solver->jacobian != NULL ) solver->jacobian ( X, M_Jac );

        Jac = dynamic_cast<MatrixEpetra*>( M_Jac.get() )->mat();
        //std::cout << "X.use_count()=" << X.use_count() << std::endl;
        //std::cout << "M_Jac.use_count()=" << M_Jac.use_count() << std::endl;
        //printf("End computeJacobian...\n");
        return true;
    }
    bool computePrecMatrix( const Epetra_Vector & x,
                            Epetra_RowMatrix & M )
    {
        //printf("End computePrecMatrix...\n");
        return true;
    }
    bool computePreconditioner( const Epetra_Vector & x,
                                Epetra_Operator & O )
    {
        //printf("End computePreconditioner...\n");
        return true;
    }
};


// SolverNonLinearTrilinos<> methods
template <typename T>
void SolverNonLinearTrilinos<T>::clear ()
{
    if ( this->initialized() )
    {
        this->M_is_initialized = false;
    }
}

template <typename T>
void SolverNonLinearTrilinos<T>::init ()
{
    if ( !this->initialized() )
    {
        this->M_is_initialized = true;


    }
}

template <typename T>
std::pair<int, typename SolverNonLinearTrilinos<T>::real_type>
SolverNonLinearTrilinos<T>::solve ( sparse_matrix_ptrtype&  jac_in,  // System Jacobian Matrix
                                    vector_ptrtype& x_in,    // Solution vector
                                    vector_ptrtype& r_in,    // Residual vector
                                    const double,              // Stopping tolerance
                                    const unsigned int )
{
    //printf("Entering solve...\n");
    MatrixEpetra* jac = dynamic_cast<MatrixEpetra *>( jac_in.get() );
    VectorEpetra<T>* x  = dynamic_cast<VectorEpetra<T>*>( x_in.get() );

    // We cast to pointers so we can be sure that they succeeded
    // by comparing the result against NULL.
    //    assert(jac != NULL); assert(jac->mat() != NULL);
    //    assert(x   != NULL); assert(x->vec()   != NULL);
    //    assert(r   != NULL); assert(r->vec()   != NULL);

    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
        Teuchos::rcp( new Teuchos::ParameterList );
    Teuchos::ParameterList& nlParams = *( nlParamsPtr.get() );

    // Set the nonlinear solver method
    nlParams.set( "Nonlinear Solver", "Line Search Based" );

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams.sublist( "Printing" );
    printParams.set( "Output Precision", 10 );
    printParams.set( "Output Processor", 0 );
    printParams.set( "Output Information",
                     NOX::Utils::OuterIteration +
                     NOX::Utils::OuterIterationStatusTest +
                     NOX::Utils::InnerIteration +
                     NOX::Utils::Parameters +
                     NOX::Utils::Details +
                     NOX::Utils::Warning );

    // start definition of nonlinear solver parameters
    // Sublist for line search
    Teuchos::ParameterList& searchParams = nlParams.sublist( "Line Search" );
    searchParams.set( "Method", "Polynomial" );

    // Sublist for direction
    Teuchos::ParameterList& dirParams = nlParams.sublist( "Direction" );
    dirParams.set( "Method", "Newton" );

    Teuchos::ParameterList& newtonParams = dirParams.sublist( "Newton" );
    newtonParams.set( "Forcing Term Method", "Constant" );

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist( "Linear Solver" );
    lsParams.set( "Aztec Solver", "GMRES" );
    lsParams.set( "Max Iterations", 800 );
    lsParams.set( "Tolerance", 1e-7 );
    lsParams.set( "Output Frequency", 50 );
    lsParams.set( "Aztec Preconditioner", "ilu" );

    // -> A : Jacobian for the first iteration
    // -> InitialGuess : first value x0
    //printf("convert vectors...\n");
    boost::shared_ptr<Epetra_Vector> InitialGuess = x->epetraVector();

    // has_ownership=false in order to let the matrix jac be destroyed by boost
    // and not by Teuchos::RCP
    Teuchos::RCP<Epetra_CrsMatrix> A =
        Teuchos::rcp( ( ( boost::shared_ptr<Epetra_CrsMatrix> )( jac->matrix() ) ).get(),false );

    //std::cout << "A.has_ownership()=" << A.has_ownership() << std::endl;

    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq =
        Teuchos::rcp( new SolverNonLinearTrilinosInterface( this ) );
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac =
        Teuchos::rcp( new SolverNonLinearTrilinosInterface( this ) );
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp( new NOX::Epetra::LinearSystemAztecOO( printParams,
                      lsParams,
                      iReq,
                      iJac,
                      A,
                      *InitialGuess ) );

    // Need a NOX::Epetra::Vector for constructor
    NOX::Epetra::Vector noxInitGuess( *InitialGuess, NOX::DeepCopy );
    Teuchos::RCP<NOX::Epetra::Group> grpPtr =
        Teuchos::rcp( new NOX::Epetra::Group( printParams,
                      iReq,
                      noxInitGuess,
                      linSys ) );

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> testNormF =
        Teuchos::rcp( new NOX::StatusTest::NormF( 1.0e-7 ) );
    Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters =
        Teuchos::rcp( new NOX::StatusTest::MaxIters( 20 ) );

    // this will be the convergence test to be used
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
        Teuchos::rcp( new NOX::StatusTest::Combo( NOX::StatusTest::Combo::OR,
                      testNormF, testMaxIters ) );

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( grpPtr, combo, nlParamsPtr );

    // Solve the nonlinesar system
    NOX::StatusTest::StatusType status = solver->solve();

    if ( NOX::StatusTest::Converged  != status )
        std::cout << "\n" << "-- NOX solver did not converged --" << "\n";

    else
        std::cout << "\n" << "-- NOX solver converged --" << "\n";

    // Print the answer
    std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
    solver->getList().print( cout );

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group & finalGroup =
        dynamic_cast<const NOX::Epetra::Group&>( solver->getSolutionGroup() );
    const Epetra_Vector & finalSolution =
        ( dynamic_cast<const NOX::Epetra::Vector&>( finalGroup.getX() ) ).getEpetraVector();

    //cout << "Computed solution : " << endl;
    //cout << finalSolution;
    x_in = boost::shared_ptr<VectorEpetra<T> > ( new VectorEpetra<T>( &finalSolution ) );

    //std::cout << "InitialGuess.use_count()=" << InitialGuess.use_count() << std::endl;
    //std::cout << "jac.use_count()=" << jac->matrix().use_count() << std::endl;
    //std::cout << "x_in.use_count()=" << x_in.use_count() << std::endl;

    return std::make_pair( 1,finalGroup.getNormF() );
}

template <typename T>
std::pair<unsigned int, typename SolverNonLinearTrilinos<T>::real_type>
SolverNonLinearTrilinos<T>::solve ( dense_matrix_type&  jac_in,  // System Jacobian Matrix
                                    dense_vector_type& x_in,    // Solution vector
                                    dense_vector_type& r_in,    // Residual vector
                                    const double,              // Stopping tolerance
                                    const unsigned int )
{

}



//------------------------------------------------------------------
// Explicit instantiations
template class SolverNonLinearTrilinos<double>;




} // namespace Feel
#endif
