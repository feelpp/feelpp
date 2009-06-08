/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
       Date: 2009-05-25

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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

#include <life/lifealg/glas.hpp>
#include <life/lifealg/solvernonlineartrilinos.hpp>

namespace Life
{
template <typename T>
bool SolverNonLinearTrilinos<T>::computeF( const Epetra_Vector & x,
                                           Epetra_Vector & f,
                                           NOX::Epetra::Interface::Required::FillType ft )
{
    boost::shared_ptr<Vector<double> > X( new VectorEpetra<double>(&x));
    boost::shared_ptr<VectorEpetra<double> > F;

    if (this->residual != NULL) this->residual (X, (boost::shared_ptr<Vector<double> >&) F );
    f=*(F->epetraVector().get());
    return true;
}
template <typename T>
bool SolverNonLinearTrilinos<T>::computeJacobian( const Epetra_Vector & x,
                                                  Epetra_Operator & Jac )
{
    boost::shared_ptr<Vector<double> > X( new VectorEpetra<double>(&x));
    boost::shared_ptr<MatrixEpetra> M_Jac;
    //boost::shared_ptr<MatrixSparse<double> > M_Jac;
    if (this->jacobian != NULL) this->jacobian (X, (boost::shared_ptr<MatrixSparse<double> >&) M_Jac );

    Jac = M_Jac->mat();
    return true;
}
template <typename T>
bool SolverNonLinearTrilinos<T>::computePrecMatrix( const Epetra_Vector & x,
                                                    Epetra_RowMatrix & M )
{
    return true;
}
template <typename T>
bool SolverNonLinearTrilinos<T>::computePreconditioner( const Epetra_Vector & x,
                                                        Epetra_Operator & O )
{
    return true;
}

// SolverNonLinearTrilinos<> methods
template <typename T>
void SolverNonLinearTrilinos<T>::clear ()
{
    if (this->initialized())
        {
            this->M_is_initialized = false;
        }
}

template <typename T>
void SolverNonLinearTrilinos<T>::init ()
{
    if (!this->initialized())
        {
            this->M_is_initialized = true;


        }
}

template <typename T>
std::pair<unsigned int, typename SolverNonLinearTrilinos<T>::real_type>
SolverNonLinearTrilinos<T>::solve ( sparse_matrix_ptrtype&  jac_in,  // System Jacobian Matrix
                                    vector_ptrtype& x_in,    // Solution vector
                                    vector_ptrtype& r_in,    // Residual vector
                                    const double,              // Stopping tolerance
                                    const unsigned int)
{
    MatrixEpetra* jac = dynamic_cast<MatrixEpetra *>( jac_in.get() );
    VectorEpetra<double>* x  = dynamic_cast<VectorEpetra<double>*>( x_in.get() );

    // We cast to pointers so we can be sure that they succeeded
    // by comparing the result against NULL.
//    assert(jac != NULL); assert(jac->mat() != NULL);
//    assert(x   != NULL); assert(x->vec()   != NULL);
//    assert(r   != NULL); assert(r->vec()   != NULL);

    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
        Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

    // Set the nonlinear solver method
    nlParams.set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
    printParams.set("Output Precision", 5);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information",
                    NOX::Utils::OuterIteration +
                    NOX::Utils::OuterIterationStatusTest +
                    NOX::Utils::InnerIteration +
                    NOX::Utils::Parameters +
                    NOX::Utils::Details +
                    NOX::Utils::Warning);

    // start definition of nonlinear solver parameters
    // Sublist for line search
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    searchParams.set("Method", "Polynomial");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    dirParams.set("Method", "Newton");

    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", 800);
    lsParams.set("Tolerance", 1e-5);
    lsParams.set("Output Frequency", 50);
    lsParams.set("Aztec Preconditioner", "ilu");

    // -> A : Jacobian for the first iteration
    // -> InitialGuess : first value x0

    boost::shared_ptr<Epetra_Vector> InitialGuess = x->epetraVector();
Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp(((boost::shared_ptr<Epetra_CrsMatrix>)(jac->matrix())).get());

    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = Teuchos::rcp(this);
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = Teuchos::rcp(this);
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                          lsParams,
                                                          iReq,
                                                          iJac,
                                                          A,
                                                          *InitialGuess));

    // Need a NOX::Epetra::Vector for constructor
    NOX::Epetra::Vector noxInitGuess(*InitialGuess, NOX::DeepCopy);
    Teuchos::RCP<NOX::Epetra::Group> grpPtr =
        Teuchos::rcp(new NOX::Epetra::Group(printParams,
                                            iReq,
                                            noxInitGuess,
                                            linSys));

  // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> testNormF =
        Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
    Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters =
        Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
    // this will be the convergence test to be used
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                              testNormF, testMaxIters));

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver->solve();

  if( NOX::StatusTest::Converged  != status )
      cout << "\n" << "-- NOX solver did not converged --" << "\n";
  else
      cout << "\n" << "-- NOX solver converged --" << "\n";

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver->getList().print(cout);

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group & finalGroup =
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector & finalSolution =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  cout << "Computed solution : " << endl;
  cout << finalSolution;

return std::make_pair(1,finalGroup.getNormF());

}

template <typename T>
std::pair<unsigned int, typename SolverNonLinearTrilinos<T>::real_type>
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
