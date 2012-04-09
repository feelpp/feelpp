/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-09-29

 Copyright (C) 2004 EPFL

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
   \file solverbase.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-09-22
*/
#ifndef SolverBase_H
#define SolverBase_H 1

class GetPot;

namespace Feel
{

/*!
  \class SolverBase
  \brief abstract base class for wrapped (linear) solvers

  @author Christoph Winkelmann
*/
class SolverBase
{
public:

    /** @name Typedefs
     */
    //@{

    typedef double value_type;
    typedef SolverBase solver_type;
    typedef Vector array_type;

    //@}

    /** @name Constructor, destructor
     */
    //@{

    //! create a new instance
    virtual static SolverBase* New() = 0;

    //! destructor
    virtual ~SolverBase() = 0;

    //@}

    /** @name Accessors
     */
    //@{

    virtual double residualNorm() const = 0;

    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix
    template <typename MatrixType>
    virtual void setMatrix( MatrixType const& newMatrix ) = 0;

    virtual void setTolerance( double newTolerance ) = 0;

    //@}

    /** @name  Methods
     */
    //@{

    /*
      solve the problem \f$ A x = b \f$

      \c A has been entered via \c setMatrix .

    */
    virtual void solve( array_type& x, array_type const& b ) = 0;

    //@}

    /*! Sets options from data file for this solver.
     *  @param dataFile GetPot object containing the options from the data file
     *  @param section section in the GetPot object containing the solver stuff
     */
    virtual void setOptionsFromGetPot( GetPot const& dataFile,
                                       std::string section ) = 0;
};

}
#endif /* SolverBase_H */
