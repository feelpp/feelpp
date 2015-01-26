/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-02

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file solvernonlinearpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-02
 */
#ifndef __SolverNonLinearPetsc_H
#define __SolverNonLinearPetsc_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvernonlinear.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>

// Petsc include files.
#if defined( FEELPP_HAS_PETSC_H )

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# include <petsc.h>
# include <petscsnes.h>
}
#else
# include <petscversion.h>
# include <petsc.h>
# include <petscsnes.h>
#endif



/**

 *
 * @author Benjamin Kirk, 2002-2005
 */

namespace Feel
{

/**
 * \class SolverNonLinearPetsc
 * \brief Petsc non linear solvers interface
 *
 * This class provides an interface to PETSc iterative solvers that is
 * compatible with the \p SolverNonLinear<> base class
 *
 * @author Christophe Prud'homme
 */
template<typename T>
class SolverNonLinearPetsc
    :
public SolverNonLinear<T>
{
    typedef SolverNonLinear<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef SolverNonLinearPetsc<T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super::vector_ptrtype vector_ptrtype;

    typedef typename super::dense_matrix_type dense_matrix_type;
    typedef typename super::dense_vector_type dense_vector_type;

    //eigen
    typedef typename super::map_dense_matrix_type map_dense_matrix_type;
    typedef typename super::map_dense_vector_type map_dense_vector_type;

    typedef DataMap datamap_type;
    typedef boost::shared_ptr<datamap_type> datamap_ptrtype;

    typedef typename super::solve_return_type solve_return_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  Constructor. Initializes Petsc data structures
     */
    SolverNonLinearPetsc( std::string const& prefix = "", WorldComm const& worldComm=Environment::worldComm() );
    SolverNonLinearPetsc( SolverNonLinearPetsc const & );

    /**
     * Destructor.
     */
    ~SolverNonLinearPetsc();

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init ();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //@}

    /** @name  Mutators
     */
    //@{
    void setReuse( int jac=1, int prec=1 );
    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Release all memory and clear data structures.
     */
    virtual void clear ();

    /**
     * Call the Petsc solver.  It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     */
    virtual solve_return_type solve ( sparse_matrix_ptrtype&,    // System Jacobian Matrix
            vector_ptrtype&,          // Solution vector
            vector_ptrtype&,          // Residual vector
            const double,        // Stopping tolerance
            const unsigned int ); // N. Iterations

    virtual std::pair<unsigned int, real_type> solve ( dense_matrix_type&,    // System Jacobian Matrix
            dense_vector_type&,          // Solution vector
            dense_vector_type&,          // Residual vector
            const double,        // Stopping tolerance
            const unsigned int ); // N. Iterations

    //use eigen
    virtual std::pair<unsigned int, real_type> solve ( map_dense_matrix_type&,    // System Jacobian Matrix
            map_dense_vector_type&,          // Solution vector
            map_dense_vector_type&,          // Residual vector
            const double,        // Stopping tolerance
            const unsigned int ); // N. Iterations



    //@}
    datamap_type const& mapRow() const
    {
        return *M_mapRow;
    }
    datamap_type const& mapCol() const
    {
        return *M_mapCol;
    }
    datamap_ptrtype const& mapRowPtr() const
    {
        return M_mapRow;
    }
    datamap_ptrtype const& mapColPtr() const
    {
        return M_mapCol;
    }

    void setMapRow( datamap_ptrtype const& d )
    {
        M_mapRow=d;
    }
    void setMapCol( datamap_ptrtype const& d )
    {
        M_mapCol=d;
    }
protected:
    void check( int err ) { CHKERRABORT( this->worldComm().globalComm(), err ); }
private:
    /**
     * Tells PETSC to use the user-specified solver stored in
     * \p _solver_type
     */
    void setPetscNlSolverType ();

    /**
     * Tells PETSC to use the user-specified solver stored in
     * \p _solver_type
     */
    void setPetscKspSolverType ();

    /**
     * Tells PETSC to use the user-specified preconditioner stored in
     * \p _preconditioner_type
     */
    void setPetscPreconditionerType ();


    /**
     * Nonlinear solver context
     */
    SNES M_snes;

    uint16_type M_prec_mat_structure;

    /**
     * Preconditioner context
     */
    PC M_pc;

    /**
     * Krylov subspace context
     */
    KSP M_ksp;


    datamap_ptrtype M_mapRow,M_mapCol;

};

template <typename T>
inline
SolverNonLinearPetsc<T>::SolverNonLinearPetsc( std::string const& prefix, WorldComm const& worldComm )
:
    super( prefix,worldComm ),
    M_mapRow(new datamap_type(worldComm)),
    M_mapCol(new datamap_type(worldComm))
{}



template <typename T>
inline
SolverNonLinearPetsc<T>::~SolverNonLinearPetsc ()
{
    this->clear ();
}

} // Feel

#endif /* FEELPP_HAS_PETSC */
#endif /* __SolverNonLinearPetsc_H */
