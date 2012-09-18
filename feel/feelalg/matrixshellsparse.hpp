/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-12-28

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)

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
   \file matrixshellsparse.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-12-28
 */
#ifndef __MatrixShellSparse_H
#define __MatrixShellSparse_H 1

#include <feel/feelalg/matrixshell.hpp>
#include <feel/feelalg/matrixsparse.hpp>

namespace Feel
{
/**
 *\class MatrixShellSparse
 *\brief Allow all sparse matrices to be shell matrices
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T>
class MatrixShellSparse : public MatrixShell<T>
{
    typedef MatrixShell<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;
    typedef typename super::vector_type vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixShellSparse( sparse_matrix_ptrtype m ) : M_m( m ) {}
    ~MatrixShellSparse() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    virtual size_type size1 () const
    {
        return M_m->size1();
    }

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    virtual size_type size2 () const
    {
        return M_m->size2();
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! copies the diagonal of the matrix into \p v.
    void diagonal( vector_type& v ) const;

    //! Multiplies the matrix with arg and stores the result in dest.
    void multVector( vector_type const& in, vector_type& out ) const;

    //! Multiplies the matrix with arg and adds the result to dest.
    void multAddVector( vector_type const& in, vector_type& out ) const;

    //@}



protected:

private:
    sparse_matrix_ptrtype M_m;
};
}
#endif /* __MatrixShellSparse_H */
