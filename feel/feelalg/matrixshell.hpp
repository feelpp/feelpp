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
   \file matrixshell.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-12-28
 */
#ifndef __MatrixShell_H
#define __MatrixShell_H 1

#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/vector.hpp>

namespace Feel
{
/**
 * \class MatrixShell
 * \brief matrices that define its action against a vector
 *
 * Generic shell matrix, i.e. a matrix that does not define anything
 * but its action on a vector. This class contains pure virtual
 * members that must be overloaded in derived classes.
 *
 * @author Christophe Prud'homme
 * @see MatrixSparse
 */
template <typename T>
class MatrixShell
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;
    typedef Vector<T> vector_type;
    typedef boost::shared_ptr<Vector<T> > vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixShell() {}
    virtual ~MatrixShell() {}

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
    virtual size_type size1 () const = 0;

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    virtual size_type size2 () const = 0;


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! copies the diagonal of the matrix into \p v.
    virtual void diagonal( vector_type& v ) const = 0;

    void diagonal ( boost::shared_ptr<Vector<T> >& dest ) const
    {
        this->diagonal( *dest );
    }

    //! Multiplies the matrix with in and stores the result in out.
    virtual void multVector( vector_type const& in, vector_type& out ) const = 0;

    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    void multVector ( const boost::shared_ptr<Vector<T> >& arg,
                      boost::shared_ptr<Vector<T> >& dest ) const
    {
        this->multVector( *arg, *dest );
    }

    //! Multiplies the matrix with in and adds the result to out.
    virtual void multAddVector( vector_type const& in, vector_type& out ) const = 0;

    void multAddVector ( const boost::shared_ptr<Vector<T> >& arg,
                         boost::shared_ptr<Vector<T> >& dest ) const
    {
        this->multAddVector( *arg, *dest );
    }
    //@}



protected:

private:

};
} // Feel
#endif /* __MatrixShell_H */
