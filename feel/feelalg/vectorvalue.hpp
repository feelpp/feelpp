/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-13

  Copyright (C) 2005,2006 EPFL

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
   \file vectorgmm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-13
 */
#ifndef __VectorValue_H
#define __VectorValue_H 1

#include <set>

#include <boost/numeric/ublas/vector.hpp>


namespace Feel
{
/*!
 * \class VectorValue
 * \brief interface to vector
 *
 * \code
 * VectorValue<T> m;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T>
class VectorValue
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;

    typedef value_type vector_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    VectorValue( value_type acc = value_type( 0 ) )
        :
        M_vec( acc )
    {}
    VectorValue( VectorValue const & m )
        :
        M_vec( m.M_vec )
    {}

    ~VectorValue()
    {
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * \return the value
     */
    value_type& operator()( size_type /*i*/ )
    {
        return M_vec;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns \p m, the row-dimension of
     * the vector where the marix is \f$ M \times N \f$.
     */
    unsigned int size () const
    {
        return 1;
    }

    /**
     * return row_start, the index of the first
     * vector row stored on this processor
     */
    unsigned int rowStart () const
    {
        return 0;
    }

    /**
     * return row_stop, the index of the last
     * vector row (+1) stored on this processor
     */
    unsigned int rowStop () const
    {
        return 0;
    }

    /**
     * \return true if vector is initialized/usable, false otherwise
     */
    bool isInitialized() const
    {
        return true;
    }

    /**
     * \c close the gmm vector, that will copy the content of write
     * optimized vector into a read optimized vector
     */
    void close () const;


    /**
     * see if vector has been closed
     * and fully assembled yet
     */
    bool closed() const
    {
        return true;
    }


    /**
     * Returns the read optimized gmm vector.
     */
    vector_type const& vec () const
    {
        return M_vec;
    }

    /**
     * Returns the read optimized gmm vector.
     */
    vector_type & vec ()
    {
        return M_vec;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear ()
    {
        M_vec = 0;
    }

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ()
    {
        M_vec = 0;
    }

    void zero ( size_type /*start1*/, size_type /*stop1*/ )
    {
        M_vec = 0;
    }

    /**
     * Add \p value to the value already accumulated
     */
    void add ( const unsigned int /*i*/,
               const value_type value )
    {
        M_vec += value;
    }

    /**
     * set to \p value
     */
    void set ( const unsigned int /*i*/,
               const value_type value )
    {
        M_vec = value;
    }



    /**
     * Print the contents of the vector in Matlab's
     * sparse vector forvec. Optionally prints the
     * vector to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL", bool renumber = false ) const;


    //@}



protected:

private:

    /**
     * the gmm sparse vector data structure
     */
    mutable vector_type M_vec;

};

template<typename T>
void
VectorValue<T>::close() const
{
}

template<typename T>
void
VectorValue<T>::printMatlab( const std::string /*filename*/, bool renumber ) const
{
}

} // Feel
#endif /* __VectorValue_H */
