/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-22

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
   \file matrixtriplet.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-22
 */
#ifndef __MatrixTriplet_H
#define __MatrixTriplet_H 1

namespace Feel
{
/*!
  \class MatrixTriplet
  \brief brief description

  @author Christophe Prud'homme
  @see
  @version $Id: devel.el,v 1.1.1.1 2001/05/23 21:11:14 prudhomm Exp $
*/
template<typename T>
class MatrixTriplet
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixTriplet( int nr, int nc,
                   std::vector<int> const& _Ti, std::vector<int> const& _Tj, std::vector<double> const& _Tx )
        :
        M_nr( nr ),
        M_nc( nc ),
        M_Ti( _Ti ),
        M_Tj( _Tj ),
        M_Tx( _Tx )
    {}

    MatrixTriplet( MatrixTriplet const & mt )
        :
        M_nr( mt.M_nr ),
        M_nc( mt.M_nc ),
        M_Ti ( mt.M_Ti ),
        M_Tj ( mt.M_Tj ),
        M_Tx ( mt.M_Tx )
    {}

    ~MatrixTriplet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    int nrows() const
    {
        return M_nr;
    }
    int ncols() const
    {
        return M_nc;
    }
    int nz() const
    {
        return M_Ti.size();
    }
    int const* Ti() const
    {
        return &M_Ti[0];
    }
    int const* Tj() const
    {
        return &M_Tj[0];
    }
    value_type const* Tx() const
    {
        return &M_Tx[0];
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    int M_nr;
    int M_nc;
    std::vector<int> M_Ti;
    std::vector<int> M_Tj;
    std::vector<value_type> M_Tx;
};
}
#endif /* __MatrixTriplet_H */
