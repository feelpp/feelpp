/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-09-30

  Copyright (C) 2012 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file boundarycondition.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-09-30
 */
#ifndef __BoundaryCondition_H
#define __BoundaryCondition_H 1

namespace Feel
{
namespace detail
{
/**
 * \class BoundaryCondition
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class BoundaryCondition
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    BoundaryCondition();
    //! copy constructor
    BoundaryCondition( BoundaryCondition const & );
    //! destructor
    ~BoundaryCondition();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    BoundaryCondition& operator=( BoundaryCondition const & o)
        {
            if (this != &o )
            {
            }
            return *this;
        }
    //@}

    /** @name Accessors
     */
    //@{


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

};
} // detail

class BoundaryCondition : public detail::BoundaryCondition
{
public:
    BOOST_PARAMETER_CONSTRUCTOR(
        BoundaryCondition, (detail::BoundaryCondition), tag,
        (required
         (marker,*))
        (optional
         (bc,*)))
};



} // Feel
#endif /* __BoundaryCondition_H */
