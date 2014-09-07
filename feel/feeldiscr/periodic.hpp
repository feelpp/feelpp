/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-09-30

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file periodic.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-09-30
 */
#ifndef __Periodic_H
#define __Periodic_H 1

#include <feel/feelalg/glas.hpp>


namespace Feel
{
namespace detail
{
struct periodicity_base {};
}

/**
 * \class Periodic
 * \brief Periodic class holder
 *
 the periodic condition is set both as template parameter for the
 functionspace and constructor parameter. In the former the tags that
 are periodic/linked together are set as well as the fact that we have
 a periodic condition. In the former we set the translation Trans that
 allow to link the dofs of tag2 and tag1.

 This data structure is then passed to the Dof class as well that will
 actually do the real work to treat the periodic condition. Its job is
 to set the same dof identifier for the dof on Tag2 as the Dof on Tag1
 with respect to the translation Trans. The 'periodic condition' check
 must be done for all entities of an element (vertex,edge,face,volume)
 even for the volume (think P0 discontinuous). During the check the
 dof points coordinates of Tag1 are put into a data structure DS as
 well as the corresponding dof identifier for rapid localisation then
 the dof points are translation by Trans and look for in DS. One and
 only one point must be found, then the dof identifier of the original
 point in Tag2 is set to the one in Tag1.

 Issues:

 - the dof are constructed by entity (vertex,edge,face,volume) so
 for each entity one must check that the dof points belongs to
 either Tag1 or Tag2. This is not yet clear how to do this.

 - the data structure DS must be efficient to find the proper point in
 Tag1 list, see sedgewick for reference.

 - the implementation of the periodic logic could be done either in
 dof or in periodic, with the former we have access to all information
 a priori but Dof get some added weight, with the latter we need to
 pass extra information and we make sure that the Dof stays more or
 less as it is. The other issue is that we must implement the same
 interface for all related periodic conditions.


 * @author Christophe Prud'homme
 * @see
 */
template<typename T = double >
class Periodic : public Feel::detail::periodicity_base
{
public:

    /** @name Constants
     */
    //@{

    typedef typename node<T>::type node_type;

    //@}

    /** @name Constants
     */
    //@{

    static const bool is_periodic = true;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Periodic() : M_tag1(invalid_uint16_type_value), M_tag2(invalid_uint16_type_value), M_trans() {}
    Periodic( uint16_type tag1, uint16_type tag2, node_type const& trans ) : M_tag1( tag1 ), M_tag2( tag2 ), M_trans( trans ) {}
    Periodic( Periodic const & p ) : M_tag1( p.M_tag1 ), M_tag2( p.M_tag2 ), M_trans( p.M_trans ) {}
    ~Periodic() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //! return whether the condition is periodic or not
    static bool isPeriodic()
        {
            return is_periodic;
        }

    //! return the translation condition that should be applied on Tag2
    node_type const& translation()
        {
            return M_trans;
        }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    uint16_type tag1() const { return M_tag1; }
    uint16_type tag2() const { return M_tag2; }

    //@}



protected:

private:
    uint16_type M_tag1;
    uint16_type M_tag2;
    node_type M_trans;
};
/**
 * \class NoPeriodicity
 * \brief NoPeriodicity class holder
 *
 * This class the default parameter for function spaces without
 * periodic boundary conditions
 *
 *
 * @author Christophe Prud'homme
 * @see
 */
class NoPeriodicity : public Feel::detail::periodicity_base
{
public:

    /** @name Constants
     */
    //@{

    static const bool is_periodic = false;
    //static const uint16_type tag1 = invalid_uint16_type_value;
    //static const uint16_type tag2 = invalid_uint16_type_value;

    typedef node<double>::type node_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    NoPeriodicity() {}
    /* for compatibility reasons with Periodic( int, int, node_type) */
    NoPeriodicity( uint16_type tag1, uint16_type tag2, node_type const& trans ) {}

    //@}

    /** @name Accessors
     */
    //@{

    //! return whether the condition is periodic or not
    static bool isPeriodic()
        {
            return is_periodic;
        }

    //! return the translation condition that should be applied on Tag2
    node_type translation()
        {
            return node_type();
        }

    uint16_type tag1() const { return invalid_uint16_type_value; }
    uint16_type tag2() const { return invalid_uint16_type_value; }

    //@}
};

}
#endif /* __Periodic_H */
