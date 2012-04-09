/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-04-20

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file dof.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-04-20
 */
#ifndef __Dof_H
#define __Dof_H 1

#include <boost/tuple/tuple.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>

namespace Feel
{
/**
 * \class Dof
 * \brief class that represents a degree of freedom
 *
 * @author Christophe Prud'homme
 * @see
 */
class Dof : public boost::tuple<size_type, int16_type, uint16_type, bool, Marker1>
{
    typedef boost::tuple<size_type, int16_type, uint16_type, bool, Marker1> super;
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
    Dof()
        :
        super()
    {
    }

    /**
     *
     */
    Dof( size_type _index, int16_type _sign, uint16_type _entity, bool _location, size_type _marker  )
        :
        super( boost::make_tuple( _index, _sign, _entity, _location, _marker ) )
    {}

    //! copy constructor
    Dof( Dof const & dof )
        :
        super( dof )
    {}

    //! destructor
    ~Dof()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Dof& operator=( Dof const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    /// \return the global index
    size_type index() const
    {
        return this->get<0>();
    }

    /// \return the sign
    int16_type sign() const
    {
        return this->get<1>();
    }

    /// \return the entity type (0: vertex, 1:edge, 2:face, 3:volume)
    uint16_type entity() const
    {
        return this->get<2>();
    }

    /// \return the location
    bool isOnBoundary() const
    {
        return this->get<3>();
    }

    /// \return the marker
    Marker1 marker() const
    {
        return this->get<4>();
    }

    /**
     * \return the coordinates
     */
    ublas::vector<double> const& coords() const
    {
        return M_coords;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * dof coordinates
     */
    void
    setCoordinates( ublas::vector<double> const& coords )
    {
        M_coords = coords;
    }
    //@}


    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    ublas::vector<double> M_coords;
};

typedef multi_index::multi_index_container<
Dof,
multi_index::indexed_by<

// sort by less<int> on index()
multi_index::ordered_unique<multi_index::const_mem_fun<Dof,
size_type,
&Dof::index> >,

// sort by less<int> on entity()
multi_index::ordered_non_unique<multi_index::tag<detail::by_entity>,
multi_index::const_mem_fun<Dof,
uint16_type,
&Dof::entity> >,
// sort by less<int> on location()
multi_index::ordered_non_unique<multi_index::tag<detail::by_location>,
multi_index::const_mem_fun<Dof,
bool,
&Dof::isOnBoundary> >,
// sort by less<int> on marker()
multi_index::ordered_non_unique<multi_index::tag<detail::by_marker>,
multi_index::const_mem_fun<Dof,
Marker1,
&Dof::marker> >
> > dof_container_type;


} // Feel
#endif /* __Dof_H */
