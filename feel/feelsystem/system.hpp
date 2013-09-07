/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file system.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef FEELPP_SYSTEM_HPP
#define FEELPP_SYSTEM_HPP 1

namespace Feel
{
/**
 * \class System
 *  \brief System of PDE associated to a function space
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename SpaceType>
class System
{
public:


    /** @name Typedefs
     */
    //@{

    static const uint16_type Dim = SpaceType::nDim;
    typedef System<SpaceType> system_type;

    typedef typename SpaceType::value_type value_type;

    typedef SpaceType functionspace_type;
    typedef boost::shared_ptr<SpaceType> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    System( functionspace_ptrtype const& Xh, po::variables_map const& vm ) :  M_Xh( Xh ), M_vm( vm ) {}
    System( System const & s ) : M_Xh( s.M_Xh ), M_vm( s.M_vm )  {}
    virtual ~System() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    System& operator=( System const& s )
    {
        if ( this != &s )
        {
            M_Xh = s.M_Xh;
            M_vm = s.M_vm;
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    //! \return the variables map
    po::variables_map const& vm() const
    {
        return M_vm;
    }

    //! \return the function space
    functionspace_ptrtype const& functionSpace() const
    {
        return M_Xh;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the variables map
    void setVm( po::variables_map const& vm )
    {
        M_vm = vm;
    }

    //! set the function space
    void setFunctionSpace( functionspace_ptrtype const& Xh )
    {
        M_Xh = Xh;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Assemble the system
     */
    virtual void assemble() = 0;

    /**
     * solve the system and retrieve the solution in \arg u an
     * element of the function space Xh
     */
    virtual void solve( element_type& u ) = 0;

    //@}



protected:

private:

    po::variables_map M_vm;

    functionspace_ptrtype M_Xh;

};
} // Feel
#endif /* FEELPP_SYSTEM_HPP */
