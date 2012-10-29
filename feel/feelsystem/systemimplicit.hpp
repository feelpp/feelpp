/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file systemimplicit.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */

#ifndef __SystemImplicit_H
#define __SystemImplicit_H 1

#include <feel/feelsystem/systemexplicit.hpp>

namespace Feel
{
/**
 * \class SystemImplicit
 * \brief describes an implicit system
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename SpaceType>
class SystemImplicit : SystemExplicit<SpaceType>
{
    typedef SystemExplicit<SpaceType> super;
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    typedef SystemImplicit<SpaceType> system_type;

    typedef typename super::value_type value_type;
    typedef typename super::functionspace_type functionspace_type;
    typedef typename super::functionspace_type functionspace_ptrtype;
    typedef typename super::element_type element_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    /* vector */
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    SystemImplicit( functionspace_ptrtype const& Xh, po::variables_map const& vm );

    //! copy constructor
    SystemImplicit( SystemImplicit const & si );

    //! destructor
    ~SystemImplicit() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    SystemImplicit& operator=( SystemImplicit const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );
            M_backend = o.M_backend;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! set the backend
    backend_ptrtype& backend()
    {
        return M_backend;
    }

    //! \return the backend
    backend_ptrtype const& backend() const
    {
        return M_backend;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the backend
    void backend( backend_ptrtype const& b )
    {
        M_backend = b;
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

    backend_ptrtype M_backend;


private:

};
template<typename SpaceType>
SystemImplicit<SpaceType>::SystemImplicit( functionspace_ptrtype const& Xh,
        po::variables_map const& vm )
    :
    super( Xh, vm ),
    M_backend( backend_type::build( vm ) )
{}
template<typename SpaceType>
SystemImplicit<SpaceType>::SystemImplicit( SystemImplicit const& sil )
    :
    super( sil ),
    M_backend( sil.M_backend )
{}


} // Feel
#endif /* __SystemImplicit_H */
