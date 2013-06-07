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
   \file systemimplicitlinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __SystemImplicitLinear_H
#define __SystemImplicitLinear_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelsystem/systemimplicit.hpp>

namespace Feel
{
/**
 * \class SystemImplicitLinear
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename SpaceType>
class SystemImplicitLinear : public SystemImplicit<SpaceType>
{
    typedef SystemImplicit<SpaceType> super;
public:


    /** @name Constants
     */
    //@{

    //@}

    /** @name Typedefs
     */
    //@{

    typedef SystemImplicitLinear<SpaceType> system_type;

    typedef typename super::value_type value_type;
    typedef typename super::functionspace_type functionspace_type;
    typedef typename super::functionspace_type functionspace_ptrtype;
    typedef typename super::element_type element_type;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef typename super::sparse_matrix_type sparse_matrix_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super::vector_type vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    SystemImplicitLinear( functionspace_ptrtype const& Xh, po::variables_map const& vm );
    SystemImplicitLinear( SystemImplicitLinear const & sil );
    ~SystemImplicitLinear() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    SystemImplicitLinear& operator=( SystemImplicitLinear const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );


            M_lhs = o.M_lhs;
            M_rhs = o.M_rhs;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return left hand side (lhs) matrix
    sparse_matrix_ptrtype const& lhs() const
    {
        return M_lhs;
    }

    //! \return left hand side (lhs) matrix
    sparse_matrix_ptrtype& lhs()
    {
        return M_lhs;
    }

    //! \return right hand side (rhs) vector
    vector_ptrtype const& rhs() const
    {
        return M_rhs;
    }

    //! \return right hand side (rhs) vector
    vector_ptrtype& rhs()
    {
        return M_rhs;
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
     * solve the linear system
     */
    void solve( element_type& u )
    {
        vector_ptrtype U( this->backend->newVector( _test=this->functionSpace() ) );
        this->backend->solve( M_lhs, M_lhs, U, M_rhs );
        u = *U;
    }

    //@}



protected:

    sparse_matrix_ptrtype M_lhs;
    vector_ptrtype M_rhs;

private:



};
template<typename SpaceType>
SystemImplicitLinear<SpaceType>::SystemImplicitLinear( functionspace_ptrtype const& Xh,
        po::variables_map const& vm )
    :
    super( Xh, vm ),
    M_lhs( this->backend()->newMatrix( Xh, Xh ) ),
    M_rhs( this->backend()->newVector( Xh ) )
{}
template<typename SpaceType>
SystemImplicitLinear<SpaceType>::SystemImplicitLinear( SystemImplicitLinear const& sil )
    :
    super( sil ),
    M_lhs( sil.M_lhs ),
    M_rhs( sil.M_rhs )
{}


}
#endif /* __SystemImplicitLinear_H */
