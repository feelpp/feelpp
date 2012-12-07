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
   \file systemimplicitnonlinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __SystemImplicitNonLinear_H
#define __SystemImplicitNonLinear_H 1

#include <feel/feelsystem/systemimplicit.hpp>


namespace Feel
{
/**
 * \class SystemImplicitNonLinear
 * \brief Describes nonlinear implicit systems
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename SpaceType>
class SystemImplicitNonLinear : public SystemImplicit<SpaceType>
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

    typedef SystemImplicitNonLinear<SpaceType> self_type;
    typedef SystemImplicitNonLinear<SpaceType> system_type;

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

    //! default constructor
    SystemImplicitNonLinear( functionspace_ptrtype const& Xh, po::variables_map const& vm );
    //! copy constructor
    SystemImplicitNonLinear( SystemImplicitNonLinear const & );
    //! destructor
    ~SystemImplicitNonLinear();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    SystemImplicitNonLinear& operator=( SystemImplicitNonLinear const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );
            M_J = o.M_J;
            M_R = o.M_R;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return the jacobian matrix
    sparse_matrix_ptrtype const& jacobian() const
    {
        return M_J;
    }

    //! \return the jacobian matrix
    sparse_matrix_ptrtype& jacobian()
    {
        return M_J;
    }

    //! \return the residual
    vector_ptrtype const& residual() const
    {
        return M_R;
    }

    //! \return the residual
    vector_ptrtype & residual()
    {
        return M_R;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! function that update the residual at each non linear iteration
    virtual void updateResidual( const vector_ptrtype& X, vector_ptrtype& R ) = 0;

    //! function that update the jacobian at each non linear iteration
    virtual void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J ) = 0;

    //! default implementation that solves the nonlinear system
    void solve( element_type& u )
    {
        vector_ptrtype U( M_backend->newVector( this->functionSpace() ) );
        *U = u;

        this->updateResidual( U, M_R );
        this->updateJacobian( U, M_J );

        M_backend->nlSolve( M_J, U, M_R, 1e-10, 10 );
        u = *U;
    }
    //@}



protected:

    sparse_matrix_ptrtype M_J;
    vector_ptrtype M_R;

private:

};
template<typename SpaceType>
SystemImplicitNonLinear<SpaceType>::SystemImplicitNonLinear( functionspace_ptrtype const& Xh,
        po::variables_map const& vm )
    :
    super( Xh, vm ),
    M_J( M_backend->newMatrix( Xh, Xh ) ),
    M_R( M_backend->newVector( Xh ) )

{
    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

}
template<typename SpaceType>
SystemImplicitNonLinear<SpaceType>::SystemImplicitNonLinear( SystemImplicitNonLinear const& sil )
    :
    super( sil ),
    M_J( sil.M_J ),
    M_R( sil.M_R )

{}


}
#endif /* __SystemImplicitNonLinear_H */
