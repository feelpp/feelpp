/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-11-22

  Copyright (C) 2013 Université de Strasbourg

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
   \file drivencavity.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-11-22
 */
#ifndef FEELPP_DRIVENCAVITY_HPP_H
#define FEELPP_DRIVENCAVITY_HPP_H 1

#include <feel/feel.hpp>

namespace Feel
{
template<int Dim>
class DrivenCavity : public Application
{
    typedef Application super;
public:


    typedef double value_type;

    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/

    typedef Lagrange<2, Vectorial> basis_u_type;
    typedef Lagrange<1, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    //#if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
    //#else
    //typedef bases<basis_u_type,basis_p_type> basis_type;
    //#endif


    typedef bases<basis_u_type> basis_type_U;
    typedef FunctionSpace<mesh_type, basis_type_U> space_type_U;
    typedef boost::shared_ptr<space_type_U> space_ptrtype_U;

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef typename space_type::element_type element_type;
    //# endmarker3 #

    /* export */
    typedef Exporter<mesh_type> export_type;

    DrivenCavity( );

    // init mesh and space
    FEELPP_DONT_INLINE void init();

    FEELPP_DONT_INLINE void run();

    FEELPP_DONT_INLINE void exportResults( element_type const& U );

    FEELPP_DONT_INLINE void Jacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    FEELPP_DONT_INLINE void Residual(const vector_ptrtype& X, vector_ptrtype& R);

private:

    double meshSize;

    double Re;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Vh;

    boost::shared_ptr<export_type> exporter;
}; // DrivenCavity
} // namespace Feel
#endif /* FEELPP_DRIVENCAVITY_HPP_H */
