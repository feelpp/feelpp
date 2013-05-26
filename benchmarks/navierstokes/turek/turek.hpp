/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-10

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)

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
   \file turek.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-10
 */
#ifndef _TUREK_HPP_
#define _TUREK_HPP_ 1

#include <iostream>
#include <fstream>
#include <iomanip>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feeldiscr/bdf2.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <data.hpp>

namespace Feel
{
/**
 * Turek Model
 *
 */
template<int Dim, int Order, int GeoOrder>
class Turek   : public Data
{
    typedef Data super;
public:
#define Entity Simplex
    /**
     * Typedefs  and Constants
     */
    static const uint16_type imOrder = 2*Order;
    static const uint16_type uOrder = Order;
    static const uint16_type pOrder = Order-1;

    typedef Turek<Dim,Order,GeoOrder> self_type;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim, GeoOrder> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    /* velocity */
    typedef fusion::vector<Lagrange<Order, Vectorial> > velocity_basis_type;
    typedef FunctionSpace<mesh_type, velocity_basis_type, value_type> velocity_functionspace_type;
    typedef boost::shared_ptr<velocity_functionspace_type> velocity_functionspace_ptrtype;
    typedef typename velocity_functionspace_type::element_type velocity_element_type;
    typedef boost::shared_ptr<velocity_element_type> velocity_element_ptrtype;

    /* pressure */
    typedef fusion::vector<Lagrange<Order-1, Scalar> > pressure_basis_type;
    typedef FunctionSpace<mesh_type, pressure_basis_type, value_type> pressure_functionspace_type;
    typedef boost::shared_ptr<pressure_functionspace_type> pressure_functionspace_ptrtype;
    typedef typename pressure_functionspace_type::element_type pressure_element_type;
    typedef boost::shared_ptr<pressure_element_type> pressure_element_ptrtype;
    /* fluid */
    typedef fusion::vector<Lagrange<Order, Vectorial>,
            Lagrange<Order-1, Scalar> > fluid_basis_type;

    typedef FunctionSpace<mesh_type, fluid_basis_type, value_type> fluid_functionspace_type;
    typedef boost::shared_ptr<fluid_functionspace_type> fluid_functionspace_ptrtype;
    typedef typename fluid_functionspace_type::element_type fluid_element_type;
    typedef boost::shared_ptr<fluid_element_type> fluid_element_ptrtype;
    typedef typename fluid_element_type::template sub_element<0>::type fluid_element_0_type;
    typedef typename fluid_element_type::template sub_element<1>::type fluid_element_1_type;

    /* Operators */
    typedef OperatorLinear<fluid_functionspace_type, fluid_functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef OperatorLinear<velocity_functionspace_type, velocity_functionspace_type> op_vector_type;
    typedef boost::shared_ptr<op_vector_type> op_vector_ptrtype;
    typedef OperatorLinear<pressure_functionspace_type, pressure_functionspace_type> op_scalar_type;
    typedef boost::shared_ptr<op_scalar_type> op_scalar_ptrtype;
    typedef FsFunctionalLinear<fluid_functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;
    typedef OperatorLagrangeP1<velocity_functionspace_type> velocity_oplagp1_type;
    typedef boost::shared_ptr<velocity_oplagp1_type> velocity_oplagp1_ptrtype;

    typedef OperatorLagrangeP1<pressure_functionspace_type> pressure_oplagp1_type;
    typedef boost::shared_ptr<pressure_oplagp1_type> pressure_oplagp1_ptrtype;

    typedef FunctionSpace<typename velocity_oplagp1_type::image_mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /* time */
    typedef Bdf<fluid_functionspace_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    /* export */
    typedef Exporter<typename velocity_oplagp1_type::image_mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * constructor: Xh space and some space functions are initialized
     */
    Turek( po::variables_map const& vm );

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh();


    void run();

    double normL2Div( fluid_element_type& U ) const;

    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J );

    void updateLinearOperatorsBdf1( fluid_element_type& U );
    void updateLinearOperatorsBdf2( fluid_element_type& U );

private:

    void initLinearOperators();


    /**
     * solve the system
     */
    void nlsolve( sparse_matrix_ptrtype& D, fluid_element_type& u, vector_ptrtype& F );

    void solve( sparse_matrix_ptrtype& D, fluid_element_type& u, vector_ptrtype& F );
    void solve( sparse_matrix_ptrtype& D, velocity_element_type& u, vector_ptrtype& F );
    void solve( sparse_matrix_ptrtype& D, pressure_element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, fluid_element_type& u );


private:

    backend_ptrtype M_backend;
    backend_ptrtype M_backend_symm_v;
    backend_ptrtype M_backend_symm_s;

    fluid_functionspace_ptrtype M_Xh;
    fluid_element_ptrtype Un1;
    fluid_element_ptrtype Un;

    op_vector_ptrtype M_mass_v;
    op_scalar_ptrtype M_mass_s;
    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;
    velocity_oplagp1_ptrtype M_velocity_oplagp1;
    pressure_oplagp1_ptrtype M_pressure_oplagp1;
    funlin_ptrtype M_residual;
    funlin_ptrtype M_stokes_rhs;

    export_ptrtype exporter;

    double time;
    int M_time_order_var;
    bdf_ptrtype M_bdf;

    std::ofstream M_data;
}; // Turek


} // Feel

#endif
