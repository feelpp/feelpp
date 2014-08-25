/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file stvenant_kirchhoff.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>

#include <structurebase.hpp>

namespace Feel
{
template<typename A, uint16_type i>
class Tagged : public A
{
public:
    static const uint16_type TAG = i;

};

using namespace Feel::vf;
/**
 * StVenant_Kirchhoff Model
 *
 */
template<int Dim, int Order>
class StVenantKirchhoff : public StructureBase
{
    typedef  StructureBase super;
public:
#define Entity Simplex
    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef StVenantKirchhoff<Dim,Order> self_type;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim, 1,Dim> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<fem::Lagrange<Dim, 0, Scalar, Discontinuous> > > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    typedef Tagged<fem::Lagrange<Dim, Order, Vectorial, Continuous, double, Simplex>, 0> basis_u_type;
    typedef Tagged<fem::Lagrange<Dim, Order, Vectorial, Continuous, double, Simplex>, 1> basis_v_type;
    typedef mpl::vector<basis_u_type, basis_v_type> basis_type;


    typedef FunctionSpace<mesh_type, basis_type, value_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename element_type::template sub_element<0>::type element_u_type;
    typedef typename element_type::template sub_element<1>::type element_v_type;

    typedef FunctionSpace<mesh_type, mpl::vector<basis_u_type>, value_type> displacement_functionspace_type;
    typedef boost::shared_ptr<displacement_functionspace_type> displacement_functionspace_ptrtype;
    typedef typename displacement_functionspace_type::element_type displacement_element_type;
    typedef boost::shared_ptr<displacement_element_type> displacement_element_ptrtype;


    typedef OperatorLinear<functionspace_type,functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    typedef OperatorLinear<displacement_functionspace_type,displacement_functionspace_type> opdisplacement_type;
    typedef boost::shared_ptr<opdisplacement_type> opdisplacement_ptrtype;
    typedef FsFunctionalLinear<displacement_functionspace_type> fundisplacement_type;
    typedef boost::shared_ptr<fundisplacement_type> fundisplacement_ptrtype;

    typedef OperatorLagrangeP1<displacement_functionspace_type> displacement_oplagp1_type;
    typedef boost::shared_ptr<displacement_oplagp1_type> displacement_oplagp1_ptrtype;

    /* time */
    typedef Bdf<functionspace_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef typename export_type::timeset_type timeset_type;

    StVenantKirchhoff( po::variables_map const& vm );

    //! load the mesh
    mesh_ptrtype loadMesh();

    //! run
    void run();


    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J );

    //! init the linear static elasticity problem
    void initElastoStaticProblem();

    // ! init the linear part of the Stvenant_Kirchhoff model
    void initLinearPart();

private:



    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );

    void linearSolve( sparse_matrix_ptrtype& D, displacement_element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& u, displacement_element_type& d );

private:

    backend_ptrtype M_backend;

    double meshSize;

    functionspace_ptrtype M_Xh;
    element_ptrtype un2;
    element_ptrtype un1;
    element_ptrtype un;


    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;
    funlin_ptrtype M_residual;

    displacement_functionspace_ptrtype M_Dh;
    // linear static alasticity bilinear form
    opdisplacement_ptrtype M_opelas;

    // linear static alasticity linear form
    fundisplacement_ptrtype M_lfelas;

    displacement_oplagp1_ptrtype M_displacement_oplagp1;

    export_ptrtype exporter;
    typename export_type::timeset_ptrtype timeSet;

    double E;
    double sigma;
    double mu;
    double lambda;
    double density;
    double gravity;



    double time;
    bdf_ptrtype M_bdf;
}; // StVenantKirchhoff

} // Feel




