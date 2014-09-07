/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Samuel Quinodoz
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-02-25

  Copyright (C) 2007 Samuel Quinodoz
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
#ifndef __Convection_H
#define __Convection_H 1

/**
   \file convection.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Samuel Quinodoz
   \date 2009-02-25
 */
#include <feel/options.hpp>
//#include <feel/feelcore/applicationxml.hpp>
#include <feel/feelcore/application.hpp>

// (non)linear algebra backend
#include <feel/feelalg/backend.hpp>

// quadrature rules
#include <feel/feelpoly/im.hpp>

// function space
#include <feel/feeldiscr/functionspace.hpp>

// linear operators
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

// exporter
#include <feel/feelfilters/exporter.hpp>



// use the Feel namespace
using namespace Feel;
using namespace Feel::vf;

#if !defined( CONVECTION_DIM )
#define CONVECTION_DIM 2
#endif
#if !defined( CONVECTION_ORDER_U )
#define CONVECTION_ORDER_U 2
#endif
#if !defined( CONVECTION_ORDER_P )
#define CONVECTION_ORDER_P 1
#endif
#if !defined( CONVECTION_ORDER_T )
#define CONVECTION_ORDER_T 2
#endif
#if !defined( CRB_SOLVER )
#define CRB_SOLVER 0
#endif

/**
 * \class Convection
 * The class derives from the Application class
 * the template arguments are :
 * \tparam Order_s velocity polynomial order
 * \tparam Order_t temperature polynomial order
 * \tparam Order_p pressure polynomial order
 */
//template< int Order_s, int Order_p, int Order_t >
class Convection : public Application
{
    typedef Application super;
public:

    //typedef Convection<Order_s, Order_p, Order_t> self_type;
    static const int Order_s = CONVECTION_ORDER_U;
    static const int Order_p = CONVECTION_ORDER_P;
    static const int Order_t = CONVECTION_ORDER_T;
    typedef Convection self_type;

    // Definitions pour mesh
    typedef Simplex<CONVECTION_DIM> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    // space and associated elements definitions
    typedef Lagrange<Order_s, Vectorial,Continuous,PointSetFekete> basis_u_type; // velocity space
    typedef Lagrange<Order_p, Scalar,Continuous,PointSetFekete> basis_p_type; // pressure space
    typedef Lagrange<Order_t, Scalar,Continuous,PointSetFekete> basis_t_type; // temperature space

#if defined( FEELPP_USE_LM )
    typedef Lagrange<0, Scalar> basis_l_type; // multipliers for pressure space
    typedef bases< basis_u_type , basis_p_type , basis_t_type,basis_l_type> basis_type;
#else
    typedef bases< basis_u_type , basis_p_type , basis_t_type> basis_type;
#endif

    //! numerical type is double
    typedef double value_type;

    typedef FunctionSpace<mesh_type, basis_t_type> t_space_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type:: sub_element<0>::type element_0_type;
    typedef typename element_type:: sub_element<1>::type element_1_type;
    typedef typename element_type:: sub_element<2>::type element_2_type;
#if defined( FEELPP_USE_LM )
    typedef typename element_type:: sub_element<3>::type element_3_type;
#endif

    typedef OperatorLinear<space_type,space_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<space_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    // Definition pour les exportations
    typedef Exporter<mesh_type> export_type;

    // Constructeur
    Convection( int argc , char** argv , AboutData const& , po::options_description const& );

    // generate the mesh
    Feel::gmsh_ptrtype createMesh();


    // Definition de la procedure pour faire tourner le code
    void run();

    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );


    // Definition de la procedure pour resoudre le systeme lineaire
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );

    // Definition de la procedure pour exporter les solutions
    void exportResults( boost::format, element_type& U, double t );
    void exportResults( element_type& U, int i );

private:
    void initLinearOperator( sparse_matrix_ptrtype& L );
    void initLinearOperator2( sparse_matrix_ptrtype& L );
    void updateJacobian1( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateJacobian2( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
private:

    backend_ptrtype M_backend;

    space_ptrtype Xh;
    boost::shared_ptr<OperatorLagrangeP1<typename space_type::sub_functionspace<2>::type::element_type> > P1h;

    oplin_ptrtype M_oplin;
    funlin_ptrtype M_lf;

    sparse_matrix_ptrtype M_L;
    sparse_matrix_ptrtype M_D;
    vector_ptrtype F;

    //pas de temps
    //value_type dt;

    // Exporters
    boost::shared_ptr<export_type> exporter;

    // Timers
    std::map<std::string,std::pair<boost::timer,double> > timers;

    std::vector <double> Grashofs;
    double M_current_Grashofs;
    double M_current_Prandtl;
};
#endif /* __Convection_H */
