/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-04-07

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file test_rbspace.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-07
*/

#define BOOST_TEST_MODULE test_rbspace
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feeldiscr/reducedbasisspace.hpp>

/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testrbspace( "RBSpace test options" );
    testrbspace.add_options()
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testrbspace.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_rbspace" ,
                     "test_rbspace" ,
                     "0.2",
                     "nD(n=1,2,3) test context of functionspace",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}



template<int Dim, int Order>
class Model:
    public Simget,
    public boost::enable_shared_from_this< Model<Dim,Order> >
{

public :

    typedef Mesh<Simplex<Dim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type,bases<Lagrange<Order> >, Periodicity <NoPeriodicity> > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    Model()
        :
        Simget()
    {
    }

    void run()
    {
        auto mesh=unitHypercube<Dim>();
        Xh = Pch<Order>( mesh );

        auto RbSpace = RbSpacePch<Order>( this->shared_from_this() , mesh );

        auto basis_x = vf::project( Xh , elements(mesh), Px() );
        auto basis_y = vf::project( Xh , elements(mesh), Py() );
        RbSpace->addPrimalBasisElement( basis_x );
        RbSpace->addPrimalBasisElement( basis_y );

        int rbspace_size = RbSpace->size();

        LOG( INFO ) << " rbspace_size : "<<rbspace_size;

        // FEM context and points
        auto ctxfem = Xh->context();
        std::vector< node_type > vec_t;
        node_type t1(Dim), t2(Dim), t3(Dim);
        /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2;
        /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8;
        /*x*/ t3(0)=0.75;   /*y*/ t3(1)=0.9;
        ctxfem.add( t1 );
        ctxfem.add( t2 );
        ctxfem.add( t3 );

        auto ctxrb = RbSpace->context();
        ctxrb.add( t1 );
        ctxrb.add( t2 );
        ctxrb.add( t3 );

        ctxrb.update();

        auto u = RbSpace->element();
        auto u_ptr = RbSpace->elementPtr();


        /*
         * evaluation at specific points
         */

        int x=0,y=1;

        //test with u = (1 0)
        u.setCoefficient( 0 , 1 );
        auto u_fem = u.expansion();
        auto fem_evaluations = evaluateFromContext( _context=ctxfem , _expr=idv(u_fem) );
        auto rb_evaluations = u.evaluate( ctxrb );
        Eigen::VectorXd true_values( 3 );
        true_values(0)=t1(x); true_values(1)=t2(x); true_values(2)=t3(x);
        double norm_fem_evaluations = fem_evaluations.norm();
        double norm_rb_evaluations = rb_evaluations.norm();
        double true_norm=true_values.norm();

        BOOST_CHECK_SMALL( math::abs(norm_fem_evaluations-true_norm), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_fem_evaluations-norm_rb_evaluations), 1e-14 );

        auto rb_evaluate_from_contextrb = evaluateFromContext( _context=ctxrb , _expr=idv(u) );
        auto rb_evaluate_from_contextfem = evaluateFromContext( _context=ctxfem , _expr=idv(u) );
        double norm_rb_evaluations_from_ctxrb = rb_evaluate_from_contextrb.norm();
        double norm_rb_evaluations_from_ctxfem = rb_evaluate_from_contextfem.norm();

        BOOST_CHECK_SMALL( math::abs(norm_rb_evaluations_from_ctxrb-true_norm), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_rb_evaluations_from_ctxrb-norm_rb_evaluations_from_ctxfem), 1e-14 );

        LOG( INFO ) << "rb unknown : \n"<<u;
        LOG( INFO ) << " rb_evaluations : \n"<<rb_evaluations;
        LOG( INFO ) << " rb evaluate from contextrb :\n"<<rb_evaluate_from_contextrb;
        LOG( INFO ) << " rb evaluate from contextfem :\n"<<rb_evaluate_from_contextfem;


        //grad
        auto grad_fem = evaluateFromContext( _context=ctxfem , _expr=gradv(u_fem) );
        auto grad_rb  = evaluateFromContext( _context=ctxrb , _expr=gradv(u) );
        Eigen::VectorXd true_values_grad( 6 );
        /*du/dx*/true_values_grad( 0 ) = 1; /*du/dy*/true_values_grad( 1 ) = 0;
        /*du/dx*/true_values_grad( 2 ) = 1; /*du/dy*/true_values_grad( 3 ) = 0;
        /*du/dx*/true_values_grad( 4 ) = 1; /*du/dy*/true_values_grad( 5 ) = 0;
        double norm_grad_fem = grad_fem.norm();
        double norm_grad_rb = grad_rb.norm();
        double norm_true_grad = true_values_grad.norm();
        BOOST_CHECK_SMALL( math::abs(norm_grad_fem-norm_true_grad), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_grad_rb-norm_true_grad), 1e-14 );

        //dx
        Eigen::VectorXd true_values_dx( 3 );
        /*du/dx*/true_values_dx( 0 ) = 1;
        /*du/dx*/true_values_dx( 1 ) = 1;
        /*du/dx*/true_values_dx( 2 ) = 1;
        auto dx_fem=evaluateFromContext( _context=ctxfem , _expr=dxv(u_fem) );
        auto dx_rb=evaluateFromContext( _context=ctxrb , _expr=dxv(u) );
        double norm_dx_fem = dx_fem.norm();
        double norm_dx_rb = dx_rb.norm();
        double norm_true_dx = true_values_dx.norm();
        BOOST_CHECK_SMALL( math::abs(norm_dx_fem-norm_true_dx), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_dx_rb-norm_true_dx), 1e-14 );
        LOG( INFO ) << " rb dx from contextrb :\n"<<dx_rb;

        //dy
        Eigen::VectorXd true_values_dy( 3 );
        /*du/dy*/true_values_dy( 0 ) = 0;
        /*du/dy*/true_values_dy( 2 ) = 0;
        /*du/dy*/true_values_dy( 4 ) = 0;
        auto dy_fem=evaluateFromContext( _context=ctxfem , _expr=dyv(u_fem) );
        auto dy_rb=evaluateFromContext( _context=ctxrb , _expr=dyv(u) );
        double norm_dy_fem = dy_fem.norm();
        double norm_dy_rb = dy_rb.norm();
        double norm_true_dy = true_values_dy.norm();
        BOOST_CHECK_SMALL( math::abs(norm_dy_fem-norm_true_dy), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_dy_rb-norm_true_dy), 1e-14 );
        LOG( INFO ) << " rb dy from contextrb :\n"<<dy_rb;



        //test with u = ( 0 1 )
        u.setCoefficient( 0 , 0 );
        u.setCoefficient( 1 , 1 );
        u_fem = u.expansion();
        LOG( INFO ) << "call evaluate from context fem idv(ufem)";
        fem_evaluations = evaluateFromContext( _context=ctxfem , _expr=idv(u_fem) );
        rb_evaluations = u.evaluate( ctxrb );
        true_values(0)=t1(1); true_values(1)=t2(1); true_values(2)=t3(1);
        norm_fem_evaluations = fem_evaluations.norm();
        norm_rb_evaluations = rb_evaluations.norm();
        true_norm=true_values.norm();


        BOOST_CHECK_SMALL( math::abs(norm_fem_evaluations-true_norm), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_fem_evaluations-norm_rb_evaluations), 1e-14 );

        LOG( INFO ) << "call evaluate from context rb idv(u)";
        rb_evaluate_from_contextrb = evaluateFromContext( _context=ctxrb , _expr=idv(u) );
        norm_rb_evaluations_from_ctxrb = rb_evaluate_from_contextrb.norm();
        BOOST_CHECK_SMALL( math::abs(norm_rb_evaluations_from_ctxrb-true_norm), 1e-14 );

        LOG( INFO ) << "rb unknown : \n"<<u;
        LOG( INFO ) << " rb_evaluations : \n"<<rb_evaluations;
        LOG( INFO ) << " rb evaluate from contextrb :\n"<<rb_evaluate_from_contextrb;

        //grad
        grad_fem = evaluateFromContext( _context=ctxfem , _expr=gradv(u_fem) );
        grad_rb  = evaluateFromContext( _context=ctxrb , _expr=gradv(u) );
        /*du/dx*/true_values_grad( 0 ) = 0; /*du/dy*/true_values_grad( 1 ) = 1;
        /*du/dx*/true_values_grad( 2 ) = 0; /*du/dy*/true_values_grad( 3 ) = 1;
        /*du/dx*/true_values_grad( 4 ) = 0; /*du/dy*/true_values_grad( 5 ) = 1;
        norm_grad_fem = grad_fem.norm();
        norm_grad_rb = grad_rb.norm();
        norm_true_grad = true_values_grad.norm();
        BOOST_CHECK_SMALL( math::abs(norm_grad_fem-norm_true_grad), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(norm_grad_rb-norm_true_grad), 1e-14 );
        LOG( INFO ) << " rb grad from contextrb :\n"<<grad_rb;

        /*
         * test with lambda expression
         * integrate over the domain
         */
        auto LambdaIntegrate = integrate( _range=elements(mesh), _expr=_e1 );
        auto IntegrateFem = integrate( _range=elements(mesh), _expr=idv(u_fem) );
        auto IntegrateRb = integrate( _range=elements(mesh), _expr=idv(u) );
        double value_integrate_fem = IntegrateFem.evaluate()(0,0);
        double value_integrate_rb = IntegrateRb.evaluate()(0,0);
        double value_lambda_integrate_fem = LambdaIntegrate( idv( u_fem ) ).evaluate()(0,0);
        double value_lambda_integrate_rb = LambdaIntegrate( idv( u ) ).evaluate()(0,0);
        double true_integrate = integrate( _range=elements(mesh) , _expr=Py() ).evaluate()(0,0);
        LOG( INFO ) << "value_lambda_integrate_fem : "<<value_lambda_integrate_fem ;
        LOG( INFO ) << "value_lambda_integrate_rb : "<<value_lambda_integrate_rb ;
        LOG( INFO ) << "true value : "<<true_integrate;
        BOOST_CHECK_SMALL( math::abs(value_lambda_integrate_rb-value_lambda_integrate_fem), 1e-14 );
        BOOST_CHECK_SMALL( math::abs(value_lambda_integrate_rb-true_integrate), 1e-14 );
        LOG( INFO ) << " rb grad from contextrb :\n"<<grad_rb;

    }

    space_ptrtype functionSpace() { return Xh; }

private :
    space_ptrtype Xh;


};


/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( rbspace )

BOOST_AUTO_TEST_CASE( test_1 )
{
    boost::shared_ptr<Model<2,1> > model ( new Model<2,1>() );
    model->run();
}

BOOST_AUTO_TEST_SUITE_END()
