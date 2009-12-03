
/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-04-14

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
   \file stvenant_kirchhoff.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-04-14
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifediscr/bdf.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>




inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description stvenant_kirchhoffoptions("StVenant Kirchhoff solid model options");
    stvenant_kirchhoffoptions.add_options()
        ("dt", Life::po::value<double>()->default_value( 1 ), "time step value")
        ("ft", Life::po::value<double>()->default_value( 1 ), "final time value")
        ("omega", Life::po::value<double>()->default_value( 2 ), "frequency")
        ("lambda", Life::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Stvenant_Kirchhoff problem")

        ("order", Life::po::value<int>()->default_value( 2 ), "order of time discretisation")
        ("diff", Life::po::value<double>()->default_value( 1 ), "diffusion parameter")
        ("penal", Life::po::value<double>()->default_value( 10 ), "penalisation parameter")
        ("penalbc", Life::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("bctype", Life::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("export", "export results(ensight, data file(1D)")
        ("export-mesh-only", "export mesh only in ensight format")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return stvenant_kirchhoffoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "stvenant_kirchhoff" ,
                           "stvenant_kirchhoff" ,
                           "0.1",
                           "nD(n=1,2,3) Stvenant_Kirchhoff model",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
template<typename A, uint16_type i>
class Tagged : public A
{
public:
    static const uint16_type TAG = i;

};
#define MIXED 0

using namespace Life::vf;
/**
 * StVenant_Kirchhoff Model
 *
 */
template<int Dim, int Order>
class StVenantKirchhoff
    :
    public Application
{
    typedef Application super;
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
    typedef Entity<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    typedef Tagged<Lagrange<Order, Vectorial>, 0> basis_u_type;
    typedef Tagged<Lagrange<Order, Vectorial>, 1> basis_v_type;
#if MIXED
    typedef mpl::vector<basis_u_type, basis_v_type> basis_type;
#else
    typedef mpl::vector<basis_u_type> basis_type;
#endif

    typedef FunctionSpace<mesh_type, basis_type, value_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
#if MIXED
    typedef typename element_type::template sub_element<0>::type element_u_type;
    typedef typename element_type::template sub_element<1>::type element_v_type;
#endif

    typedef OperatorLinear<functionspace_type,functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    /* time */
    typedef Bdf<functionspace_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    StVenantKirchhoff( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        M_lambda( this->vm()["lambda"].template as<double>() ),
        M_Xh(),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        dt( this->vm()["dt"].template as<double>() ),
        ft( this->vm()["ft"].template as<double>() ),
        omega( this->vm()["omega"].template as<double>() )
    {
        if ( this->vm().count( "help" ) )
            {
                std::cout << this->optionsDescription() << "\n";
                return;
            }




        this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                                % this->about().appName()
                                % entity_type::name()
                                % Order
                                % this->vm()["hsize"].template as<double>()
                            );

        /**
         * Physical data
         */
        M_time_order = this->vm()["order"].template as<int>();
        E = 21*1e5;
        sigma = 0.28;
        mu = E/(2*(1+sigma));
        lambda = E*sigma/((1+sigma)*(1-2*sigma));
        density = 1;
        gravity = -density*0.05;

        Log() << "[data] dt=" << dt << "\n";
        Log() << "[data] ft=" << ft << "\n";

        mesh_ptrtype mesh = createMesh( meshSize );

        M_Xh = functionspace_ptrtype( functionspace_type::New( mesh ) );
        un2 = element_ptrtype( new element_type( M_Xh, "un2" ) );
        un1 = element_ptrtype( new element_type( M_Xh, "un1" ) );
        un = element_ptrtype( new element_type( M_Xh, "un" ) );

    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize, double ymin = 0, double ymax = 1 );

    /**
     * run the convergence test
     */
    void run();


    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J);

private:



    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double M_lambda;

    functionspace_ptrtype M_Xh;
    element_ptrtype un2;
    element_ptrtype un1;
    element_ptrtype un;

    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;
    funlin_ptrtype M_residual;

    export_ptrtype exporter;

    double E;
    double sigma;
    double mu;
    double lambda;
    double density;
    double gravity;

    double dt;
    double ft;
    double omega;


    double time;
    double M_time_order;
    bdf_ptrtype M_bdf;
}; // StVenantKirchhoff

template<int Dim, int Order>
typename StVenantKirchhoff<Dim,Order>::mesh_ptrtype
StVenantKirchhoff<Dim,Order>::createMesh( double meshSize, double ymin, double ymax )
{
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    GmshTensorizedDomain<Dim,1,Dim,Entity> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( 0, 20 ) );
    td.setY( std::make_pair( -1, 1 ) );

    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // StVenantKirchhoff::createMesh


template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    Log() << "[updateResidual] start\n";
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();

    mesh_ptrtype mesh = M_Xh->mesh();
#if MIXED
    element_type U( M_Xh, "U" );
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();
    U = *X;
#else
    element_type u( M_Xh, "U" );
    element_type v( M_Xh, "V" );
    u = *X;
#endif

    AUTO( g, constant(0.0) );
    AUTO( defv, 0.5*( gradv(u)+trans(gradv(u)) ) );
    AUTO( def, 0.5*( grad(v)+trans(grad(v)) ) );
    AUTO( Id, (mat<Dim,Dim>( cst(1), cst(0), cst(0), cst(1.) )) );
    //std::cout << "u = " << u << "\n";

    AUTO( eta, 0.1*Px()*( Px() -5 )*(Px()-2.5)*sin( omega*M_PI*cst_ref(time)  ) );

    *M_residual =
        integrate( elements( mesh ), _Q<imOrder>(),
                   .5*mu*(trace( (gradv(u)*trans(gradv(u)))*grad(v) ) )+
                   .25*lambda*trace(gradv(u)*trans(gradv(u)))*div(v) -
                   trans(gravity*oneY())*id(v) );

#if 1
    // force applied at the bottom
    *M_residual +=
        integrate( markedfaces( mesh, 2 ), _Q<imOrder>(),
                   -trans(eta*oneY())*id(v) );
#endif

#if MIXED
    *M_residual +=
        integrate( elements( mesh ), _Q<imOrder>(),
                   - density*trans(idv( M_bdf->derivate( M_time_order, dt ).template element<1>() ) ) *id(v)
                   //-density*trans(2*idv(un->template element<0>())-idv(un1->template element<0>())) *id(v) /(dt*dt)
                   );
    *M_residual +=
        integrate( elements( mesh ), _Q<imOrder>(),
                   + trans(idv( u ))*id(vv)*M_bdf->derivateCoefficient( M_time_order, dt )
                   - trans(idv( M_bdf->derivate( M_time_order, dt ).template element<0>() ) )*id(vv)
                   );
    FsFunctionalLinear<functionspace_type> flin( M_Xh, M_backend );
    M_oplin->apply( U, flin );
#else
    *M_residual +=
        integrate( elements( mesh ), _Q<imOrder>(),
                   -density*trans(2*idv(*un)-idv(*un1)) *id(v) /(dt*dt)
                   );

    FsFunctionalLinear<functionspace_type> flin( M_Xh, M_backend );
    M_oplin->apply( u, flin );
#endif





    M_residual->add( flin );
    M_residual->close();
    *R = M_residual->container();
    Log() << "[updateResidual] done in " << ti.elapsed() << "s\n";
                   }
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    Log() << "[updateJacobian] start\n";
    static bool is_init = false;
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
#if MIXED
    element_type U( M_Xh, "U" );
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();

    U = *X;
#else
    element_type u( M_Xh, "U" );
    element_type v( M_Xh, "V" );
    u = *X;
#endif
    if ( is_init == false )
        {
            *M_jac = integrate( elements( mesh ), _Q<imOrder>(),
                                .5*mu*(trace( (gradv(u)*trans(gradt(u)))*grad(v) ) )+
                                .25*lambda*trace(gradv(u)*trans(gradt(u)))*div(v)
                                );

            is_init = true;
        }
    else
        {
            M_jac->matPtr()->zero();
            *M_jac += integrate( elements( mesh ), _Q<imOrder>(),
                                 .5*mu*(trace( (gradv(u)*trans(gradt(u)))*grad(v) ) )+
                                 .25*lambda*trace(gradv(u)*trans(gradt(u)))*div(v) );
        }
    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();
    Log() << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J)
{
}

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::run()
{
    using namespace Life::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type U( M_Xh, "U" );
#if MIXED
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();
#else
    element_type u( M_Xh, "U" );
    element_type v( M_Xh, "V" );
#endif


    M_bdf = bdf_ptrtype( new bdf_type( M_Xh ) );


    value_type penalisation = this->vm()["penal"].template as<value_type>();
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    int bctype = this->vm()["bctype"].template as<int>();
    value_type order = this->vm()["order"].template as<int>();


    Log() << "lambda = " << lambda << "\n"
          << "mu     = " << mu << "\n"
          << "gravity= " << gravity << "\n";

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    AUTO( deft, 0.5*( gradt(u)+trans(gradt(u)) ) );
    AUTO( def, 0.5*( grad(v)+trans(grad(v)) ) );
    AUTO( Id, (mat<Dim,Dim>( cst(1), cst(0), cst(0), cst(1.) )) );
    *M_oplin =
        integrate( elements(mesh), _Q<imOrder>(),
                   //density*trans(idt(uu))*id(v)*M_bdf->derivateCoefficient( M_time_order, dt ) +
                   density*trans(idt(u))*id(v)/(dt*dt)+
                   lambda*divt(u)*div(v)  +
                   2*mu*trace(trans(deft)*def)
#if MIXED
                   //
                   //+ trans(idt( u ))*id(vv)*M_bdf->derivateCoefficient( M_time_order, dt )
                   - trans(idt(uu))*id(vv)
#endif
                   );

    *M_oplin +=
        integrate( markedfaces(mesh,1), _Q<imOrder>(),
                   - trans((2*mu*deft+lambda*trace(deft)*Id )*N())*id(v)
                   - trans((2*mu*def+lambda*trace(def)*Id )*N())*idt(u)
                   + penalisation_bc*trans(idt(u))*id(v)/hFace() );

    *M_oplin +=
        integrate( markedfaces(mesh,3), _Q<imOrder>(),
                   - trans((2*mu*deft+lambda*trace(deft)*Id )*N())*id(v)
                   - trans((2*mu*def+lambda*trace(def)*Id )*N())*idt(u)
                   + penalisation_bc*trans(idt(u))*id(v)/hFace() );

    M_oplin->close();

    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );


    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    //u = project( M_Xh, elements(mesh), constant(0.)*one() );
    U.zero();

    un->zero();
    un1->zero();

    vector_ptrtype Un( M_backend->newVector( U.functionSpace() ) );

    vector_ptrtype R( M_backend->newVector( U.functionSpace() ) );
    sparse_matrix_ptrtype J;

    M_bdf->initialize( U );

    boost::timer ttotal;
    int iterations = 0;
    for( time = dt, iterations = 0; time < ft; time +=dt, ++iterations )
        {
            boost::timer ti;
            Log() << "============================================================\n";
            Log() << "time: " << time << "s, iteration: " << iterations << "\n";

            *Un = U;

            this->updateResidual( Un, R );
            this->updateJacobian( Un, J );

            solve( J, U, R );

            exportResults( time, U );

            *un1 = *un;
            *un = U;

            M_bdf->shiftRight( U );

            Log() << "time spent in iteration :  " << ti.elapsed() << "s\n";
        }

    Log() << "total time spent :  " << ttotal.elapsed() << "s\n";
    Log() << "total number of iterations :  " << iterations << "\n";


} // StVenantKirchhoff::run

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::solve( sparse_matrix_ptrtype& D,
                                      element_type& u,
                                      vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    u = *U;


} // StVenantKirchhoff::solve


template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::exportResults( double time, element_type& U)
{

    Log() << "exportResults starts\n";

    exporter->step(time)->setMesh( U.functionSpace()->mesh() );
    exporter->step(time)->add( "pid",
                   regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );

#if MIXED
    exporter->step(time)->add( "displ", U.template element<0>() );
    exporter->step(time)->add( "veloc", U.template element<1>() );
#else
    exporter->step(time)->add( "displ", U );
#endif

    exporter->save();
    exporter->step(time)->setState( STEP_ON_DISK );
} // StVenantKirchhoff::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 3;
    typedef Life::StVenantKirchhoff<nDim, nOrder> solid_type;

    /* define and run application */
    solid_type solid( argc, argv, makeAbout(), makeOptions() );

    solid.run();
}





