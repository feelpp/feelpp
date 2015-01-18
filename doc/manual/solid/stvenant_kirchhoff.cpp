/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-04-14

  Copyright (C) 2008-2012 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2012-2014 Feel++ Consortium

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
#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stvenant_kirchhoffoptions( "StVenant Kirchhoff solid model options" );
    stvenant_kirchhoffoptions.add_options()
    ( "omega", Feel::po::value<double>()->default_value( 2 ), "frequency" )
    ( "lambda", Feel::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Stvenant_Kirchhoff problem" )

    ( "order", Feel::po::value<int>()->default_value( 2 ), "order of time discretisation" )
    ( "diff", Feel::po::value<double>()->default_value( 1 ), "diffusion parameter" )
    ( "penal", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter" )
    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )
    ( "bctype", Feel::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ;
    return stvenant_kirchhoffoptions.add( Feel::feel_options() );
}


namespace Feel
{
/**
   \page StVenantKirchhoff The StVenant_Kirchhoff Model
   \author Christophe Prud'homme

   This is the documentation for the StVenant-Kirchhoff model
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

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    typedef Lagrange<Order, Vectorial> basis_u_type;
    typedef Lagrange<Order, Vectorial> basis_v_type;
    typedef bases<basis_u_type> basis_type;


    typedef FunctionSpace<mesh_type, basis_type, value_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

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

    StVenantKirchhoff()
        :
        super(),
        M_lambda( doption(_name="lambda") ),
        M_Xh(),
        dt( doption(_name="bdf.time-step") ),
        ft( doption(_name="bdf.time-final") ),
        omega( doption(_name="omega") )
        {
            this->changeRepository( boost::format( "doc/manual/solid/%1%/%2%/P%3%/h_%4%/" )
                                    % this->about().appName()
                                    % entity_type::name()
                                    % Order
                                    % doption(_name="gmsh.hsize")
                );

            /**
             * Physical data
             */
            M_time_order = ioption(_name="order");
            E = 21*1e5;
            sigma = 0.28;
            mu = E/( 2*( 1+sigma ) );
            lambda = E*sigma/( ( 1+sigma )*( 1-2*sigma ) );
            density = 1;
            gravity = -density*0.05;

            LOG(INFO) << "[data] dt=" << dt << "\n";
            LOG(INFO) << "[data] ft=" << ft << "\n";

            mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                                _desc=domain( _name=( boost::format( "beam-%1%" ) % Dim ).str() ,
                                                              _shape="hypercube",
                                                              _usenames=true,
                                                              _xmin=0., _xmax=20,
                                                              _ymin=-1., _ymax=1.,
                                                              _h=doption("gmsh.hsize") ) );

            M_Xh = functionspace_ptrtype( functionspace_type::New( mesh ) );
            un2 = element_ptrtype( new element_type( M_Xh, "un2" ) );
            un1 = element_ptrtype( new element_type( M_Xh, "un1" ) );
            un = element_ptrtype( new element_type( M_Xh, "un" ) );

        }

    /**
     * run the convergence test
     */
    void run();


    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J );

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& u );

private:

    double M_lambda;

    functionspace_ptrtype M_Xh;
    element_ptrtype un2;
    element_ptrtype un1;
    element_ptrtype un;

    oplin_ptrtype M_oplin;

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
void
StVenantKirchhoff<Dim, Order>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    LOG(INFO) << "[updateResidual] start\n";
    value_type penalisation_bc = option(_name="penalbc").template as<value_type>();

    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "U" );
    element_type v( M_Xh, "V" );
    u = *X;

    auto g = constant( 0.0 );
    auto defv = sym( gradv( u ) );
    auto def=  sym( grad( u ) );
    auto Id = eye<Dim>();
    //std::cout << "u = " << u << "\n";

    auto eta = 0.1*Px()*( Px() -5 )*( Px()-2.5 )*sin( omega*M_PI*cst_ref( time )  );

    form1( _test=M_Xh, _vector=R ) =
        integrate( elements( mesh ),
                   .5*mu*( trace( ( gradv( u )*trans( gradv( u ) ) )*grad( v ) ) )+
                   .25*lambda*trace( gradv( u )*trans( gradv( u ) ) )*div( v ) -
                   trans( gravity*oneY() )*id( v ) );

    // force applied at the bottom
    form1( _test=M_Xh, _vector=R ) +=
        integrate( markedfaces( mesh, 2 ),
                   -trans( eta*oneY() )*id( v ) );
    form1( _test=M_Xh, _vector=R ) +=
        integrate( elements( mesh ),
                   -density*trans( 2*idv( *un )-idv( *un1 ) ) *id( v ) /( dt*dt )
            );

    M_oplin->apply( u, v );

    R->add( 1., v );
    LOG(INFO) << "residual norm 2 = " << R->l2Norm() << "\n";
    LOG(INFO) << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J )
{
    boost::timer ti;
    LOG(INFO) << "[updateJacobian] start\n";
    static bool is_init = false;
    value_type penalisation_bc = option(_name="penalbc").template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "U" );
    element_type v( M_Xh, "V" );
    u = *X;

    if ( !J ) J= backend()->newMatrix( M_Xh, M_Xh );

    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) =
        integrate( elements( mesh ),
                   .5*mu*( trace( ( gradv( u )*trans( gradt( u ) ) )*grad( v ) ) )+
                   .25*lambda*trace( gradv( u )*trans( gradt( u ) ) )*div( v ) );
    J->addMatrix( 1.0, M_oplin->mat() );
    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J )
{
}

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::run()
{
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type U( M_Xh, "U" );
    element_type u( M_Xh, "U" );
    element_type v( M_Xh, "V" );


    M_bdf = bdf( _space=M_Xh );


    value_type penalisation = option(_name="penal").template as<value_type>();
    value_type penalisation_bc = option(_name="penalbc").template as<value_type>();
    int bctype = ioption(_name="bctype");
    value_type order = ioption(_name="order");


    LOG(INFO) << "lambda = " << lambda << "\n"
              << "mu     = " << mu << "\n"
              << "gravity= " << gravity << "\n";

    M_oplin = opLinear( _domainSpace=M_Xh, _imageSpace=M_Xh );
    auto deft = sym( gradt( u ) );
    auto def = sym( grad( v ) );
    auto Id = eye<Dim>();
    *M_oplin =
        integrate( elements( mesh ),
                   //density*trans(idt(uu))*id(v)*M_bdf->derivateCoefficient( M_time_order, dt ) +
                   density*trans( idt( u ) )*id( v )/( dt*dt )+
                   lambda*divt( u )*div( v )  +
                   2*mu*trace( trans( deft )*def )
            );

    *M_oplin +=
        integrate( markedfaces( mesh,1 ),
                   - trans( ( 2*mu*deft+lambda*trace( deft )*Id )*N() )*id( v )
                   - trans( ( 2*mu*def+lambda*trace( def )*Id )*N() )*idt( u )
                   + penalisation_bc*trans( idt( u ) )*id( v )/hFace() );

    *M_oplin +=
        integrate( markedfaces( mesh,3 ),
                   - trans( ( 2*mu*deft+lambda*trace( deft )*Id )*N() )*id( v )
                   - trans( ( 2*mu*def+lambda*trace( def )*Id )*N() )*idt( u )
                   + penalisation_bc*trans( idt( u ) )*id( v )/hFace() );

    backend()->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    backend()->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    U.zero();

    un->zero();
    un1->zero();

    M_bdf->initialize( U );


    boost::timer ttotal;
    int iterations = 0;

    for ( time = dt, iterations = 0; time < ft; time +=dt, ++iterations )
    {
        boost::timer ti;
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "time: " << time << "s, iteration: " << iterations << "\n";

        backend()->nlSolve( _solution=U );

        exportResults( time, U );

        *un1 = *un;
        *un = U;

        M_bdf->shiftRight( U );

        LOG(INFO) << "time spent in iteration :  " << ti.elapsed() << "s\n";
    }

    LOG(INFO) << "total time spent :  " << ttotal.elapsed() << "s\n";
    LOG(INFO) << "total number of iterations :  " << iterations << "\n";


} // StVenantKirchhoff::run



template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::exportResults( double time, element_type& U )
{

    LOG(INFO) << "exportResults starts\n";

    exporter->step( time )->setMesh( U.functionSpace()->mesh() );

    exporter->step( time )->add( "displ", U );

    exporter->save();
} // StVenantKirchhoff::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="stvenantkirchhoff",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );


    typedef Feel::StVenantKirchhoff<2, 3> solid_type;

    /* define and run application */
    solid_type solid;

    solid.run();
}
