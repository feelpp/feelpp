/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-09

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
   \file nonlinearpow.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-04-14
 */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/applicationxml.hpp>
#include <feel/feelcore/xmlparser.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelvf/vf.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description nonlinearpowoptions( "Nonlinearpow problem options" );
    nonlinearpowoptions.add_options()
    ( "lambda", Feel::po::value<double>()->default_value( 2.0 ), "power of u" )

    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )

    ( "export", "export results(ensight, data file(1D)" )
    ( "export-mesh-only", "export mesh only in ensight format" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return nonlinearpowoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "nonlinearpow" ,
                           "nonlinearpow" ,
                           "0.2",
                           "nD(n=1,2,3) NonLinearUL problem",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{
using namespace Feel::vf;
/**
 * Nonlinearpow Problem
 *
 * solve \f$ -\Delta u + u^\lambda = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
template<int Dim,
         int Order = 1,
         template<uint16_type,uint16_type,uint16_type> class Entity = Simplex>
class NonLinearPow
    :
public ApplicationXML
{
    typedef ApplicationXML super;
public:

    // -- TYPEDEFS --
    typedef NonLinearPow<Dim,Order, Entity> self_type;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;

    typedef OperatorLinear<functionspace_type,functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    NonLinearPow( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();


    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );

private:

    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u , element_type& ue );

private:

    backend_ptrtype M_backend;

    double meshSize;
    int M_lambda;
    value_type M_penalisation_bc;
    functionspace_ptrtype M_Xh;
    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;
    funlin_ptrtype M_residual;

    export_ptrtype exporter;
}; // NonLinearPow

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
NonLinearPow<Dim,Order,Entity>::NonLinearPow( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( doption("hsize") ),
    M_lambda( doption("lambda") ),
    M_penalisation_bc( this->vm()["penalbc"].template as<value_type>() ),
    M_Xh(),
    exporter()
{
    Parameter h( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.04:0.08:0.2" );
    this->
    addParameter( Parameter( _name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str() ) )
    .addParameter( Parameter( _name="order",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Order  ).c_str() ) )
    .addParameter( Parameter( _name="lambda",_type=DISC_ATTR,_latex="\\lambda",_values="1,2,3" ) )
    .addParameter( h );

    std::vector<Parameter> depend;
    std::vector<std::string> funcs;
    depend.push_back( h );
    std::ostringstream oss;
    oss << "h**" << boost::lexical_cast<std::string>( Order + 1  ) ;
    funcs.push_back( oss.str() );
    oss.str( "" );
    std::vector<std::string> funcs2;
    oss << "h**" << boost::lexical_cast<std::string>( Order ) ;
    funcs2.push_back( oss.str() );

    this->
    addOutput( Output( _name="norm_L2",_latex="\\left\\| u \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs ) );
    //.addOutput( Output(_name="norm_H1",_latex="\\left\\| u \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs2) );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="square",
                                                _shape="hypercube",
                                                _dim=Dim,
                                                _h=meshSize ) );
    M_Xh = functionspace_ptrtype( functionspace_type::New( mesh ) );

    exporter = export_ptrtype( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );
}


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    LOG(INFO) << "[updateResidual] start\n";

    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    AUTO( g, ( Px()*Px()+Py()*Py() ) );
    u = *X;

    //*M_residual =
    form1( M_Xh, _vector=R ) =
        //integrate( elements( mesh ), _Q<4*Order>(), -(+ gradv(u)*trans(grad(v)) + pow(idv(u),M_lambda)*id(v) - (-4+pow(g,M_lambda))*id(v) ) ) +
        integrate( elements( mesh ), _Q<4*Order>(), ( + gradv( u )*trans( grad( v ) ) + idv( u )*id( v ) - ( constant( -4. )+g )*id( v ) ) )+
        integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                   ( - id( v )*( gradv( u )*N() )
                     //- idv(u)*(grad(v)*N())
                     + M_penalisation_bc*( idv( u )-g )*id( v )/hFace() )
                 );

    //M_residual->close();
    R->close();

    //if ( M_jac->matPtr()->closed() )
    //*M_jac += on( boundaryfaces(mesh), u, M_residual->container(), g, ON_PENALISATION );
    //*R = M_residual->container();

    u = *X;
    v = vf::project( M_Xh, elements( mesh ), g );
    double i1 = integrate( elements( mesh ), _Q<4*Order>(), ( + gradv( u )*trans( gradv( v ) ) + idv( u )*idv( v ) - ( constant( -4. )+g )*idv( v ) ) ).evaluate()( 0, 0 );
    double i2 = integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                           - idv( v )*( gradv( u )*N() )
                           //- idv(u)*(grad(v)*N())
                         ).evaluate()( 0 , 0 );
    double i3 = integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                           + M_penalisation_bc*( idv( u )-g )*idv( v )/hFace() ).evaluate()( 0 , 0 );

    double i4 = integrate( elements( mesh ), _Q<4*Order>(), ( + gradv( u )*trans( gradv( v ) ) + idv( u )*idv( v ) - ( constant( -4. )+g )*idv( v ) ) ).evaluate()( 0,0 )+
                integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                           - idv( v )*( gradv( u )*N() )
                           //- idv(u)*(grad(v)*N())
                           + M_penalisation_bc*( idv( u )-g )*idv( v )/hFace() ).evaluate()( 0 , 0 );
    vector_ptrtype V( M_backend->newVector( u.functionSpace() ) );
    *V = v;
    std::cout << "Residual  1 = " << i1 << "\n"
              << "Residual  2 = " << i2 << "\n"
              << "Residual  3 = " << i3 << "\n"
              << "sum Residuals   = " << i4 << "\n"
              << "Residual = " << M_backend->dot( V, R ) << "\n";

    vector_ptrtype R1( M_backend->newVector( u.functionSpace() ) );
    form1( M_Xh, _vector=R1 ) =
        //integrate( elements( mesh ), _Q<4*Order>(), -(+ gradv(u)*trans(grad(v)) + pow(idv(u),M_lambda)*id(v) - (-4+pow(g,M_lambda))*id(v) ) ) +
        integrate( elements( mesh ), _Q<4*Order>(), ( + gradv( u )*trans( grad( v ) ) + idv( u )*id( v ) - ( constant( -4. )+g )*id( v ) ) );
    R1->close();

    vector_ptrtype R2( M_backend->newVector( u.functionSpace() ) );
    form1( M_Xh, _vector=R2 ) =
        integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                   - id( v )*( gradv( u )*N() )
                   //- idv(u)*(grad(v)*N())
                   //+ M_penalisation_bc*print("tt",(idv(u)-g)*id(v)/hFace())
                 );
    R2->close();
    std::cout<< "Residual r1 = " << M_backend->dot( V, R1 ) << "\n";
    std::cout<< "Residual r2 = " << M_backend->dot( V, R2 ) << "\n";
    //M_residual->close();
    R->close();
    LOG(INFO) << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J )
{
    boost::timer ti;
    LOG(INFO) << "[updateJacobian] start\n";
    static bool is_init = false;
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    u = *X;
    AUTO( g, ( Px()*Px()+Py()*Py() ) );
    //u = project( M_Xh, elements(mesh), g );
#if 1

    //*M_jac = integrate( elements( mesh ), _Q<2*Order+2>(), gradt(u)*trans(grad(v))+M_lambda*pow(idv(u),M_lambda-1)*idt(u)*id(v) )
    *M_jac =
        integrate( elements( mesh ), _Q<2*Order+2>(), gradt( u )*trans( grad( v ) ) + idt( u )*id( v ) )+
        integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                   ( - id( v )*( gradt( u )*N() )
                     //- idt(u)*(grad(v)*N())
                     + M_penalisation_bc*idt( u )*id( v )/hFace() ) );

    M_jac->close();
#else

    if ( is_init == false )
    {
        *M_jac = integrate( elements( mesh ), _Q<2*Order+2>(), M_lambda*pow( idv( u ),M_lambda-1 )*idt( u )*id( v ) );
        is_init = true;
    }

    else
    {
        M_jac->matPtr()->zero();
        *M_jac += integrate( elements( mesh ), _Q<2*Order+2>(), M_lambda*pow( idv( u ),M_lambda-1 )*idt( u )*id( v ) );
    }

    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
#endif
    J = M_jac->matPtr();
    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::run()
{

    this->addParameterValue( Dim )
    .addParameterValue( Order )
    .addParameterValue( doption("lambda") )
    .addParameterValue( doption("hsize") );

    if ( this->preProcessing() == RUN_EXIT ) return;

    using namespace Feel::vf;

    boost::timer t1;

    mesh_ptrtype mesh = M_Xh->mesh();

    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    element_type ue( M_Xh, "ue" );



    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    *M_oplin =
        integrate( elements( mesh ), _Q<2*( Order-1 )>(), gradt( u )*trans( grad( v ) ) ) +
        integrate( boundaryfaces( mesh ), _Q<2*Order>(),
                   ( - id( v )*( gradt( u )*N() )
                     - idt( u )*( grad( v )*N() )
                     + M_penalisation_bc*idt( u )*id( v )/hFace() ) );
    M_oplin->close();

    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );

    AUTO( u_exact, ( Px()*Px()+Py()*Py() ) );

    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    //u = vf::project( M_Xh, elements(mesh), constant(0.) );
    u = vf::project( M_Xh, elements( mesh ), u_exact );
    ue = vf::project( M_Xh, elements( mesh ), u_exact );

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;

    vector_ptrtype R( M_backend->newVector( u.functionSpace() ) );
    this->updateResidual( U, R );

    sparse_matrix_ptrtype J;
    this->updateJacobian( U, J );


    solve( J, u, R );

    *U = u;
    this->updateResidual( U, R );
    std::cout << "R( u ) = " << M_backend->dot( U, R ) << "\n";

    LOG(INFO) << "solution computed in " << t1.elapsed() << "s\n";
    LOG(INFO) << "R( u ) = " << M_backend->dot( U, R ) << "\n";

    t1.restart();

    double L2error2 = integrate( elements( mesh ), _Q<2*Order>(),
                                 ( idv( u )-u_exact )*( idv( u )-u_exact ) ).evaluate()( 0, 0 );
    double L2error = math::sqrt( L2error2 );

    std::cout << "||error||_L2=" << L2error << "\n";
    LOG(INFO) << "||error||_L2=" << L2error << "\n";
    LOG(INFO) << "L2 norm computed in " << t1.elapsed() << "s\n";

    exportResults( u , ue );

    this->addOutputValue( L2error );
    this->postProcessing();

} // NonLinearPow::run

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    u = *U;


} // NonLinearPow::solve


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::exportResults( element_type& U , element_type& UE )
{
    LOG(INFO) << "exportResults starts\n";
    exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

    if ( !this->vm().count( "export-mesh-only" ) )
    {
        exporter->step( 0 )->add( "u", U );
        exporter->step( 0 )->add( "ue", UE );
    }

    exporter->save();
} // NonLinearPow::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    typedef Feel::NonLinearPow<nDim, nOrder> nonlinearpow_app_type;

    /* instantiate application */
    nonlinearpow_app_type nonlinearpow( argc, argv, makeAbout(), makeOptions() );

    /* run application */
    nonlinearpow.run();
}





