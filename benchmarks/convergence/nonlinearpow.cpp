/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-04-14
 */
#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/applicationxml.hpp>
#include <life/lifecore/xmlparser.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/exporter.hpp>

#include <life/lifevf/vf.hpp>

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description nonlinearpowoptions("Nonlinearpow problem options");
    nonlinearpowoptions.add_options()
        ("lambda", Life::po::value<double>()->default_value( 2.0 ), "power of u")

        ("penalbc", Life::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")

        ("export", "export results(ensight, data file(1D)")
        ("export-mesh-only", "export mesh only in ensight format")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return nonlinearpowoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "nonlinearpow" ,
                           "nonlinearpow" ,
                           "0.2",
                           "nD(n=1,2,3) NonLinearUL problem",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
using namespace Life::vf;
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
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J);

private:

    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u , element_type& ue);

private:

    backend_ptrtype M_backend;

    double meshSize;
    int M_lambda;

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
    meshSize( this->vm()["hsize"].template as<double>() ),
    M_lambda( this->vm()["lambda"].template as<double>() ),
    M_Xh(),
    exporter()
{
    Parameter h(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.04:0.08:0.2" );
    this->
        addParameter( Parameter(_name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str()) )
        .addParameter( Parameter(_name="order",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Order  ).c_str()) )
        .addParameter( Parameter(_name="lambda",_type=DISC_ATTR,_latex="\\lambda",_values="1,2,3") )
        .addParameter( h );

    std::vector<Parameter> depend;
    std::vector<std::string> funcs;
    depend.push_back(h);
    std::ostringstream oss;
    oss << "h**" << boost::lexical_cast<std::string>( Order + 1  ) ;
    funcs.push_back(oss.str());
    oss.str("");
    std::vector<std::string> funcs2;
    oss << "h**" << boost::lexical_cast<std::string>( Order ) ;
    funcs2.push_back(oss.str());

    this->
        addOutput( Output(_name="norm_L2",_latex="\\left\\| u \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs) );
        //.addOutput( Output(_name="norm_H1",_latex="\\left\\| u \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs2) );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain(_name="square",
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
    Log() << "[updateResidual] start\n";
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    AUTO( g, (Px()*Px()+Py()*Py()) );
    u = *X;
    //u = project( M_Xh, elements(mesh), g );
    *M_residual =
        //integrate( elements( mesh ), _Q<4*Order>(), -(+ gradv(u)*trans(grad(v)) + pow(idv(u),M_lambda)*id(v) - (-4+pow(g,M_lambda))*id(v) ) ) +
        integrate( elements( mesh ), _Q<4*Order>(), -(+ gradv(u)*trans(grad(v)) + idv(u)*id(v) - (-4+g)*id(v)) ) +
        integrate( boundaryfaces(mesh), _Q<2*Order>(),
                   -( - id(v)*(gradv(u)*N())
                      //- idv(u)*(grad(v)*N())
                     + penalisation_bc*idv(u)*id(v)/hFace())+
                   g*( - grad(v)*N() + penalisation_bc*id(v)/hFace() )
                   );

    M_residual->close();
    *R = M_residual->container();

    Log() << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    Log() << "[updateJacobian] start\n";
    static bool is_init = false;
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    u = *X;
    AUTO( g, (Px()*Px()+Py()*Py()) );
    //u = project( M_Xh, elements(mesh), g );
#if 1

    //*M_jac = integrate( elements( mesh ), _Q<2*Order+2>(), gradt(u)*trans(grad(v))+M_lambda*pow(idv(u),M_lambda-1)*idt(u)*id(v) )
    *M_jac = integrate( elements( mesh ), _Q<2*Order+2>(), gradt(u)*trans(grad(v)) + idt(u)*id(v) )
        + integrate( boundaryfaces(mesh), _Q<2*Order>(),
                     ( - id(v)*(gradt(u)*N())
                       - idt(u)*(grad(v)*N())
                       + penalisation_bc*idt(u)*id(v)/hFace() ) );
    M_jac->close();
#else
    if ( is_init == false )
        {
            *M_jac = integrate( elements( mesh ), _Q<2*Order+2>(), M_lambda*pow(idv(u),M_lambda-1)*idt(u)*id(v) );
            is_init = true;
        }
    else
        {
            M_jac->matPtr()->zero();
            *M_jac += integrate( elements( mesh ), _Q<2*Order+2>(), M_lambda*pow(idv(u),M_lambda-1)*idt(u)*id(v) );
        }
    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
#endif
    J = M_jac->matPtr();
    Log() << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::run()
{

    this->addParameterValue( Dim )
        .addParameterValue( Order )
        .addParameterValue( this->vm()["lambda"].template as<double>() )
        .addParameterValue( this->vm()["hsize"].template as<double>() );

    if (this->preProcessing() == RUN_EXIT) return;

    using namespace Life::vf;

    boost::timer t1;

    mesh_ptrtype mesh = M_Xh->mesh();

    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    element_type ue( M_Xh, "ue" );


    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    *M_oplin =
        integrate( elements( mesh ), _Q<2*(Order-1)>(), gradt(u)*trans(grad(v)) ) +
        integrate( boundaryfaces(mesh), _Q<2*Order>(),
                   ( - id(v)*(gradt(u)*N())
                     - idt(u)*(grad(v)*N())
                     + penalisation_bc*idt(u)*id(v)/hFace() ) );
    M_oplin->close();

    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );

    AUTO( u_exact, (Px()*Px()+Py()*Py()) );

    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    //u = project( M_Xh, elements(mesh), constant(0.) );
    u = project( M_Xh, elements(mesh), u_exact );
    ue = project( M_Xh, elements(mesh), u_exact );

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

    Log() << "solution computed in " << t1.elapsed() << "s\n";
    Log() << "R( u ) = " << M_backend->dot( U, R ) << "\n";

    t1.restart();

    double L2error2 = integrate( elements(mesh), _Q<2*Order>(),
                                 pow(idv(u)-u_exact,2) ).evaluate()( 0, 0 );
    double L2error = math::sqrt( L2error2 );

    std::cout << "||error||_L2=" << L2error << "\n";
    Log() << "||error||_L2=" << L2error << "\n";
    Log() << "L2 norm computed in " << t1.elapsed() << "s\n";

    exportResults( u , ue);

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
NonLinearPow<Dim, Order, Entity>::exportResults( element_type& U , element_type& UE)
{
    Log() << "exportResults starts\n";
    exporter->step(0)->setMesh( U.functionSpace()->mesh() );
    if ( !this->vm().count( "export-mesh-only" ) )
        {
            exporter->step(0)->add( "u", U );
            exporter->step(0)->add( "ue", UE );
        }
    exporter->save();
} // NonLinearPow::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    typedef Life::NonLinearPow<nDim, nOrder> nonlinearpow_app_type;

    /* instantiate application */
    nonlinearpow_app_type nonlinearpow( argc, argv, makeAbout(), makeOptions() );

    /* run application */
    nonlinearpow.run();
}





