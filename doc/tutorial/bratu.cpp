/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-09

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
   \file bratu.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-04-14
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>

#include <life/lifevf/vf.hpp>

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description bratuoptions("Bratu problem options");
    bratuoptions.add_options()
        ("lambda", Life::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Bratu problem")

        ("penalbc", Life::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")

        ("export", "export results(ensight, data file(1D)")
        ("export-mesh-only", "export mesh only in ensight format")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return bratuoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "bratu" ,
                           "bratu" ,
                           "0.1",
                           "nD(n=1,2,3) Bratu problem",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
using namespace Life::vf;
/**
 * Bratu Problem
 *
 * solve \f$ -\Delta u + \lambda \exp(u) = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
template<int Dim,
         int Order = 1,
         template<uint16_type,uint16_type,uint16_type> class Entity = Simplex>
class Bratu
    :
        public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 4*Order;

    typedef Bratu<Dim,Order, Entity> self_type;

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
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /* number of dofs per element */
    static const uint16_type nLocalDof = boost::remove_reference<typename fusion::result_of::at<basis_type,mpl::int_<0> >::type>::type::nLocalDof;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;

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
    Bratu( int argc, char** argv, AboutData const& ad, po::options_description const& od );

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

private:

    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double M_lambda;

    functionspace_ptrtype M_Xh;
    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;
    funlin_ptrtype M_residual;

    export_ptrtype exporter;
}; // Bratu

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
Bratu<Dim,Order,Entity>::Bratu( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    M_lambda( this->vm()["lambda"].template as<double>() ),
    M_Xh(),
    exporter()
{

    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    this->changeRepository( boost::format( "doc/tutorial/%1%/%2%/P%3%/h_%4%/lambda_%5%" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            % this->vm()["lambda"].template as<double>()
                            );

    mesh_ptrtype mesh = createMesh( meshSize );
    M_Xh = functionspace_ptrtype( functionspace_type::New( mesh ) );

    exporter = export_ptrtype( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );
}
template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Bratu<Dim,Order,Entity>::mesh_ptrtype
Bratu<Dim,Order,Entity>::createMesh( double meshSize, double ymin, double ymax )
{
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
    td.setCharacteristicLength( meshSize );
    if ( Dim >=2 )
        td.setY( std::make_pair( ymin, ymax ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // Bratu::createMesh


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Bratu<Dim, Order, Entity>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    Log() << "[updateResidual] start\n";
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    u = *X;
    AUTO( g, constant(0.0) );

    *M_residual =
        integrate( elements( mesh ), _Q<2*Order>(), + gradv(u)*trans(grad(v)) + M_lambda*exp(idv(u))*id(v) ) +
        integrate( boundaryfaces(mesh), _Q<2*Order>(),
                   ( - trans(id(v))*(gradv(u)*N())
                     - trans(idv(u))*(grad(v)*N())
                     + penalisation_bc*trans(idv(u))*id(v)/hFace())-
                   g*( - grad(v)*N() + penalisation_bc*id(v)/hFace() )
                   );

    M_residual->close();
    *R = M_residual->container();
    Log() << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Bratu<Dim, Order, Entity>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    Log() << "[updateJacobian] start\n";
    static bool is_init = false;
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    u = *X;
    if ( is_init == false )
        {
            *M_jac = integrate( elements( mesh ), _Q<2*Order>(), M_lambda*(exp(idv(u)))*idt(u)*id(v) );
            is_init = true;
        }
    else
        {
            M_jac->matPtr()->zero();
            *M_jac += integrate( elements( mesh ), _Q<2*Order>(), M_lambda*(exp(idv(u)))*idt(u)*id(v) );
        }
    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();
    Log() << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Bratu<Dim, Order, Entity>::run()
{
    using namespace Life::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    *M_oplin = integrate( elements( mesh ), _Q<2*Order>(), gradt(u)*trans(grad(v)) );
    *M_oplin = integrate( boundaryfaces(mesh), _Q<2*Order>(),
                          ( - trans(id(v))*(gradt(u)*N())
                            - trans(idt(u))*(grad(v)*N())
                            + penalisation_bc*trans(idt(u))*id(v)/hFace()) );
    M_oplin->close();

    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );



    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    u = project( M_Xh, elements(mesh), constant(0.) );

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

    exportResults( u );

} // Bratu::run

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Bratu<Dim, Order, Entity>::solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    u = *U;


} // Bratu::solve


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Bratu<Dim, Order, Entity>::exportResults( element_type& U )
{
    Log() << "exportResults starts\n";
    exporter->step(0)->setMesh( U.functionSpace()->mesh() );
    if ( !this->vm().count( "export-mesh-only" ) )
        {
            exporter->step(0)->add( "u", U );
        }
    exporter->save();
} // Bratu::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    typedef Life::Bratu<nDim, nOrder> bratu_app_type;

    /* instantiate application */
    bratu_app_type bratu( argc, argv, makeAbout(), makeOptions() );

    /* run application */
    bratu.run();
}





