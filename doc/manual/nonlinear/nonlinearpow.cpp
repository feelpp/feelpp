/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

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
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

#include <feel/feelvf/vf.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description nonlinearpowoptions("Nonlinearpow problem options");
    nonlinearpowoptions.add_options()
        ("lambda", Feel::po::value<int>()->default_value( 2 ), "power of u")

        ("penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")
        ("hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")

        ("export", "export results(ensight, data file(1D)")
        ("export-mesh-only", "export mesh only in ensight format")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return nonlinearpowoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "nonlinearpow" ,
                           "nonlinearpow" ,
                           "0.1",
                           "nD(n=1,2,3) NonLinearUL problem",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
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
        public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 4*Order;

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
    void exportResults( element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;
    int M_lambda;

    functionspace_ptrtype M_Xh;
    oplin_ptrtype M_oplin;
    export_ptrtype exporter;
}; // NonLinearPow

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
NonLinearPow<Dim,Order,Entity>::NonLinearPow( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    M_lambda( this->vm()["lambda"].template as<int>() ),
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
                            % this->vm()["lambda"].template as<int>()
                            );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=domain( _name= (boost::format( "%1%-%2%-%3%" ) % "hypercube" % Dim % 1).str() ,
                                                      _shape="hypercube",
                                                      _h=meshSize ),
                                        _partitions=this->comm().size()  );
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

    u = *X;
    auto g = constant(0.0);

    form1( _test=M_Xh, _vector=R ) =
        integrate( elements( mesh ),  + gradv(u)*trans(grad(v)) + pow(idv(u),M_lambda)*id(v) +id(v), _Q<2*Order>() ) +
        integrate( boundaryfaces(mesh),
                   ( - trans(id(v))*(gradv(u)*N())
                     - trans(idv(u))*(grad(v)*N())
                     + penalisation_bc*trans(idv(u))*id(v)/hFace())-
                   g*( - grad(v)*N() + penalisation_bc*id(v)/hFace() )
                   );
    R->close();
    Log() << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    Log() << "[updateJacobian] start\n";
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    u = *X;
    if ( !J ) J = M_backend->newMatrix( M_Xh, M_Xh );
    form2(_test=M_Xh, _trial=M_Xh, _matrix=J )  =
        integrate( elements( mesh ),  M_lambda*pow(idv(u),M_lambda-1)*idt(u)*id(v), _Q<2*Order>() );
    J->addMatrix( 1.0, M_oplin->mat() );
    Log() << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::run()
{
    using namespace Feel::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    *M_oplin =
        integrate( elements( mesh ), gradt(u)*trans(grad(v)) ) +
        integrate( boundaryfaces(mesh),
                   ( - trans(id(v))*(gradt(u)*N())
                     - trans(idt(u))*(grad(v)*N())
                     + penalisation_bc*trans(idt(u))*id(v)/hFace()) );
    M_oplin->close();

    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    u = vf::project( M_Xh, elements(mesh), constant(0.) );

    M_backend->nlSolve( _solution=u );

    exportResults( u );

} // NonLinearPow::run


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
NonLinearPow<Dim, Order, Entity>::exportResults( element_type& U )
{
    Log() << "exportResults starts\n";
    exporter->step(0)->setMesh( U.functionSpace()->mesh() );
    if ( !this->vm().count( "export-mesh-only" ) )
        {
            exporter->step(0)->addRegions();
            exporter->step(0)->add( "u", U );
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





