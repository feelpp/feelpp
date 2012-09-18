/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \file bratu.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-04-14
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>




inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Bratu problem options" );
    bratuoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 1 ), "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ), "final time value" )
    ( "lambda", Feel::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Bratu problem" )

    ( "order", Feel::po::value<int>()->default_value( 2 ), "order of time discretisation" )
    ( "diff", Feel::po::value<double>()->default_value( 1 ), "diffusion parameter" )
    ( "penal", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter" )
    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-mesh-only", "export mesh only in ensight format" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return bratuoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "bratu" ,
                           "bratu" ,
                           "0.1",
                           "nD(n=1,2,3) Bratu problem on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{
using namespace Feel::vf;
/**
 * Bratu Problem
 *
 * solve \f$ -\Delta u + \lambda \exp(u) = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
template<int Dim,
         int Order,
         typename Cont,
         template<uint16_type,uint16_type,uint16_type> class Entity,
         template<uint16_type> class FType>
class Bratu
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 4*Order;

    typedef Bratu<Dim,Order, Cont, Entity, FType> self_type;

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

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >,Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    template<typename Conti = Cont>
    struct space
    {
        /*basis*/
#if 0
        typedef typename mpl::if_<mpl::bool_<Conti::is_continuous>,
                mpl::identity<fusion::vector<Lagrange<Order, FType> > >,
                mpl::identity<fusion::vector<OrthonormalPolynomialSet<Order, FType> > > >::type::type basis_type;
#else
        typedef fusion::vector<Lagrange<Order, FType> > basis_type;

#endif

        /*space*/
        typedef FunctionSpace<mesh_type, basis_type, Conti, value_type> type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef typename type::element_type element_type;
        typedef typename element_type::template sub_element<0>::type element_0_type;
        typedef typename element_type::template sub_element<1>::type element_1_type;
    };
    typedef typename space<Cont>::type functionspace_type;
    typedef typename space<Cont>::ptrtype functionspace_ptrtype;
    typedef typename space<Cont>::element_type element_type;

    typedef OperatorLinear<functionspace_type,functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    Bratu( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        M_lambda( this->vm()["lambda"].template as<double>() ),
        M_Xh(),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
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

        mesh_ptrtype mesh = createMesh( meshSize );

        M_Xh = functionspace_ptrtype( space<Cont>::type::New( mesh ) );
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
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J );

private:



    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( double time,
                        f1_type& u,
                        f2_type& v,
                        f3_type& e );

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

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
typename Bratu<Dim,Order,Cont,Entity,FType>::mesh_ptrtype
Bratu<Dim,Order,Cont,Entity,FType>::createMesh( double meshSize, double ymin, double ymax )
{
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    GmshHypercubeDomain td( entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,entity_type::is_hypercube );
    td.setCharacteristicLength( meshSize );

    if ( Dim >=2 )
        td.setY( std::make_pair( ymin, ymax ) );

    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // Bratu::createMesh


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    Debug() << "[updateResidual] start\n";
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    im_type im;
    u = *X;
    AUTO( g, constant( 0.0 ) );
    //std::cout << "u = " << u << "\n";
    *M_residual =
        integrate( elements( mesh ), im, + gradv( u )*trans( grad( v ) ) + M_lambda*exp( idv( u ) )*id( v ) ) +
        integrate( boundaryfaces( mesh ), im,
                   //integrate( markedfaces(mesh,1), im,
                   ( - trans( id( v ) )*( gradv( u )*N() )
                     - trans( idv( u ) )*( grad( v )*N() )
                     + penalisation_bc*trans( idv( u ) )*id( v )/hFace() )-
                   g*( - grad( v )*N() + penalisation_bc*id( v )/hFace() )
                 );

    M_residual->close();
    *R = M_residual->container();
    Debug() << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J )
{
    boost::timer ti;
    Debug() << "[updateJacobian] start\n";
    static bool is_init = false;
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    u = *X;
    im_type im;

    if ( is_init == false )
    {
        *M_jac = integrate( elements( mesh ), im, M_lambda*( exp( idv( u ) ) )*idt( u )*id( v ) );
        is_init = true;
    }

    else
    {
        M_jac->matPtr()->zero();
        *M_jac += integrate( elements( mesh ), im, M_lambda*( exp( idv( u ) ) )*idt( u )*id( v ) );
    }

    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();
    Debug() << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J )
{
}

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::run()
{
    using namespace Feel::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    im_type im;

    value_type penalisation = this->vm()["penal"].template as<value_type>();
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    int bctype = this->vm()["bctype"].template as<int>();
    value_type dt = this->vm()["dt"].template as<value_type>();
    value_type ft = this->vm()["ft"].template as<value_type>();
    value_type order = this->vm()["order"].template as<int>();

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    *M_oplin =
        integrate( elements( mesh ), im, gradt( u )*trans( grad( v ) ) ) +
        integrate( boundaryfaces( mesh ), im,
                   ( - trans( id( v ) )*( gradt( u )*N() )
                     - trans( idt( u ) )*( grad( v )*N() )
                     + penalisation_bc*trans( idt( u ) )*id( v )/hFace() ) );
    M_oplin->close();

    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );



    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    u = vf::project( M_Xh, elements( mesh ), constant( 0. ) );

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

    exportResults( 0.0, u, u, u );

} // Bratu::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::solve( sparse_matrix_ptrtype& D,
        element_type& u,
        vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    u = *U;


} // Bratu::solve


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
template<typename f1_type, typename f2_type, typename f3_type>
void
Bratu<Dim, Order, Cont, Entity,FType>::exportResults( double time,
        f1_type& U,
        f2_type& V,
        f3_type& E )
{

    LOG(INFO) << "exportResults starts\n";

    exporter->step( time )->setMesh( U.functionSpace()->mesh() );

    //exporter->step(time)->setMesh( this->createMesh( meshSize/2, 0.5, 1 ) );
    //exporter->step(time)->setMesh( this->createMesh( meshSize/Order, 0, 1 ) );
    //exporter->step(time)->setMesh( this->createMesh( meshSize, 0, 1 ) );
    if ( !this->vm().count( "export-mesh-only" ) )
    {
        exporter->step( time )->add( "pid",
                                     regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );


        exporter->step( time )->add( "u", U );
        exporter->step( time )->add( "v", V );
        exporter->step( time )->add( "e", E );
    }

    exporter->save();


    if ( Dim == 1 )
    {
        std::ostringstream fname_u;
        fname_u << "u-" << Application::processId() << "." << boost::format( "%.2f" ) % time << ".dat";
        std::ofstream ofs3( fname_u.str().c_str() );
        typename mesh_type::element_iterator it = U.functionSpace()->mesh()->beginElementWithProcessId( Application::processId() );
        typename mesh_type::element_iterator en = U.functionSpace()->mesh()->endElementWithProcessId( Application::processId() );

        if ( !U.areGlobalValuesUpdated() )
            U.updateGlobalValues();

        for ( ; it!=en; ++it )
        {
            for ( size_type i = 0; i < space<Cont>::type::basis_type::nLocalDof; ++i )
            {
                size_type dof0 = boost::get<0>( U.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                ofs3 << std::setw( 5 ) << it->id() << " "
                     << std::setw( 5 ) << i << " "
                     << std::setw( 5 ) << dof0 << " "
                     << std::setw( 15 ) << U.globalValue( dof0 ) << " ";
                value_type a = it->point( 0 ).node()[0];
                value_type b = it->point( 1 ).node()[0];

                if ( i == 0 )
                    ofs3 << a;

                else if ( i == 1 )
                    ofs3 <<  b;

                else
                    ofs3 <<  a + ( i-1 )*( b-a )/( space<Continuous>::type::basis_type::nLocalDof-1 );

                ofs3 << "\n";

            }
        }

        ofs3.close();

        std::ostringstream fname_v;
        fname_v << "values-" << Application::processId() << "." << boost::format( "%.2f" ) % time << ".dat";
        std::ofstream ofs2( fname_v.str().c_str() );
        it = V.functionSpace()->mesh()->beginElementWithProcessId( Application::processId() );
        en = V.functionSpace()->mesh()->endElementWithProcessId(  Application::processId() );

        if ( !V.areGlobalValuesUpdated() ) V.updateGlobalValues();

        if ( !E.areGlobalValuesUpdated() ) E.updateGlobalValues();

        for ( ; it!=en; ++it )
        {
            for ( size_type i = 0; i < space<Continuous>::type::basis_type::nLocalDof; ++i )
            {
                size_type dof0 = boost::get<0>( V.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                ofs2 << std::setw( 5 ) << it->id() << " "
                     << std::setw( 5 ) << i << " "
                     << std::setw( 5 ) << dof0 << " "
                     << std::setw( 15 ) << V.globalValue( dof0 ) << " "
                     << std::setw( 15 ) << E.globalValue( dof0 ) << " ";
                value_type a = it->point( 0 ).node()[0];
                value_type b = it->point( 1 ).node()[0];

                if ( i == 0 )
                    ofs2 << a;

                else if ( i == 1 )
                    ofs2 <<  b;

                else
                    ofs2 <<  a + ( i-1 )*( b-a )/( space<Continuous>::type::basis_type::nLocalDof-1 );

                ofs2 << "\n";

            }
        }

    }


} // Bratu::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    typedef Continuous MyContinuity;
    //typedef Discontinuous MyContinuity;
    //typedef Feel::Bratu<nDim, nOrder, MyContinuity, Hypercube, Scalar> bratu_type;

    typedef Feel::Bratu<nDim, nOrder, MyContinuity, Simplex, Scalar> bratu_type;

    /* define and run application */
    bratu_type bratu( argc, argv, makeAbout(), makeOptions() );

    bratu.run();
}





