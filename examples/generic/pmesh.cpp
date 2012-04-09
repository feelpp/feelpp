/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-11-23

  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-11-23
 */
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/options.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
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
    Feel::po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
    ( "diff", Feel::po::value<double>()->default_value( 1 ), "diffusion parameter" )
    ( "penal", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter" )
    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )
    ( "f", Feel::po::value<double>()->default_value( 1 ), "forcing term" )
    ( "g", Feel::po::value<double>()->default_value( 0 ), "boundary term" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-mesh-only", "export mesh only in ensight format" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return laplacianoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "pmesh" ,
                           "pmesh" ,
                           "0.2",
                           "nD(n=1,2,3) Laplacian on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006, 2007, 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;

}


namespace Feel
{
using namespace vf;
template<int Dim>struct ExactSolution {};
template<>
struct ExactSolution<1>
{
    typedef __typeof__( sin( M_PI*Px() ) ) type;
    typedef __typeof__( M_PI*M_PI*sin( M_PI*Px() ) ) laplacian_type;
};

template<>
struct ExactSolution<2>
{
    typedef __typeof__( sin( M_PI*Px() )*cos( M_PI*Py() ) ) type;
    typedef __typeof__( 2*M_PI*M_PI*sin( M_PI*Px() )*cos( M_PI*Py() ) ) laplacian_type;
};

template<>
struct ExactSolution<3>
{
    typedef __typeof__( sin( M_PI*Px() )*cos( M_PI*Py() )*cos( M_PI*Pz() ) ) type;
    typedef __typeof__( 3*M_PI*M_PI*sin( M_PI*Px() )*cos( M_PI*Py() )*cos( M_PI*Pz() ) ) laplacian_type;
};
/**
 * Laplacian Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim,
         int Order,
         typename Cont,
         template<uint16_type,uint16_type,uint16_type> class Entity,
         template<uint16_type> class FType>
class Laplacian
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    template<typename Conti = Cont>
    struct space
    {
        /*basis*/
        typedef fusion::vector<Lagrange<Order, FType> > basis_type;
#if 0
        typedef typename mpl::if_<mpl::bool_<Conti::is_continuous>,
                mpl::identity<fusion::vector<Lagrange<Order, FType> > >,
                mpl::identity<fusion::vector<OrthonormalPolynomialSet<Order, FType> > > >::type::type basis_type;
#endif
        /*space*/
        typedef FunctionSpace<mesh_type, basis_type, Conti, value_type> type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef typename type::element_type element_type;
        typedef typename element_type::template sub_element<0>::type element_0_type;
        typedef typename element_type::template sub_element<1>::type element_1_type;
    };


    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers(),
        stats()
    {
        Log() << "[Laplacian] hsize = " << meshSize << "\n";
        Log() << "[Laplacian] export = " << this->vm().count( "export" ) << "\n";

    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize, double ymin = 0, double ymax = 1 );

    /**
     * alias for run()
     */
    void operator()()
    {
        run();
    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve the system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solve( Mat& D, Vec1& u, Vec2& F, bool is_sym );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( f1_type& u,
                        f2_type& v,
                        f3_type& e );

private:

    double meshSize;

    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // Laplacian

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
typename Laplacian<Dim,Order,Cont,Entity,FType>::mesh_ptr_type
Laplacian<Dim,Order,Cont,Entity,FType>::createMesh( double meshSize, double ymin, double ymax )
{
    timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );
    //mesh->setRenumber( false );

    GmshHypercubeDomain td( entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,entity_type::is_hypercube );
    td.setCharacteristicLength( meshSize );
    td.setY( std::make_pair( ymin, ymax ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    mesh->partition();
    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // Laplacian::createMesh


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Laplacian<Dim, Order, Cont, Entity, FType>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;


    this->changeRepository( boost::format( "%1%/%2%/P%3%/%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                          );
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( meshSize );
    //stats["nelt"] = mesh->elements().size();

    boost::shared_ptr<p0_space_type> P0h( new p0_space_type( mesh ) );
    typename p0_space_type::element_type pid( P0h, "pid" );
    pid = regionProcess( P0h );
    typename mesh_type::element_iterator it = mesh->beginElementWithProcessId( Application::processId() );
    typename mesh_type::element_iterator en = mesh->endElementWithProcessId( Application::processId() );
    Log() << "first local index = " << pid.firstLocalIndex() << "\n";
    Log() << "last local index = " << pid.lastLocalIndex() << "\n";

    for ( ; it != en; ++it )
    {
        size_type dof0 = boost::get<0>( P0h->dof()->localToGlobal( it->id(), 0 ) );

        Log() << "element " << it->id() << " pid = " << it->processId() << " dof = " << dof0 << " value= " << pid( dof0 ) << "\n";
        //FEELPP_ASSERT( pid( dof0
    }

    if ( this->vm().count( "export" ) )
    {
        exporter->step( 1. )->setMesh( mesh );
#if 1
        exporter->step( 1. )->add( "pid", pid );

#endif
        exporter->save();
    }


    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    Log() << "defining Xh\n";
    typename space<Cont>::ptrtype Xh = space<Cont>::type::New( mesh );
    Log() << "defining Xh done\n";
    //Xh->dof()->showMe();
    Log() << "defining u\n";
    typename space<Cont>::element_type u( Xh, "u" );
    Log() << "defining u done\n";
    typename space<Cont>::element_type v( Xh, "v" );
    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();
#if 0

    if ( this->vm().count( "export-mesh-only" ) )
        this->exportResults( u, u, u );

    /*
     * a quadrature rule for numerical integration
     */
    im_type im;

    value_type penalisation = this->vm()["penal"].template as<value_type>();
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    int bctype = this->vm()["bctype"].template as<int>();

    value_type pi = M_PI;
    //__typeof__( sin(pi*Px()) ) g= sin(pi*Px());
    //__typeof__( pi*pi*sin(pi*Px()) ) f = pi*pi*sin(pi*Px());
    AUTO( g, val( sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) ) );
    AUTO( f, val( Dim*pi*pi*sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) ) );

    vector_ptrtype F( backend_type::newVector( Xh ) );
    timers["assembly"].first.restart();


    form( Xh, *F )  = integrate( elements( *mesh ), im, trans( f )*id( v ) );

    if ( bctype == 1 || !Cont::is_continuous )
        form( Xh, *F, false ) += integrate( boundaryfaces( *mesh ), im, trans( g )*( - grad( v )*N() + penalisation_bc*id( v )/hFace() ) );

    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F.m" );

    timers["assembly"].second = timers["assembly"].first.elapsed();

    backend_ptrtype b( backend_type::build( this->vm() ) );
    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( b->newMatrix( Xh, Xh ) );

    value_type diff = this->vm()["diff"].template as<double>();

    timers["assembly"].first.restart();

    form( Xh, Xh, *D ) = integrate( elements( *mesh ), im, ( diff*gradt( u )*trans( grad( v ) ) ) );

    if ( !Cont::is_continuous )
        form( Xh, Xh, *D, false ) +=integrate( internalfaces( *mesh ), im,
                                               // - {grad(u)} . [v]
                                               -averaget( gradt( u ) )*jump( id( v ) )
                                               // - [u] . {grad(v)}
                                               -average( grad( v ) )*jumpt( idt( u ) )
                                               // penal*[u] . [v]/h_face
                                               + penalisation* ( trans( jumpt( idt( u ) ) )*jump( id( v ) ) )/hFace()
                                             );

    if ( bctype == 1 || !Cont::is_continuous )
        form( Xh, Xh, *D, false ) += integrate( boundaryfaces( *mesh ), im,
                                                ( - trans( id( v ) )*( gradt( u )*N() )
                                                        - trans( idt( u ) )*( grad( v )*N() )
                                                        + penalisation_bc*trans( idt( u ) )*id( v )/hFace() ) );

    else if ( bctype == 0 )
        form( Xh, Xh, *D, false ) += on( boundaryfaces( *mesh ), u, *F, g );

    D->close();
    timers["assembly"].second += timers["assembly"].first.elapsed();

    if ( this->vm().count( "export-matlab" ) )
        D->printMatlab( "D" );

    this->solve( *D, u, *F, ( bctype == 1 || !Cont::is_continuous ) );

    typename space<Continuous>::ptrtype Xch = space<Continuous>::type::New( mesh );
    typename space<Continuous>::element_type uEx( Xch, "uEx" );
    sparse_matrix_ptrtype M( b->newMatrix( Xch, Xch ) );
    form( Xch, Xch, *M ) = integrate( elements( *mesh ), im, trans( idt( uEx ) )*id( uEx ) );
    M->close();
    vector_ptrtype L( b->newVector( Xch ) );
    form( Xch, *L ) = integrate( elements( *mesh ), im, trans( g )*id( uEx ) );
    this->solve( *M, uEx, *L, true );

    std::cout << "||error||_0 = " << math::sqrt( integrate( elements( *mesh ), im, val( ( idv( u )-g )^2 ) ).evaluate()( 0,0 ) ) << "\n";
    std::cout << "||error||_0 = " << math::sqrt( integrate( elements( *mesh ), im, val( trans( idv( u )-g )*( idv( u )-g ) ) ).evaluate()( 0,0 ) ) << "\n";
    std::cout << "||error||_0 = " << math::sqrt( integrate( elements( *mesh ), im, trans( idv( u )-idv( uEx ) )*( idv( u )-idv( uEx ) ) ).evaluate()( 0,0 ) ) << "\n";

    if ( Cont::is_continuous )
        this->exportResults( u, u, uEx );

    else
    {

        form( Xch, *L ) = integrate( elements( *mesh ), im, trans( idv( u ) )*id( uEx ) );
        typename space<Continuous>::element_type uc( Xch, "uc" );
        this->solve( *M, uc, *L, true );
        this->exportResults( u, uc, uEx );
    }

    std::cout << "[timer] run(): init (" << mesh->numElements() << " Elems): " << timers["init"].second << "\n";
    std::cout << "[timer] run(): assembly (" << Xh->dof()->nDof() << " DOFs): " << timers["assembly"].second << "\n";
#endif
} // Laplacian::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
template<typename Mat, typename Vec1, typename Vec2>
void
Laplacian<Dim, Order, Cont, Entity, FType>::solve( Mat& D,
        Vec1& u,
        Vec2& F,
        bool is_sym  )
{
    timers["solver"].first.restart();

    backend_ptrtype b( backend_type::build( this->vm() ) );
    //backend->set_symmetric( is_sym );

    vector_ptrtype U( b->newVector( u.functionSpace() ) );
    b->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    Log() << "[timer] solve: " << timers["solver"].second << "\n";
} // Laplacian::solve


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
template<typename f1_type, typename f2_type, typename f3_type>
void
Laplacian<Dim, Order, Cont, Entity,FType>::exportResults( f1_type& U,
        f2_type& V,
        f3_type& E )
{
    timers["export"].first.restart();

    if ( this->vm().count( "export" ) )
    {
        Log() << "exportResults starts\n";
        exporter->step( 1. )->setMesh( U.functionSpace()->mesh() );

        //exporter->step(1.)->setMesh( this->createMesh( meshSize/2, 0.5, 1 ) );
        //exporter->step(1.)->setMesh( this->createMesh( meshSize, 0, 1 ) );
        if ( !this->vm().count( "export-mesh-only" ) )
        {
            exporter->step( 1. )->add( "pid",
                                       regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );


            exporter->step( 1. )->add( "u", U );
            exporter->step( 1. )->add( "v", V );
            exporter->step( 1. )->add( "e", E );
        }

        exporter->save();


        if ( Dim == 1 )
        {
            std::ostringstream fname_u;
            fname_u << "u.dat";
            std::ofstream ofs3( fname_u.str().c_str() );
            typename mesh_type::element_iterator it = U.functionSpace()->mesh()->beginElement();
            typename mesh_type::element_iterator en = U.functionSpace()->mesh()->endElement();

            for ( ; it!=en; ++it )
            {
                for ( size_type i = 0; i < U.nbLocalDof(); ++i )
                {
                    size_type dof0 = boost::get<0>( U.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                    ofs3 << std::setw( 5 ) << it->id() << " "
                         << std::setw( 5 ) << i << " "
                         << std::setw( 5 ) << dof0 << " "
                         << std::setw( 15 ) << U( dof0 ) << " ";
                    value_type a = it->point( 0 ).node()[0];
                    value_type b = it->point( 1 ).node()[0];

                    if ( i == 0 )
                        ofs3 << a;

                    else if ( i == 1 )
                        ofs3 <<  b;

                    else
                        ofs3 <<  a + ( i-1 )*( b-a )/( V.nbLocalDof()-1 );

                    ofs3 << "\n";

                }
            }

            ofs3.close();

            std::ostringstream fname_v;
            fname_v << "values.dat";
            std::ofstream ofs2( fname_v.str().c_str() );
            it = V.functionSpace()->mesh()->beginElement();
            en = V.functionSpace()->mesh()->endElement();

            for ( ; it!=en; ++it )
            {
                for ( size_type i = 0; i < V.nbLocalDof(); ++i )
                {
                    size_type dof0 = boost::get<0>( V.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                    ofs2 << std::setw( 5 ) << it->id() << " "
                         << std::setw( 5 ) << i << " "
                         << std::setw( 5 ) << dof0 << " "
                         << std::setw( 15 ) << V( dof0 ) << " "
                         << std::setw( 15 ) << E( dof0 ) << " ";
                    value_type a = it->point( 0 ).node()[0];
                    value_type b = it->point( 1 ).node()[0];

                    if ( i == 0 )
                        ofs2 << a;

                    else if ( i == 1 )
                        ofs2 <<  b;

                    else
                        ofs2 <<  a + ( i-1 )*( b-a )/( V.nbLocalDof()-1 );

                    ofs2 << "\n";

                }
            }

        }
    }

    timers["export"].second = timers["export"].first.elapsed();
    Log() << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Laplacian::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 1;
    typedef Continuous MyContinuity;
    //typedef Discontinuous MyContinuity;
    //typedef Feel::Laplacian<nDim, nOrder, MyContinuity, Hypercube, Scalar> laplacian_type;

    typedef Feel::Laplacian<nDim, nOrder, MyContinuity, Simplex, Scalar> laplacian_type;

    /* define and run application */
    laplacian_type laplacian( argc, argv, makeAbout(), makeOptions() );

    laplacian.run();
}




