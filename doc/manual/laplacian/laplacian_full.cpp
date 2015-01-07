/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-23

  Copyright (C) 2006-2009 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-11-23
 */
#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain" )
    ( "dt", Feel::po::value<double>()->default_value( 1 ), "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ), "final time value" )

    ( "diff", Feel::po::value<double>()->default_value( 1 ), "diffusion parameter" )
    ( "penal", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter" )
    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
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
    Feel::AboutData about( "laplacian" ,
                           "laplacian" ,
                           "0.2",
                           "nD(n=1,2,3) Laplacian on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006, 2007 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
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

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    template<typename Conti = Cont>
    struct space
    {
        /*basis*/
        typedef bases<Lagrange<Order, FType> > basis_type;


#if 0
        typedef typename mpl::if_<mpl::bool_<Conti::is_continuous>,
                mpl::identity<bases<Lagrange<Order, FType> > >,
                mpl::identity<bases<OrthonormalPolynomialSet<Order, FType> > > >::type::type basis_type;
#endif
        /*space*/
        typedef FunctionSpace<mesh_type, basis_type, Conti, value_type> type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef typename type::element_type element_type;
        //typedef typename element_type::template sub_element<0>::type element_0_type;
        //typedef typename element_type::template sub_element<1>::type element_1_type;
    };

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    Laplacian()
        :
        super(),
        meshSize( doption("hsize") ),
        shape( soption("shape") ),
        b(  backend_type::build( soption("backend") ) ),
        bc( backend_type::build( soption("backend") ) ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers(),
        stats()
    {
    }

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
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( double time,
                        f1_type& u,
                        f2_type& v,
                        f3_type& e );

private:

    double meshSize;
    std::string shape;
    backend_ptrtype b,bc;
    export_ptrtype exporter;
    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // Laplacian


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


    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % doption("hsize")
                          );
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                _usenames=true,
                                                _convex=entity_type::type(),
                                                _shape=shape,
                                                _dim=Dim,
                                                _h=meshSize ) );
    stats["nelt"] = mesh->elements().size();

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    auto Xh = space<Cont>::type::New( mesh );
    //Xh->dof()->showMe();
    auto u = Xh->element( "u" );
    auto v = Xh->element( "v" );
    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();

    if ( this->vm().count( "export-mesh-only" ) )
        this->exportResults( 0., u, u, u );

    value_type penalisation = this->vm()["penal"].template as<value_type>();
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    int bctype = ioption("bctype");

    double t = 0;
    auto g = val( exp( -cst_ref( t ) )*sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) );
    auto f = ( -1*( t > 0 )+pi*pi*Dim )*g;

    timers["assembly"].second = timers["assembly"].first.elapsed();

    //
    // Compute domain and boundary size
    //
    double local_domain_size = integrate( elements( mesh ),  constant( 1.0 ) ).evaluate()( 0,0 );
    double global_domain_size;
    mpi::all_reduce( Application::comm(), local_domain_size, global_domain_size, std::plus<double>() );
    LOG(INFO) << "int_Omega = " << global_domain_size << "[ " << local_domain_size << " ]\n";
    double local_boundary_size = integrate( boundaryfaces( mesh ),  constant( 1.0 ) ).evaluate()( 0,0 );
    double global_boundary_size;
    mpi::all_reduce( Application::comm(), local_boundary_size, global_boundary_size, std::plus<double>() );
    LOG(INFO) << "int_Omega = " << global_boundary_size << "[ " << local_boundary_size << " ]\n";

    /*
     * Construction of the left hand side
     */
    auto D  = b->newMatrix( Xh, Xh );
    auto Mt = b->newMatrix( Xh, Xh );

    value_type diff = doption("diff");
    value_type dt   = doption("dt"  );
    value_type ft   = doption("ft"  );



    timers["assembly"].first.restart();

    size_type pattern = ( Cont::is_continuous?Pattern::COUPLED:Pattern::COUPLED|Pattern::EXTENDED );
    form2( Xh, Xh, _matrix=Mt, _init=true, _pattern=pattern ) =
        integrate( elements( mesh ),
                   idt( u )*id( v )/dt );
    Mt->close();
    form2( Xh, Xh, _matrix=D, _init=true, _pattern=pattern ) =
        integrate( elements( mesh ),
                   diff*gradt( u )*trans( grad( v ) )
                 );
    D->close();
    vector_ptrtype Un( b->newVector( u.functionSpace() ) );
    vector_ptrtype Vn( b->newVector( u.functionSpace() ) );
    u = vf::project( Xh, elements( mesh ), g );
    v = vf::project( Xh, elements( mesh ), constant( 1.0 ) );

    std::cout << "int Px() = "  << integrate( elements( mesh ),  Px() ).evaluate() << "\n";
    std::cout << "int idv(1) = "  << integrate( elements( mesh ),  idv( v ) ).evaluate() << "\n";
    v = vf::project( Xh, elements( mesh ), Px() );
    std::cout << "int idv(x) = "  << integrate( elements( mesh ),  idv( v ) ).evaluate() << "\n";
    std::cout << "int gradv(x) = "  << integrate( elements( mesh ),  gradv( v ) ).evaluate() << "\n";
    v = vf::project( Xh, elements( mesh ), Px()*Px() );
    std::cout << "int Px()*Px() = "  << integrate( elements( mesh ),  Px()*Px() ).evaluate() << "\n";
    std::cout << "int idv(x) = "  << integrate( elements( mesh ),  idv( v ) ).evaluate() << "\n";
    std::cout << "int gradv(x) = "  << integrate( elements( mesh ),  gradv( v ) ).evaluate() << "\n";

    for ( size_type i = 0; i < u.size(); ++i )
    {
        Un->set( i, u( i ) );
        Vn->set( i, v( i ) );
    }

    Un->close();
    Vn->close();
    std::cout << "||u||_energy = " << D->energy( Vn, Un ) << "\n"
              << "||u||_energy_2 = " << integrate( elements( mesh ),  gradv( u )*trans( gradv( v ) ) ).evaluate() << "\n"
              << "||u||_energy_2' = " << integrate( boundaryfaces( mesh ),  ( gradv( u )*N() )*idv( v ) ).evaluate() << "\n"
              << "||u||_energy_3 = " << integrate( elements( mesh ),  -trace( hessv( u ) )*idv( v ) ).evaluate() << "\n";

    form2( Xh, Xh, _matrix=D, _init=true, _pattern=pattern ) =
        integrate( elements( mesh ),
                   diff*gradt( u )*trans( grad( v ) )
                   //-trace(hesst(u))*id(v)
                 );

    if ( !Cont::is_continuous )
        form2( Xh, Xh, _matrix=D ) +=integrate( internalfaces( mesh ),
                                        // - {grad(u)} . [v]
                                        -averaget( gradt( u ) )*jump( id( v ) )
                                        // - [u] . {grad(v)}
                                        -average( grad( v ) )*jumpt( idt( u ) )
                                        // penal*[u] . [v]/h_face
                                        + penalisation* ( trans( jumpt( idt( u ) ) )*jump( id( v ) ) )/hFace()
                                      );

    if ( bctype == 1 || !Cont::is_continuous )
    {
        form2( Xh, Xh, _matrix=D ) += integrate( boundaryfaces( mesh ),
                                         ( - trans( id( v ) )*( gradt( u )*N() )
                                           - trans( idt( u ) )*( grad( v )*N() )
                                           + penalisation_bc*trans( idt( u ) )*id( v )/hFace() ) );
    }

    D->close();

    if ( this->vm().count( "export-matlab" ) )
        D->printMatlab( "D.m" );

    //
    // Construct L2-projection operator
    //
    auto Xch = space<Continuous>::type::New( mesh );
    auto uEx = Xch->element( "uEx" );
    auto M = bc->newMatrix( Xch, Xch );
    form2( Xch, Xch, _matrix=M, _init=true ) = integrate( elements( mesh ),  trans( idt( uEx ) )*id( uEx ) );
    M->close();
    auto L = bc->newVector( Xch );

    //
    // initially u = g(0,x,y,z)
    //
    u.zero();
    //
    //  construct right hand side
    //
    vector_ptrtype Ft( b->newVector( Xh ) );
    timers["assembly"].first.restart();
    form1( _test=Xh, _vector=Ft, _init=true )  = integrate( elements( mesh ),
            trans( f )*id( v ) );

    if ( bctype == 1 || !Cont::is_continuous )
        form1( Xh, _vector=Ft ) += integrate( boundaryfaces( mesh ),
                                      trans( g )*( - grad( v )*N() + penalisation_bc*id( v )/hFace() ) );

    Ft->close();

    if ( this->vm().count( "export-matlab" ) )
        Ft->printMatlab( "F.m" );

    timers["assembly"].second += timers["assembly"].first.elapsed();

    form1( Xch, _vector=L, _init=true ) = integrate( elements( mesh ),  trans( g )*id( uEx ) );
    L->close();

    if ( this->vm().count( "export-matlab" ) )
        L->printMatlab( "L.m" );

    // compute PDE solution
    b->solve( _matrix=D, _solution=u, _rhs=Ft );

    // compute L2 projection of exact solution
    bc->solve( _matrix=M, _solution=uEx, _rhs=L );



    // compute local error wrt the exact solution
    double error = integrate( elements( mesh ),  val( trans( idv( u )-g )*( idv( u )-g ) ) ).evaluate()( 0,0 );
    double errorex = integrate( elements( mesh ),  trans( idv( u )-idv( uEx ) )*( idv( u )-idv( uEx ) ) ).evaluate()( 0,0 );
    double local_jump = integrate( internalfaces( mesh ),  trans( jumpv( idv( u ) ) )*( jumpv( idv( u ) ) ) ).evaluate()( 0,0 );

    // compute global error : suming all contributions from
    // all processors
    double global_error = 0;
    mpi::all_reduce( Application::comm(), error, global_error, std::plus<double>() );
    double global_errorex = 0;
    mpi::all_reduce( Application::comm(), errorex, global_errorex, std::plus<double>() );
    double global_jump = 0;
    mpi::all_reduce( Application::comm(), local_jump, global_jump, std::plus<double>() );

    LOG(INFO) << "||[u]||= " << math::sqrt( global_jump )  << "\n";
    //LOG(INFO) << "local ||error||_0 = " << math::sqrt( error ) << "\n";
    LOG(INFO) << "global ||error||_0 = " << math::sqrt( global_error ) << "\n";
    //LOG(INFO) << "local ||error(l2proj)||_0 = " << math::sqrt( errorex ) << "\n";
    LOG(INFO) << "global ||errorex(l2proj)||_0 = " << math::sqrt( global_errorex ) << "\n";

    // exporting the results using the exporter set by the
    // command line: the default being the 'ensight' exporter
    if ( Cont::is_continuous )
        this->exportResults( t, u, u, uEx );

    else
    {

        form1( Xch, _vector=L, _init=true ) = integrate( elements( mesh ),  trans( idv( u ) )*id( uEx ) );
        typename space<Continuous>::element_type uc( Xch, "uc" );
        bc->solve( _matrix=M, _solution=uc, _rhs=L );
        this->exportResults( t, u, uc, uEx );
    }

    D->addMatrix( 1.0, *Mt );

    //
    // Time loop
    //
    for ( t = dt; t < ft; t += dt )
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "time : " << t << "s dt: " << dt << " ft: "<< ft << "s\n";

        //
        //  construct right hand side
        //
        timers["assembly"].first.restart();
        form1( _test=Xh, _vector=Ft, _init=true )  = integrate( elements( mesh ),
                trans( f )*id( v ) + idv( u )*id( v )/dt );

        if ( bctype == 1 || !Cont::is_continuous )
            form1( Xh, _vector=Ft ) += integrate( boundaryfaces( mesh ),
                                          trans( g )*( - grad( v )*N() + penalisation_bc*id( v )/hFace() ) );

        Ft->close();
        timers["assembly"].second += timers["assembly"].first.elapsed();

        if ( this->vm().count( "export-matlab" ) )
        {
            std::ostringstream ostr;
            ostr.precision( 3 );
            ostr << "Ft-" << t;
            Ft->printMatlab( ostr.str() );
        }

        b->solve( _matrix=D, _solution=u, _rhs=Ft );

        if ( this->vm().count( "export-matlab" ) )
        {
            std::ostringstream ostr;
            ostr.precision( 3 );
            ostr << "u-" << t;
            u.printMatlab( ostr.str() );
        }

        // find the L2-projection (with a continuous expansion) of
        // the exact solution
        form1( Xch, _vector=L, _init=true ) = integrate( elements( mesh ),  trans( g )*id( uEx ) );
        bc->solve( _matrix=M, _solution=uEx, _rhs=L );

        // compute local error wrt the exact solution
        double error = integrate( elements( mesh ),  val( trans( idv( u )-g )*( idv( u )-g ) ) ).evaluate()( 0,0 );
        double errorex = integrate( elements( mesh ),  trans( idv( u )-idv( uEx ) )*( idv( u )-idv( uEx ) ) ).evaluate()( 0,0 );

        // compute global error : suming all contributions from
        // all processors
        double global_error = 0;
        mpi::all_reduce( Application::comm(), error, global_error, std::plus<double>() );
        double global_errorex = 0;
        mpi::all_reduce( Application::comm(), errorex, global_errorex, std::plus<double>() );

        LOG(INFO) << "local ||error||_0 = " << math::sqrt( error ) << "\n";
        LOG(INFO) << "global ||error||_0 = " << math::sqrt( global_error ) << "\n";
        LOG(INFO) << "local ||error(l2proj)||_0 = " << math::sqrt( errorex ) << "\n";
        LOG(INFO) << "global ||errorex(l2proj)||_0 = " << math::sqrt( global_errorex ) << "\n";

        // compute the value of the field u at the middle of the
        // domain
        node_type __n( Dim );

        __n[0]=0.5;

        if ( Dim >= 2 )
            __n[1]=0.5;

        if ( Dim >= 3 )
            __n[2]=0.5;

        LOG(INFO) << "u value at (0.5,0.5,0.5)= " << u( __n ) << "\n";
        LOG(INFO) << "e value at (0.5,0.5,0.5)= " << uEx( __n ) << "\n";


        // exporting the results using the exporter set by the
        // command line: the default being the 'ensight' exporter
        if ( Cont::is_continuous )
            this->exportResults( t, u, u, uEx );

        else
        {

            form1( Xch, _vector=L, _init=true ) = integrate( elements( mesh ),  trans( idv( u ) )*id( uEx ) );
            typename space<Continuous>::element_type uc( Xch, "uc" );
            bc->solve( _matrix=M, _solution=uc, _rhs=L );
            this->exportResults( t, u, uc, uEx );
        }

        LOG(INFO) << "[timer] run(): init (" << mesh->numElements() << " Elems): " << timers["init"].second << "\n";
        LOG(INFO) << "[timer] run(): assembly (" << Xh->dof()->nDof() << " DOFs): " << timers["assembly"].second << "\n";
    } // time loop

} // Laplacian::run


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
template<typename f1_type, typename f2_type, typename f3_type>
void
Laplacian<Dim, Order, Cont, Entity,FType>::exportResults( double time,
        f1_type& U,
        f2_type& V,
        f3_type& E )
{
    timers["export"].first.restart();

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
        fname_u << "u-" << Application::processId() << ".dat";
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
        fname_v << "values-" << Application::processId() << ".dat";
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

    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Laplacian::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    //typedef Continuous MyContinuity;
    typedef Discontinuous MyContinuity;
    typedef Feel::Laplacian<nDim, nOrder, MyContinuity, Hypercube, Scalar> laplacian_type;

    //typedef Feel::Laplacian<nDim, nOrder, MyContinuity, Simplex, Scalar> laplacian_type;

    /* define and run application */
    laplacian_type laplacian;

    laplacian.run();
}




