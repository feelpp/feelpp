/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009-2011 Université Joseph Fourier (Grenoble I)

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
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#if !defined( __FEELPP_BENCH_LAPLACIAN_HPP)
#define __FEELPP_BENCH_LAPLACIAN_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/simget.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>



#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * \class Laplacian class
 * \brief solves the laplacian equations
 *
 */
template<int Dim,
         typename BasisU,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Laplacian
    :
public Simget
{
    typedef Simget super;
public:

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #
    typedef BasisU basis_u_type;
    typedef bases<basis_u_type> basis_type;
    //# endmarker1 #

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Laplacian( std::string const& basis_name )
        :
        super(),
        M_backend(),
        M_basis_name( basis_name ),
        M_exporter()
    {
        mu = option(_name="mu").template as<value_type>();
        penalbc = option(_name="bccoeff").template as<value_type>();
    }


    std::string name() const
    {
        return M_basis_name;
    }

    /**
     * run the convergence test
     */
    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N )
    {
        run();
    }

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    std::string M_basis_name;
    double mu;
    double penalbc;

    boost::shared_ptr<export_type> M_exporter;
}; // Laplacian


template<int Dim, typename BasisU, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, BasisU, Entity>::run()
{
    using namespace Feel::vf;

    int nparts = Environment::worldComm().size();
    bool prepare = boption(_name="benchmark.prepare");
    if ( prepare )
        nparts = ioption(_name="benchmark.partitions");

    this->changeRepository( boost::format( "%1%/%2%/%3%/h_%4%/l_%5%/parts_%6%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % M_basis_name
                            % meshSizeInit()
                            % level()
                            % nparts );


    //! init backend
    M_backend = backend_type::build(soption("backend"));


    boost::mpi::timer t;

    auto mesh = loadMesh( _mesh=new mesh_type, _h=meshSizeInit()/std::pow(2,level()-1), _partitions=nparts );

    M_stats.put( "t.init.mesh",t.elapsed() );
    t.restart();
    if ( prepare ) return;
    /*
     * The function space and some associate elements are then defined
     */
    //# marker4 #

    space_ptrtype Xh = space_type::New( mesh );

    auto u = Xh->element();
    auto v = Xh->element();

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << meshSize() << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";
    LOG(INFO) << "[mesh]   number of elements: " << Xh->mesh()->numElements() << "\n";
    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";

    size_type gnelts=0;
    mpi::all_reduce( this->comm(), Xh->mesh()->numElements() , gnelts, [] ( size_type x, size_type y ) {return x + y;} );

    M_stats.put( "n.space.nelts", gnelts );
    M_stats.put( "n.space.nlocalelts",Xh->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nDof() );
    M_stats.put( "n.space.nlocaldof",Xh->nLocalDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    double penalbc = option(_name="bccoeff").template as<value_type>();
    double mu = option(_name="mu").template as<value_type>();

    //
    // simulation setup: rhs, Dirichlet and exact solution
    //
    std::string dim_str =  boost::str( boost::format( "%1%D" ) % Dim );
    auto vars=symbols<Dim>();
    auto u_exact_g = parse( soption(_name="exact",_prefix=dim_str), vars );
    LOG(INFO) << "u_exact=" << u_exact_g;
    auto u_exact = expr<3>( u_exact_g, vars, "u_exact" );
	auto f_g = -mu*laplacian( u_exact_g, vars );
    LOG(INFO) << "f=" << f_g;
    auto f = expr<3>( f_g, vars, "f" );
    auto grad_exact_g = grad( u_exact_g, vars );
    LOG(INFO) << "grad_exact = " << grad_exact_g;
    auto grad_exact = expr<1,Dim,3>( grad_exact_g, vars, "grad_exact" );

    boost::mpi::timer subt;
    // right hand side
    auto F = M_backend->newVector( Xh );
    form1( Xh, _vector=F, _init=true );
    M_stats.put( "t.init.vector",t.elapsed() );
    LOG(INFO) << "  -- time for vector init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;
    form1( Xh, _vector=F ) = integrate( elements( mesh ), trans( f )*id( v ) );
    M_stats.put( "t.assembly.vector.source",subt.elapsed() );
    subt.restart();

    if ( this->vm()[ "bctype" ].template as<int>() == 1  )
    {
        form1( Xh, _vector=F ) += integrate( _range=boundaryfaces( mesh ), _expr=-mu*trans( u_exact )*dn( u ) );
        M_stats.put( "t.assembly.vector.dirichlet1",subt.elapsed() );
        subt.restart();
        form1( Xh, _vector=F ) += integrate( _range=boundaryfaces( mesh ), _expr=mu*penalbc*inner( u_exact,id( v ) )/hFace() );
        M_stats.put( "t.assembly.vector.dirichlet2",subt.elapsed() );
        LOG(INFO) << "   o time for rhs weak dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    M_stats.put( "t.assembly.vector.total",t.elapsed() );
    LOG(INFO) << "  -- time vector global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    /*
     * Construction of the left hand side
     */
    t.restart();
    auto D = M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, _matrix=D, _init=true );
    M_stats.put( "t.init.matrix",t.elapsed() );
    LOG(INFO) << "  -- time for matrix init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    subt.restart();
    t.restart();

    form2( Xh, Xh, _matrix=D ) =integrate( _range=elements( mesh ),_expr=mu*( gradt( u )*trans(grad( v ) )) );
    M_stats.put( "t.assembly.matrix.diffusion",subt.elapsed() );
    LOG(INFO) << "   o time for diffusion terms: " << subt.elapsed() << "\n";
    subt.restart();

    if ( this->vm()[ "bctype" ].template as<int>() == 1  )
    {
        form2( Xh, Xh, _matrix=D )+=integrate( _range=boundaryfaces( mesh ),_expr=mu*(-trans( dnt( u ) )*id( v )-trans( dn( u ) )*idt( v ) ) );
        M_stats.put( "t.assembly.matrix.dirichlet1",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, _matrix=D )+=integrate( _range=boundaryfaces( mesh ),_expr=mu*penalbc*inner( idt( u ),id( v ) )/hFace() );
        M_stats.put( "t.assembly.matrix.dirichlet2",subt.elapsed() );
        subt.restart();
        LOG(INFO) << "   o time for weak dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    //# endmarker7 #
    D->close();
    F->close();

    if ( this->vm()[ "bctype" ].template as<int>() == 0  )
    {
        form2( Xh, Xh, _matrix=D ) += on( _range=boundaryfaces( mesh ), _element=u, _rhs=F, _expr=u_exact );
        M_stats.put( "t.assembly.matrix.dirichlet",subt.elapsed() );
        LOG(INFO) << "   o time for strong dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    M_stats.put( "t.assembly.matrix.total",t.elapsed() );
    LOG(INFO) << " -- time matrix global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;



    t.restart();

    if ( !this->vm().count( "no-solve" ) )
    {
        auto r = M_backend->solve( _matrix=D, _solution=u, _rhs=F );
        M_stats.put( "d.solver.bool.converged",r.template get<0>() );
        M_stats.put( "d.solver.int.nit",r.template get<1>() );
        M_stats.put( "d.solver.double.residual",r.template get<2>() );
    }

    M_stats.put( "t.solver.total",t.elapsed() );
    LOG(INFO) << " -- time for solver : "<<t.elapsed()<<" seconds \n";
    t.restart();


    double meas = integrate( _range=elements( mesh ), _expr=constant( 1.0 ) ).evaluate()( 0, 0 );
    LOG(INFO) << "[laplacian] measure(Omega)=" << meas << " (should be equal to 4)\n";
    std::cout << "[laplacian] measure(Omega)=" << meas << " (should be equal to 4)\n";

    double mean_u = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 )/meas;
    LOG(INFO) << "[laplacian] mean(u)=" << mean_u << "\n";
    std::cout << "[laplacian] mean(u)=" << mean_u << "\n";

    double mean_uexact = integrate( elements( mesh ), u_exact ).evaluate()( 0, 0 )/meas;
    std::cout << "[laplacian] mean(uexact)=" << mean_uexact << "\n";
    M_stats.put( "t.integrate.mean",t.elapsed() );
    size_type nnz = 0 ;
    auto nNz = D->graph()->nNz() ;

    for ( auto iter = nNz.begin(); iter!=nNz.end(); ++iter )
        nnz += ( *iter ) ;


    size_type gnnz=0;
    mpi::all_reduce( this->comm(), nnz, gnnz, [] ( size_type x, size_type y ) {return x + y;} );
    LOG(INFO) << "[laplacian] matrix NNZ local:"<< nnz << " global: "  << gnnz << "\n";
    M_stats.put( "n.matrix.nnz",gnnz );
    M_stats.put( "n.matrix.nlocalnz",nnz );

    v.zero();
    v = vf::project( Xh, elements(mesh), u_exact );
    double interp_errorL2 = integrate( _range=elements( mesh ), _expr=trans( idv( v )-u_exact )*( idv( v )-u_exact ), _quad=_Q<8>() ).evaluate()( 0, 0 );

    t.restart();
    double u_errorL2 = integrate( _range=elements( mesh ), _expr=trans( idv( u )-u_exact )*( idv( u )-u_exact ) ).evaluate()( 0, 0 );
    double u_exactL2 = integrate( _range=elements( mesh ), _expr=trans( u_exact )*( u_exact ) ).evaluate()( 0, 0 );
    M_stats.put( "t.integrate.l2norm",t.elapsed() );
    t.restart();
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";
    LOG(INFO) << "||u_error||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";
    M_stats.put( "e.l2.i",math::sqrt( interp_errorL2 ) );
    M_stats.put( "e.l2.u",math::sqrt( u_errorL2 ) );
    M_stats.put( "e.l2.u.n",math::sqrt( u_errorL2/u_exactL2 ) );

    double u_errorsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=trace( ( gradv( u )-grad_exact )*trans( gradv( u )-grad_exact ) ) ).evaluate()( 0, 0 );
    double u_exactsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=trace( ( grad_exact )*trans( grad_exact ) ) ).evaluate()( 0, 0 );
    double u_error_H1 = math::sqrt( ( u_errorL2+u_errorsemiH1 )/( u_exactL2+u_exactsemiH1 ) );

    M_stats.put( "t.integrate.h1norm",t.elapsed() );
    t.restart();
    LOG(INFO) << "||u_error||_1= " << u_error_H1 << "\n";
    M_stats.put( "e.h1.u",  math::sqrt( u_errorsemiH1 + u_errorL2 ));
    M_stats.put( "e.h1.u.n",u_error_H1 );
    M_stats.put( "e.semih1.u.n", math::sqrt( u_errorsemiH1/u_exactsemiH1 ) );
    M_stats.put( "e.semih1.u", math::sqrt( u_errorsemiH1 ) );

    // compute flux through the boundary
    // first using normal derivative
    double flux1_exact = integrate( _range=boundaryfaces(mesh), _expr=-mu*grad_exact*N(), _quad=_Q<10>()).evaluate()(0,0);
    LOG(INFO) << "flux exact : "  << std::scientific << std::setprecision(5) << flux1_exact;
    double flux1_approx = integrate( _range=boundaryfaces(mesh), _expr=-mu*gradv(u)*N()).evaluate()(0,0);
    LOG(INFO) << "flux approx : "  << std::scientific << std::setprecision(5) << flux1_approx;
    //double flux1_error = math::abs(flux1_exact-flux1_approx)/math::abs(flux1_exact);
    double flux1_error = math::abs(flux1_exact-flux1_approx);
    LOG(INFO) << "relative error flux = " <<   flux1_error;
    M_stats.put( "e.flux.dn", flux1_error );

    double flux1_exact_part = integrate( _range=markedfaces(mesh, "Dirichlet"), _expr=-mu*grad_exact*N(), _quad=_Q<10>()).evaluate()(0,0);
    LOG(INFO) << "flux exact part Dirichlet : "  << std::scientific << std::setprecision(5) << flux1_exact_part;
    double flux1_approx_part = integrate( _range=markedfaces(mesh, "Dirichlet"), _expr=-mu*gradv(u)*N()).evaluate()(0,0);
    LOG(INFO) << "flux approx part : "  << std::scientific << std::setprecision(5) << flux1_approx_part;
    //double flux1_error = math::abs(flux1_exact-flux1_approx)/math::abs(flux1_exact);
    double flux1_error_part = math::abs(flux1_exact_part-flux1_approx_part);
    LOG(INFO) << "relative error flux part = " <<   flux1_error_part;
    M_stats.put( "e.part.dn", flux1_error_part );

    // second using residual from equation without dirichlet conditions
    v.zero();
    v = vf::project( Xh, boundaryfaces(mesh), cst(1.) );
    double flux2_approx = integrate( _range=boundaryelements(mesh), _expr=f*idv(v)-mu*gradv(u)*trans(gradv(v))).evaluate()(0,0);
    double flux2_error = math::abs(flux1_exact-flux2_approx);
    LOG(INFO) << "residual approx : "  << std::scientific << std::setprecision(5) << flux2_approx;
    LOG(INFO) << "relative error flux 2 = " <<   flux2_error;
    M_stats.put( "e.flux.residual", flux2_error );
    double flux2_exact = integrate( _range=boundaryelements(mesh), _expr=f*idv(v)-mu*grad_exact*trans(gradv(v)),_quad=_Q<10>()).evaluate()(0,0);
    double flux3_error = math::abs(flux2_exact-flux2_approx);
    LOG(INFO) << "relative error flux 3 = " <<   flux3_error;
    M_stats.put( "e.flux.residual2", flux3_error );


    // second using residual from equation without dirichlet conditions on a part of the boundary
    v.zero();
    v = vf::project( Xh, markedfaces(mesh,"Dirichlet"), cst(1.) );
    double flux2_approx_part = integrate( _range=boundaryelements(mesh), _expr=f*idv(v)-mu*gradv(u)*trans(gradv(v))).evaluate()(0,0);
    flux2_approx_part += integrate( _range=markedfaces(mesh, "Neumann"), _expr=mu*gradv(u)*N()*idv(v)).evaluate()(0,0);
    double flux2_error_part = math::abs(flux1_exact_part-flux2_approx_part);
    LOG(INFO) << "residual approx part : "  << std::scientific << std::setprecision(5) << flux2_approx_part;
    LOG(INFO) << "relative error flux 2 part = " <<   flux2_error_part;
    M_stats.put( "e.part.residual", flux2_error_part );
    double flux2_exact_part = integrate( _range=boundaryelements(mesh), _expr=f*idv(v)-mu*grad_exact*trans(gradv(v)),_quad=_Q<10>()).evaluate()(0,0);
    flux2_exact_part += integrate( _range=markedfaces(mesh, "Neumann"), _expr=mu*grad_exact*N()*idv(v), _quad=_Q<10>()).evaluate()(0,0);
    double flux3_error_part = math::abs(flux2_exact_part-flux2_approx_part);
    LOG(INFO) << "relative error flux 3 part = " <<   flux3_error_part;
    M_stats.put( "e.part.residual2", flux3_error_part );
    M_stats.put( "d.part.fluxexact1", flux1_exact_part );
    M_stats.put( "d.part.fluxexact2", flux2_exact_part );
    M_stats.put( "d.part.flux", flux1_approx_part );
    M_stats.put( "d.part.residual", flux2_approx_part );


    // export results
    v = vf::project( Xh, elements( Xh->mesh() ), u_exact );
    M_stats.put( "t.export.projection",t.elapsed() );
    t.restart();
    this->exportResults( u, v );
    M_stats.put( "t.export.total",t.elapsed() );
    t.restart();

} // Laplacian::run


template<int Dim, typename BasisU, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, BasisU, Entity>::exportResults( element_type& u, element_type& v )
{
    M_exporter =  exporter( _mesh=u.functionSpace()->mesh() );
    if ( M_exporter->doExport() )
    {
        //M_exporter->step( 0 )->setMesh( u.functionSpace()->mesh() );
        M_exporter->step( 0 )->add( "u", u );
        M_exporter->step( 0 )->add( "u_exact", v );
        M_exporter->save();
    }
} // Laplacian::export
} // Feel

#endif // __FEELPP_BENCH_LAPLACIAN_HPP
