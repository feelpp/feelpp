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
   \file stokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#if !defined( __FEELPP_BENCH_STOKES_HPP)
#define __FEELPP_BENCH_STOKES_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/simget.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>



#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

void printStats( std::ostream& out, std::vector<std::map<std::string,boost::any> > & stats );



namespace Feel
{
/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
template<int Dim,
         typename BasisU,
         typename BasisP,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Stokes
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
    typedef BasisP basis_p_type;
    typedef Lagrange<0,Scalar, Continuous> basis_l_type;
    typedef bases<basis_u_type,basis_p_type,basis_l_type> basis_type;
    //# endmarker1 #

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Stokes( std::string const& basis_name )
        :
        super(),
        M_backend(),
        M_basis_name( basis_name ),
        M_exp()
    {
        mu = option(_name="mu").template as<value_type>();
        penalbc = option(_name=prefixvm( name(),"bccoeff" )).template as<value_type>();
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

    boost::shared_ptr<export_type> M_exp;
}; // Stokes


template<int Dim, typename BasisU, typename BasisP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, BasisU, BasisP, Entity>::run()
{
    using namespace Feel::vf;

    int nparts = Environment::worldComm().size();
    bool prepare = option(_name="benchmark.prepare").template as<bool>();
    if ( prepare )
        nparts = option(_name="benchmark.partitions").template as<int>();


    this->changeRepository( boost::format( "perf/%1%/%2%/%3%/h_%4%/l_%5%/parts_%6%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % M_basis_name
                            % meshSizeInit()
                            % level()
                            % nparts );

    //! init backend
    M_backend = backend(_rebuild=true);


    boost::mpi::timer t;
    auto mesh = loadMesh( _mesh=new mesh_type, _refine=level() );

    M_stats.put( "t.init.mesh",t.elapsed() );
    t.restart();
    /*
     * The function space and some associate elements are then defined
     */
    //# marker4 #

    space_ptrtype Xh = space_type::New( mesh );

    auto U = Xh->element( "u" );
    auto V = Xh->element( "v" );
    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();
    auto nu = U.template element<2>();


    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << M_meshSize << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";
    LOG(INFO) << "[mesh]   number of elements: " << Xh->mesh()->numElements() << "\n";
    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    M_stats.put( "h",M_meshSize );
    M_stats.put( "n.space.nelts",Xh->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nLocalDof() );
    M_stats.put( "n.space.ndof.u",Xh->template functionSpace<0>()->nDof() );
    M_stats.put( "n.space.ndof.p",Xh->template functionSpace<1>()->nDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
    t.restart() ;



    //# marker5 #
    auto deft = gradt( u );
    auto def = grad( v );
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+mu*dnt( u );

    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+mu*dn( u );
    //# endmarker6 #

    double betacoeff = option(_name="beta").template as<value_type>();
    bool add_convection = ( math::abs( betacoeff  ) > 1e-10 );
    double mu = option(_name="mu").template as<value_type>();

    auto r=1.;
    auto L=5.;
    auto pin = 10.;
    auto pout = 0.;
    symbol x("x"),y("y"),z("z");
    auto u1 = (1-y*y/r*r)*(pin-pout)/(2*mu*L);
    auto u2 = 0.;
    matrix u_exact_g = matrix(Dim,1);
    u_exact_g = u1,u2;
    auto p_exact_g = (pout-pin)*x/L+pin;


    auto u_exact = expr<2,1,2>( u_exact_g, {x,y} );
    auto p_exact = expr( p_exact_g, {x,y} );
	auto f_g = -mu*laplacian( u_exact_g, {x,y} ) + grad( p_exact_g, {x,y} ).transpose();
    auto f = expr<2,1,2>( f_g, {x,y} );
    LOG(INFO) << "f = " << f_g << "\n";
    auto beta=u_exact;


    auto gradu_exact_g = grad( u_exact_g, {x,y} );
    auto divu_exact_g = div( u_exact_g, {x,y} );
    LOG(INFO) << "gradu_exact_g = " << gradu_exact_g;
    LOG(INFO) << "divu_exact_g = " << divu_exact_g;
    auto gradu_exact = expr<2,2,2>( gradu_exact_g, {x,y} );
    auto divu_exact = expr<1,1,2>( divu_exact_g, {x,y} );
    auto convection=gradu_exact*beta;

    boost::mpi::timer subt;
    // right hand side
    auto rhs = form1( _test=Xh );
    M_stats.put( "t.init.rhs",t.elapsed() );
    LOG(INFO) << "  -- time for vector init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;
    rhs = integrate( elements( mesh ), trans( f )*id( v ) );
    M_stats.put( "t.assembly.rhs.source",subt.elapsed() );
    subt.restart();


    if ( add_convection )
    {
        rhs += integrate( elements( mesh ), trans( convection )*id( v ) );
        M_stats.put( "t.assembly.rhs.convection",subt.elapsed() );
        subt.restart();
    }


    if ( option(_name= prefixvm( name(),"bctype" ) ).template as<int>() == 1  )
    {
        rhs += integrate( _range=boundaryfaces( mesh ), _expr=-trans( u_exact )*SigmaN );
        M_stats.put( "t.assembly.rhs.dirichletup",subt.elapsed() );
        subt.restart();
        rhs += integrate( _range=boundaryfaces( mesh ), _expr=penalbc*inner( u_exact,id( v ) )/hFace() );
        //form1( Xh, F ) += integrate( _range=boundaryfaces(mesh), _expr=penalbc*max(betacoeff,mu/hFace())*(trans(id(v))*N())*N());
        M_stats.put( "t.assembly.rhs.dirichletp",subt.elapsed() );
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
    int pattern = Pattern::DEFAULT;
    int patternsym = Pattern::DEFAULT;

    auto a = form2( _test=Xh, _trial=Xh );
    M_stats.put( "t.init.lhs",t.elapsed() );
    LOG(INFO) << "  -- time for matrix init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    subt.restart();
    t.restart();

    a =integrate( _range=elements( mesh ),_expr=mu*( inner( gradt( u ),grad( v ) ) ) );
    M_stats.put( "t.assembly.lhs.diffusion",subt.elapsed() );
    LOG(INFO) << "   o time for diffusion terms: " << subt.elapsed() << "\n";
    subt.restart();

    if ( add_convection )
    {
        a += integrate( _range=elements( mesh ),_expr=trans( gradt( u )*beta )*id( v ) );
        M_stats.put( "t.assembly.lhs.convection",subt.elapsed() );
        LOG(INFO) << "   o time for convection terms: " << subt.elapsed() << "\n";
        subt.restart();
    }
    a +=integrate( _range=elements( mesh ),_expr=-div( v )*idt( p ) );
    a+=integrate( _range=elements( mesh ),_expr=divt( u )*id( q ) );
    M_stats.put( "t.assembly.lhs.up",subt.elapsed() );
    LOG(INFO) << "   o time for velocity/pressure terms: " << subt.elapsed() << "\n";
    subt.restart();

    a+=integrate( _range=elements( mesh ),_expr=id( p )*idt( nu ) );
    a+=integrate( _range=elements( mesh ),_expr=idt( p )*id( nu ) );
    M_stats.put( "t.assembly.lhs.pl",subt.elapsed() );
    LOG(INFO) << "   o time for pressure/multiplier terms: " << subt.elapsed() << "\n";
    subt.restart();

    if ( option(_name= prefixvm( name(),"bctype" ) ).template as<int>() == 1  )
    {
        //form2( Xh, Xh, D, _pattern=patternsym )+=integrate( _range=boundaryfaces(mesh),_expr=-trans(SigmaNt)*id(v) );
        a+=integrate( _range=boundaryfaces( mesh ),_expr=-trans( SigmaNt )*id( v ) );
        M_stats.put( "t.assembly.lhs.dirichlet1",subt.elapsed() );
        subt.restart();
        a+=integrate( _range=boundaryfaces( mesh ),_expr=-trans( SigmaN )*idt( v ) );
        M_stats.put( "t.assembly.lhs.dirichlet2",subt.elapsed() );
        subt.restart();

        a+=integrate( _range=boundaryfaces( mesh ),_expr=+penalbc*inner( idt( u ),id( v ) )/hFace() );
        M_stats.put( "t.assembly.lhs.dirichlet_u*u",subt.elapsed() );
        subt.restart();
        LOG(INFO) << "   o time for weak dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    M_stats.put( "t.assembly.lhs.total",t.elapsed() );
    LOG(INFO) << " -- time lhs global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    if ( option(_name=prefixvm( name(),"bctype" ) ).template as<int>() == 0  )
    {
        a += on( _range=boundaryfaces( mesh ), _element=u, _rhs=rhs.vectorPtr(), _expr=u_exact );
        M_stats.put( "t.assembly.lhs.dirichlet",subt.elapsed() );
        LOG(INFO) << "   o time for strong dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }



    t.restart();

    if ( !this->vm().count( "no-solve" ) )
    {
        auto r = a.solve( _rhs=rhs, _solution=U );
        M_stats.put( "d.solver.bool.converged",r.template get<0>() );
        M_stats.put( "d.solver.int.nit",r.template get<1>() );
        M_stats.put( "d.solver.double.residual",r.template get<2>() );
    }

    M_stats.put( "t.solver.total",t.elapsed() );
    LOG(INFO) << " -- time for solver : "<<t.elapsed()<<" seconds \n";


    double meas = mean( _range=elements( mesh ), _expr=constant( 1.0 ) )(0,0);
    LOG(INFO) << "[stokes] measure(Omega)=" << meas << " (should be equal to 4)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 4)\n";

    double mean_p = mean( elements( mesh ), idv( p ))(0,0);
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    // get the zero mean pressure
    p.add( - mean_p );
    mean_p = mean( elements( mesh ), idv( p ) )(0,0);
    LOG(INFO) << "[stokes] mean(p-mean(p))=" << mean_p << "\n";
    double mean_pexact = mean( elements( mesh ), p_exact )(0,0);
    size_type nnz = 0 ;
    auto nNz = a.matrixPtr()->graph()->nNz() ;

    for ( auto iter = nNz.begin(); iter!=nNz.end(); ++iter )
        nnz += ( *iter ) ;

    LOG(INFO) << "[stokes] matrix NNZ "<< nnz << "\n";
    M_stats.put( "n.matrix.nnz",nnz );
    double u_errorL2 = normL2( _range=elements( mesh ), _expr=trans( idv( u )-u_exact ) );
    double u_exactL2 = normL2( _range=elements( mesh ), _expr=trans( u_exact )*( u_exact ) );
    LOG(INFO) << "||u_error_rel||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";
    LOG(INFO) << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";
    M_stats.put( "e.l2.u",math::sqrt( u_errorL2 ) );
    M_stats.put( "e.l2.urel",math::sqrt( u_errorL2/u_exactL2 ) );

    double u_errorsemiH1 = normL2( _range=elements( mesh ), _expr=trace( ( gradv( u )-gradu_exact ) ) );
    double u_exactsemiH1 = normL2( _range=elements( mesh ), _expr=trace( ( gradu_exact )*trans( gradu_exact ) ) );
    double u_error_rel_H1 = math::sqrt( ( u_errorL2+u_errorsemiH1 )/( u_exactL2+u_exactsemiH1 ) );
    double u_error_H1 = math::sqrt( u_errorL2+u_errorsemiH1 );
    LOG(INFO) << "||u_error||_1= " << u_error_H1 << "\n";
    LOG(INFO) << "||u_error_rel||_1= " << u_error_rel_H1 << "\n";
    M_stats.put( "e.h1.u",u_error_H1 );
    M_stats.put( "e.h1.urel",u_error_rel_H1 );

    double p_errorL2 = normL2( _range=elements( mesh ), _expr=( idv( p )-( p_exact-mean_pexact ) ) );
    double p_exactL2 = normL2( _range=elements( mesh ), _expr=( p_exact-mean_pexact ) );
    LOG(INFO) << "||p_error_rel||_2 = " << math::sqrt( p_errorL2/p_exactL2 ) << "\n";
    LOG(INFO) << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";
    M_stats.put( "e.l2.prel",math::sqrt( p_errorL2/p_exactL2 ) );
    M_stats.put( "e.l2.p",math::sqrt( p_errorL2 ) );

    double mean_div_u = mean( elements( mesh ), divv( u ) )(0,0);
    LOG(INFO) << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double mean_div_uexact = mean( elements( mesh ), divu_exact )(0,0);
    LOG(INFO) << "[stokes] mean_div(uexact)=" << mean_div_uexact << "\n";

    double div_u_error_L2 = normL2( elements( mesh ), divv( u ) );
    M_stats.put( "e.l2.div_u",math::sqrt( div_u_error_L2 ) );
    LOG(INFO) << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );

    this->exportResults( U, V );

} // Stokes::run


template<int Dim, typename BasisU, typename BasisP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, BasisU, BasisP, Entity>::exportResults( element_type& U, element_type& V )
{
    M_exp =  exporter( _mesh=U.functionSpace()->mesh() );
    if ( M_exp->doExport() )
    {
        M_exp->step( 0 )->setMesh( U.functionSpace()->mesh() );
        M_exp->step( 0 )->add( "u", U.template element<0>() );
        M_exp->step( 0 )->add( "p", U.template element<1>() );
        M_exp->step( 0 )->add( "u_exact", V.template element<0>() );
        M_exp->step( 0 )->add( "p_exact", V.template element<1>() );
        M_exp->save();
    }
} // Stokes::export
} // Feel

#endif // __FEELPP_BENCH_STOKES_HPP
