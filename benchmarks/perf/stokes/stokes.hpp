/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
    typedef bases<basis_u_type,basis_p_type> basis_type;
    //# endmarker1 #

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Stokes( std::string const& basis_name,
            po::variables_map const& vm, AboutData const& ad )
        :
        super( vm, ad ),
        M_backend(),
        M_basis_name( basis_name ),
        exporter()
    {
        mu = this->vm()["mu"].template as<value_type>();
        penalbc = this->vm()[prefixvm( name(),"bccoeff" )].template as<value_type>();
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

    boost::shared_ptr<export_type> exporter;
}; // Stokes


template<int Dim, typename BasisU, typename BasisP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, BasisU, BasisP, Entity>::run()
{
    using namespace Feel::vf;

    if ( this->vm().count( "nochdir" ) == false )
    {
        this->changeRepository( boost::format( "perf/%1%/%2%/%3%/h_%4%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % M_basis_name
                                % meshSize() );
    }

    //! init backend
    M_backend = backend_type::build( this->vm() );
    exporter =  boost::shared_ptr<export_type>( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );

    boost::timer t;
#if defined(KOVASZNAY)
    double xmin = -0.5, xmax=1.5;
    double ymin =  0, ymax=2;
#else
    double xmin = -1, xmax=1;
    double ymin = -1, ymax=1;
#endif
    double shear = this->vm()["shear"].template as<value_type>();
    bool recombine = this->vm()["recombine"].template as<bool>();
    /*
     * First we create the mesh, in the case of quads we wish to have
     * non-regular meshes to ensure that we don't have some super-convergence
     * associated to the mesh. Hence we use recombine=true to recombine
     * triangles generated by a Delaunay algorithm into quads
     */
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % Dim % 1 ).str() ,
                                        _shape="hypercube",
                                        _usenames=true,
                                        _convex=( ( !recombine )&&convex_type::is_hypercube )?"Hypercube":"Simplex",
                                        _recombine=( recombine&&convex_type::is_hypercube ), // generate quads which are not regular
                                        _dim=Dim,
                                        _h=M_meshSize,
                                        _shear=shear,
                                        _xmin=xmin,_xmax=xmax,
                                        _ymin=ymin,_ymax=ymax ) );

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

    Log() << "Data Summary:\n";
    Log() << "   hsize = " << M_meshSize << "\n";
    Log() << "  export = " << this->vm().count( "export" ) << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";
    Log() << "[mesh]   number of elements: " << Xh->mesh()->numElements() << "\n";
    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    M_stats.put( "h",M_meshSize );
    M_stats.put( "n.space.nelts",Xh->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nLocalDof() );
    M_stats.put( "n.space.ndof.u",Xh->template functionSpace<0>()->nDof() );
    M_stats.put( "n.space.ndof.p",Xh->template functionSpace<1>()->nDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    Log() << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
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

    double betacoeff = this->vm()["beta"].template as<value_type>();
    bool add_convection = ( math::abs( betacoeff  ) > 1e-10 );
    double mu = this->vm()["mu"].template as<value_type>();

#if FEELPP_SOLUTION_1
    // u exact solution
    //auto u_exact = vec(cos(Px())*cos(Py()), sin(Px())*sin(Py()));
    auto u_exact = val( -exp( Px() )*( Py()*cos( Py() )+sin( Py() ) )*unitX()+ exp( Px() )*Py()*sin( Py() )*unitY()+ Pz()*unitZ() );

    auto du_dx = val( -exp( Px() )*( Py()*cos( Py() )+sin( Py() ) ) );
    auto du_dy = val( - exp( Px() )*( 2*cos( Py() )-Py()*sin( Py() ) ) );
    auto dv_dx = val( exp( Px() )*Py()*sin( Py() ) );
    auto dv_dy = val( exp( Px() )*( sin( Py() )+Py()*cos( Py() ) ) );
    auto grad_exact = val( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = val( du_dx + dv_dy );
    auto beta = vec( cst( betacoeff ),cst( betacoeff ) );
    auto convection = val( grad_exact*beta );

    // this is the exact solution which has zero mean : the mean of
    // cos(x)*sin(y) is sin(1)*(1-cos(1))) on [0,1]^2
    //auto p_exact = cos(Px())*sin(Py())-(sin(1.0)*(1.-cos(1.0)));
    auto p_exact = val( 2.0*exp( Px() )*sin( Py() ) );

    // f is such that f = \Delta u_exact + \nabla p_exact
    //auto f = vec( (2*cos(Px())*cos(Py())-sin(Px())*sin(Py())),
    //              2*sin(Px())*sin(Py())+cos(Px())*cos(Py()) );
    auto f = val( -2*exp( Px() )*( mu-1. )*vec( sin( Py() ),cos( Py() ) ) );
#endif
#if FEELPP_SOLUTION_KOVASNAY
    //
    // the Kovasznay flow (2D)
    //
    double pi = M_PI;
    double lambda = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*pi*pi );
    // total stress tensor (test)
    auto u1 = 1. - exp( lambda * Px() ) * cos( 2.*pi*Py() );
    auto u2 = ( lambda/( 2.*pi ) ) * exp( lambda * Px() ) * sin( 2.*pi*Py() );
    auto u_exact = val( vec( u1,u2 ) );

    auto du_dx = ( -lambda*exp( lambda * Px() )*cos( 2.*pi*Py() ) );
    auto du_dy = ( 2*pi*exp( lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dx = ( ( lambda*lambda/( 2*pi ) )*exp( lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dy = ( lambda*exp( lambda * Px() )*cos( 2.*pi*Py() ) );
    auto grad_exact = val( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = val( du_dx + dv_dy );

    auto beta = vec( cst( betacoeff ),cst( betacoeff ) );
    auto convection = val( grad_exact*beta );

    auto p_exact = val( ( -exp( 2.*lambda*Px() ) )/2.0-0.125*( exp( -1.0*lambda )-1.0*exp( 3.0*lambda ) )/lambda );

    //auto f1 = (exp( lambda * Px() )*((lambda*lambda - 4.*pi*pi)*mu*cos(2.*pi*Py()) - lambda*exp( lambda * Px() )));
    auto f1 = ( -mu*( -lambda*lambda*exp( lambda*Px() )*cos( 2.0*pi*Py() )+4.0*exp( lambda*Px() )*cos( 2.0*pi*Py() )*pi*pi )-lambda*exp( 2.0*lambda*Px() ) );

    //auto f2 = (exp( lambda * Px() )*mu*(lambda/(2.*pi))*sin(2.*pi*Py())*(-lambda*lambda +4*pi*pi));
    auto f2 = ( -mu*( lambda*lambda*lambda*exp( lambda*Px() )*sin( 2.0*pi*Py() )/pi/2.0-2.0*lambda*exp( lambda*Px() )*sin( 2.0*pi*Py() )*pi ) );

    auto f = val( vec( f1,f2 ) ); //+ convection;

    //double pmean = integrate( elements(mesh), p_exact ).evaluate()( 0, 0 )/mesh->measure();
    double pmean = -0.125*( math::exp( -1.0*lambda )-1.0*math::exp( 3.0*lambda ) )/lambda;

#endif
#if FEELPP_SOLUTION_ETHIERSTEINMANN
    //
    // the EthierSteinmann flow (3D)
    //
    double pi = M_PI;
    double lambda = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*pi*pi );
    // total stress tensor (test)
    auto u1 = 1. - exp( lambda * Px() ) * cos( 2.*pi*Py() );
    auto u2 = ( lambda/( 2.*pi ) ) * exp( lambda * Px() ) * sin( 2.*pi*Py() );
    auto u_exact = val( vec( u1,u2 ) );

    auto du_dx = ( -lambda*exp( lambda * Px() )*cos( 2.*pi*Py() ) );
    auto du_dy = ( 2*pi*exp( lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dx = ( ( lambda*lambda/( 2*pi ) )*exp( lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dy = ( lambda*exp( lambda * Px() )*cos( 2.*pi*Py() ) );
    auto grad_exact = val( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = val( du_dx + dv_dy );

    auto beta = vec( cst( betacoeff ),cst( betacoeff ) );
    auto convection = val( grad_exact*beta );

    auto p_exact = val( ( -exp( 2.*lambda*Px() ) )/2.0-0.125*( exp( -1.0*lambda )-1.0*exp( 3.0*lambda ) )/lambda );

    //auto f1 = (exp( lambda * Px() )*((lambda*lambda - 4.*pi*pi)*mu*cos(2.*pi*Py()) - lambda*exp( lambda * Px() )));
    auto f1 = ( -mu*( -lambda*lambda*exp( lambda*Px() )*cos( 2.0*pi*Py() )+4.0*exp( lambda*Px() )*cos( 2.0*pi*Py() )*pi*pi )-lambda*exp( 2.0*lambda*Px() ) );

    //auto f2 = (exp( lambda * Px() )*mu*(lambda/(2.*pi))*sin(2.*pi*Py())*(-lambda*lambda +4*pi*pi));
    auto f2 = ( -mu*( lambda*lambda*lambda*exp( lambda*Px() )*sin( 2.0*pi*Py() )/pi/2.0-2.0*lambda*exp( lambda*Px() )*sin( 2.0*pi*Py() )*pi ) );

    auto f = val( vec( f1,f2 ) ); //+ convection;

    //double pmean = integrate( elements(mesh), p_exact ).evaluate()( 0, 0 )/mesh->measure();
    double pmean = -0.125*( math::exp( -1.0*lambda )-1.0*math::exp( 3.0*lambda ) )/lambda;

#endif

    boost::timer subt;
    // right hand side
    auto F = M_backend->newVector( Xh );
    form1( Xh, F, _init=true );
    M_stats.put( "t.init.vector",t.elapsed() );
    Log() << "  -- time for vector init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;
    form1( Xh, F ) = integrate( elements( mesh ), trans( f )*id( v ) );
    M_stats.put( "t.assembly.vector.source",subt.elapsed() );
    subt.restart();

    if ( add_convection )
    {
        form1( Xh, F ) += integrate( elements( mesh ), trans( convection )*id( v ) );
        M_stats.put( "t.assembly.vector.convection",subt.elapsed() );
        subt.restart();
    }


    if ( this->vm()[ prefixvm( name(),"bctype" ) ].template as<int>() == 1  )
    {
        form1( Xh, F ) += integrate( _range=boundaryfaces( mesh ), _expr=-trans( u_exact )*SigmaN );
        M_stats.put( "t.assembly.vector.dirichletup",subt.elapsed() );
        subt.restart();
        form1( Xh, F ) += integrate( _range=boundaryfaces( mesh ), _expr=penalbc*inner( u_exact,id( v ) )/hFace() );
        //form1( Xh, F ) += integrate( _range=boundaryfaces(mesh), _expr=penalbc*max(betacoeff,mu/hFace())*(trans(id(v))*N())*N());
        M_stats.put( "t.assembly.vector.dirichletp",subt.elapsed() );
        Log() << "   o time for rhs weak dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    M_stats.put( "t.assembly.vector.total",t.elapsed() );
    Log() << "  -- time vector global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    /*
     * Construction of the left hand side
     */
    size_type pattern = Pattern::COUPLED;
    size_type patternsym = Pattern::COUPLED;

    if ( this->vm()[ "faster" ].template as<int>() == 1 )
    {
        pattern = Pattern::COUPLED;
        patternsym = Pattern::COUPLED|Pattern::SYMMETRIC;
    }

    if ( this->vm()[ "faster" ].template as<int>() == 2 )
    {
        pattern = Pattern::DEFAULT;
        patternsym = Pattern::DEFAULT;
    }

    if ( this->vm()[ "faster" ].template as<int>() == 3 )
    {
        pattern = Pattern::DEFAULT;
        patternsym = Pattern::DEFAULT|Pattern::SYMMETRIC;
    }

    //# marker7 #
    t.restart();
    auto D = M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, D, _init=true );
    M_stats.put( "t.init.matrix",t.elapsed() );
    Log() << "  -- time for matrix init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    subt.restart();
#if 0

    if (  this->vm().count( "extra-terms" ) )
    {
        //form2( Xh, Xh, D ) =integrate( _range=elements(mesh),_expr=mu*(inner(dxt(u),dx(v))+inner(dyt(u),dy(v))));
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=trans( gradt( u )*idv( u ) )*id( v ) );
        M_stats.put( "t.assembly.matrix.convection1",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D, _pattern=pattern ) =integrate( _range=elements( mesh ),_expr=trans( gradt( u )*idv( u ) )*id( v ) );
        M_stats.put( "t.assembly.matrix.convection11",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=inner( gradt( u )*idv( u ),id( v ) ) );
        M_stats.put( "t.assembly.matrix.convection2",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=trans( dxt( u )*idv( u )( 0 )+dyt( u )*idv( u )( 1 )+dzt( u )*idv( u )( 2 ) )*id( v ) );
        M_stats.put( "t.assembly.matrix.convection3",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=mu*( trans( dxt( u ) )*dx( v )+trans( dyt( u ) )*dy( v )+trans( dzt( u ) )*dz( v ) ) );
        M_stats.put( "t.assembly.matrix.diffusion0",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=mu*trace( trans( gradt( u ) )*grad( v ) ) );
        M_stats.put( "t.assembly.matrix.diffusion1",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=inner( mu*gradt( u ),grad( v ) ) );
        M_stats.put( "t.assembly.matrix.diffusion2",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=inner( mu*dxt( u ),dx( v ) )+inner( mu*dyt( u ),dy( v ) )+inner( mu*dzt( u ),dz( v ) ) );
        M_stats.put( "t.assembly.matrix.diffusion3",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D ) =integrate( _range=elements( mesh ),_expr=mu*( inner( dxt( u ),dx( v ) )+inner( dyt( u ),dy( v ) )+inner( dzt( u ),dz( v ) ) ) );
        M_stats.put( "t.assembly.matrix.diffusion4",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, D, _pattern=pattern ) =integrate( _range=elements( mesh ),_expr=mu*( inner( gradt( u ),grad( v ) ) ) );
        M_stats.put( "t.assembly.matrix.diffusionfast",subt.elapsed() );
        subt.restart();

    }

#endif
    t.restart();

    form2( Xh, Xh, D, _pattern=patternsym ) =integrate( _range=elements( mesh ),_expr=mu*( inner( gradt( u ),grad( v ) ) ) );
    M_stats.put( "t.assembly.matrix.diffusion",subt.elapsed() );
    Log() << "   o time for diffusion terms: " << subt.elapsed() << "\n";
    subt.restart();

    if ( add_convection )
    {
        form2( Xh, Xh, D, _pattern=pattern ) +=integrate( _range=elements( mesh ),_expr=trans( gradt( u )*beta )*id( v ) );
        M_stats.put( "t.assembly.matrix.convection",subt.elapsed() );
        Log() << "   o time for convection terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    form2( Xh, Xh, D )+=integrate( _range=elements( mesh ),_expr=-div( v )*idt( p ) );
    form2( Xh, Xh, D )+=integrate( _range=elements( mesh ),_expr=divt( u )*id( q ) );
    M_stats.put( "t.assembly.matrix.up",subt.elapsed() );
    Log() << "   o time for velocity/pressure terms: " << subt.elapsed() << "\n";
    subt.restart();

    if ( this->vm()[ prefixvm( name(),"bctype" ) ].template as<int>() == 1  )
    {
        if (  !this->vm().count( "extra-terms" ) )
        {
            //form2( Xh, Xh, D, _pattern=patternsym )+=integrate( _range=boundaryfaces(mesh),_expr=-trans(SigmaNt)*id(v) );
            form2( Xh, Xh, D, _pattern=patternsym )+=integrate( _range=boundaryfaces( mesh ),_expr=-trans( SigmaNt )*id( v ) );
            M_stats.put( "t.assembly.matrix.dirichlet1",subt.elapsed() );
            subt.restart();
            form2( Xh, Xh, D, _pattern=patternsym )+=integrate( _range=boundaryfaces( mesh ),_expr=-trans( SigmaN )*idt( v ) );
            M_stats.put( "t.assembly.matrix.dirichlet2",subt.elapsed() );
            subt.restart();
        }

        else
        {
#if 0
            ProfilerStart( "/tmp/bfaces" );
            form2( Xh, Xh, D )+=integrate( _range=boundaryfaces( mesh ),_expr=idt( p )*( trans( N() )*id( v ) ) );
            ProfilerStop();
            M_stats.put( "t.assembly.matrix.dirichlet_pn*v",subt.elapsed() );
            subt.restart();
            form2( Xh, Xh, D )+=integrate( _range=boundaryfaces( mesh ),_expr=-trans( dnt( u ) )*id( v ) );
            M_stats.put( "t.assembly.matrix.dirichlet_dnt(u)*v",subt.elapsed() );
            subt.restart();
            //form2( Xh, Xh, D )+=integrate( _range=boundaryfaces(mesh),_expr=-trans(SigmaN)*idt(u) );
            form2( Xh, Xh, D )+=integrate( _range=boundaryfaces( mesh ),_expr=id( p )*( trans( N() )*idt( u ) ) );
            //form2( Xh, Xh, D )+=integrate( _range=boundaryfaces(mesh),_expr=(idt(u)(0)*Nx()+idt(u)(1)*Ny()+idt(u)(2)*Nz())*id(p), _verbose=true );
            M_stats.put( "t.assembly.matrix.dirichlet_pN*u",subt.elapsed() );
            subt.restart();
            form2( Xh, Xh, D )+=integrate( _range=boundaryfaces( mesh ),_expr=-trans( dn( u ) )*idt( u ) );
            M_stats.put( "t.assembly.matrix.dirichlet_dn(v)*u",subt.elapsed() );
            subt.restart();
#endif
        }

        form2( Xh, Xh, D, _pattern=patternsym )+=integrate( _range=boundaryfaces( mesh ),_expr=+penalbc*inner( idt( u ),id( v ) )/hFace() );
        //form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), +penalbc*(trans(idt(u))*N())*(trans(id(v))*N())*max(betacoeff,mu/hFace()) );
        M_stats.put( "t.assembly.matrix.dirichlet_u*u",subt.elapsed() );
        subt.restart();
        Log() << "   o time for weak dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    //# endmarker7 #
    D->close();
    F->close();
    M_stats.put( "t.assembly.matrix.total",t.elapsed() );
    Log() << " -- time matrix global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    if ( this->vm()[ prefixvm( name(),"bctype" ) ].template as<int>() == 0  )
    {
        form2( Xh, Xh, D ) += on( _range=boundaryfaces( mesh ), _element=u, _rhs=F, _expr=u_exact );
        M_stats.put( "t.assembly.matrix.dirichlet",subt.elapsed() );
        Log() << "   o time for strong dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }



    t.restart();

    if ( !this->vm().count( "no-solve" ) )
    {
        auto r = M_backend->solve( _matrix=D, _solution=U, _rhs=F, _constant_null_space=true );
        M_stats.put( "d.solver.bool.converged",r.template get<0>() );
        M_stats.put( "d.solver.int.nit",r.template get<1>() );
        M_stats.put( "d.solver.double.residual",r.template get<2>() );
    }

    M_stats.put( "t.solver.total",t.elapsed() );
    Log() << " -- time for solver : "<<t.elapsed()<<" seconds \n";


    double meas = integrate( _range=elements( mesh ), _expr=constant( 1.0 ) ).evaluate()( 0, 0 );
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 4)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 4)\n";

    double mean_p = integrate( elements( mesh ), idv( p ) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    // get the zero mean pressure
    p.add( - mean_p );
    mean_p = integrate( elements( mesh ), idv( p ) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p-mean(p))=" << mean_p << "\n";
    std::cout << "[stokes] mean(p-mean(p))=" << mean_p << "\n";
    double mean_pexact = integrate( elements( mesh ), p_exact ).evaluate()( 0, 0 )/meas;
    std::cout << "[stokes] mean(pexact)=" << mean_pexact << "\n";
    size_type nnz = 0 ;
    auto nNz = D->graph()->nNz() ;

    for ( auto iter = nNz.begin(); iter!=nNz.end(); ++iter )
        nnz += ( *iter ) ;

    Log() << "[stokes] matrix NNZ "<< nnz << "\n";
    M_stats.put( "n.matrix.nnz",nnz );
    double u_errorL2 = integrate( _range=elements( mesh ), _expr=trans( idv( u )-u_exact )*( idv( u )-u_exact ) ).evaluate()( 0, 0 );
    double u_exactL2 = integrate( _range=elements( mesh ), _expr=trans( u_exact )*( u_exact ) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";;
    Log() << "||u_error||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";;
    M_stats.put( "e.l2.u",math::sqrt( u_errorL2/u_exactL2 ) );

    double u_errorsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=trace( ( gradv( u )-grad_exact )*trans( gradv( u )-grad_exact ) ) ).evaluate()( 0, 0 );
    double u_exactsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=trace( ( grad_exact )*trans( grad_exact ) ) ).evaluate()( 0, 0 );
    double u_error_H1 = math::sqrt( ( u_errorL2+u_errorsemiH1 )/( u_exactL2+u_exactsemiH1 ) );
    std::cout << "||u_error||_1= " << u_error_H1 << "\n";
    M_stats.put( "e.h1.u",u_error_H1 );

    double p_errorL2 = integrate( _range=elements( mesh ), _expr=( idv( p )-( p_exact-mean_pexact ) )*( idv( p )-( p_exact-mean_pexact ) ) ).evaluate()( 0, 0 );
    double p_exactL2 = integrate( _range=elements( mesh ), _expr=( p_exact-mean_pexact )*( p_exact-mean_pexact ) ).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2/p_exactL2 ) << "\n";;
    Log() << "||p_error||_2 = " << math::sqrt( p_errorL2/p_exactL2 ) << "\n";;
    M_stats.put( "e.l2.p",math::sqrt( p_errorL2/p_exactL2 ) );
    Log() << "[stokes] solve for D done\n";

    double mean_div_u = integrate( elements( mesh ), divv( u ) ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double mean_div_uexact = integrate( elements( mesh ), div_exact ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(uexact)=" << mean_div_uexact << "\n";
    std::cout << "[stokes] mean_div(uexact)=" << mean_div_uexact << "\n";

    double div_u_error_L2 = integrate( elements( mesh ), divv( u )*divv( u ) ).evaluate()( 0, 0 );
    M_stats.put( "e.l2.div_u",math::sqrt( div_u_error_L2 ) );
    Log() << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";
    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );

    this->exportResults( U, V );

} // Stokes::run


template<int Dim, typename BasisU, typename BasisP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, BasisU, BasisP, Entity>::exportResults( element_type& U, element_type& V )
{
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", U.template element<0>() );
        exporter->step( 0 )->add( "p", U.template element<1>() );
        exporter->step( 0 )->add( "u_exact", V.template element<0>() );
        exporter->step( 0 )->add( "p_exact", V.template element<1>() );
        exporter->save();
    }
} // Stokes::export
} // Feel

#endif // __FEELPP_BENCH_STOKES_HPP

