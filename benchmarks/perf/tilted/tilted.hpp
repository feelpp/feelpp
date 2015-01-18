/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009-2012 Université Joseph Fourier (Grenoble I)

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
   \file tilted.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#if !defined( __FEELPP_BENCH_TILTED_HPP)
#define __FEELPP_BENCH_TILTED_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/simget.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/projector.hpp>


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * \class Tilted class
 * \brief solves the tilted equations
 *
 */
template<int Dim,
         typename BasisU,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Tilted
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

    typedef FunctionSpace<mesh_type, Lagrange<0,Scalar,Discontinuous> > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Tilted( std::string const& basis_name,
               po::variables_map const& vm, AboutData const& ad )
        :
        super(),
        M_backend(),
        M_basis_name( basis_name ),
        exporter()
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
    void exportResults( element_type& u, element_type& v, p0_element_type& k, element_type& vbdy );

private:

    backend_ptrtype M_backend;
    std::string M_basis_name;
    double mu;
    double penalbc;

    boost::shared_ptr<export_type> exporter;
}; // Tilted


template<int Dim, typename BasisU, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Tilted<Dim, BasisU, Entity>::run()
{
    using namespace Feel::vf;

    if ( Environment::vm().count( "nochdir" ) == false )
    {
        this->changeRepository( boost::format( "perf/%1%/%2%/%3%/h_%4%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % M_basis_name
                                % meshSize() );
    }

    //! init backend
    M_backend = backend_type::build( soption("backend"));
    exporter =  boost::shared_ptr<export_type>( Exporter<mesh_type>::New( Environment::vm(),
                                                                          (boost::format( "%1%-%2%" ) % this->about().appName() % this->level() ).str() ) );

    boost::timer t;
    double shear = option(_name="shear").template as<value_type>();
    bool recombine = boption(_name="recombine");
    /*
     * First we create the mesh, in the case of quads we wish to have
     * non-regular meshes to ensure that we don't have some super-convergence
     * associated to the mesh. Hence we use recombine=true to recombine
     * triangles generated by a Delaunay algorithm into quads
     */
#if 1
    auto mesh = loadMesh( _mesh=new mesh_type, _filename=(boost::format("tilted%1%.msh") % this->level()).str() );
#else
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % 2 ).str() ,
                                              _usenames=true,
                                              _shape="hypercube",
                                              _convex="hypercube",
                                              _h=this->meshSize(),
                                              _xmin=0,
                                              _ymin=0 ),
                                _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                _partitions=this->comm().size()  );
#endif
    M_stats.put( "t.init.mesh",t.elapsed() );
    t.restart();
    /*
     * The function space and some associate elements are then defined
     */
    //# marker4 #

    space_ptrtype Xh = space_type::New( mesh );
    p0_space_ptrtype P0h = p0_space_type::New( mesh );

    auto u = Xh->element();
    auto v = Xh->element();

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << M_meshSize << "\n";
    LOG(INFO) << "  export = " << Environment::vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";
    LOG(INFO) << "[mesh]   number of elements: " << Xh->mesh()->numElements() << "\n";
    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";

    M_stats.put( "h",M_meshSize );
    M_stats.put( "n.space.nelts",Xh->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nLocalDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    double penalbc = option(_name="bccoeff").template as<value_type>();
    double mu = option(_name="mu").template as<value_type>();

    double kappa = option(_name="kappa").template as<value_type>();
    double pi = constants::pi<double>();

#if defined(FEELPP_SOLUTION_1)
    double alpha=(3./pi)*math::atan(math::sqrt(1.+2./kappa));
    double beta=math::cos(alpha*pi/3.)/math::cos(2*alpha*pi/3.);

    auto r = sqrt((Px()-0.5)*(Px()-0.5)+(Py()-0.5)*(Py()-0.5));
    auto theta1=acos( (Px()-0.5)/r );
    auto theta =
        chi((Py()-0.5)>=0.)*theta1 +
        chi((Py()-0.5)<0)*(2*pi-theta1);

    auto ralpha=pow(r,alpha);
    auto u_exact =
        chi(0<=theta && theta<2.*pi/3.)*ralpha*cos(alpha*(theta-pi/3.)) +
        chi(2.*pi/3.<=theta && theta<=2.*pi)*ralpha*beta*cos(alpha*(4.*pi/3.-theta));
    auto f = cst(0.);
    auto k =
        chi(0<=theta && theta<2.*pi/3.)*cst(kappa) + // k1
        chi(2.*pi/3.<=theta && theta<=2.*pi)* cst(1.0); // k2

    auto ralphamo = pow(r,alpha-1);
    auto u_r=trans(vec((Px()-0.5)/r,(Py()-0.5)/r));
    auto u_theta=trans(vec(-(Py()-0.5)/r,(Px()-0.5)/r));
    auto grad_u_exact =
        chi(0<=theta && theta<2.*pi/3.)*alpha*ralphamo*(-sin(alpha*(theta-pi/3))*u_theta
                                                        +cos(alpha*(theta-pi/3))*u_r)+
        chi(2.*pi/3.<=theta && theta<=2.*pi)*alpha*ralphamo*beta*(sin(alpha*(4*pi/3-theta))*u_theta
                                                                  +cos(alpha*(4*pi/3-theta))*u_r);

#endif
#if defined(FEELPP_SOLUTION_2)
    double alpha=(3./pi)*math::atan(math::sqrt(1.+2.*kappa));
    double beta=1./(2.*math::cos(alpha*pi/3.));

    auto r = sqrt((Px()-0.5)*(Px()-0.5)+(Py()-0.5)*(Py()-0.5));
    auto theta1=acos( (Px()-0.5)/r );
    auto theta =
        chi((Py()-0.5)>=0.)*theta1 +
        chi((Py()-0.5)<0)*(2*pi-theta1);

    auto ralpha=pow(r,alpha);
    auto u_exact =
        chi(0<=theta && theta<2.*pi/3.)*ralpha*sin(alpha*(theta-pi/3.)) +
        chi(2.*pi/3.<=theta && theta<=2.*pi)*ralpha*beta*sin(alpha*(4.*pi/3.-theta));
    auto f = cst(0.);
    auto k =
        chi(0<=theta && theta<2.*pi/3.)*cst(kappa) + // k1
        chi(2.*pi/3.<=theta && theta<=2.*pi)* cst(1.0); // k2

    auto ralphamo = pow(r,alpha-1);
    auto u_r=trans(vec((Px()-0.5)/r,(Py()-0.5)/r));
    auto u_theta=trans(vec(-(Py()-0.5)/r,(Px()-0.5)/r));
    auto grad_u_exact =
        chi(0<=theta && theta<2.*pi/3.)*alpha*ralphamo*(cos(alpha*(theta-pi/3))*u_theta
                                                        +sin(alpha*(theta-pi/3))*u_r)+
        chi(2.*pi/3.<=theta && theta<=2.*pi)*alpha*ralphamo*beta*(-cos(alpha*(4*pi/3-theta))*u_theta
                                                                  +sin(alpha*(4*pi/3-theta))*u_r);
#endif

#if defined(FEELPP_SOLUTION_3)
    double alpha=(6./pi)*math::atan(1./math::sqrt(1.+2.*kappa));
    double beta=math::cos(alpha*pi/3.)*(math::sin(alpha*pi/6.));

    auto r = sqrt((Px()-0.5)*(Px()-0.5)+(Py()-0.5)*(Py()-0.5));
    auto theta1=acos( (Px()-0.5)/r );
    auto theta =
        chi((Py()-0.5)>=0.)*theta1 +
        chi((Py()-0.5)<0)*(2*pi-theta1);

    auto ralpha=pow(r,alpha);
    auto u_exact =
        chi((0<=theta) && (theta<2.*pi/3.))*ralpha*cos(alpha*(theta-pi/3.)) +
        chi((2.*pi/3.<=theta) && (theta<=2.*pi))*ralpha*beta*sin(alpha*(5.*pi/6.-theta));
    auto f = cst(0.);
    auto k =
        chi((0<=theta) && (theta<2.*pi/3.))*cst(kappa) + // k1
        chi((2.*pi/3.<=theta) && (theta<=2.*pi))* cst(1.0); // k2

#endif

    auto proj = opProjection(_domainSpace=Xh, _imageSpace=Xh);
    auto projP0 = opProjection(_domainSpace=P0h, _imageSpace=P0h);
    v =  proj->project( _expr=u_exact, _quad=_Q<5>() );
    auto tke = projP0->project( _expr=k, _quad=_Q<0>() );
    auto ke = vf::project( _space=P0h, _expr=chi((idv(tke) >= 0.55) )*1+chi(idv(tke) <0.45)*0.1 );
    auto vbdy = proj->project( _expr=u_exact, _range=boundaryfaces(mesh), _quad=_Q<5>() );

    {
        double u_errorL2 = integrate( _range=elements( mesh ), _expr=(idv( v )-u_exact )*( idv( v )-u_exact), _quad=_Q<10>() ).evaluate()( 0, 0 );
        double u_exactL2 = integrate( _range=elements( mesh ), _expr=u_exact*u_exact, _quad=_Q<10>() ).evaluate()( 0, 0 );
        std::cout << "||v_error||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";;
        M_stats.put( "e.l2.v",math::sqrt( u_errorL2/u_exactL2 ) );
    }

    boost::timer subt;
    // right hand side
    auto F = M_backend->newVector( Xh );
    form1( Xh, _vector=F, _init=true );
    M_stats.put( "t.init.vector",t.elapsed() );
    LOG(INFO) << "  -- time for vector init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;
    form1( Xh, _vector=F ) = integrate( elements( mesh ), f*id( v ) );
    M_stats.put( "t.assembly.vector.source",subt.elapsed() );
    subt.restart();

    if ( ioption(_name="bctype") == 1  )
    {
        form1( Xh, _vector=F ) += integrate( _range=boundaryfaces( mesh ), _expr=-idv(v)*k*grad( u )*N() );
        M_stats.put( "t.assembly.vector.dirichlet1",subt.elapsed() );
        subt.restart();
        form1( Xh, _vector=F ) += integrate( _range=boundaryfaces( mesh ), _expr=k*penalbc*idv(v)*id( v )/hFace() );
        //form1( Xh, _vector=F ) += integrate( _range=boundaryfaces(mesh), _expr=penalbc*max(betacoeff,mu/hFace())*(trans(id(v))*N())*N());
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
    size_type pattern = Pattern::COUPLED;
    size_type patternsym = Pattern::COUPLED;

    if ( ioption(_name="faster") == 1 )
    {
        pattern = Pattern::COUPLED;
        patternsym = Pattern::COUPLED|Pattern::SYMMETRIC;
    }

    if ( ioption(_name="faster") == 2 )
    {
        pattern = Pattern::DEFAULT;
        patternsym = Pattern::DEFAULT;
    }

    if ( ioption(_name="faster") == 3 )
    {
        pattern = Pattern::DEFAULT;
        patternsym = Pattern::DEFAULT|Pattern::SYMMETRIC;
    }

    //# marker7 #
    t.restart();
    auto D = M_backend->newMatrix( Xh, Xh );
    M_stats.put( "t.init.matrix",t.elapsed() );
    LOG(INFO) << "  -- time for matrix init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    subt.restart();
    t.restart();

    //form2( _trial=Xh, _test=Xh, _matrix=D ) =integrate( _range=elements( mesh ),_expr=k*( gradt( u )*trans(grad( v ) ) ) );
    form2( _trial=Xh, _test=Xh, _matrix=D ) =integrate( _range=elements( mesh ),_expr=k*( gradt( u )*trans(grad( v ) ) ), _quad=_Q<4>() );
    M_stats.put( "t.assembly.matrix.diffusion",subt.elapsed() );
    LOG(INFO) << "   o time for diffusion terms: " << subt.elapsed() << "\n";
    subt.restart();

    if ( ioption(_name="bctype") == 1  )
    {
        form2( Xh, Xh, _matrix=D )+=integrate( _range=boundaryfaces( mesh ),_expr=-k*( gradt( u )*N() )*id( v )-k*trans( grad( u )*N() )*idt( v ) );
        M_stats.put( "t.assembly.matrix.dirichlet1",subt.elapsed() );
        subt.restart();
        form2( Xh, Xh, _matrix=D )+=integrate( _range=boundaryfaces( mesh ),_expr=+k*penalbc*idt( u )*id( v )/hFace() );
        M_stats.put( "t.assembly.matrix.dirichlet2",subt.elapsed() );
        subt.restart();
        LOG(INFO) << "   o time for weak dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    //# endmarker7 #
    D->close();
    F->close();

    if ( ioption(_name="bctype") == 0  )
    {
        form2( Xh, Xh, _matrix=D ) += on( _range=boundaryfaces( mesh ), _element=u, _rhs=F, _expr=u_exact);
        M_stats.put( "t.assembly.matrix.dirichlet",subt.elapsed() );
        LOG(INFO) << "   o time for strong dirichlet terms: " << subt.elapsed() << "\n";
        subt.restart();
    }

    M_stats.put( "t.assembly.matrix.total",t.elapsed() );
    LOG(INFO) << " -- time matrix global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    t.restart();

    if ( !Environment::vm().count( "no-solve" ) )
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
    LOG(INFO) << "[tilted] measure(Omega)=" << meas << " (should be equal to 4)\n";
    std::cout << "[tilted] measure(Omega)=" << meas << " (should be equal to 4)\n";

    double mean_u = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 )/meas;
    LOG(INFO) << "[tilted] mean(u)=" << mean_u << "\n";
    std::cout << "[tilted] mean(u)=" << mean_u << "\n";

    // get the zero mean pressure
    //u.add( - mean_u );
    mean_u = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 )/meas;
    LOG(INFO) << "[tilted] mean(u-mean(u))=" << mean_u << "\n";
    std::cout << "[tilted] mean(u-mean(u))=" << mean_u << "\n";
    double mean_uexact = integrate( elements( mesh ), idv(v) ).evaluate()( 0, 0 )/meas;
    std::cout << "[tilted] mean(uexact)=" << mean_uexact << "\n";
    M_stats.put( "t.integrate.mean",t.elapsed() );
    size_type nnz = 0 ;
    auto nNz = D->graph()->nNz() ;

    for ( auto iter = nNz.begin(); iter!=nNz.end(); ++iter )
        nnz += ( *iter ) ;

    LOG(INFO) << "[tilted] matrix NNZ "<< nnz << "\n";
    M_stats.put( "n.matrix.nnz",nnz );

    t.restart();
    double u_errorL2 = integrate( _range=elements( mesh ), _expr=(idv( u )-u_exact )*( idv( u )-u_exact), _quad=_Q<10>() ).evaluate()( 0, 0 );
    double u_exactL2 = integrate( _range=elements( mesh ), _expr=u_exact*u_exact, _quad=_Q<10>() ).evaluate()( 0, 0 );
    M_stats.put( "t.integrate.l2norm",t.elapsed() );
    t.restart();
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;
    //std::cout << "||u_exact||_2 = " << math::sqrt( u_exactL2 ) << "\n";;
    //std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2/u_exactL2 ) << "\n";;
    LOG(INFO) << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;
    M_stats.put( "e.l2.u",math::sqrt( u_errorL2 ) );

#if defined(FEELPP_SOLUTION_1) || defined(FEELPP_SOLUTION_2)
    double u_errorsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=( gradv( u )-grad_u_exact )*trans( gradv( u )-grad_u_exact ), _quad=_Q<10>() ).evaluate()( 0, 0 );
    double u_exactsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=( grad_u_exact )*trans( grad_u_exact ), _quad=_Q<10>() ).evaluate()( 0, 0 );
    //double u_error_H1 = math::sqrt( ( u_errorL2+u_errorsemiH1 )/( u_exactL2+u_exactsemiH1 ) );
    double u_error_H1 = math::sqrt( ( u_errorL2+u_errorsemiH1 ));

    M_stats.put( "t.integrate.h1norm",t.elapsed() );
    t.restart();
    std::cout << "||u_error||_1= " << u_error_H1 << "\n";
    M_stats.put( "e.h1.u",u_error_H1 );
#endif

    M_stats.put( "t.export.projection",t.elapsed() );
    t.restart();
    this->exportResults( u, v, ke, vbdy );
    M_stats.put( "t.export.total",t.elapsed() );
    t.restart();

} // Tilted::run


template<int Dim, typename BasisU, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Tilted<Dim, BasisU, Entity>::exportResults( element_type& u, element_type& v, p0_element_type& k, element_type& vbdy )
{
    if ( exporter->doExport() )
    {
#if defined(FEELPP_SOLUTION_1)
        exporter->addPath( boost::format("heterogeneous1/%1%") % u.functionSpace()->basisName()  );
#endif /* FEELPP_SOLUTION_1 */
#if defined(FEELPP_SOLUTION_2)
        exporter->addPath( boost::format("heterogeneous2/%1%") % u.functionSpace()->basisName()  );
#endif /* FEELPP_SOLUTION_1 */
#if defined(FEELPP_SOLUTION_3)
        exporter->addPath( boost::format("heterogeneous3/%1%") % u.functionSpace()->basisName()  );
#endif /* FEELPP_SOLUTION_1 */

        exporter->step( 0 )->setMesh( u.functionSpace()->mesh() );
        exporter->step( 0 )->add( (boost::format("u_%1%")%this->level()).str(), u );
        exporter->step( 0 )->add( (boost::format("u_exact_%1%")%this->level()).str(), v );
        exporter->step( 0 )->add( (boost::format("k_%1%")%this->level()).str(), (boost::format("k_%1%")%this->level()).str(), k, mpl::bool_<false>(), mpl::bool_<false>() );
        exporter->step( 0 )->add( (boost::format("u_exact_bdy_%1%")%this->level()).str(), vbdy );
        exporter->save();
    }
} // Tilted::export
} // Feel

#endif // __FEELPP_BENCH_TILTED_HPP

