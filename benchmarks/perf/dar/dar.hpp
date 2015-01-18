/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-02-05

  Copyright (C) 2012 University Joseph Fourier

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
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * Generate the Ring mesh
 */
gmsh_ptrtype createRing( int Dim, int Order, double meshSize, std::string const& convex );

/**
 * DAR Solver using discontinous approximation spaces
 *
 * solve \f$ -\beta\cdot\nabla u + \mu u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma_{in}\f$
 */
template<int Dim,
         int Order,
         typename Cont,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class DAR
    :
public Simget
{
    typedef Simget super;
public:

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Entity<Dim,Order,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    template<typename Conti = Cont, template<uint16_type D> class FieldType = Scalar>
    struct space
    {
        typedef bases<Lagrange<Order, FieldType,Conti> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type> type;
    };

    /* export */
    typedef Exporter<mesh_type,Order> export_type;

    DAR( po::variables_map const& vm, AboutData const& ad )
        :
        super(),
        bcCoeff( doption(_name="bccoeff") ),
        geomap( ( GeomapStrategyType )ioption(_name="geomap") )
    {
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
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( f1_type& u,
                        f2_type& v,
                        f3_type& e );

private:


    double bcCoeff;

    GeomapStrategyType geomap;


}; // DAR

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
void
DAR<Dim, Order, Cont, Entity>::run()
{
    using namespace Feel::vf;

    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % Environment::about().appName()
                            % entity_type::name()
                            % Order
                            % meshSizeInit()
                          );
    value_type penalisation = option(_name="penal").template as<value_type>();
    int bctype = ioption(_name="bctype");

    double beta_x = option(_name="bx").template as<value_type>();
    double beta_y = option(_name="by").template as<value_type>();
    value_type mu = option(_name="mu").template as<value_type>();
    value_type stiff = option(_name="stiff").template as<value_type>();
    bool ring = boption(_name="ring");

    std::cout << "[DAR] hsize = " << meshSizeInit() << "\n";
    std::cout << "[DAR] bx = " << beta_x << "\n";
    std::cout << "[DAR] by = " << beta_y << "\n";
    std::cout << "[DAR] mu = " << mu << "\n";
    std::cout << "[DAR] stiff = " << stiff << "\n";
    std::cout << "[DAR] ring = " << ring << "\n";
    std::cout << "[DAR] bccoeff = " << bcCoeff << "\n";
    std::cout << "[DAR] bctype = " << bctype << "\n";

    auto exporter =  boost::shared_ptr<export_type>( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );

    boost::timer t;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh;

    if ( ring )
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=createRing( Dim,Order,meshSizeInit(),entity_type::name() ) );

    else
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=domain( _name=( boost::format( "hypercube-%1%" ) % Dim ).str() ,
                                             _order=Order,
                                             _shape="hypercube",
                                             _dim=Dim,
                                             _h=meshSize() ) );

    M_stats.put( "t.init.mesh",t.elapsed() );
    t.restart();

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space<Cont>::type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    M_stats.put( "h",M_meshSize );
    M_stats.put( "n.space.nelts",Xh->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nLocalDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    LOG(INFO) << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    auto r = sqrt( Px()*Px()+Py()*Py() );
    //auto r = sqrt(trans(P())*P());
    auto beta = beta_x*( ring*Py()/r + !ring )*oneX()+beta_y*( -ring*Px()/r + !ring )*oneY();
    //auto beta = (ones<Dim,1>());
    auto beta_N = ( trans( N() )*beta );
    auto beta_abs = abs( beta_N );
    auto beta_minus = constant( 0.5 )*( beta_abs-beta_N );
    auto beta_plus = constant( 0.5 )*( beta_abs+beta_N );
    auto g = ( ( !ring*exp( -mu*Px() )*atan( ( Py()-	0.5 )/stiff ) + ring*exp( -mu*r*acos( Py()/r ) )*atan( ( r-0.5 )/stiff ) ) );
    auto f = ( constant( 0.0 ) );

    auto F = backend()->newVector( Xh );
    M_stats.put( "t.init.vector",t.elapsed() );
    LOG(INFO) << "  -- time for vector init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;
    form1( _test=Xh, _vector=F, _init=true )  =
        integrate( elements( mesh ), trans( f )*id( v ) );

    if ( bctype == 1 || !Cont::is_continuous )
    {
        form1( _test=Xh, _vector=F ) +=
            integrate( _range=boundaryfaces( mesh ), _expr=trans( beta_minus*g )*id( v ),_geomap=geomap );
    }

    M_stats.put( "t.assembly.vector.total",t.elapsed() );
    LOG(INFO) << "  -- time vector global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;


    auto D = backend()->newMatrix( _test=Xh, _trial=Xh, _pattern=Pattern::COUPLED|Pattern::EXTENDED );
    M_stats.put( "t.init.matrix",t.elapsed() );
    LOG(INFO) << "  -- time for matrix init done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements( mesh ), _quad=_Q<2*Order>(),
                   // -(u,beta*grad(v))+(mu*u,v)-(u,div(beta)*v)
                   _expr=-trans( idt( u ) )*( grad( v )*beta ) + mu*trans( idt( u ) )*id( v ),
                   //- idt(u)*id(v)*(dx(beta_x)+dy(beta_y))
                   _geomap=geomap
                 );

    if ( !Cont::is_continuous )
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=internalfaces( mesh ), _quad=_Q<2*Order>(),
                       // {beta u} . [v]
                       //( trans(averaget(trans(beta)*idt(u))) * jump(trans(id(v))) )
                       _expr=( averaget( trans( beta )*idt( u ) ) * jump( id( v ) ) )
                             // penal*[u] . [v]
                             + penalisation*beta_abs*( trans( jumpt( trans( idt( u ) ) ) )*jump( trans( id( v ) ) ) ),
                       _geomap=geomap );
    }

    else // continuous case: stabilization by interior penalty
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( internalfaces( mesh ),
                       // penal*[grad(u)] . [grad(v)]
                       + penalisation*beta_abs*hFace()*hFace()*( trans( jumpt( gradt( u ) ) )*jump( grad( v ) ) ) );
    }

    if ( bctype == 1 || !Cont::is_continuous )
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=boundaryfaces( mesh ), _expr=beta_plus*trans( idt( u ) )*id( v ),_geomap=geomap );
    }

    else if ( bctype == 0 )
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            on( _range=boundaryfaces( mesh ), _element=u, _rhs=F, _expr=g );
    }

    M_stats.put( "t.assembly.matrix.total",t.elapsed() );
    LOG(INFO) << " -- time matrix global assembly done in "<<t.elapsed()<<" seconds \n";
    t.restart() ;

    if ( this->vm().count( "export-matlab" ) )
    {
        F->printMatlab( "F.m" );
        D->printMatlab( "D.m" );
    }

    backend(_rebuild=true )->solve( _matrix=D, _solution=u, _rhs=F );
    M_stats.put( "t.solver.total",t.elapsed() );
    LOG(INFO) << " -- time for solver : "<<t.elapsed()<<" seconds \n";
    t.restart();

    double c = integrate( internalfaces( mesh ), trans( jumpv( idv( u ) ) )*jumpv( idv( u ) )  ).evaluate()( 0, 0 );
    M_stats.put( "t.integrate.jump",t.elapsed() );
    t.restart();
    double error = integrate( elements( mesh ), trans( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0, 0 );
    M_stats.put( "t.integrate.l2norm",t.elapsed() );
    t.restart();

    LOG(INFO) << "||error||_0 =" << error << "\n";
    std::cout << "c =" << c << "\n";
    std::cout << "||error||_0 =" << error << "\n";
    M_stats.put( "e.l2.error",math::sqrt( error ) );


    //this->exportResults( u, g, beta );
    auto Xch = space<Continuous>::type::New( mesh );
    auto uEx = Xch->element();
    auto uC = Xch->element();
    auto Xvch = space<Continuous,Vectorial>::type::New( mesh );
    auto betaC = Xvch->element();

    auto L2Proj = projector( Xch, Xch );
    uEx = L2Proj->project( g );
    uC = L2Proj->project( vf::idv( u ) );
    auto L2Projv = projector( Xvch, Xvch );
    betaC = L2Projv->project( ( beta ) );

    exporter->step( 0 )->setMesh( u.functionSpace()->mesh() );
    exporter->step( 0 )->add( "u", u );
    exporter->step( 0 )->add( "uC", uC );
    exporter->step( 0 )->add( "exact", uEx );
    exporter->step( 0 )->add( "beta", betaC );
    exporter->save();
} // DAR::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename f1_type, typename f2_type, typename f3_type>
void
DAR<Dim, Order, Cont, Entity>::exportResults( f1_type& U,
        f2_type& E,
        f3_type& beta )
{
#if 0
    auto Xch = space<Continuous>::type::New( mesh );
    auto uEx = Xch->element();
    auto uC = Xch->element();
    auto Xvch = space<Continuous,Vectorial>::type::New( mesh );
    auto betaC = Xvch->element();

    auto L2Proj = projector( Xch, Xch );
    uEx = L2Proj->project( E );
    uC = L2Proj->project( vf::idv( U ) );
    auto L2Projv = projector( Xvch, Xvch );
    betaC = L2Projv->project( trans( beta ) );

    exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 0 )->add( "u", U );
    exporter->step( 0 )->add( "uC", uC );
    exporter->step( 0 )->add( "exact", uEx );
    exporter->step( 0 )->add( "beta", betaC );
    exporter->save();
#endif

} // DAR::export
} // Feel




