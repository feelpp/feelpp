/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)

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
#include <feel/feel.hpp>




namespace Feel
{

/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
template<int nDim, int uOrder, int geoOrder=1>
class Stokes
    :
public Simget
{
    typedef Simget super;
public:


    typedef double value_type;

    /*mesh*/
    typedef Mesh<Simplex<nDim,geoOrder> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef Mesh<Simplex<nDim-1,geoOrder,nDim>> mesh_lag_type;
    typedef boost::shared_ptr<mesh_lag_type> mesh_lag_ptrtype;

    /*basis*/
    typedef Lagrange<uOrder, Vectorial> basis_u_type;
    typedef Lagrange<uOrder-1, Scalar> basis_p_type;

    typedef Lagrange<uOrder, Scalar> basis_l_type;

    // use lagrange multipliers to ensure zero mean pressure
#if defined( FEELPP_USE_LM )
#if defined (DIM2)
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
#else
    typedef bases<basis_u_type,basis_p_type, basis_l_type, basis_l_type> basis_type;
#endif
#else
    typedef bases<basis_u_type,basis_p_type> basis_type;
#endif
    //# endmarker1 #

    /*space*/
    //# marker2 #
#if defined( FEELPP_USE_LM )
#if defined(DIM2)
    typedef FunctionSpace<meshes<mesh_type,mesh_type, mesh_lag_type>, basis_type> space_type;
#else
    typedef FunctionSpace<meshes<mesh_type,mesh_type, mesh_lag_type, mesh_lag_type>, basis_type> space_type;
#endif
#else
    typedef FunctionSpace<mesh_type, basis_type> space_type;
#endif
    typedef boost::shared_ptr<space_type> space_ptrtype;

    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef typename space_type::element_type element_type;
    //# endmarker3 #


    FEELPP_DONT_INLINE
    Stokes( std::string const& n );

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

    std::string name() const
    {
        return config_name;
    }

private:


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename ExprUExact, typename ExprPExact>
    FEELPP_DONT_INLINE
    void exportResults( ExprUExact uexact, ExprPExact pexact,
                        element_type& u, element_type& v );

private:

    std::string config_name;
    double meshSize;
    double L;
    double p_in,p_out;
    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    mesh_lag_ptrtype mesh_lag;
    space_ptrtype Xh;
}; // Stokes



template<int nDim, int uOrder, int geoOrder>
Stokes<nDim,uOrder,geoOrder>::Stokes(std::string const& n )
    :
    super(),
    config_name( n ),
    meshSize( doption("hsize") ),
    L( this->vm()["L"].template as<value_type>() ),
    p_in( this->vm()["p_in"].template as<value_type>() ),
    p_out( this->vm()["p_out"].template as<value_type>() ),
    mu( this->vm()["mu"].template as<value_type>() ),
    penalbc( this->vm()["bccoeff"].template as<value_type>() )
{

}

template<int nDim, int uOrder, int geoOrder>
void
Stokes<nDim,uOrder,geoOrder>::init()
{
    std::string geoname = fs::path( soption( _name="geofile" ) ).stem().string();
    Environment::changeRepository( boost::format( "benchmarks/%1%/%2%/%3%D/P%4%P%5%G%6%/h_%7%/l_%8%" )
                                   % this->about().appName() % geoname
                                   % nDim
                                   % basis_u_type::nOrder % basis_p_type::nOrder % geoOrder
                                   % meshSizeInit()
                                   % level());



#if defined( DIM2 )
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo(_filename=soption(_name="geofile"),
                                     _h=meshSizeInit()),
                           _refine=level());
#elif defined(DIM3)
    if ( geoname == "straighttube" )
    {
        GeoTool::Node Centre(0,0,0);
        GeoTool::Node Rayon( 1);
        GeoTool::Node Dir(1,0,0);
        GeoTool::Node Lg(5,0,0);
        GeoTool::Cylindre C( meshSizeInit(),"Cyl",Centre,Dir,Rayon,Lg);
        C.setMarker(_type="surface",_name="inlet",_marker1=true);
        C.setMarker(_type="surface",_name="outlet",_marker2=true);
        C.setMarker(_type="surface",_name="wall",_marker3=true);
        C.setMarker(_type="volume",_name="Omega",_markerAll=true);

        mesh = C.createMesh(_mesh= new mesh_type,
                            _name="straighttube",
                            _partitions=Environment::worldComm().localSize(),
                            _worldcomm=Environment::worldComm(),
                            _refine=level() );

    }
    else
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=geo(_filename=soption(_name="geofile"),
                                         _h=meshSizeInit()),
                               _refine=level());
    }

#endif

#if defined (FEELPP_USE_LM)
    mesh_lag = merge( mesh->trace( markedfaces(mesh,"inlet") ),
                      mesh->trace( markedfaces(mesh,"outlet") ) );
#if defined( DIM2 )
    auto meshv = fusion::make_vector(mesh,mesh,mesh_lag);
#else
    auto meshv = fusion::make_vector(mesh,mesh,mesh_lag,mesh_lag);
#endif
    Xh = space_type::New( _mesh=meshv );
#else
    Xh = space_type::New( mesh );
#endif
}
template<int nDim, int uOrder, int geoOrder>
void
Stokes<nDim,uOrder,geoOrder>::run()
{
    mpi::timer chrono;
    this->init();
    LOG(INFO) << "chrono init: " << chrono.elapsed() << "\n";

    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
#if defined( FEELPP_USE_LM )
    auto lambda1 = U.template element<2>();
    auto nu1 = V.template element<2>();
#if defined(DIM3)
    auto lambda2 = U.template element<3>();
    auto nu2 = V.template element<3>();
#endif
#endif
    //# endmarker4 #

    M_stats.put( "h",M_meshSize );
    M_stats.put( "n.space.nelts",Xh->template functionSpace<0>()->mesh()->numElements() );
    M_stats.put( "n.space.ndof",Xh->nDof() );
    M_stats.put( "n.space.ndof.u",Xh->template functionSpace<0>()->nDof() );
    M_stats.put( "n.space.ndof.p",Xh->template functionSpace<1>()->nDof() );

    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]      number of dof(U): " << Xh-> template functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh-> template functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]      number of dof(P): " << Xh-> template functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh-> template functionSpace<1>()->nLocalDof()  << "\n";
#if defined( FEELPP_USE_LM )
    LOG(INFO) << "[dof]      number of dof(L): " << Xh-> template functionSpace<2>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(L): " << Xh-> template functionSpace<2>()->nLocalDof()  << "\n";
#endif
    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << meshSize << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";google::FlushLogFiles(google::GLOG_INFO);

    auto deft = sym(gradt( u ));
    auto def = sym(grad( v ));

    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+2*mu*deft*N();

    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+2*mu*def*N();

    std::string inlet("inlet"), outlet("outlet"), wall("wall");

    chrono.restart();
    // right hand side
    auto stokes_rhs = form1( _test=Xh );
    stokes_rhs += integrate( markedfaces( mesh,inlet ),  -trans(p_in*N())*id( v ) );
    stokes_rhs += integrate( markedfaces( mesh,outlet ), -trans(p_out*N())*id( v ) );
    LOG(INFO) << "chrono lhs: " << chrono.elapsed() << "\n";
    LOG(INFO) << "[stokes] vector local assembly done\n";google::FlushLogFiles(google::GLOG_INFO);

    // left hand side
    auto stokes = form2( _test=Xh, _trial=Xh );

    stokes += integrate( elements( mesh ), 2*mu*trace(trans( deft)*def ) );
    LOG(INFO) << "chrono mu*inner(deft,def): " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) - divt( u )*id( q ) );
    LOG(INFO) << "chrono (u,p): " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    chrono.restart();
    stokes +=integrate( _range=markedfaces( mesh, "wall" ), _expr=- trans(SigmaNt)*id(v) - trans(SigmaN)*idt(v) + doption(_name="bccoeff")*trans(idt(v))*id(v)/hFace());
    auto t = vec(Ny(),-Nx());

#if defined( FEELPP_USE_LM )
#if defined( DIM2 )
    stokes +=integrate( markedfaces( mesh,inlet ), -trans(id( v ))*t*idt( lambda1 ));
    stokes +=integrate( markedfaces( mesh,inlet ), -trans( idt( u ))*t*id( nu1 ) );
    stokes +=integrate( markedfaces( mesh,outlet ), -trans(id( v ))*t*idt( lambda1 ) );
    stokes +=integrate( markedfaces( mesh,outlet ), -trans( idt( u ))*t*id( nu1 ) );
    LOG(INFO) << "chrono (lambda,p): " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    chrono.restart();
    // set lagrange multipliers to 0 on boundary weakly
    stokes += integrate( boundaryfaces(mesh_lag), doption(_name="bccoefflag")*idt(nu1)*id(nu1)/hFace() );
    LOG(INFO) << "chrono (lambda,lambda): " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);

#elif defined( DIM3 )
    auto alpha = 1./sqrt(1-Nz()*Nz());
    auto C = alpha*mat<3,2>( cst(0.), Ny(), cst(0.), -Nx(), cst(1.), cst(0.) );
    auto lagt=vec(idt(lambda1),idt(lambda2));
    auto lag=vec(id(lambda1),id(lambda2));
    auto Clag = alpha*vec( id(lambda2)*Ny(), -id(lambda2)*Nx(), id(lambda1) );
    auto Clagt = alpha*vec( idt(lambda2)*Ny(), -idt(lambda2)*Nx(), idt(lambda1) );
    //stokes +=integrate( markedfaces( mesh,inlet ), -trans(id(v))*(C*lagt));
    for( auto bdy : { inlet, outlet } )
    {
        stokes +=integrate( markedfaces( mesh, bdy ),
                            -trans(cross(id(v),N()))*(Clagt) );
        stokes +=integrate( markedfaces( mesh, bdy ),
                            -trans( cross(idt( u ),N()))*(Clag) );
    }
    //stokes +=integrate( markedfaces( mesh,outlet ),-trans(cross(id(v),N()))*(C*lagt) -trans( cross(idt( u ),N()))*(C*lag) );
    //stokes += integrate( markedfaces( mesh,inlet ), doption("eps")*trans(idt(lambda))*id( nu ) );
    //stokes += integrate( markedfaces( mesh,outlet ), doption("eps")*trans(idt(lambda))*id( nu ) );
    // set lagrange multipliers to 0 on boundary weakly
    stokes += integrate( boundaryfaces(mesh_lag), doption(_name="bccoefflag")*(idt(nu1)*id(nu1)+idt(nu2)*id(nu2))/hFace() );
#endif // DIM3


    LOG(INFO) << "chrono (lambda,p): " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);

    chrono.restart();
#endif


    //stokes +=integrate( markedfaces( mesh,wall ), -inner( SigmaNt,id( v ) ) );
    //stokes +=integrate( markedfaces( mesh,wall ), -inner( SigmaN,idt( u ) ) );
    //stokes +=integrate( markedfaces( mesh,wall ), +penalbc*inner( idt( u ),id( v ) )/hFace() );

#if! defined( FEELPP_USE_LM )
    auto cTau = doption("gamma-tau")*mu/hFace();
    stokes +=integrate( markedfaces( mesh,inlet ), cTau*trans(cross(idt(v),N()))*cross(id(v),N()) );
    stokes +=integrate( markedfaces( mesh,outlet ), cTau*trans(cross(idt(v),N()))*cross(id(v),N()) );
#endif
    LOG(INFO) << "chrono bc: " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    chrono.restart();
    //# endmarker7 #
    //stokes+=on(_range=markedfaces(mesh,"wall"), _element=u, _rhs=stokes_rhs, _expr=zero<nDim,1>());
    //stokes+=on(_range=boundaryfaces(mesh), _element=u, _rhs=stokes_rhs, _expr=zero<nDim,1>());

    chrono.restart();
    auto retsolver = stokes.solve( _solution=U, _rhs=stokes_rhs, _rebuild=true );

    M_stats.put( "d.solver.bool.converged",retsolver.template get<0>() );
    M_stats.put( "d.solver.int.nit",retsolver.template get<1>() );
    M_stats.put( "d.solver.double.residual",retsolver.template get<2>() );

#if defined( DIM2 )
    LOG(INFO) << "int_outlet T = " << integrate( markedfaces( mesh,outlet ), t ).evaluate() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    LOG(INFO) << "int_inlet T = " << integrate( markedfaces( mesh,inlet ), t ).evaluate() << "\n";
    LOG(INFO) << "int_outlet u.T = " << integrate( markedfaces( mesh,outlet ), trans(idv(u))*t ).evaluate() << "\n";
    LOG(INFO) << "int_inlet  u.T = " << integrate( markedfaces( mesh,inlet ), trans(idv(u))*t ).evaluate() << "\n";
    LOG(INFO) << "int_outlet x u.T = " << integrate( markedfaces( mesh,outlet ), Px()*trans(idv(u))*t ).evaluate() << "\n";
    LOG(INFO) << "int_inlet x u.T = " << integrate( markedfaces( mesh,inlet ), Px()*trans(idv(u))*t ).evaluate() << "\n";
    LOG(INFO) << "int_outlet x^2 u.T = " << integrate( markedfaces( mesh,outlet ), Px()*Px()*trans(idv(u))*t ).evaluate() << "\n";
    LOG(INFO) << "int_inlet x^2 u.T = " << integrate( markedfaces( mesh,inlet ), Px()*Px()*trans(idv(u))*t ).evaluate() << "\n";google::FlushLogFiles(google::GLOG_INFO);
#elif defined( DIM3 )
    LOG(INFO) << "int_outlet u.T = " << integrate( markedfaces( mesh,outlet ), cross(idv(u),N()) ).evaluate() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    LOG(INFO) << "int_inlet  u.T = " << integrate( markedfaces( mesh,inlet ), cross(idv(u),N()) ).evaluate() << "\n";
    LOG(INFO) << "int_outlet x u.T = " << integrate( markedfaces( mesh,outlet ), Px()*cross(idv(u),N()) ).evaluate() << "\n";
    LOG(INFO) << "int_inlet x u.T = " << integrate( markedfaces( mesh,inlet ), Px()*cross(idv(u),N()) ).evaluate() << "\n";
    LOG(INFO) << "int_outlet x^2 u.T = " << integrate( markedfaces( mesh,outlet ), Px()*Px()*cross(idv(u),N() )).evaluate() << "\n";
    LOG(INFO) << "int_inlet x^2 u.T = " << integrate( markedfaces( mesh,inlet ), Px()*Px()*cross(idv(u),N())).evaluate() << "\n";google::FlushLogFiles(google::GLOG_INFO);
#endif
    LOG(INFO) << "chrono solver: " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);

    auto r=1;
    auto L=5;
    auto P_inlet=p_in;
    auto P_outlet=p_out;
#if defined(DIM2)
    auto u_exact=vec((1-(Py()*Py())/(r*r))*(P_inlet-P_outlet)/(2*mu*L),cst(0.) );

    auto p_exact=(P_outlet-P_inlet)*Px()/L + P_inlet;
#elif defined(DIM3)
    auto u_exact=vec(  (P_inlet-P_outlet)*r*r*(1-(Py()*Py()+Pz()*Pz())/(r*r))/(4*mu*L) , cst(0.) , cst(0.) );

    auto p_exact=(-Px()*(P_inlet-P_outlet)/L + P_inlet);
#endif
    this->exportResults( u_exact, p_exact, U, V );
    LOG(INFO) << "chrono export: " << chrono.elapsed() << "\n";google::FlushLogFiles(google::GLOG_INFO);

} // Stokes::run


template<int nDim, int uOrder, int geoOrder>
template<typename ExprUExact, typename ExprPExact>
void
Stokes<nDim,uOrder,geoOrder>::exportResults( ExprUExact u_exact, ExprPExact p_exact,
                                             element_type& U, element_type& V )
{
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto v = V.template element<0>();
    auto q = V.template element<1>();
#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
#endif

    auto Uh = Pchv<5>( mesh );

    auto u_exact_proj=project(_space=Uh,_range=elements(mesh),_expr=u_exact);

    double u_errorL2 = normL2( _range=elements( u.mesh() ), _expr=( idv( u )-u_exact ) );
    LOG(INFO) << "||u_error||_2 = " << u_errorL2 << "\n";;
    M_stats.put( "e.l2.u", u_errorL2 );

    double u_errorsemiH1 = integrate( _range=elements( mesh ),
                                      _expr=trace( ( gradv( u )-gradv(u_exact_proj) )*trans( gradv( u )-gradv(u_exact_proj) ) ) ).evaluate()( 0, 0 );
    double u_error_H1 = math::sqrt( u_errorL2*u_errorL2+u_errorsemiH1 );
    std::cout << "||u_error||_1= " << u_error_H1 << "\n";
    M_stats.put( "e.h1.u",u_error_H1 );
    M_stats.put( "e.semih1.u", math::sqrt(u_errorsemiH1) );

    double mean_p = mean( _range=elements( u.mesh() ), _expr=idv( p ) )(0,0);
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";

    double p_errorL2 = normL2( _range=elements( u.mesh() ), _expr=( idv( p ) - p_exact ) );
    LOG(INFO) << "||p_error||_2 = " << p_errorL2 << "\n";
    M_stats.put( "e.l2.p", p_errorL2 );

    double mean_div_u = mean( _range=elements( u.mesh() ), _expr=divv( u ) )(0,0);
    LOG(INFO) << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = normL2( _range=elements( u.mesh() ), _expr=divv( u ) );
    LOG(INFO) << "[stokes] ||div(u)||_2=" << div_u_error_L2 << "\n";
    M_stats.put( "e.l2.div", div_u_error_L2 );

    auto deff = sym(gradv(u));
    auto SigmaNN =(-idv(p)*N()+2*mu*deff*N());
    auto s = integrate( markedfaces( mesh,"inlet" ) , SigmaNN,_quad=_Q<12>()).evaluate();
    std::cout << "s inlet="  << s.norm() << "\n";
    std::cout << "error s inlet="  << math::abs(pi-s.norm()) << "\n";
    M_stats.put( "e.output.Fin", math::abs(pi-s.norm()) );
    M_stats.put( "d.output.Fin", s.norm() );

    v = vf::project( u.functionSpace(), elements( u.mesh() ), u_exact );
    q = vf::project( p.functionSpace(), elements( p.mesh() ), p_exact );

#if defined( FEELPP_USE_LM )
    auto exporter1d1 = exporter( _mesh=lambda.mesh(), _name=this->about().appName()+"_lambda" );
    exporter1d1->add( "lambda", lambda );
    exporter1d1->save();
#endif

    auto myexporter = exporter( _mesh=mesh, _name=this->about().appName() );
    myexporter->add( "u", u );
    myexporter->add( "p", p );
    myexporter->add( "u_exact", v );
    myexporter->add( "p_exact", q );
    myexporter->save();

} // Stokes::export
} // Feel
