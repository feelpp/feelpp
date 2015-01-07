/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \file stokes_kovasnay_curve.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#ifndef __FEELPP_BENCH_KOVASNAY_CURVE_HPP
#define __FEELPP_BENCH_KOVASNAY_CURVe_HPP 1

#include <boost/any.hpp>
#include <boost/utility.hpp>

#include <feel/feel.hpp>


#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>



#include <feel/feelmesh/elements.hpp>
#include <feel/feel.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
/*inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesKovasnayCurveOptions( "Stokes_Kovasnay_Curve options" );
    stokesKovasnayCurveOptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return stokesKovasnayCurveOptions.add( Feel::feel_options() ) ;
}
*/

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
/*inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "stokes_kovasnay_curve" ,
                           "stokes_kovasnay_curve" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2012 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;

}

*/

namespace Feel
{
using namespace Feel::vf;

/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
template<int POrder=1,int GeoOrder=1>
class Stokes_Kovasnay_Curve
    :
public Simget
{
    typedef Simget super;
public:


    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<2,GeoOrder> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/

    typedef Lagrange<POrder+1, Vectorial> basis_u_type;
    typedef Lagrange<POrder, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;

    // use lagrange multipliers to ensure zero mean pressure
    //#if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
    //#else
    //typedef bases<basis_u_type,basis_p_type> basis_type;
    //#endif

    typedef bases<basis_u_type> basis_type_U;
    typedef FunctionSpace<mesh_type, basis_type_U> space_type_U;
    typedef boost::shared_ptr<space_type_U> space_ptrtype_U;

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef typename space_type::element_type element_type;
    //# endmarker3 #

    /* export */
    typedef Exporter<mesh_type> export_type;

    // FEELPP_DONT_INLINE
    //Stokes_Kovasnay_Curve(std::string const& basis_name );

    // init mesh and space
    //FEELPP_DONT_INLINE
    void init();


    std::string name() const
    {
        return M_basis_name;
    }

   Stokes_Kovasnay_Curve( std::string const& basis_name)
        :
        super( ),
        M_basis_name( basis_name ),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        mu( this->vm()["mu"].template as<value_type>() ),
        penalbc( this->vm()["bccoeff"].template as<value_type>() ),
        exporter2( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
 
    }
    /**
     * run the convergence test
     */
    //FEELPP_DONT_INLINE
    void run();

private:


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename ExprUExact, typename ExprPExact>
    //FEELPP_DONT_INLINE
    void exportResults( ExprUExact uexact, ExprPExact pexact,
                        element_type& u, element_type& v);

private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string M_basis_name;
    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    space_ptrtype_U Uh;
    sparse_matrix_ptrtype M,D,A;
    vector_ptrtype F;

    boost::shared_ptr<export_type> exporter2;
}; // Stokes

    //template<int POrder, int GeoOrder>




template<int POrder, int GeoOrder>
void
Stokes_Kovasnay_Curve<POrder,GeoOrder>::init()
{
    /* if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
        }*/

    if ( this->vm().count( "nochdir" ) == false )
        this->changeRepository( boost::format( "benchmarks/%1%/%2%/Part%6%/P%3%P%4%/h_%5%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % basis_u_type::nOrder % basis_p_type::nOrder
                                % doption("hsize")
                                % Environment::numberOfProcessors() );


    //---------------------------------------------------------------------------------------
   mesh = createGMSHMesh( _mesh=new mesh_type,
                          _desc=domain( _name="kovaznay",
                                        _usenames=false,
                                        _shape="hypercube",
                                        _h=meshSize,
                                        _xmin=-0.5, _xmax=1,
                                        _ymin=-0.5, _ymax=1.5 ) );
    auto Wh = Pchv<POrder>( mesh );
    auto uu = Wh->element();
    auto vv = Wh->element();

    auto FF = M_backend->newVector( Wh );
    auto DD =  M_backend->newMatrix( Wh, Wh );
  
    auto ll = form1( _test=Wh,  _vector=FF  );

    auto aa = form2( _trial=Wh, _test=Wh, _matrix=DD );
    aa = integrate(_range=elements(mesh),
                  _expr=trace(gradt(uu)*trans(grad(vv))) );
    aa+=on(_range=markedfaces(mesh,1), _rhs=ll, _element=uu,
          _expr=zero<2,1>() );
    aa+=on(_range=markedfaces(mesh,3), _rhs=ll, _element=uu,
          _expr=zero<2,1>() );
    aa+=on(_range=markedfaces(mesh,4), _rhs=ll, _element=uu,
          _expr=zero<2,1>() );
    aa+=on(_range=markedfaces(mesh,2), _rhs=ll, _element=uu,
          _expr=vec(cst(0.),0.08*(Px()+0.5)*(Px()-1)*(Px()*Px()-1)));

    aa.solve(_rhs=ll,_solution=uu);
  

    //auto m1 = lagrangeP1(_space=Wh)->mesh();
    /*//-----------------------------------------
    auto XhVisu = Pchv<1>(m1);
    auto opIVisu = opInterpolation(_domainSpace=Vh,
                                   _imageSpace=XhVisu,
                                   _type=InterpolationNonConforme(false,true,false) );
    auto uVisu = opIVisu->operator()(u);
    auto e = exporter( _mesh=m1, _name="initial_visu" );
    e->step(0)->setMesh( m1 );
    e->step(0)->add( "u", uVisu );
    e->save();

    meshMove( m1, uVisu );

    auto e1 = exporter( _mesh=m1, _name="moved_visu" );
    e1->step(0)->setMesh( m1  );
    e1->step(0)->add( "u_visu", uVisu );
    e1->save();*/

    //-----------------------
	/*auto e2 = exporter( _mesh=m1, _name="initial" );
    e2->step(0)->setMesh( m1 );
    e2->step(0)->add( "uu", uu );
    e2->save();*/

    meshMove( mesh, uu );
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";
    LOG(INFO) << "number of dof in Wh: " << Wh->nDof() << "\n";
    std::cout << "number of dof in Wh: " << Wh->nDof() << "\n";
    /* auto Vh4 = Pchv<4>( mesh );
    /* auto Vh4 = Pchv<4>( mesh );
    auto m2 = lagrangeP1(_space=Vh4)->mesh();
    auto e3 = exporter( _mesh=m2, _name="moved" );
    e3->step(0)->setMesh( m2  );
    e3->step(0)->add( "u", u );
    e3->save();*/
    //-------------------


    mesh=straightenMesh(_mesh=mesh);

    /* //********************** Rectangle ***************************************
    mesh = createGMSHMesh( _mesh=new mesh_type,
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                            _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % convex_type().dimension() % 1 ).str() ,
                                          _shape="hypercube",
                                          _dim=convex_type().dimension(),
                                          _h=meshSize ) );

     mesh->addMarkerName( "inlet", 1, 1 );
     mesh->addMarkerName( "outlet", 3, 1 );
     mesh->addMarkerName( "wall1", 2, 1 );
     mesh->addMarkerName( "wall2", 4, 1 );
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";
    */




    Xh = space_type::New( mesh );
    Uh=space_type_U::New(mesh);
    F = M_backend->newVector( Xh );
    D =  M_backend->newMatrix( Xh, Xh );
    A =  M_backend->newMatrix( Xh, Xh );
}
template<int POrder, int GeoOrder>
void
Stokes_Kovasnay_Curve<POrder,GeoOrder>::run()
{
    this->init();

    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    //#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    //#endif
    //# endmarker4 #

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "hsize = " << meshSize << "\n";
    M_stats.put( "hsize", meshSize );
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";




    //# marker5 #
    auto deft = gradt( u );//sym
    auto def = grad( v );//sym
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+mu*deft*N();//*2
    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+mu*def*N();//*2
    //# endmarker6 #


    auto r=1;
    auto L=5;
    auto P_inlet=1;
    auto P_outlet=0.;

    //*********************  solution exacte kovasnay 2D *****************
    //auto mu = 0.035;
    double lambdaa = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*pi*pi );
    // total stress tensor (test)
    auto u1 = 1. - exp( lambdaa * Px() ) * cos( 2.*pi*Py() );
    auto u2 = ( lambdaa/( 2.*pi ) ) * exp( lambdaa * Px() ) * sin( 2.*pi*Py() );
    auto u_exact = val( vec( u1,u2 ) );

    auto du_dx = ( -lambdaa*exp( lambdaa * Px() )*cos( 2.*pi*Py() ) );
    auto du_dy = ( 2*pi*exp( lambdaa * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dx = ( ( lambdaa*lambdaa/( 2*pi ) )*exp( lambdaa * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dy = ( lambdaa*exp( lambdaa * Px() )*cos( 2.*pi*Py() ) );
    auto grad_exact = val( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = val( du_dx + dv_dy );

    //auto p_exact = val( ( -exp( 2.*lambdaa*Px() ) )/2.0-0.125*( exp( -1.0*lambdaa )-1.0*exp( 3.0*lambdaa ) )/lambdaa );
    auto p_exact = val( ( -0.5*exp( 2.*lambdaa*Px() ) ) );

    //auto f1 = (exp( lambdaa * Px() )*((lambdaa*lambdaa - 4.*pi*pi)*mu*cos(2.*pi*Py()) - lambdaa*exp( lambdaa * Px() )));
    auto f1 = ( -mu*( -lambdaa*lambdaa*exp( lambdaa*Px() )*cos( 2.0*pi*Py() )+4.0*exp( lambdaa*Px() )*cos( 2.0*pi*Py() )*pi*pi )-lambdaa*exp( 2.0*lambdaa*Px() ) );

    //auto f2 = (exp( lambdaa * Px() )*mu*sin(2.*pi*Py())*(-lambdaa*lambdaa +4*pi*pi));
    auto f2 = ( -mu*( lambdaa*lambdaa*lambdaa*exp( lambdaa*Px() )*sin( 2.0*pi*Py() )/pi/2.0-2.0*lambdaa*exp( lambdaa*Px() )*sin( 2.0*pi*Py() )*pi ) );

    auto f = val( vec( f1,f2 ) ); //+ convection;

    //double pmean = integrate( elements(mesh), p_exact ).evaluate()( 0, 0 )/mesh->measure();
    //double pmean = -0.125*( math::exp( -1.0*lambda )-1.0*math::exp( 3.0*lambda ) )/lambda;
    //*************************************************************************************

    std::cout<< "number of dof in Xh: " << Xh->nDof() << "\n";
    auto taille=5.5+1.503061665;
    auto mean_p_exact =integrate( elements( mesh ),  p_exact,_quad=_Q<20>() ).evaluate()(0,0) /mesh->measure();//taille;
    std::cout << "[stokes] mean(p_exact)=" << mean_p_exact << "\n";

    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( elements( mesh ),inner( f,id( v ) ),_quad=_Q<20>());

   //#if defined( FEELPP_USE_LM )
    stokes_rhs += integrate( elements( mesh ),id(nu)* p_exact,_quad=_Q<20>() );   //add the mean of the exact pressure
    //#endif
    //************
    LOG(INFO) << "[stokes] vector local assembly done\n";

    // Construction of the left hand side
    //# marker7 #
    auto stokes = form2( _test=Xh, _trial=Xh, _matrix=D );
    boost::timer chrono;
    stokes += integrate( elements( mesh ), mu*inner( deft,def ) ,_quad=_Q<20>());//*2
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ) ,_quad=_Q<20>());
    //#if defined( FEELPP_USE_LM )
    stokes +=integrate( elements( mesh ), id( q )*idt( lambda ) + idt( p )*id( nu ),_quad=_Q<20>());
    chrono.restart();
    //#endif

    stokes+=on( _range=boundaryfaces(mesh), _element=u,_rhs=stokes_rhs,
                _expr=u_exact );

    M_backend->solve( _matrix=D, _solution=U, _rhs=F );
    // stokes.solve( _rhs=stokes_rhs, _solution=U);


    size_type nnz = 0 ;
    auto nNz = D->graph()->nNz() ;
    for ( auto iter = nNz.begin(); iter!=nNz.end(); ++iter )
        nnz += ( *iter ) ;
    size_type gnnz=0;
    mpi::all_reduce( this->comm(), nnz, gnnz, [] ( size_type x, size_type y ) {return x + y;} );
    LOG(INFO) << "[nnz]       number of local nnz: " << nnz << " global: "  << gnnz << "\n";
    LOG(INFO) << "[dof]             number of dof in Xh: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]        number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]          number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof]     number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]          number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof]     number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    std::cout << "[nnz]       number of local nnz: " << nnz << " global: "  << gnnz << "\n";
    std::cout << "[dof]             number of dof: " << Xh->nDof() << "\n";
    std::cout << "[dof]        number of dof/proc: " << Xh->nLocalDof() << "\n";
    std::cout << "[dof]          number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    std::cout << "[dof]     number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    std::cout << "[dof]          number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    std::cout << "[dof]     number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    M_stats.put("nnz" , nnz);
    M_stats.put("dof ",Xh);


    this->exportResults( u_exact, p_exact, U, V );
} // Stokes::run


template<int POrder, int GeoOrder>
template<typename ExprUExact, typename ExprPExact>
void
Stokes_Kovasnay_Curve<POrder,GeoOrder>::exportResults( ExprUExact u_exact, ExprPExact p_exact,
                                                       element_type& U, element_type& V )
{
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto v = V.template element<0>();
    auto q = V.template element<1>();
    //#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    LOG(INFO) << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";
    std::cout << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";
    //#endif

    auto u_exact_proj=vf::project(Uh,elements(mesh),u_exact);

    double u_errorL2 = integrate( elements( u.mesh() ), trans( idv( u )-u_exact )*( idv( u )-u_exact ) ,_quad=_Q<20>()).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";
    LOG(INFO) <<"||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";
    M_stats.put("u.l2 ", math::sqrt( u_errorL2 ));

    double meas = integrate( elements( u.mesh() ), cst( 1.0 ),_quad=_Q<20>()).evaluate()( 0, 0 );
#if (STOKESPRESSMESHTYPE ==2)
    LOG(INFO) << "[stokes] measure(Omega)=" << meas << " (should be equal to "<< 4*math::atan(1.)*5 << ")\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to "<< 4*math::atan(1.)*5 << ")\n";
#elif (STOKESPRESSMESHTYPE == 1)
    LOG(INFO) << "[stokes] measure(Omega)=" << meas << " (should be equal to  1)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to  1)\n";
#endif

    double mean_p = integrate( elements( u.mesh() ), idv( p ), _quad=_Q<20>()).evaluate()( 0, 0 )/meas;
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";

    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    // double p_errorL2 = integrate( elements( u.mesh() ), ( idv( p )+mean_p_exact - p_exact )*( idv( p )+mean_p_exact-p_exact ), _quad=_Q<QuadOrder>() _quad=_Q<20>()).evaluate()( 0, 0 );
    double p_errorL2 = integrate( elements( u.mesh() ), ( idv( p ) - p_exact )*( idv( p )-p_exact ),_quad=_Q<20>()).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";
    LOG(INFO) <<"||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";
    M_stats.put("p.l2 ", math::sqrt( p_errorL2 ));
    LOG(INFO) << "[stokes] solve for D done\n";



    double u_errorH1 = integrate( elements( u.mesh() ),  trans( idv( u )-u_exact )*( idv( u )-u_exact ),_quad=_Q<20>()).evaluate()( 0, 0 ) +  integrate( elements( u.mesh() ),  trans( gradv( u ) -  gradv ( u_exact_proj ) )*( gradv( u ) -  gradv ( u_exact_proj )), _quad=_Q<20>()).evaluate()( 0, 0 );
    double H1_u=math::sqrt( u_errorH1 );
    std::cout << "||u_errorH1||_2 = " <<H1_u<< "\n";
    LOG(INFO) << "||u_errorH1||_2 = " <<H1_u<< "\n";
    M_stats.put("u.H1 ", H1_u);

    double mean_div_u = integrate( elements( u.mesh() ), divv( u ) ,_quad=_Q<20>()).evaluate()( 0, 0 );
    LOG(INFO) << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements( u.mesh() ), divv( u )*divv( u ) ,_quad=_Q<20>()).evaluate()( 0, 0 );
    LOG(INFO) << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";


    //#if (STOKESPRESSMESHTYPE ==2)
    auto Du= gradv(u);//sym
    auto Dv= gradv(v);//sym
    auto SigmaNN =-idv(p)*vf::N()+mu*Du*vf::N();//*2

    //**************  F ******************
    /* auto FappCur = integrate(markedfaces( mesh,2) , inner(SigmaNN,idv(v)), _quad=_Q<20>()).evaluate();
    std::cout.precision(17);
    std::cout << "FappCur = "<<FappCur << "\n" ;
    std::cout << "|Fex-Fappcur| = "<< math::abs(-0.6709310448-FappCur(0,0)) << "\n" ;
    // LOG(INFO) << "FappCur = "<<math::abs(2.486163775- FappCur) << "\n" ;
    */

    auto pI = integrate(markedfaces( mesh,2) , -idv(p)*vf::N(), _quad=_Q<20>()).evaluate();
    std::cout.precision(17);
    std::cout << "pI1 = "<<pI(0,0) << "\n" ;
    std::cout << "pI2 = "<<pI(1,0) << "\n" ;
    auto gradient = integrate(markedfaces( mesh,2) , mu*Du*vf::N(), _quad=_Q<20>()).evaluate();
    std::cout.precision(17);
    std::cout << "mu*Gradu.n1 = "<<gradient(0,0) << "\n" ;
    std::cout << "mu*Gradu.n2 = "<<gradient(1,0) << "\n" ;

    auto SigmaN = integrate(markedfaces( mesh,2) , SigmaNN, _quad=_Q<20>()).evaluate();
    std::cout.precision(17);
    std::cout << " SigmaN1 = "<< SigmaN(0,0) << "\n" ;
    std::cout << " SigmaN2 = "<< SigmaN(1,0) << "\n" ;
    std::cout << "||Fex-Fappcur||_2 = "<< math::sqrt((0.08217721411-SigmaN(0,0))*(0.08217721411-SigmaN(0,0))+(-0.7531082589-SigmaN(1,0))*(-0.7531082589-SigmaN(1,0))) << "\n" ;

    M_stats.put("Fapp ",  math::sqrt((0.08217721411-SigmaN(0,0))*(0.08217721411-SigmaN(0,0))+(-0.7531082589-SigmaN(1,0))*(-0.7531082589-SigmaN(1,0))));


    //**************  Somme des integrales  ******************
    auto v1=vec(cst(1.), cst(0.));
    v=vf::project(Xh->template functionSpace<0>(), markedfaces(mesh, 2), v1 );

    double lambda2 = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*pi*pi );
    auto f11 =  ( -mu*( -lambda2*lambda2*exp( lambda2*Px() )*cos( 2.0*pi*Py() )+4.0*exp( lambda2*Px() )*cos( 2.0*pi*Py() )*pi*pi )-lambda2*exp( 2.0*lambda2*Px() ) );

    auto f22 = ( -mu*( lambda2*lambda2*lambda2*exp( lambda2*Px() )*sin( 2.0*pi*Py() )/pi/2.0-2.0*lambda2*exp( lambda2*Px() )*sin( 2.0*pi*Py() )*pi ) );

    auto ff = val( vec( f11,f22 ) ); //+ convection;

    auto sum1 =integrate( elements( mesh ),inner( -ff,idv( v ) ), _quad=_Q<20>()).evaluate();
    sum1 +=integrate( elements( mesh ),mu*inner( Du,Dv ) - divv( v )*idv( p), _quad=_Q<20>()).evaluate();//*2
    sum1 +=integrate( markedfaces( mesh,1 ),inner(idv(p)*vf::N()-mu*Du*vf::N(),idv(v)), _quad=_Q<20>()).evaluate();//*2
    sum1 +=integrate( markedfaces( mesh,3 ),inner(idv(p)*vf::N()-mu*Du*vf::N(),idv(v)),_quad=_Q<20>()).evaluate();//*2
    std::cout.precision(17);
    std::cout << "Sum1 = "<< sum1 << "\n" ;


    auto v2=vec(cst(0.), cst(1.));
    v=vf::project(Xh->template functionSpace<0>(), markedfaces(mesh, 2), v2 );
    auto sum2 =integrate( elements( mesh ),inner( -ff,idv( v ) ), _quad=_Q<20>()).evaluate();
    sum2 +=integrate( elements( mesh ),mu*inner( Du,Dv ) - divv( v )*idv( p), _quad=_Q<20>()).evaluate();//*2
    sum2 +=integrate( markedfaces( mesh,1 ),inner(idv(p)*vf::N()-mu*Du*vf::N(),idv(v)), _quad=_Q<20>()).evaluate();//*2
    sum2 +=integrate( markedfaces( mesh,3 ),inner(idv(p)*vf::N()-mu*Du*vf::N(),idv(v)),_quad=_Q<20>()).evaluate();//*2
    std::cout.precision(17);
    std::cout << "Sum2 = "<< sum2 << "\n" ;

    std::cout << "||Fex-R||_2 = "<< math::sqrt((0.08217721411-sum1(0,0))*(0.08217721411-sum1(0,0))+(-0.7531082589-sum2(0,0))*(-0.7531082589-sum2(0,0))) << "\n" ;

    M_stats.put("Sum ", math::sqrt((0.08217721411-sum1(0,0))*(0.08217721411-sum1(0,0))+(-0.7531082589-sum2(0,0))*(-0.7531082589-sum2(0,0))));


    auto u_ex = vf::project( u.functionSpace(), elements( u.mesh() ), u_exact );
    auto p_ex = vf::project( p.functionSpace(), elements( p.mesh() ), p_exact );


    if ( exporter2->doExport() )
    {
        /*auto m11 = lagrangeP1(_space=Xh)->mesh();
        auto XhVisu = Pchv<1>(m11);
        auto opIVisu = opInterpolation(_domainSpace=Xh,
                                       _imageSpace=XhVisu,
                                       _type=InterpolationNonConforme(false,true,false) );
        auto uVisu = opIVisu->operator()(u);
        auto pVisu = opIVisu->operator()(p);
        auto u_exactVisu = opIVisu->operator()(u_exact);
        auto p_exactVisu = opIVisu->operator()(p_exact);

        auto e22 = exporter2( _mesh=m11, _name="initial_visu" );
        e22->step(0)->setMesh( m11 );
        e22->step(0)->add( "u", uVisu );
        e22->step(0)->add( "p", pVisu );
        e22->step(0)->add( "u_exact", u_exactVisu );
        e22->step(0)->add( "p_exact", p_exactVisu );

        e22->save();

        /* meshMove( m1, uVisu );
        auto e1 = exporter( _mesh=m1, _name="moved_visu" );
        e1->step(0)->setMesh( m1  );
        e1->step(0)->add( "u_visu", uVisu );
        e1->save();*/

        exporter2->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter2->step( 0 )->addRegions();
        exporter2->step( 0 )->add( "v", V.template element<0>() );
        exporter2->step( 0 )->add( "u", U.template element<0>() );
        exporter2->step( 0 )->add( "u", U.template element<0>() );
        exporter2->step( 0 )->add( "p", U.template element<1>() );

        exporter2->step( 0 )->add( "u_exact", u_ex );
        exporter2->step( 0 )->add( "p_exact", p_ex );
        exporter2->save();
        }

} // Stokes::export
} // Feel
#endif// __FEELPP_BENCH_KOVASNAY_CURVe_HPP
