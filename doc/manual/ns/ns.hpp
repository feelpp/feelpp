/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-09-30

  Copyright (C) 2012 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file ns.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-09-30
 */

#include <feel/feel.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesoptions( "Stokes options" );
    stokesoptions.add_options()
        ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
        ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
        ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
        ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
        ( "geofile", Feel::po::value<std::string>()->default_value( "backstep.geo" ), "geometry file name" )
        ( "export-matlab", "export matrix and vectors in matlab" )
        ;
    return stokesoptions.add( Feel::feel_options() ) ;
}




namespace Feel
{
/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
class NavierStokes
    :
public Simget
{
    typedef Simget super;
public:


    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Simplex<FEELPP_NS_DIM> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #,
    typedef Lagrange<FEELPP_NS_ORDER+1, Vectorial> basis_u_type;
    typedef Lagrange<FEELPP_NS_ORDER, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    // use lagrange multipliers to ensure zero mean pressure
#if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
#else
    typedef bases<basis_u_type,basis_p_type> basis_type;
#endif
    //# endmarker1 #

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef space_type::element_type element_type;
    //# endmarker3 #

    typedef BDF<space_type> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    NavierStokes();

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

private:


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    FEELPP_DONT_INLINE
    void exportResults( double time, element_type& u );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Xh;

    boost::shared_ptr<export_type> exporter;
}; // NavierStokes


NavierStokes::NavierStokes()
    :
    super(),
    M_backend( backend_type::build( soption("backend") ) ),
    meshSize( this->vm()["hsize"].as<double>() ),
    mu( this->vm()["mu"].as<value_type>() ),
    penalbc( this->vm()["bccoeff"].as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

void
NavierStokes::init()
{
    Environment::changeRepository( boost::format( "doc/manual/ns/%1%/%2%/P%3%P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % convex_type::name()
                                   % basis_u_type::nOrder % basis_p_type::nOrder
                                   % this->vm()["hsize"].as<double>() );


    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename=this->vm()["geofile"].as<std::string>()
                                      _dim=FEELPP_NS_DIM,
                                      _h=meshSize ) );


    Xh = space_type::New( mesh );

}
void
NavierStokes::run()
{
    this->init();

    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "u" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "p" );
#if defined( FEELPP_USE_LM )
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
#endif
    //# endmarker4 #

    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]      number of dof(U): " << Xh->functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh->functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]      number of dof(P): " << Xh->functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh->functionSpace<1>()->nLocalDof()  << "\n";

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << meshSize << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";




    //# marker5 #
    auto deft = gradt( u )+trans(gradt(u));
    auto def = grad( v )+trans(grad(v));
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+mu*deft*N();

    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+mu*def*N();
    //# endmarker6 #

    auto F = M_backend->newVector( Xh );
    auto D =  M_backend->newMatrix( Xh, Xh );

    // right hand side
    auto ns_rhs = form1( _test=Xh, _vector=F );


    LOG(INFO) << "[navier-stokes] vector local assembly done\n";

    // construction of the BDF
    auto bdfns=bdf(_space=Xh);

    /*
     * Construction of the left hand side
     */

    auto navierstokes = form2( _test=Xh, _trial=Xh, _matrix=D );
    mpi::timer chrono;
    navierstokes += integrate( elements( mesh ), mu*inner( deft,def )+ trans(idt( u ))*id( v )*bdfns->polyDerivCoefficient( 0 ) );
    LOG(INFO) << "mu*inner(deft,def)+(bdf(u),v): " << chrono.elapsed() << "\n";
    chrono.restart();
    navierstokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ) );
    LOG(INFO) << "(u,p): " << chrono.elapsed() << "\n";
    chrono.restart();
#if defined( FEELPP_USE_LM )
    navierstokes +=integrate( elements( mesh ), id( q )*idt( lambda ) + idt( p )*id( nu ) );
    LOG(INFO) << "(lambda,p): " << chrono.elapsed() << "\n";
    chrono.restart();
#endif

    std::for_each( inflow_conditions.begin(), inflow_conditions.end(),
                   [&]( BoundaryCondition const&  bc )
                   {
                       // right hand side
                       ns_rhs += integrate( markedfaces( mesh, bc.marker() ), inner( idf(&bc,BoundaryCondition::operator()),-SigmaN+penalbc*id( v )/hFace() ) );

                       navierstokes +=integrate( boundaryfaces( mesh ), -inner( SigmaNt,id( v ) ) );
                       navierstokes +=integrate( boundaryfaces( mesh ), -inner( SigmaN,idt( u ) ) );
                       navierstokes +=integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
                   });
    std::for_each( wall_conditions.begin(), wall_conditions.end(),
                   [&]( BoundaryCondition const&  bc )
                   {
                       navierstokes +=integrate( boundaryfaces( mesh ), -inner( SigmaNt,id( v ) ) );
                       navierstokes +=integrate( boundaryfaces( mesh ), -inner( SigmaN,idt( u ) ) );
                       navierstokes +=integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
                   });
    std::for_each( outflow_conditions.begin(), outflow_conditions.end(),
                   [&]( BoundaryCondition const&  bc )
                   {
                       ns_rhs += integrate( markedfaces( mesh, bc.marker() ), inner( idf(&bc,BoundaryCondition::operator()),N() ) );
                   });

    LOG(INFO) << "bc: " << chrono.elapsed() << "\n";
    chrono.restart();

    u = vf::project( _space=Xh->functionSpace<0>(), _expr=cst(0.) );
    p = vf::project( _space=Xh->functionSpace<1>(), _expr=cst(0.) );

    M_bdf->initialize( U );

    for( bdfns->start(); bdfns->isFinished(); bdfns->next() )
    {
        // add time dependent terms
        auto bdf_poly = bdfns->polyDeriv();
        form1( _test=Xh, _vector=Ft ) =
            integrate( _range=elements(mesh), _expr=trans(idv( bdf_poly ))*id( v ) );
        // add convective terms
        form1( _test=Xh, _vector=Ft ) +=
            integrate( _range=elements(mesh), _expr=trans(gradv(u)*idv( u ))*id(v) );
        // add contrib from time independent terms
        Ft->add( 1., F );

        // add time stepping terms from BDF to right hand side
        backend()->solve( _matrix=D, _solution=U, _rhs=Ft );

        this->exportResults( bdfns->time(), U );
    }





} // NavierNavierstokes::run



void
Navierstokes::exportResults( double time, element_type& U )
{
    auto u = U.element<0>();
    auto p = U.element<1>();

    auto v = V.element<0>();
    auto q = V.element<1>();
#if defined( FEELPP_USE_LM )
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
    LOG(INFO) << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";

#endif

    double div_u_error_L2 = normL2( _range=elements( u.mesh() ), _expr=divv( u ) );
    LOG(INFO) << "[navierstokes] ||div(u)||_2=" << div_u_error_L2 << "\n";

    if ( exporter->doExport() )
    {
        exporter->step( time )->setMesh( U.functionSpace()->mesh() );
        exporter->step( time )->addRegions();
        exporter->step( time )->add( {"u","p","l"}, U );
        exporter->save();
    }

} // NavierStokes::export
} // Feel




