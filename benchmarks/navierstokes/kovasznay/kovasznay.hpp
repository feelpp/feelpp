/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009-2011 Christophe Prud'homme
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
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>



#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>
#include <fstream>
#include <sstream>

#include <feel/feelcore/applicationxml.hpp>
#include <feel/feelcore/xmlparser.hpp>

#define WITH_LAGRANGE_MULTIPLIERS 1

Feel::po::options_description makeOptions();
Feel::AboutData makeAbout();


namespace Feel
{
using namespace vf;
/**
 * \class Kovasznay class
 * \brief solves the stokes equations
 *
 */
template<int _OrderU,
         int _OrderP = _OrderU-1,
         typename Entity = Simplex<2,1> >
class Kovasznay
    :
public ApplicationXML
{
    typedef ApplicationXML super;
public:

    static const uint16_type Dim  = 2;
    static const uint16_type OrderU = _OrderU;
    static const uint16_type OrderP = _OrderP;
    static const bool is_equal_order = ( OrderU==OrderP );

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef Lagrange<OrderU, Vectorial> basis_u_type;
    typedef Lagrange<OrderP, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
#if defined(WITH_LAGRANGE_MULTIPLIERS)
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
#else
    typedef bases<basis_u_type,basis_p_type> basis_type;
#endif
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    BOOST_MPL_ASSERT( ( boost::is_same<typename space_type::bases_list, basis_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<0> >::type, typename basis_u_type::template ChangeTag<0>::type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<1> >::type, typename basis_p_type::template ChangeTag<1>::type> ) );
    //BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<2> >::type, basis_l_type> ) );
    typedef boost::shared_ptr<space_type> space_ptrtype;

    /* functions */
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Kovasznay( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    /**
     * run the convergence test
     */
    void run();

private:

    void buildLhs();
    void buildRhs();

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    double meshSize;
    double M_stabP;
    double M_stabD;
    space_ptrtype Xh;
    element_type U, V;
    sparse_matrix_ptrtype D;
    vector_ptrtype F;
    double mu;
    double M_lambda;
    bool M_weak_dirichlet;
    double penalbc;
    double M_beta;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;
    mesh_ptrtype mesh;
private:

    template<typename element_1_type, typename StabExpr>
    void addPressureStabilisation( element_1_type& p, element_1_type& q, StabExpr& stabexpr );

    template<typename element_0_type, typename StabExpr>
    void addDivergenceStabilisation( element_0_type& u, element_0_type& v, StabExpr& stabexpr );

}; // Kovasznay

template<int _OrderU, int _OrderP, typename Entity>
Kovasznay<_OrderU, _OrderP, Entity>::Kovasznay( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( doption("hsize") ),
    M_stabP( doption("penalisation") ),
    M_stabD( doption("penalisation") ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{
    M_weak_dirichlet = this->vm().count( "weak" );
    mu = this->vm()["mu"].template as<value_type>();
    Parameter h;

    switch ( OrderU )
    {
    case 1:
        h = Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.02:0.025:0.05" );
        break;

    case 2:
        h = Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.02:0.025:0.1" );
        break;

    case 3:
        h = Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.035:0.025:0.2" );
        break;
    }

    this->
    //addParameter( Parameter(_name="mu",_type=CONT_ATTR,_latex="\\mu", _values=boost::lexical_cast<std::string>( mu ).c_str()))
    addParameter( Parameter( _name="convex",_type=DISC_ATTR,_values= boost::lexical_cast<std::string>( convex_type::is_simplex  ).c_str() ) )
    .addParameter( Parameter( _name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str() ) )
    .addParameter( Parameter( _name="orderU",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderU  ).c_str() ) )
    .addParameter( Parameter( _name="orderP",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderP  ).c_str() ) )
    .addParameter( h );

    std::vector<Parameter> depend;
    std::vector<std::string> funcs;
    depend.push_back( h );
    std::ostringstream oss;

    if ( OrderP == OrderU )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP  );

    else if ( OrderP == OrderU-1 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderU  );

    else if ( OrderP == OrderU-2 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP+1  );

    funcs.push_back( oss.str() );
    std::vector<std::string> divfuncs;
    oss.str( "" );

    if ( OrderP == OrderU )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP-1./2.  );

    else if ( OrderP == OrderU-1 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderU-1./2.  );

    else if ( OrderP == OrderU-2 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP+1-1./2.  );

    divfuncs.push_back( oss.str() );


    oss.str( "" );
    std::vector<std::string> funcs2;
    oss << "h**" << boost::lexical_cast<std::string>( OrderP+1 ) ;
    funcs2.push_back( oss.str() );
    oss.str( "" );
    std::vector<std::string> funcs3;

    if ( OrderP == OrderU-2 )
        oss << "h**" << boost::lexical_cast<std::string>( OrderU ) ;

    else
        oss << "h**" << boost::lexical_cast<std::string>( OrderU+1 ) ;

    funcs3.push_back( oss.str() );

    this->
    addOutput( Output( _name="norm_L2_u",_latex="\\left\\| u \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs3 ) )
    .addOutput( Output( _name="norm_H1_u",_latex="\\left\\| u \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs ) )
    .addOutput( Output( _name="norm_L2_divu",_latex="\\left\\| \\nabla \\cdot u \\right\\|_{L^2}",_dependencies=depend,_funcs=divfuncs ) )
    .addOutput( Output( _name="norm_L2_p",_latex="\\left\\| p \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs2 ) );


    M_lambda = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*M_PI*M_PI );
    penalbc = this->vm()["bccoeff"].template as<value_type>();
    M_beta = this->vm()["beta"].template as<value_type>();

    this->//addParameterValue( mu )
    addParameterValue( convex_type::is_simplex )
    .addParameterValue( Dim )
    .addParameterValue( OrderU )
    .addParameterValue( OrderP )
    .addParameterValue( doption("hsize") );

    boost::timer t;

    /*
     * First we create the mesh : a square [0,1]x[0,1] with characteristic
     * length = meshSize
     */
    LOG(INFO) << "creating mesh with hsize=" << meshSize << "\n";
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name="square",
                                         _shape="hypercube",
                                         _usenames=true,
                                         _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                         _dim=Dim,
                                         _h=meshSize,
                                         _xmin=-0.5,_xmax=1.,
                                         _ymin=-0.5,_ymax=1.5 ) );

    LOG(INFO) << "mesh created in " << t.elapsed() << "s\n";
    t.restart();

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );

    U = Xh->element( "U" );
    V = Xh->element( "V" );

    // project in V the exact solution (v,q) are just views over V
    double pi = M_PI;
    auto u1 = val( 1. - exp( M_lambda * Px() ) * cos( 2.*pi*Py() ) );
    auto u2 = val( ( M_lambda/( 2.*pi ) ) * exp( M_lambda * Px() ) * sin( 2.*pi*Py() ) );
    auto u_exact = vec( u1,u2 );
    auto p_exact = val( ( 1-exp( 2.*M_lambda*Px() ) )/2.0 );

    auto v = V.template element<0>();
    auto q = V.template element<1>();
    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << meshSize << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";
    LOG(INFO) << "Xh and elements created in " << t.elapsed() << "s\n";
    t.restart();
    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    F = vector_ptrtype( M_backend->newVector( Xh ) );
    D = sparse_matrix_ptrtype(  M_backend->newMatrix( Xh, Xh ) );
    LOG(INFO) << "D and F allocated in " << t.elapsed() << "s\n";
    t.restart();
}


template<int _OrderU, int _OrderP, typename Entity>
template<typename element_1_type, typename PressureStabExpr>
void
Kovasznay<_OrderU, _OrderP, Entity>::addPressureStabilisation( element_1_type& p,
        element_1_type& q,
        PressureStabExpr& p_stabexpr )
{

    if ( is_equal_order && boption("stab-p") )
    {
        boost::timer t;
        LOG(INFO) << "[assembly] add stabilisation terms for equal order approximation ( orderU="
              << OrderU << ", orderP=" << OrderP << " )\n";
        size_type pattern = Pattern::COUPLED|Pattern::EXTENDED;
        form2( Xh, Xh, D, _pattern=pattern )  +=
            integrate( internalfaces( mesh ),
                       ( p_stabexpr )*( trans( jumpt( gradt( p ) ) )*jump( grad( q ) ) ),
                       _Q<2*OrderP+2>() );
        LOG(INFO) << "[assembly] form2 D equal order stabilisation terms in " << t.elapsed() << "s\n";
        t.restart();
    }

}

template<int _OrderU, int _OrderP, typename Entity>
template<typename element_0_type, typename DivStabExpr>
void
Kovasznay<_OrderU, _OrderP, Entity>::addDivergenceStabilisation( element_0_type& u,
        element_0_type& v,
        DivStabExpr& d_stabexpr )
{
#if 0

    if ( boption("stab-div") )
    {
        boost::timer t;
        LOG(INFO) << "[assembly] add stabilisation terms for divergence ( orderU="
              << OrderU << ", orderP=" << OrderP << " )\n";
        size_type pattern = Pattern::COUPLED|Pattern::EXTENDED;
        form2( Xh, Xh, D, _pattern=pattern )  +=
            integrate( internalfaces( mesh ),
                       ( d_stabexpr )*( trans( jumpt( divt( u ) ) )*jump( div( v ) ) ) );
        LOG(INFO) << "[assembly] form2 D divergence stabilisation terms in " << t.elapsed() << "s\n";
        t.restart();
    }

#endif
}

template<int _OrderU, int _OrderP, typename Entity>
void
Kovasznay<_OrderU, _OrderP, Entity>::buildRhs()
{
    boost::timer t;
    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();
#if defined(WITH_LAGRANGE_MULTIPLIERS)
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
#endif

    auto deft = .5*( gradt( u )+trans( gradt( u ) ) );
    auto def = .5*( grad( v )+trans( grad( v ) ) );

    // total stress tensor (trial)
    auto SigmaNt = ( -idt( p )*N()+2*mu*deft*N() );
    auto SigmaN = ( -id( p )*N()+2*mu*def*N() );

    //
    // the Kovasznay flow (2D)
    //
    // total stress tensor (test)
    double pi = M_PI;
    auto u1 = val( 1. - exp( M_lambda * Px() ) * cos( 2.*pi*Py() ) );
    auto u2 = val( ( M_lambda/( 2.*pi ) ) * exp( M_lambda * Px() ) * sin( 2.*pi*Py() ) );
    auto u_exact = vec( u1,u2 );

    auto du_dx = val( -M_lambda*exp( M_lambda * Px() )*cos( 2.*pi*Py() ) );
    auto du_dy = val( 2*pi*exp( M_lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dx = val( ( M_lambda*M_lambda/( 2*pi ) )*exp( M_lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dy = val( M_lambda*exp( M_lambda * Px() )*cos( 2.*pi*Py() ) );
    auto grad_exact = ( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = du_dx + dv_dy;

    auto beta = vec( cst( M_beta ),cst( M_beta ) );
    auto convection = grad_exact*beta;

    auto p_exact = val( ( 1-exp( 2.*M_lambda*Px() ) )/2.0 );

    auto f1 = val( exp( M_lambda * Px() )*( ( M_lambda*M_lambda - 4.*pi*pi )*mu*cos( 2.*pi*Py() ) - M_lambda*exp( M_lambda * Px() ) ) );
    auto f2 = val( exp( M_lambda * Px() )*mu*( M_lambda/( 2.*pi ) )*sin( 2.*pi*Py() )*( -M_lambda*M_lambda +4*pi*pi ) );

    auto f = vec( f1,f2 )+ convection;

    double pmean = integrate( elements( mesh ), p_exact ).evaluate()( 0, 0 )/mesh->measure();

    form1( Xh, F, _init=true ) =integrate( elements( mesh ), trans( f )*id( v ) );
#if defined(WITH_LAGRANGE_MULTIPLIERS)
    // impose pmean mean value for the solution through the lagrange multipliers
    form1( Xh, F ) += integrate( elements( mesh ), pmean*id( nu ) );
#endif

    if ( M_weak_dirichlet )
    {
        form1( Xh, F )+= integrate( boundaryfaces( mesh ),
                                    trans( u_exact )*( -SigmaN+
                                            penalbc*id( v )/hFace() +
                                            penalbc*( trans( id( v ) )*N() )*N()*
                                            max( M_beta,mu/hFace() ) ) );
        //max(sqrt(trans(beta)*beta),mu/hFace()) ) );
    }

    LOG(INFO) << "[stokes] vector local assembly done\n";
    LOG(INFO) << "form1 F created in " << t.elapsed() << "s\n";
    t.restart();

}
template<int _OrderU, int _OrderP, typename Entity>
void
Kovasznay<_OrderU, _OrderP, Entity>::buildLhs()
{
    boost::timer t;

    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();
#if defined(WITH_LAGRANGE_MULTIPLIERS)
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
#endif


    auto deft = .5*( gradt( u )+trans( gradt( u ) ) );
    auto def = .5*( grad( v )+trans( grad( v ) ) );

    // total stress tensor (trial)
    auto SigmaNt = ( -idt( p )*N()+2*mu*deft*N() );
    auto SigmaN = ( -id( p )*N()+2*mu*def*N() );

    double pi = M_PI;
    auto u1 = val( 1. - exp( M_lambda * Px() ) * cos( 2.*pi*Py() ) );
    auto u2 = val( ( M_lambda/( 2.*pi ) ) * exp( M_lambda * Px() ) * sin( 2.*pi*Py() ) );
    auto u_exact = vec( u1,u2 );

    auto du_dx = val( -M_lambda*exp( M_lambda * Px() )*cos( 2.*pi*Py() ) );
    auto du_dy = val( 2*pi*exp( M_lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dx = val( ( M_lambda*M_lambda/( 2*pi ) )*exp( M_lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dy = val( M_lambda*exp( M_lambda * Px() )*cos( 2.*pi*Py() ) );
    auto grad_exact = ( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = du_dx + dv_dy;

    auto beta = vec( cst( M_beta ),cst( M_beta ) );
    auto convection = grad_exact*beta;

    auto p_exact = val( ( 1-exp( 2.*M_lambda*Px() ) )/2.0 );

    auto f1 = val( exp( M_lambda * Px() )*( ( M_lambda*M_lambda - 4.*pi*pi )*mu*cos( 2.*pi*Py() ) - M_lambda*exp( M_lambda * Px() ) ) );
    auto f2 = val( exp( M_lambda * Px() )*mu*( M_lambda/( 2.*pi ) )*sin( 2.*pi*Py() )*( -M_lambda*M_lambda +4*pi*pi ) );

    auto f = vec( f1,f2 )+ convection;

    size_type pattern = Pattern::COUPLED;

    if ( ( is_equal_order &&
            boption("stab-p") ) ||
            boption("stab-div") )
        pattern |= Pattern::EXTENDED;

    Feel::Context graph( pattern );
    LOG(INFO) << "[stokes] test : " << ( graph.test ( Pattern::DEFAULT ) || graph.test ( Pattern::EXTENDED ) ) << "\n";
    LOG(INFO) << "[stokes]  : graph.test ( Pattern::DEFAULT )=" <<  graph.test ( Pattern::DEFAULT ) << "\n";
    LOG(INFO) << "[stokes]  : graph.test ( Pattern::COUPLED )=" <<  graph.test ( Pattern::COUPLED ) << "\n";
    LOG(INFO) << "[stokes]  : graph.test ( Pattern::EXTENDED)=" <<  graph.test ( Pattern::EXTENDED ) << "\n";
    LOG(INFO) << "[assembly] add diffusion terms\n";
    form2( Xh, Xh, D, _init=true, _pattern=pattern );
    LOG(INFO) << "[assembly] form2 D init in " << t.elapsed() << "s\n";
    t.restart();
    form2( Xh, Xh, D )+= integrate( elements( mesh ), 2*mu*trace( deft*trans( def ) ) + trans( gradt( u )*beta )*id( v ) );
    LOG(INFO) << "[assembly] form2 D convection and viscous terms in " << t.elapsed() << "s\n";
    t.restart();
    LOG(INFO) << "[assembly] add velocity/pressure terms\n";
    form2( Xh, Xh, D )+=integrate( elements( mesh ),- div( v )*idt( p ) + divt( u )*id( q ) );
    LOG(INFO) << "[assembly] form2 D velocity/pressure terms in " << t.elapsed() << "s\n";
    t.restart();
    LOG(INFO) << "[assembly] add lagrange multipliers terms for zero mean pressure\n";
#if defined(WITH_LAGRANGE_MULTIPLIERS)
    form2( Xh, Xh, D )+=integrate( elements( mesh ), id( q )*idt( lambda ) + idt( p )*id( nu ) );
    LOG(INFO) << "[assembly] form2 D pressure/multipliers terms in " << t.elapsed() << "s\n";
    t.restart();
#endif

    if ( M_weak_dirichlet )
    {
        LOG(INFO) << "[assembly] add terms for weak Dirichlet condition handling\n";
        form2( Xh, Xh, D )+=integrate( boundaryfaces( mesh ), -trans( SigmaNt )*id( v ) );
        form2( Xh, Xh, D )+=integrate( boundaryfaces( mesh ), -trans( SigmaN )*idt( u ) );
        form2( Xh, Xh, D )+=integrate( boundaryfaces( mesh ), +penalbc*trans( idt( u ) )*id( v )/hFace() );
        form2( Xh, Xh, D )+=integrate( boundaryfaces( mesh ), +penalbc*( trans( idt( u ) )*N() )*( trans( id( v ) )*N() )*max( M_beta,mu/hFace() ) );
        LOG(INFO) << "[assembly] form2 D boundary terms in " << t.elapsed() << "s\n";
        t.restart();
    }

    if ( math::abs( M_beta ) < 1e-10 )
    {
        double pterm = math::pow( double( OrderU ), 7./2. );
        auto p_stabexpr = M_stabP*hFace()*hFace()*hFace()/( pterm );
        this->addPressureStabilisation( p, q, p_stabexpr  );
    }

    else
    {
        double pterm = math::pow( double( OrderU ), 4. );
        auto p_stabexpr = M_stabP*hFace()*hFace()*hFace()/max( mu*pterm,hFace()*M_beta );
        this->addPressureStabilisation( p, q, p_stabexpr  );
    }

    //auto d_stabexpr = M_stabD*hFace()*hFace()/math::pow(double(OrderU), 7./2.);
    //this->addDivergenceStabilisation( u, v, d_stabexpr  );

    LOG(INFO) << "[stokes] matrix local assembly done\n";

    D->close();
    F->close();

    if ( !M_weak_dirichlet )
    {
        form2( Xh, Xh, D ) += on( boundaryfaces( mesh ), u, F, u_exact );
    }

    LOG(INFO) << "[stokes] vector/matrix global assembly done\n";
    LOG(INFO) << "form2 D created in " << t.elapsed() << "s\n";
    t.restart();

}
template<int _OrderU, int _OrderP, typename Entity>
void
Kovasznay<_OrderU, _OrderP, Entity>::run()
{
    if ( this->preProcessing() == RUN_EXIT ) return;

    using namespace Feel::vf;

    this->buildRhs();
    this->buildLhs();

    if ( this->vm().count( "export-matlab" ) )
    {
        boost::timer t;
        D->printMatlab( "D.m" );
        F->printMatlab( "F.m" );
        LOG(INFO) << "system saved in " << t.elapsed() << "s\n";
        t.restart();
    }

    boost::timer t;
    M_backend->solve( _matrix=D, _solution=U, _rhs=F, _rtolerance=1e-14 );

    LOG(INFO) << "system solved in " << t.elapsed() << "s\n";
    t.restart();

    this->exportResults( U, V );
    LOG(INFO) << "postprocessing done in " << t.elapsed() << "s\n";
    t.restart();

    this->postProcessing();

} // Kovasznay::run

template<int _OrderU, int _OrderP, typename Entity>
void
Kovasznay<_OrderU, _OrderP, Entity>::exportResults( element_type& U, element_type& V )
{
    const double pi = M_PI;
    auto u1 = val( 1. - exp( M_lambda * Px() ) * cos( 2.*pi*Py() ) );
    auto u2 = val( ( M_lambda/( 2.*pi ) ) * exp( M_lambda * Px() ) * sin( 2.*pi*Py() ) );
    auto u_exact = vec( u1,u2 );

    auto du_dx = val( -M_lambda*exp( M_lambda * Px() )*cos( 2.*pi*Py() ) );
    auto du_dy = val( 2*pi*exp( M_lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dx = val( ( M_lambda*M_lambda/( 2*pi ) )*exp( M_lambda * Px() )*sin( 2.*pi*Py() ) );
    auto dv_dy = val( M_lambda*exp( M_lambda * Px() )*cos( 2.*pi*Py() ) );
    auto grad_exact = ( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) );
    auto div_exact = du_dx + dv_dy;

    auto p_exact = val( ( 1-exp( 2.*M_lambda*Px() ) )/2.0 );

    auto u = U.template element<0>();
    auto p = U.template element<1>();
    double u_error_L2_2 = integrate( elements( mesh ),
                                     trans( idv( u )-u_exact )*( idv( u )-u_exact ) ).evaluate()( 0, 0 );
    double uex_L2_2 = integrate( elements( mesh ), trans( u_exact )*( u_exact ) ).evaluate()( 0, 0 );
    double u_errorL2 = math::sqrt( u_error_L2_2/uex_L2_2 );
    std::cout << "||u_error||_0/||uex||_0 = " << u_errorL2 << "\n";;


    double u_errorsemiH1 = integrate( elements( mesh ),
                                      trace( ( gradv( u )-grad_exact )*trans( gradv( u )-grad_exact ) ) ).evaluate()( 0, 0 );
    double u_error_H1 = math::sqrt( u_error_L2_2+u_errorsemiH1 );
    double uex_semiH1_2 = integrate( elements( mesh ), trace( ( grad_exact )*trans( grad_exact ) ) ).evaluate()( 0, 0 );
    double uex_H1 = math::sqrt( uex_L2_2+uex_semiH1_2 );
    double u_errorH1 = u_error_H1/uex_H1;
    std::cout << "||u_error||_1/||uex||_1 = " << u_errorH1 << "\n";


    double mean_p = integrate( elements( mesh ), idv( p ) ).evaluate()( 0, 0 )/mesh->measure();
    double mean_pexact = integrate( elements( mesh ), p_exact ).evaluate()( 0, 0 )/mesh->measure();
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";
    LOG(INFO) << "[stokes] mean(p_exact)=" << mean_pexact << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p_exact)=" << mean_pexact << "\n";


    LOG(INFO) << "[stokes] measure(Omega)=" << mesh->measure() << " (should be equal to 3)\n";
    std::cout << "[stokes] measure(Omega)=" << mesh->measure() << " (should be equal to 3)\n";

    double p_errorL2_2 = integrate( elements( mesh ),
                                    ( ( idv( p )-mean_p )-( p_exact-mean_pexact ) )*( ( idv( p )-mean_p )-( p_exact-mean_pexact ) ) ).evaluate()( 0, 0 );
    double pex_L2 = integrate( elements( mesh ), ( p_exact-mean_pexact )*( p_exact-mean_pexact ) ).evaluate()( 0, 0 );
    double p_errorL2 = math::sqrt( p_errorL2_2/pex_L2 );
    std::cout << "||p_error||_0/||pex||_0 = " <<  p_errorL2 << "\n";;

    LOG(INFO) << "[stokes] solve for D done\n";



    double mean_div_u = integrate( elements( mesh ), divv( u ) ).evaluate()( 0, 0 );
    LOG(INFO) << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_errorL2_2 = integrate( elements( mesh ), divv( u )*divv( u ) ).evaluate()( 0, 0 );
    double uex_div_L2 = integrate( elements( mesh ), div_exact ).evaluate()( 0, 0 );
    double uex_n_L2 = integrate( boundaryfaces( mesh ), trans( u_exact )*N() ).evaluate()( 0, 0 );
    std::cout << "[stokes] ||div(uexact)||=" << uex_div_L2 << "\n";
    std::cout << "[stokes] ||uexact,n||=" << uex_n_L2 << "\n";
    double div_u_errorL2 = math::sqrt( div_u_errorL2_2 );
    LOG(INFO) << "[stokes] ||div(u)||_2=" << div_u_errorL2 << "\n";
    std::cout << "[stokes] ||div(u)||=" << div_u_errorL2 << "\n";

    this->addOutputValue( u_errorL2 )
    .addOutputValue( u_errorH1 )
    .addOutputValue( div_u_errorL2 )
    .addOutputValue( p_errorL2 );

    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", U.template element<0>() );
        exporter->step( 0 )->add( "p", U.template element<1>() );
        exporter->step( 0 )->add( "u_exact", V.template element<0>() );
        exporter->step( 0 )->add( "p_exact", V.template element<1>() );
        exporter->save();
    }
} // Kovasznay::export

template<int _OrderU,int _OrderP, typename Entity> const uint16_type Kovasznay<_OrderU, _OrderP, Entity>::Dim;
template<int _OrderU,int _OrderP, typename Entity> const uint16_type Kovasznay<_OrderU, _OrderP, Entity>::OrderU;
template<int _OrderU,int _OrderP, typename Entity> const uint16_type Kovasznay<_OrderU, _OrderP, Entity>::OrderP;
} // Feel
