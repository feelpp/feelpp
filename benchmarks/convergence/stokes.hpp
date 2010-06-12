/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>

#include <life/lifepoly/im.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>



#include <life/lifemesh/elements.hpp>

#include <life/lifevf/vf.hpp>
#include <fstream>
#include <sstream>

#include <life/lifecore/applicationxml.hpp>
#include <life/lifecore/xmlparser.hpp>


/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Life::Application subclass.
 *
 * \return the list of options
 */
inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description stokesoptions("Stokes options");
    stokesoptions.add_options()
        ("penal", Life::po::value<double>()->default_value( 0.5 ), "penalisation parameter")
        ("mu", Life::po::value<double>()->default_value( 1.0/40. ), "viscosity coefficient (default from Sherwin/Karnyadakis book)")
        ("beta", Life::po::value<double>()->default_value( 0.0 ), "convection coefficient")
        ("hsize", Life::po::value<double>()->default_value( 0.1 ), "first h value to start convergence")
        ("bctype", Life::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("bccoeff", Life::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
        ("penalisation", Life::po::value<double>()->default_value( 1e-3 ), "penalisation parameter for equal order approximation")
        ("stab", Life::po::value<bool>()->default_value( true ), "0 = no stabilisation for equal order approx., 1 = stabilisation for equal order approx.")

        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return stokesoptions.add( Life::life_options() ) ;
}


/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Life::Application subclass.
 *
 * \return some data about the application.
 */
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "stokes" ,
                           "stokes" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2009 Universite de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
   return about;

}


namespace Life
{
using namespace vf;
/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
template<int Dim,
         int _OrderU,
         int _OrderP = _OrderU-1,
         template<uint16_type,uint16_type,uint16_type> class Entity = Simplex>
class Stokes
    :
        public ApplicationXML
{
    typedef ApplicationXML super;
public:

    static const uint16_type OrderU = _OrderU;
    static const uint16_type OrderP = _OrderP;
    static const bool is_equal_order = (OrderU==OrderP);

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef Lagrange<OrderU, Vectorial> basis_u_type;
    typedef Lagrange<OrderP, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    BOOST_MPL_ASSERT( ( boost::is_same<typename space_type::bases_list, basis_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<0> >::type, basis_u_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<1> >::type, basis_p_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<2> >::type, basis_l_type> ) );
    typedef boost::shared_ptr<space_type> space_ptrtype;

    /* functions */
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;
    typedef typename element_type::template sub_element<2>::type element_2_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    /**
     * run the convergence test
     */
    void run();

private:


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    double meshSize;
    double M_stabP;
    space_ptrtype Xh;
    sparse_matrix_ptrtype D;
    double mu;
    double M_lambda;
    double penalbc;
    double M_beta;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;
    mesh_ptrtype mesh;
private:

    void addStabilisation(element_1_type& p, element_1_type& q );

}; // Stokes

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
Stokes<Dim, _OrderU, _OrderP, Entity>::Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    M_stabP( this->vm()["penalisation"].template as<double>()/math::pow(double(OrderU), 7./2.) ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

    mu = this->vm()["mu"].template as<value_type>();
    Parameter h;
    switch( OrderU )
    {
    case 1:
        h = Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.02:0.025:0.05" );
        break;
    case 2:
        h = Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.02:0.025:0.1" );
        break;

    case 3:
        h = Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.035:0.025:0.2" );
        break;
    }
    this->
        //addParameter( Parameter(_name="mu",_type=CONT_ATTR,_latex="\\mu", _values=boost::lexical_cast<std::string>( mu ).c_str()))
        addParameter( Parameter(_name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str()) )
        .addParameter( Parameter(_name="orderU",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderU  ).c_str()) )
        .addParameter( Parameter(_name="orderP",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderP  ).c_str()) )
        .addParameter( h );

    std::vector<Parameter> depend;
    std::vector<std::string> funcs;
    depend.push_back(h);
    std::ostringstream oss;
    if ( OrderP == OrderU )
        oss << "h**" << boost::lexical_cast<std::string>( OrderP  );
    else
        oss << "h**" << boost::lexical_cast<std::string>( OrderP+1  );
    funcs.push_back(oss.str());
    oss.str("");
    std::vector<std::string> funcs2;
    oss << "h**" << boost::lexical_cast<std::string>( OrderP+1 ) ;
    funcs2.push_back(oss.str());

    this->
        addOutput( Output(_name="norm_H1_u",_latex="\\left\\| u \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs) )
        .addOutput( Output(_name="norm_L2_p",_latex="\\left\\| p \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs2) );


    M_lambda = 1./(2.*mu) - math::sqrt( 1./(4.*mu*mu) + 4.*M_PI*M_PI);
    penalbc = this->vm()["bccoeff"].template as<value_type>();
    M_beta = this->vm()["beta"].template as<value_type>();
}

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::addStabilisation( element_1_type& p,
                                                         element_1_type& q )
{

    if ( is_equal_order && this->vm()["stab"].template as<bool>() )
    {
        boost::timer t;
        Log() << "[assembly] add stabilisation terms for equal order approximation ( orderU="
              << OrderU << ", orderP=" << OrderP << " )\n";
        size_type pattern = DOF_PATTERN_COUPLED|DOF_PATTERN_NEIGHBOR;
        form2( Xh, Xh, D, _pattern=pattern )  +=
            integrate( internalfaces(mesh),
                       M_stabP*hFace()*hFace()*(trans(jumpt(gradt(p)))*jump(grad(q))) );
        Log() << "[assembly] form2 D stabilisation terms in " << t.elapsed() << "s\n"; t.restart();
    }

}

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::run()
{
    this->//addParameterValue( mu )
        addParameterValue( Dim )
        .addParameterValue( OrderU )
        .addParameterValue( OrderP )
        .addParameterValue( this->vm()["hsize"].template as<double>() );

    if (this->preProcessing() == RUN_EXIT) return;

    using namespace Life::vf;

    boost::timer t;

    /*
     * First we create the mesh : a square [0,1]x[0,1] with characteristic
     * length = meshSize
     */
    Log() << "creating mesh with hsize=" << meshSize << "\n";
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name="square",
                                         _shape="hypercube",
                                         _dim=Dim,
                                         _h=meshSize,
                                         _xmin=-0.5,_xmax=1.,
                                         _ymin=-0.5,_ymax=1.5 ) );

    Log() << "mesh created in " << t.elapsed() << "s\n"; t.restart();

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type E( Xh, "u" );
    element_type V( Xh, "v" );
    element_0_type ue = E.template element<0>();
    element_1_type pe = E.template element<1>();
    element_2_type le = E.template element<2>();
    element_0_type u = U.template element<0>();
    element_0_type v = V.template element<0>();
    element_1_type p = U.template element<1>();
    element_1_type q = V.template element<1>();
    element_2_type lambda = U.template element<2>();
    element_2_type nu = V.template element<2>();

    Log() << "Data Summary:\n";
    Log() << "   hsize = " << meshSize << "\n";
    Log() << "  export = " << this->vm().count("export") << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";
    Log() << "functionspace and elements created in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

    vector_ptrtype F( M_backend->newVector( Xh ) );

    auto deft = gradt(u);
    auto def = grad(v);

    // total stress tensor (trial)
    auto SigmaNt = (-idt(p)*N()+mu*deft*N());
    auto SigmaN = (-id(p)*N()+mu*def*N());
#if 1 // the Kovasznay flow (2D)
    // total stress tensor (test)

    double pi = 4*math::atan(1.0);
    AUTO( u1, val(1. - exp( M_lambda * Px() ) * cos(2.*pi*Py())) );
    AUTO( u2, val((M_lambda/(2.*pi)) * exp( M_lambda * Px() ) * sin(2.*pi*Py())) );

    AUTO( u_exact, u1*oneX() + u2*oneY());

    AUTO( du_dx, val(-M_lambda*exp( M_lambda * Px() )*cos(2.*pi*Py())));
    AUTO( du_dy, val(2*pi*exp( M_lambda * Px() )*sin(2.*pi*Py())));
    AUTO( dv_dx, val((M_lambda*M_lambda/(2*pi))*exp( M_lambda * Px() )*sin(2.*pi*Py())));
    AUTO( dv_dy, val(M_lambda*exp( M_lambda * Px() )*cos(2.*pi*Py())));

    AUTO( grad_exact, (mat<2,2>(du_dx, du_dy, dv_dx, dv_dy)) );

    AUTO( beta, M_beta*(oneX() + oneY()) );

    AUTO( convection, grad_exact*beta);

    AUTO( p_exact_kov, val((1-exp(2.*M_lambda*Px()))/2.0) );
    double pmeas = integrate( elements(mesh), constant(1.) ).evaluate()( 0, 0 );
    double pmean = integrate( elements(mesh), p_exact_kov ).evaluate()( 0, 0 )/pmeas;
    AUTO( p_exact, p_exact_kov-pmean );

    AUTO( f1, val(exp( M_lambda * Px() )*((M_lambda*M_lambda - 4.*pi*pi)*mu*cos(2.*pi*Py()) - M_lambda*exp( M_lambda * Px() ))) );
    AUTO( f2, val(exp( M_lambda * Px() )*mu*(M_lambda/(2.*pi))*sin(2.*pi*Py())*(-M_lambda*M_lambda +4*pi*pi)) );

    AUTO( f, f1*oneX() + f2*oneY() + convection );
#else
    // u exact solution
    AUTO( u_exact, vec(cos(Px())*cos(Py()),sin(Px())*sin(Py()) ) );


    // this is the exact solution which has zero mean : the mean of
    // cos(x)*sin(y) is sin(1)*(1-cos(1))) on [0,1]^2
    AUTO( p_exact, cos(Px())*sin(Py())-(sin(1.0)*(1.-cos(1.0)) ) );

    // f is such that f = \Delta u_exact + \nabla p_exact
    AUTO( f, vec(2*cos(Px())*cos(Py())-sin(Px())*sin(Py()),
                 2*sin(Px())*sin(Py())+cos(Px())*cos(Py()) ) );
#endif
    v = vf::project( v.functionSpace(), elements(mesh), beta );
    ue = vf::project( ue.functionSpace(), elements(mesh), u_exact );
    Log() << "convection terms projectd in " << t.elapsed() << "s\n"; t.restart();

    // right hand side
    form1( Xh, F, _init=true ) =integrate( elements(mesh), trans(f)*id(v) );
    form1( Xh, F )+= integrate( boundaryfaces(mesh),
                                trans(u_exact)*(-SigmaN+penalbc*id(v)/hFace() ) );

    Log() << "[stokes] vector local assembly done\n";
    Log() << "form1 F created in " << t.elapsed() << "s\n"; t.restart();
    /*
     * Construction of the left hand side
     */
    D = sparse_matrix_ptrtype(  M_backend->newMatrix( Xh, Xh ) );
    size_type pattern = DOF_PATTERN_COUPLED;
    if ( is_equal_order && this->vm()["stab"].template as<bool>() )
        pattern |= DOF_PATTERN_NEIGHBOR;
    Life::Context graph( pattern );
    Log() << "[stokes] test : " << ( graph.test ( DOF_PATTERN_DEFAULT ) || graph.test ( DOF_PATTERN_NEIGHBOR ) ) << "\n";
    Log() << "[stokes]  : graph.test ( DOF_PATTERN_DEFAULT )=" <<  graph.test ( DOF_PATTERN_DEFAULT ) << "\n";
    Log() << "[stokes]  : graph.test ( DOF_PATTERN_COUPLED )=" <<  graph.test ( DOF_PATTERN_COUPLED ) << "\n";
    Log() << "[stokes]  : graph.test ( DOF_PATTERN_NEIGHBOR)=" <<  graph.test ( DOF_PATTERN_NEIGHBOR ) << "\n";
    Log() << "[assembly] add diffusion terms\n";
    form2( Xh, Xh, D, _init=true, _pattern=pattern );
    Log() << "[assembly] form2 D init in " << t.elapsed() << "s\n"; t.restart();
    form2( Xh, Xh, D )+= integrate( elements(mesh), mu*trace(deft*trans(def)) + trans(gradt(u)*idv(v))*id(v) );
    Log() << "[assembly] form2 D convection and viscous terms in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[assembly] add velocity/pressure terms\n";
    form2( Xh, Xh, D )+=integrate( elements(mesh),- div(v)*idt(p) + divt(u)*id(q) );
    Log() << "[assembly] form2 D velocity/pressure terms in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[assembly] add lagrange multipliers terms for zero mean pressure\n";
    form2( Xh, Xh, D )+=integrate( elements(mesh), id(q)*idt(lambda) + idt(p)*id(nu) );
    Log() << "[assembly] form2 D pressure/multipliers terms in " << t.elapsed() << "s\n"; t.restart();
    Log() << "[assembly] add terms for weak Dirichlet condition handling\n";
    form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), -trans(SigmaNt)*id(v) );
    form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), -trans(SigmaN)*idt(u) );
    form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), +penalbc*trans(idt(u))*id(v)/hFace() );
    Log() << "[assembly] form2 D boundary terms in " << t.elapsed() << "s\n"; t.restart();

    this->addStabilisation( p, q  );

    Log() << "[stokes] matrix local assembly done\n";

    D->close();
    F->close();
    Log() << "[stokes] vector/matrix global assembly done\n";
    Log() << "form2 D created in " << t.elapsed() << "s\n"; t.restart();

    if( this->vm().count( "export-matlab" ) )
    {
        D->printMatlab( "S.m" );
        F->printMatlab( "F.m" );
    }

    vector_ptrtype X( M_backend->newVector( U.functionSpace() ) );
    M_backend->solve( _matrix=D, _solution=X, _rhs=F, _rtolerance=1e-16 );
    U = *X;

    Log() << "system solved in " << t.elapsed() << "s\n"; t.restart();
    Log() << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";
    std::cout << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";

    v = vf::project( u.functionSpace(), elements(mesh), idv(u)-u_exact );
    double u_errorL2 = integrate( elements(mesh),
                                  trans(idv(u)-u_exact)*(idv(u)-u_exact) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;


    double u_errorsemiH1 = integrate( elements(mesh),
                                      trace((gradv(u)-grad_exact)*trans(gradv(u)-grad_exact))).evaluate()( 0, 0 );
    double u_error_H1 = math::sqrt( u_errorL2+u_errorsemiH1 );
    std::cout << "||u_error||_2 = " << u_error_H1 << "\n";


    double p_errorL2 = integrate( elements(mesh), (idv(p)-p_exact)*(idv(p)-p_exact) ).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";;

    Log() << "[stokes] solve for D done\n";

    double meas = integrate( elements(mesh), constant(1.0) ).evaluate()( 0, 0);
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";

    double mean_p = integrate( elements(mesh), idv(p) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    double mean_div_u = integrate( elements(mesh), divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements(mesh), divv(u)*divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";
    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );
    Log() << "postprocessing done in " << t.elapsed() << "s\n"; t.restart();
    this->exportResults( U, V );
    Log() << "exporting done in " << t.elapsed() << "s\n"; t.restart();


    this->addOutputValue( u_error_H1 ).addOutputValue( math::sqrt( p_errorL2 ) );
    this->postProcessing();

} // Stokes::run

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::exportResults( element_type& U, element_type& V )
{
    if ( exporter->doExport() )
    {
        exporter->step(0)->setMesh( U.functionSpace()->mesh() );
        exporter->step(0)->add( "u", U.template element<0>() );
        exporter->step(0)->add( "p", U.template element<1>() );
        exporter->step(0)->add( "u_exact", V.template element<0>() );
        exporter->step(0)->add( "p_exact", V.template element<1>() );
        exporter->save();
    }
} // Stokes::export





template<int Dim,int _OrderU,int _OrderP,template<uint16_type,uint16_type,uint16_type> class Entity>
const uint16_type Stokes<Dim, _OrderU, _OrderP, Entity>::OrderU;
template<int Dim,int _OrderU,int _OrderP,template<uint16_type,uint16_type,uint16_type> class Entity>
const uint16_type Stokes<Dim, _OrderU, _OrderP, Entity>::OrderP;
} // Life
