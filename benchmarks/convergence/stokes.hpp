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
        ("penalisation", Life::po::value<double>()->default_value( 1.0 ), "penalisation parameter for equal order approximation")
        ("stab", Life::po::value<bool>()->default_value( true ), "0 = no stabilisation for equal order approx., 1 = stabilisation for equal order approx.")

        ("export", "export results(ensight, data file(1D)")
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
#if 0
    BOOST_MPL_ASSERT_MSG(
        ((OrderU == OrderP-1) || (OrderU == OrderP-2) || (OrderU==OrderP))
        , ORDER_VELOCITY_AND_ORDER_PRESSURE_INCOMPATIBLE
        , (mpl::int_<OrderU>,mpl::int_<OrderP>)
        );
#endif // 0
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
#if 0
    typedef Product<Expansion<Continuous,basis_u_type>,
                    Expansion<Continuous,basis_p_type>,
                    Expansion<Continuous,basis_l_type> > basis_type;
#endif
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

    Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {

        Parameter h(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.04:0.08:0.2" );
        this->
            addParameter( Parameter(_name="mu",_type=CONT_ATTR,_latex="\\mu",_values="0.01:1:10") )
            .addParameter( Parameter(_name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str()) )
            .addParameter( Parameter(_name="orderU",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderU  ).c_str()) )
            .addParameter( Parameter(_name="orderP",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( OrderP  ).c_str()) )
            .addParameter( h );

        std::vector<Parameter> depend;
        std::vector<std::string> funcs;
        depend.push_back(h);
        std::ostringstream oss;
        oss << "h**" << boost::lexical_cast<std::string>( OrderU + 1  ) ;
        funcs.push_back(oss.str());
        oss.str("");
        std::vector<std::string> funcs2;
        oss << "h**" << boost::lexical_cast<std::string>( OrderU ) ;
        funcs2.push_back(oss.str());

        this->
            addOutput( Output(_name="norm_L2_u",_latex="\\left\\| u \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs) )
            .addOutput( Output(_name="norm_L2_p",_latex="\\left\\| p \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs2) );

        mu = this->vm()["mu"].template as<value_type>();
        M_lambda = 1./(2.*mu) - math::sqrt( 1./(4.*mu*mu) + 4.*M_PI*M_PI);
        penalbc = this->vm()["bccoeff"].template as<value_type>();
        M_beta = this->vm()["beta"].template as<value_type>();
    }


    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve system
     */
    void solve( sparse_matrix_ptrtype const& D, element_type& u, vector_ptrtype const& F, bool is_sym );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double M_lambda;
    double penalbc;
    double M_beta;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;
}; // Stokes

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
    Stokes<Dim, _OrderU, _OrderP, Entity>::run()
{
    this->addParameterValue( this->vm()["mu"].template as<double>() )
        .addParameterValue( Dim )
        .addParameterValue( OrderU )
        .addParameterValue( OrderP )
        .addParameterValue( this->vm()["hsize"].template as<double>() );

    if (this->preProcessing() == RUN_EXIT) return;

    using namespace Life::vf;

    /*
     * First we create the mesh : a square [0,1]x[0,1] with characteristic
     * length = meshSize
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="square",
                                                      _shape="hypercube",
                                                      _dim=Dim,
                                                      _h=meshSize,
                                                      _xmin=-0.5,_xmax=1.,
                                                      _ymin=-0.5,_ymax=1.5 ) );

    /*
     * The function space and some associate elements are then defined
     */
    space_ptrtype Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
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


    vector_ptrtype F( M_backend->newVector( Xh ) );

#if 0
    // viscous stress tensor (trial) : 0.5 ( \nabla u + \nabla u ^T )
    AUTO( deft, 0.5*(gradt(u)+trans(gradt(u)) ));
    /x/ viscous stress tensor (test) : 0.5 ( \nabla u + \nabla u ^T )
    AUTO( def, 0.5*(grad(v)+trans(grad(v))) );
#else
    AUTO( deft, gradt(u) );
    AUTO( def, grad(v) );
#endif
    // total stress tensor (trial)
    AUTO( SigmaNt, (-idt(p)*N()+mu*deft*N()) );
#if 1 // the Kovasznay flow (2D)
    // total stress tensor (test)
    AUTO( SigmaN, (-id(p)*N()+mu*def*N()) );
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
    double pmeas = integrate( elements(mesh), _Q<OrderU>(), constant(1.) ).evaluate()( 0, 0 );
    double pmean = integrate( elements(mesh), _Q<OrderU>(), p_exact_kov ).evaluate()( 0, 0 )/pmeas;
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
    // right hand side
    form1( Xh, F, _init=true )  =
        integrate( elements(mesh), _Q<OrderU+5>(), trans(f)*id(v) )+
        integrate( boundaryfaces(mesh),
                   // higher order quadrature to accurately integrate u_exact
                   _Q<OrderU+1>(),
                   trans(u_exact)*(-SigmaN+penalbc*id(v)/hFace() ) );

    Log() << "[stokes] vector local assembly done\n";

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
    size_type pattern = DOF_PATTERN_COUPLED|DOF_PATTERN_NEIGHBOR;
    Log() << "[assembly] add diffusion terms\n";
    form2( Xh, Xh, D, _init=true, _pattern=pattern )=
        integrate( elements(mesh), _Q<2*(OrderU-1)>(),
                   mu*trace(deft*trans(def)) +
                   trans(gradt(u)*beta)*id(v) );
    Log() << "[assembly] add velocity/pressure terms\n";
    form2( Xh, Xh, D )+=integrate( elements(mesh), _Q<2*(OrderU-1)>(),
                                   - div(v)*idt(p) + divt(u)*id(q) );
    Log() << "[assembly] add lagrange multipliers terms for zero mean pressure\n";
    form2( Xh, Xh, D )+=integrate( elements(mesh), _Q<OrderU-1>(),
                                   id(q)*idt(lambda) + idt(p)*id(nu) );
     Log() << "[assembly] add terms for weak Dirichlet condition handling\n";
    form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), _Q<2*OrderU>(),
                                   -trans(SigmaNt)*id(v)
                                   -trans(SigmaN)*idt(u)
                                   +penalbc*trans(idt(u))*id(v)/hFace() );

    if ( is_equal_order && this->vm()["stab"].template as<bool>() )
    {
        Log() << "[assembly] add stabilisation terms for equal order approximation ( orderU=" << OrderU << ", orderP=" << OrderP << " )\n";
        double p_term = double(OrderU);
        p_term = math::pow(p_term, 7./2.);
        AUTO( penalisation_term, constant(this->vm()["penalisation"].template as<double>())*hFace()*hFace()/p_term );
        form2( Xh, Xh, D ) += integrate( internalfaces(mesh), _Q<2*(OrderP-1)>(),
                                         penalisation_term*(trans(jumpt(gradt(p)))*jump(grad(q))) );
    }

    Log() << "[stokes] matrix local assembly done\n";
    D->close();
    F->close();
    Log() << "[stokes] vector/matrix global assembly done\n";


    this->solve( D, U, F, false );

    Log() << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";
    std::cout << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";

    double u_errorL2 = integrate( elements(mesh), _Q<2*OrderU>(), trans(idv(u)-u_exact)*(idv(u)-u_exact) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;


    double p_errorL2 = integrate( elements(mesh), _Q<2*OrderU>(), (idv(p)-p_exact)*(idv(p)-p_exact) ).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";;

    Log() << "[stokes] solve for D done\n";

    double meas = integrate( elements(mesh), _Q<0>(), constant(1.0) ).evaluate()( 0, 0);
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";

    double mean_p = integrate( elements(mesh), _Q<OrderP>(), idv(p) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    double mean_div_u = integrate( elements(mesh), _Q<OrderU-1>(), divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements(mesh), _Q<2*(OrderU-1)>(), divv(u)*divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";
    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );

    this->exportResults( U, V );

    this->addOutputValue( math::sqrt( u_errorL2 ) ).addOutputValue( math::sqrt( p_errorL2 ) );
    this->postProcessing();

    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";
} // Stokes::run

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::solve( sparse_matrix_ptrtype const& D,
                                              element_type& u,
                                              vector_ptrtype const& F,
                                              bool is_sym )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F, false );
    u = *U;
} // Stokes::solve

template<int Dim, int _OrderU, int _OrderP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, _OrderU, _OrderP, Entity>::exportResults( element_type& U, element_type& V )
{
    if ( this->vm().count("export") )
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
