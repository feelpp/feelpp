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
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>

#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>



#include <life/lifemesh/elements.hpp>

#include <life/lifevf/vf.hpp>

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
        ("f", Life::po::value<double>()->default_value( 0 ), "forcing term")
        ("mu", Life::po::value<double>()->default_value( 1.0 ), "reaction coefficient component")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("bctype", Life::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("bccoeff", Life::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
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
         int Order,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Stokes
    :
        public Application
{
    typedef Application super;
public:


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
    //# marker1 #
    typedef Lagrange<Order, Vectorial> basis_u_type;
    typedef Lagrange<Order-1, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
    //# endmarker1 #

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    BOOST_MPL_ASSERT( ( boost::is_same<typename space_type::bases_list, basis_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<0> >::type, basis_u_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<1> >::type, basis_p_type> ) );
    BOOST_MPL_ASSERT( ( boost::is_same<typename mpl::at<typename space_type::bases_list,mpl::int_<2> >::type, basis_l_type> ) );

    /* functions */
    //# marker3 #
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;
    typedef typename element_type::template sub_element<2>::type element_2_type;
    //# endmarker3 #

    /* export */
    typedef Exporter<mesh_type> export_type;

    Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
        mu = this->vm()["mu"].template as<value_type>();
        penalbc = this->vm()["bccoeff"].template as<value_type>();
    }


    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );


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
    double penalbc;

    boost::shared_ptr<export_type> exporter;
}; // Stokes

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Stokes<Dim,Order,Entity>::mesh_ptrtype
Stokes<Dim,Order,Entity>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );


    GmshTensorizedDomain<convex_type::nDim,convex_type::nOrder,convex_type::nRealDim,Entity> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( convex_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    return mesh;
} // Stokes::createMesh


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, Order, Entity>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    using namespace Life::vf;

    this->changeRepository( boost::format( "doc/tutorial/%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            );

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );

    /*
     * The function space and some associate elements are then defined
     */
    //# marker4 #
    space_ptrtype Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
    element_0_type u = U.template element<0>();
    element_0_type v = V.template element<0>();
    element_1_type p = U.template element<1>();
    element_1_type q = V.template element<1>();
    element_2_type lambda = U.template element<2>();
    element_2_type nu = V.template element<2>();
    //# endmarker4 #

    Log() << "Data Summary:\n";
    Log() << "   hsize = " << meshSize << "\n";
    Log() << "  export = " << this->vm().count("export") << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";


    vector_ptrtype F( M_backend->newVector( Xh ) );

    //# marker5 #
    AUTO( deft, gradt(u) );
    AUTO( def, grad(v) );
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    AUTO( SigmaNt, (-idt(p)*N()+mu*deft*N()) );

    // total stress tensor (test)
    AUTO( SigmaN, (-id(p)*N()+mu*def*N()) );
    //# endmarker6 #

    // u exact solution
    AUTO( u_exact, vec(cos(Px())*cos(Py()),sin(Px())*sin(Py()) ) );

    // this is the exact solution which has zero mean : the mean of
    // cos(x)*sin(y) is sin(1)*(1-cos(1))) on [0,1]^2
    AUTO( p_exact, cos(Px())*sin(Py())-(sin(1.0)*(1.-cos(1.0)) ) );

    // f is such that f = \Delta u_exact + \nabla p_exact
    AUTO( f, vec(2*cos(Px())*cos(Py())-sin(Px())*sin(Py()),
                 2*sin(Px())*sin(Py())+cos(Px())*cos(Py()) ) );

    // right hand side
    form1( Xh, F, _init=true )  =
        integrate( elements(mesh), _Q<Order+5>(), trans(f)*id(v) )+
        integrate( boundaryfaces(mesh),
                   // higher order quadrature to accurately integrate u_exact
                   _Q<3*Order>(),
                   trans(u_exact)*(-SigmaN+penalbc*id(v)/hFace() ) );

    Log() << "[stokes] vector local assembly done\n";

    /*
     * Construction of the left hand side
     */
    //# marker7 #
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh, Xh, D, _init=true )=integrate( elements(mesh), _Q<2*(Order-1)>(),
                                              mu*trace(deft*trans(def)) );
    form2( Xh, Xh, D )+=integrate( elements(mesh), _Q<2*(Order-1)>(),
                                   - div(v)*idt(p) + divt(u)*id(q) );
    form2( Xh, Xh, D )+=integrate( elements(mesh), _Q<Order-1>(),
                                   id(q)*idt(lambda) + idt(p)*id(nu) );
    form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), _Q<2*Order>(),
                                   -trans(SigmaNt)*id(v)
                                   -trans(SigmaN)*idt(u)
                                   +penalbc*trans(idt(u))*id(v)/hFace() );
    //# endmarker7 #
    Log() << "[stokes] matrix local assembly done\n";
    D->close();
    F->close();
    Log() << "[stokes] vector/matrix global assembly done\n";


    this->solve( D, U, F, false );

    Log() << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";
    std::cout << "value of the Lagrange multiplier lambda= " << lambda(0) << "\n";

    double u_errorL2 = integrate( elements(mesh), _Q<2*Order>(), trans(idv(u)-u_exact)*(idv(u)-u_exact) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;


    double p_errorL2 = integrate( elements(mesh), _Q<2*Order>(), (idv(p)-p_exact)*(idv(p)-p_exact) ).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";;

    Log() << "[stokes] solve for D done\n";

    double meas = integrate( elements(mesh), _Q<0>(), constant(1.0) ).evaluate()( 0, 0);
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";

    double mean_p = integrate( elements(mesh), _Q<Order-1>(), idv(p) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    double mean_div_u = integrate( elements(mesh), _Q<Order-1>(), divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements(mesh), _Q<2*(Order-1)>(), divv(u)*divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";
    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );

    this->exportResults( U, V );

    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";
} // Stokes::run

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, Order, Entity>::solve( sparse_matrix_ptrtype const& D,
                                   element_type& u,
                                   vector_ptrtype const& F,
                                   bool is_sym )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F, false );
    u = *U;
} // Stokes::solve

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, Order, Entity>::exportResults( element_type& U, element_type& V )
{
    exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 0 )->add( "u", U.template element<0>() );
    exporter->step( 0 )->add( "p", U.template element<1>() );
    exporter->step( 0 )->add( "u_exact", V.template element<0>() );
    exporter->step( 0 )->add( "p_exact", V.template element<1>() );
    exporter->save();
} // Stokes::export
} // Life

int
main( int argc, char** argv )
{

    using namespace Life;
    /* assertions handling */
    Life::Assert::setLog( "stokes.assert");

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;

    typedef Life::Stokes<nDim, nOrder, Simplex> stokes_type;


    /* define and run application */
    stokes_type stokes( argc, argv, makeAbout(), makeOptions() );
    stokes.run();
}





