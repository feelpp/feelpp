/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file stokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifemesh/elements.hpp>

#include <life/lifevf/vf.hpp>

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

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim, 1,Dim> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef fem::Lagrange<Dim, Order, Vectorial,
                          Continuous, double, Entity> basis_u_type;
    typedef fem::Lagrange<Dim, Order-1, Scalar,
                          Continuous, double, Entity> basis_p_type;
    typedef fem::Lagrange<Dim, 0, Scalar,
                          Continuous, double, Entity> basis_l_type;
    typedef fusion::vector<basis_u_type,
                            basis_p_type, basis_l_type> basis_type;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    /* functions */
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;
    typedef typename element_type::template sub_element<2>::type element_2_type;

    /*quadrature*/
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;


    Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "stokes" ) )
    {

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "stokes" );

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
    void exportResults( element_type& u );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double penalbc;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;
}; // Stokes

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Stokes<Dim,Order,Entity>::mesh_ptrtype
Stokes<Dim,Order,Entity>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );


    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( entity_type::name().c_str() );

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
                            % entity_type::name()
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
    space_ptrtype Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
    element_0_type u = U.template element<0>();
    element_0_type v = V.template element<0>();
    element_1_type p = U.template element<1>();
    element_1_type q = V.template element<1>();
    element_2_type lambda = U.template element<2>();
    element_2_type nu = V.template element<2>();

    /*
     * a quadrature rule for numerical integration
     */
    im_type im;


    Log() << "Data Summary:\n";
    Log() << "   hsize = " << meshSize << "\n";
    Log() << "  export = " << this->vm().count("export") << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";


    vector_ptrtype F( M_backend->newVector( Xh ) );

    // viscous stress tensor (trial) : 0.5 ( \nabla u + \nabla u ^T )
    AUTO( deft, 0.5*( gradt(u)+trans(gradt(u)) ) );
    // viscous stress tensor (test) : 0.5 ( \nabla u + \nabla u ^T )
    AUTO( def, 0.5*( grad(v)+trans(grad(v)) ) );

    // total stress tensor (trial)
    AUTO( SigmaNt, (-idt(p)*N()+2*mu*deft*N()) );

    // total stress tensor (test)
    AUTO( SigmaN, (-id(p)*N()+2*mu*def*N()) );

    // top lid velocity
    AUTO( g, oneX() );

    // right hand side
    form1( Xh, F, _init=true )  =
        integrate( markedfaces(mesh,4), im,
                   trans(g)*(-SigmaN+penalbc*id(v)/hFace() ) );

    Log() << "[stokes] vector local assembly done\n";

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh, Xh, D, _init=true )=integrate( elements(mesh), im,
                                              mu*trace(deft*trans(def)) );
    form2( Xh, Xh, D )+=integrate( elements(mesh), im,
                                   - div(v)*idt(p) + divt(u)*id(q) );
    form2( Xh, Xh, D )+=integrate( elements(mesh), im,
                                   id(q)*idt(lambda) + idt(p)*id(nu) );
    form2( Xh, Xh, D )+=integrate( boundaryfaces(mesh), im,
                                   -trans(SigmaNt)*id(v)
                                   -trans(SigmaN)*idt(u)
                                   +penalbc*trans(idt(u))*id(v)/hFace() );

    Log() << "[stokes] matrix local assembly done\n";
    D->close();
    F->close();
    Log() << "[stokes] vector/matrix global assembly done\n";

    this->solve( D, U, F, false );

    Log() << "[stokes] solve for D done\n";

    double meas = integrate( elements(mesh), im, constant(1.0) ).evaluate()( 0, 0);
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 1)\n";

    double mean_p = integrate( elements(mesh), im, idv(p) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    this->exportResults( U );

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
Stokes<Dim, Order, Entity>::exportResults( element_type& U )
{
    typename timeset_type::step_ptrtype timeStep = timeSet->step( 1.0 );
    timeStep->setMesh( U.functionSpace()->mesh() );
    timeStep->add( "u", U.template element<0>() );
    timeStep->add( "p", U.template element<1>() );
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





