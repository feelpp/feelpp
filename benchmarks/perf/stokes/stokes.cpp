/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>



#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

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
    Feel::po::options_description stokesoptions("Stokes options");
    stokesoptions.add_options()
        ("penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter")
        ("f", Feel::po::value<double>()->default_value( 0 ), "forcing term")
        ("mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component")
        ("hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence")
        ("bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return stokesoptions.add( Feel::feel_options() ) ;
}


/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "stokes" ,
                           "stokes" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2010 Universite de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
   return about;

}


namespace Feel
{
/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
template<int Dim,
         typename BasisU,
         typename BasisP,
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

    /*mesh*/
    typedef Entity<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #
    typedef BasisU basis_u_type;
    typedef BasisP basis_p_type;
    typedef bases<basis_u_type,basis_p_type> basis_type;
    //# endmarker1 #

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

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

    double mu;
    double penalbc;

    boost::shared_ptr<export_type> exporter;
}; // Stokes


template<int Dim, typename BasisU, typename BasisP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, BasisU, BasisP, Entity>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    using namespace Feel::vf;

    if ( this->vm().count( "nochdir" ) == false )
    {
        this->changeRepository( boost::format( "perf/stokes/%1%/%2%/P%3%/h_%4%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % BasisU::nOrder
                                % this->vm()["hsize"].template as<double>() );
    }

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=domain( _name= (boost::format( "%1%-%2%-%3%" ) % "hypercube" % Dim % 1).str() ,
                                                      _shape="hypercube",
                                                      _dim=Dim,
                                                      _h=meshSize,
                                                      _xmin=-1.,_xmax=1.,
                                                      _ymin=-1.,_ymax=1. ) );


    /*
     * The function space and some associate elements are then defined
     */
    //# marker4 #
    boost::timer t;
    space_ptrtype Xh = space_type::New( mesh );

    auto U = Xh->element( "u" );
    auto V = Xh->element( "v" );
    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();

    Log() << "Data Summary:\n";
    Log() << "   hsize = " << meshSize << "\n";
    Log() << "  export = " << this->vm().count("export") << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";

    Log() << "  -- time space and functions construction "<<t.elapsed()<<" seconds \n"; t.restart() ;

    auto F = M_backend->newVector( Xh );

    //# marker5 #
    auto deft = gradt(u);
    auto def = grad(v);
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    auto SigmaNt = -idt(p)*N()+mu*deft*N();

    // total stress tensor (test)
    auto SigmaN = -id(p)*N()+mu*def*N();
    //# endmarker6 #

    // u exact solution
    //auto u_exact = vec(cos(Px())*cos(Py()), sin(Px())*sin(Py()));
    auto u_exact = vec(-exp(Px())*(Py()*cos(Py())+sin(Py())), exp(Px())*Py()*sin(Py()));

    // this is the exact solution which has zero mean : the mean of
    // cos(x)*sin(y) is sin(1)*(1-cos(1))) on [0,1]^2
    //auto p_exact = cos(Px())*sin(Py())-(sin(1.0)*(1.-cos(1.0)));
    auto p_exact = 2.0*exp(Px())*sin(Py());

    // f is such that f = \Delta u_exact + \nabla p_exact
    //auto f = vec( (2*cos(Px())*cos(Py())-sin(Px())*sin(Py())),
    //              2*sin(Px())*sin(Py())+cos(Px())*cos(Py()) );
    auto f = vec(cst(0.),cst(0.)) ;

    // right hand side
    form1( Xh, F, _init=true )  = integrate( elements(mesh), trans(f)*id(v) );
    if ( this->vm()[ "bctype" ].template as<int>() == 1  )
    {
        form1( Xh, F, _init=true ) += integrate( boundaryfaces(mesh), trans(u_exact)*(-SigmaN+penalbc*id(v)/hFace() ) );
    }

    Log() << "  -- vector local assembly done in "<<t.elapsed()<<" seconds \n"; t.restart() ;

    /*
     * Construction of the left hand side
     */
    //# marker7 #
    auto D = M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, D, _init=true );
    boost::timer subt;
    form2( Xh, Xh, D ) =integrate( _range=elements(mesh),_expr=mu*(trans(dxt(u))*dx(v)+trans(dyt(u))*dy(v)));
    Log() << "   o time for diffusion terms: " << subt.elapsed() << "\n";subt.restart();

    form2( Xh, Xh, D )+=integrate( _range=elements(mesh),_expr=-div(v)*idt(p));
    form2( Xh, Xh, D )+=integrate( _range=elements(mesh),_expr=divt(u)*id(q) );
    Log() << "   o time for velocity/pressure terms: " << subt.elapsed() << "\n";subt.restart();

    if ( this->vm()[ "bctype" ].template as<int>() == 1  )
    {
        form2( Xh, Xh, D )+=integrate( _range=boundaryfaces(mesh),_expr=-trans(SigmaNt)*id(v) );
        form2( Xh, Xh, D )+=integrate( _range=boundaryfaces(mesh),_expr=-trans(SigmaN)*idt(u) );
        form2( Xh, Xh, D )+=integrate( _range=boundaryfaces(mesh),_expr=+penalbc*trans(idt(u))*id(v)/hFace() );
        Log() << "   o time for weak dirichlet terms: " << subt.elapsed() << "\n";subt.restart();
    }

    //# endmarker7 #
    D->close();
    F->close();
    if ( this->vm()[ "bctype" ].template as<int>() == 0  )
    {
        form2( Xh, Xh, D ) += on( boundaryfaces(mesh), u, F, u_exact );
        Log() << "   o time for strong dirichlet terms: " << subt.elapsed() << "\n";subt.restart();
    }
    Log() << " -- time vector/matrix global assembly done in "<<t.elapsed()<<" seconds \n"; t.restart() ;


    t.restart();

    M_backend->solve( _matrix=D, _solution=U, _rhs=F );

    Log() << " -- time for solver : "<<t.elapsed()<<" seconds \n";


    double meas = integrate( _range=elements(mesh), _expr=constant(1.0) ).evaluate()( 0, 0);
    Log() << "[stokes] measure(Omega)=" << meas << " (should be equal to 4)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to 4)\n";

    double mean_p = integrate( elements(mesh), idv(p) ).evaluate()( 0, 0 )/meas;
    Log() << "[stokes] mean(p)=" << mean_p << "\n";
    std::cout << "[stokes] mean(p)=" << mean_p << "\n";

    // get the zero mean pressure
    p.add( - mean_p );

    size_type nnz = 0 ;
    auto nNz = D->graph()->nNz() ;
    for(auto iter = nNz.begin();iter!=nNz.end();++iter)
      nnz += (*iter) ;
    Log() << "[stokes] matrix NNZ "<< nnz << "\n";

    double u_errorL2 = integrate( _range=elements(mesh), _expr=trans(idv(u)-u_exact)*(idv(u)-u_exact) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;
    Log() << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";;


    double p_errorL2 = integrate( _range=elements(mesh), _expr=(idv(p)-p_exact)*(idv(p)-p_exact) ).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";;
    Log() << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";;

    Log() << "[stokes] solve for D done\n";

    double mean_div_u = integrate( elements(mesh), divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements(mesh), divv(u)*divv(u) ).evaluate()( 0, 0 );
    Log() << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";
    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( Xh->template functionSpace<0>(), elements( Xh->mesh() ), u_exact );
    q = vf::project( Xh->template functionSpace<1>(), elements( Xh->mesh() ), p_exact );

    this->exportResults( U, V );

    Log() << "[mesh]   number of elements: " << Xh->mesh()->numElements() << "\n";
    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";
} // Stokes::run


template<int Dim, typename BasisU, typename BasisP, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Stokes<Dim, BasisU, BasisP, Entity>::exportResults( element_type& U, element_type& V )
{
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", U.template element<0>() );
        exporter->step( 0 )->add( "p", U.template element<1>() );
        exporter->step( 0 )->add( "u_exact", V.template element<0>() );
        exporter->step( 0 )->add( "p_exact", V.template element<1>() );
        exporter->save();
    }
} // Stokes::export
} // Feel

int
main( int argc, char** argv )
{

    using namespace Feel;
    /* assertions handling */
    Feel::Assert::setLog( "stokes.assert");

    const int nDim = 2;
    //typedef Feel::Stokes<nDim, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex> stokes_type;
    //typedef Feel::Stokes<nDim, Lagrange<1, Vectorial>,Lagrange<1, Scalar>, Simplex> stokes_type;
    typedef Feel::Stokes<nDim, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<0, Scalar,Discontinuous>, Simplex> stokes_type;
    //typedef Feel::Stokes<nDim, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<1, Scalar>, Simplex> stokes_type;


    /* define and run application */
    stokes_type stokes( argc, argv, makeAbout(), makeOptions() );
    stokes.run();
}
