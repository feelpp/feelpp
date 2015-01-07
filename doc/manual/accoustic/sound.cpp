/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-05

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file sound.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-05
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/solvereigen.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>

namespace Feel
{
gmsh_ptrtype createRoom( int Dim, double meshSize );

inline
po::options_description
makeOptions()
{
    po::options_description soundoptions( "Sound options" );
    soundoptions.add_options()
    ( "kc2", Feel::po::value<double>()->default_value( 1 ), "k/c parameter" )
    ( "sigma", Feel::po::value<double>()->default_value( 20 ), "shift parameter for the eigenvalue problem" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ;
    return soundoptions.add( Feel::feel_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "sound" ,
                     "sound" ,
                     "0.2",
                     "nD(n=1,2,3) acoustics in an amphitheater",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2007-2011 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 * Sound Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim, int Order>
class Sound
    :
public Simget
{
    typedef Simget super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar> >,Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Sound( std::string const& name )
        :
        super(),
        M_name( name )
    {
    }

    std::string name() const
    {
        return M_name;
    }

    /**
     * run the convergence test
     */
    void run();

private:
    std::string M_name;
}; // Sound


template<int Dim, int Order>
void
Sound<Dim, Order>::run()
{
    LOG(INFO) << "[Sound] hsize = " << this->meshSize() << "\n";
    this->changeRepository( boost::format( "%1%/%2%/P%3%/%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % doption("hsize")
                          );
    //! backend
    auto backend = backend_type::build( soption("backend") );

    //! eigen solver
    auto eigen = SolverEigen<value_type>::build( this->vm() );

    //! exporter to paraview or gmsh
    auto exporter = Exporter<mesh_type>::New( this->vm(), this->about().appName() );

    boost::timer t;

    /*
     * First we create the mesh
     */
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc = createRoom( Dim,this->meshSize() ) );

    M_stats.put( "h",this->meshSize() );
    M_stats.put( "n.space.nelts",mesh->numElements() );
    M_stats.put( "t.init.mesh",t.elapsed() );
    t.restart();

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    M_stats.put( "n.space.ndof",Xh->nLocalDof() );
    M_stats.put( "t.init.space",t.elapsed() );
    t.restart();

    auto F = backend->newVector( Xh );

    if ( Dim == 2 )
        form1( _test=Xh, _vector=F, _init=true )  = integrate( _range=markedfaces( mesh,2 ), _expr=val( Py()*( 1-Py() ) )*id( v ) );

    else
        form1( _test=Xh, _vector=F, _init=true )  = integrate( _range=markedfaces( mesh,51 ),  _expr=val( Py()*( 1-Py() )*Pz()*( 1-Pz() ) )*id( v ) );

    M_stats.put( "t.assembly.vector.total",t.elapsed() );
    t.restart();

    /*
     * Construction of the left hand side
     */
    auto D = backend->newMatrix( Xh, Xh );

    double kc2 = doption("kc2");

    form2( _test=Xh, _trial=Xh, _matrix=D,_init=true ) = integrate( _range=elements( mesh ),  _expr=( kc2*idt( u )*id( v )-gradt( u )*trans( grad( v ) ) ) );
    D->close();

    M_stats.put( "t.assembly.matrix.total",t.elapsed() );

    t.restart();

    backend->solve( _matrix=D, _solution=u, _rhs=F );

    M_stats.put( "t.solver.total",t.elapsed() );
    t.restart();

    // eigen modes
    double sigma = doption("sigma");
    auto S = backend->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=S, _init=true ) = integrate( _range=elements( mesh ),  _expr=gradt( u )*trans( grad( v ) ) );

    M_stats.put( "t.assembly.matrix.A",t.elapsed() );
    t.restart();

    auto M = backend->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M, _init=true ) = integrate( _range=elements( mesh ),  _expr=idt( u )*id( v ) );
    M_stats.put( "t.assembly.matrix.B",t.elapsed() );
    t.restart();


    int maxit = ioption("solvereigen.maxiter");
    int tol =   ioption("solvereigen.tol"    );
    int nev =   ioption("solvereigen.nev"    );
    int ncv =   ioption("solvereigen.ncv"    );

    double eigen_real, eigen_imag;

    vector_ptrtype modepetsc( backend->newVector( Xh ) );
    int nconv;
    //boost::tie( nconv, eigen_real, eigen_imag, modepetsc )  =

    SolverEigen<double>::eigenmodes_type modes;
    modes=
        eigs( _matrixA=S,
              _matrixB=M,
              // this is a generalized Hermitian Eigenvalue Problem
              _problem=( EigenProblemType )GHEP,
              _nev=nev,
              _ncv=ncv,
              _maxit=maxit,
              _tolerance=tol,
              _spectrum=( PositionOfSpectrum )ioption("solvereigen.position") );

    element_type mode( Xh, "mode" );

    if ( !modes.empty() )
    {
        LOG(INFO) << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";
        std::cout << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>()
                  << "," <<  modes.begin()->second.get<1>() << ")\n";

        mode = *modes.begin()->second.get<2>();
    }

    //    this->exportResults( u, mode );
    //    LOG(INFO) << "[timer] run(): init (" << mesh->numElements() << " Elems): " << timers["init"].second << "\n";
    //    LOG(INFO) << "[timer] run(): assembly (" << Xh->dof()->nDof() << " DOFs): " << timers["assembly"].second << "\n";
    M_stats.put( "t.eigensolver.total",t.elapsed() );
    t.restart();

    exporter->step( 0. )->setMesh( u.functionSpace()->mesh() );

    exporter->step( 0. )->add( "p", u );
    exporter->step( 0. )->add( "mode", mode );
    exporter->save();
} // Sound::run
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,_desc=makeOptions(),_about=makeAbout() );

    /* assertions handling */
    Feel::Assert::setLog( "sound.assert" );

    Application sound;

    sound.setStats( boost::assign::list_of( "n.space" )( "t.init" )( "t.assembly.vector" )( "t.assembly.matrix" )( "t.solver" )( "t.eigensolver" ) );
    sound.add( new Sound<2, 1>( "2D-P1" ) );
    sound.add( new Sound<3, 1>( "3D-P1" ) );
    sound.run();
    sound.printStats( std::cout );
}





