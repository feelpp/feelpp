/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-05

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-05
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>
#include <life/lifealg/solvereigen.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>

std::pair<std::string,std::string> createRoom( int Dim, double meshSize );



inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description soundoptions("Sound options");
    soundoptions.add_options()
        ("kc2", Life::po::value<double>()->default_value( 1 ), "k/c parameter")
        ("sigma", Life::po::value<double>()->default_value( 20 ), "shift parameter for the eigenvalue problem")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")

        ("export", "export results(ensight, data file(1D)")
        ("export-mesh-only", "export mesh only in ensight format")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return soundoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "sound" ,
                            "sound" ,
                            "0.2",
                            "nD(n=1,2,3) acoustics in an amphitheater",
                            Life::AboutData::License_GPL,
                            "Copyright (c) 2007 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
using namespace vf;
/**
 * Sound Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim,
         int Order,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Sound
    :
    public Application
{
    typedef Application super;
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
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >,Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;

    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Sound( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        eigen( SolverEigen<value_type>::build( this->vm() ) ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers(),
        stats()
    {
        Log() << "[Sound] hsize = " << meshSize << "\n";
        Log() << "[Sound] export = " << this->vm().count("export") << "\n";
    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve the system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solve( Mat& D, Vec1& u, Vec2& F, bool is_sym );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& mode0 );

private:

    backend_ptrtype M_backend;
    double meshSize;
    boost::shared_ptr<SolverEigen<value_type> > eigen;


    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // Sound

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Sound<Dim,Order,Entity>::mesh_ptr_type
Sound<Dim,Order,Entity>::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );
    //mesh->setRenumber( false );

    Gmsh gmsh;
    gmsh.setOrder( GMSH_ORDER_ONE );
    std::string mesh_name, mesh_desc;
    boost::tie( mesh_name, mesh_desc ) = ::createRoom(Dim,meshSize);
    std::string fname = gmsh.generate( mesh_name, mesh_desc );

    ImporterGmsh<mesh_type> import( fname );
    import.setVersion( "2.0" );
    //mesh->setRenumber( false );
    mesh->accept( import );

    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // Sound::createMesh


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Sound<Dim, Order, Entity>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    //    int maxIter = 10.0/meshSize;
    using namespace Life::vf;


    this->changeRepository( boost::format( "%1%/%2%/P%3%/%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            );
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( meshSize );
    stats["nelt"] = mesh->elements().size();

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    space_ptrtype Xh = space_type::New( mesh );
    //Xh->dof()->showMe();
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();

    /*
     * a quadrature rule for numerical integration
     */
    im_type im;


    vector_ptrtype F( M_backend->newVector( Xh ) );
    timers["assembly"].first.restart();


    if ( Dim == 2 )
        form( Xh, *F )  = integrate( markedfaces(mesh,2), im, val(Py()*(1-Py()))*id(v) );
    else
        form( Xh, *F )  = integrate( markedfaces(mesh,51), im, val(Py()*(1-Py())*Pz()*(1-Pz()))*id(v) );

    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F.m" );
    timers["assembly"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );

    double kc2 = this->vm()["kc2"].template as<double>();

    timers["assembly"].first.restart();

    form2( Xh, Xh, D ) = integrate( elements(mesh), im, ( kc2*idt(u)*id(v)-gradt(u)*trans(grad(v))) );
    D->close();

    timers["assembly"].second += timers["assembly"].first.elapsed();

    if ( this->vm().count( "export-matlab" ) )
        D->printMatlab( "D.m" );

    this->solve( D, u, F, true );

    // eigen modes
    double sigma = this->vm()["sigma"].template as<double>();
    sparse_matrix_ptrtype S( M_backend->newMatrix( Xh, Xh ) );
    form2( Xh, Xh, S ) = integrate( elements(mesh), im, gradt(u)*trans(grad(v)) );
    S->close();
    if ( this->vm().count( "export-matlab" ) )
        S->printMatlab( "S" );

    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );
    form2( Xh, Xh, M ) = integrate( elements(mesh), im, idt(u)*id(v));
    M->close();
    if ( this->vm().count( "export-matlab" ) )
        M->printMatlab( "M.m" );

    int maxit = this->vm()["solvereigen-maxiter"].template as<int>();
    int tol = this->vm()["solvereigen-tol"].template as<double>();

    int nev = this->vm()["solvereigen-nev"].template as<int>();

    int ncv = this->vm()["solvereigen-ncv"].template as<int>();;

    double eigen_real, eigen_imag;

    vector_ptrtype modepetsc( M_backend->newVector( Xh ) );
    int nconv;
    boost::tie( nconv, eigen_real, eigen_imag, modepetsc )  =
        eigs( _matrixA=S,
              _matrixB=M,
              _nev=nev,
              _ncv=ncv,
              _spectrum=LARGEST_MAGNITUDE );
    element_type mode( Xh, "mode" );

    Log() << "eigenvalue " << 0 << " = (" << eigen_real << "," << eigen_imag  << ")\n";
    //Log() << "eigenvalue " << 0 << " relative error = " << eigen->relativeError( 0 ) << "\n";

    mode = *modepetsc;
    this->exportResults( u, mode );

    Log() << "[timer] run(): init (" << mesh->numElements() << " Elems): " << timers["init"].second << "\n";
    Log() << "[timer] run(): assembly (" << Xh->dof()->nDof() << " DOFs): " << timers["assembly"].second << "\n";

} // Sound::run

template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename Mat, typename Vec1, typename Vec2>
void
Sound<Dim, Order, Entity>::solve( Mat& D, Vec1& u, Vec2& F, bool is_sym  )
{
    timers["solver"].first.restart();

    //backend.set_symmetric( is_sym );

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    Log() << "[timer] solve: " << timers["solver"].second << "\n";
} // Sound::solve


template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Sound<Dim, Order, Entity>::exportResults( element_type& U, element_type& Mode )
{
    timers["export"].first.restart();

    if ( this->vm().count( "export" ) )
        {
            Log() << "exportResults starts\n";
            exporter->step(1.)->setMesh( U.functionSpace()->mesh() );
            //exporter->step(1.)->setMesh( this->createMesh( meshSize/2, 0.5, 1 ) );
            //exporter->step(1.)->setMesh( this->createMesh( meshSize/Order, 0, 1 ) );
            //exporter->step(1.)->setMesh( this->createMesh( meshSize ) );
            if ( !this->vm().count( "export-mesh-only" ) )
                {
                    exporter->step(1.)->add( "p", U );
                    exporter->step(1.)->add( "mode", Mode );
                }
            exporter->save();
        }
    timers["export"].second = timers["export"].first.elapsed();
    Log() << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Sound::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 3; // 2 or 3
    const int nOrder = 1;

    typedef Life::Sound<nDim, nOrder, Simplex> sound_type;

    /* assertions handling */
    Life::Assert::setLog( "sound.assert");

    /* define and run application */
    sound_type sound( argc, argv, makeAbout(), makeOptions() );

    sound.run();
}





