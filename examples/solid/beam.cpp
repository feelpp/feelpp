/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-02-18

  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble)

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
   \file beam.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-02-18
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifepoly/polynomialset.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifemesh/meshmover.hpp>

#include <life/lifevf/vf.hpp>



inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description beamoptions("Beam options");
    beamoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("beta", Life::po::value<double>()->default_value( 1.0 ), "beta value in -Delta u + beta u = f")
        ("bccoeff", Life::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
        ("bctype", Life::po::value<int>()->default_value( 1 ), "Dirichlet condition type(0=elimination,1=penalisation, 2=weak")
        ("scale", Life::po::value<double>()->default_value( 10 ), "scale factor for mesh mover")
        ("export", "export results(ensight, data file(1D)")
        ("export-matlab", "export matrix and vectors in matlab" )

        ;
    return beamoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "beam" ,
                           "beam" ,
                           "0.1",
                           "Linear elasticity model for a beam",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2007,2008 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
template<int nDim, int nOrder>
class Beam
    :
    public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type Dim = nDim;
    static const uint16_type feOrder = nOrder;
    static const uint16_type imOrder = nOrder;

    typedef double value_type;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef fusion::vector<Lagrange<feOrder, Vectorial> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    Beam( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        beta( this->vm()["beta"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        M_bctype( this->vm()["bctype"].template as<int>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "beam" ) ),
        timers()
    {
        Log() << "[Beam] hsize = " << meshSize << "\n";
        Log() << "[Beam] beta = " << beta << "\n";
        Log() << "[Beam] bccoeff = " << bcCoeff << "\n";
        Log() << "[Beam] bctype = " <<  M_bctype << "\n";
        Log() << "[Beam] export = " << this->vm().count("export") << "\n";

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "beam" );
    }

    ~Beam()
    {
        std::map<std::string,std::pair<boost::timer,double> >::iterator it = timers.begin();
        std::map<std::string,std::pair<boost::timer,double> >::iterator en = timers.end();
        for( ; it != en; ++it )
            {
                Log() << it->first << " : " << it->second.second << " s elapsed\n";
            }
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
     * solve symmetric system
     */
    void solve( sparse_matrix_ptrtype const& D, element_type& u, vector_ptrtype const& F );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type const& u, element_type const& v );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double beta;
    double bcCoeff;
    int M_bctype;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;

    std::map<std::string,std::pair<boost::timer,double> > timers;
}; // Beam

template<int nDim, int nOrder>
typename Beam<nDim,nOrder>::mesh_ptrtype
Beam<nDim,nOrder>::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );

    GmshTensorizedDomain<Dim,1,Dim,Simplex> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( 0, 20 ) );
    td.setY( std::make_pair( -1, 1 ) );
    ImporterGmsh<mesh_type> import( td.generate( entity_type::name().c_str() ) );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    return mesh;
} // Beam::createMesh

template<int nDim, int nOrder>
void
Beam<nDim,nOrder>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    this->changeRepository( boost::format( "%1%/%2%/P%3%/%4%/" )
                           % this->about().appName()
                           % entity_type::name()
                           % nOrder
                           % this->vm()["hsize"].template as<double>()
                           );
    this->setLogs();

    using namespace Life::vf;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    space_ptrtype Xh = space_type::New( mesh );
    Xh->printInfo();

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    timers["init"].second = timers["init"].first.elapsed();

    /*
     * Data associated with the simulation
     */
    const double E = 21*1e5;
    const double sigma = 0.28;
    const double mu = E/(2*(1+sigma));
    const double lambda = E*sigma/((1+sigma)*(1-2*sigma));
    const double density = 1;
    const double gravity = -density*0.05;
    Log() << "lambda = " << lambda << "\n"
              << "mu     = " << mu << "\n"
              << "gravity= " << gravity << "\n";

    /*
     * Construction of the right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ y \f$ direction.
     */
    vector_ptrtype F( M_backend->newVector( Xh ) );
    F->zero();
    timers["assembly"].first.restart();
    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F0.m" );

    form1( Xh, F, _init=true ) =   integrate( elements(mesh), _Q<imOrder>(), trans(gravity*oneY())*id(v) );
    F->close();
    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F1.m" );
    timers["assembly"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
    timers["assembly"].first.restart();
    AUTO( deft, 0.5*( gradt(u)+trans(gradt(u)) ) );
    AUTO( def, 0.5*( grad(v)+trans(grad(v)) ) );
    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), _Q<imOrder>(),
                   lambda*divt(u)*div(v)  +
                   2*mu*trace(trans(deft)*def) );
    if ( M_bctype == 1 ) // weak Dirichlet bc
        {
            AUTO( Id, (mat<nDim,nDim>( cst(1), cst(0), cst(0), cst(1.) )) );
            form2( Xh, Xh, D ) +=
                integrate( markedfaces(mesh,1), _Q<imOrder>(),
                           - trans((2*mu*deft+lambda*trace(deft)*Id )*N())*id(v)
                           - trans((2*mu*def+lambda*trace(def)*Id )*N())*idt(u)
                           + bcCoeff*trans(idt(u))*id(v)/hFace() );
        }


    D->close();
    if ( M_bctype == 0 )
        form2( Xh, Xh, D ) += on( markedfaces(mesh,(nDim==2)?1:23), u, F, constant(0)*one() );

    if ( this->vm().count( "export-matlab" ) )
        {
            F->printMatlab( "F2.m" );
            D->printMatlab( "elas.m" );
        }
    timers["assembly"].second += timers["assembly"].first.elapsed();

    this->solve( D, u, F );


    v = project( Xh, elements(Xh->mesh()), P() );

    this->exportResults( 0, u, v );

    MeshMover<mesh_type> meshmove;
    u.vec() *= this->vm()["scale"].template as<double>();
    meshmove.apply( Xh->mesh(), u );

    element_type w( Xh, "w" );
    w = project( Xh,
                 elements(Xh->mesh()),
                 P() );

    this->exportResults( 1, u, w );

    std::cout << "||v||_2 = " << v.l2Norm() << " (P() before move)\n";
    std::cout << "||w||_2 = " << w.l2Norm() << " (P() after move)\n";


    w.add( -1.0, v );

    std::cout << "||u||_2 = " << w.l2Norm() << " (displacement u)\n";
    std::cout << "||w-v||_2 = " << w.l2Norm() << " (displacement w-v\n";

    w.add( -1.0, u );
    std::cout << "||(w-v)-u||_2 = " << w.l2Norm() << " (should be 0, ie u=w-v)\n";
} // Beam::run

template<int nDim, int nOrder>
void
Beam<nDim,nOrder>::solve( sparse_matrix_ptrtype const& D, element_type& u, vector_ptrtype const& F )
{
    timers["solver"].first.restart();

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    Log() << "[timer] solve: " << timers["solver"].second << "\n";
}
template<int nDim, int nOrder>
void
Beam<nDim,nOrder>::exportResults( double time, element_type const& u, element_type const &v  )
{
    timers["export"].first.restart();

    typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
    timeStep->setMesh( u.functionSpace()->mesh() );
    timeStep->add( "displ", u );
    timeStep->add( "P", v );
    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
} // Beam::export
} // Life


int
main( int argc, char** argv )
{
    const int nDim = 2;
    const int nOrder = 3;

    Life::Beam<nDim,nOrder> beam( argc, argv, makeAbout(), makeOptions());
    beam.run();
}




