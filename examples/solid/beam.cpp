/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-02-18

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble)

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
   \file beam.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-02-18
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feelmesh/meshmover.hpp>

#include <feel/feelvf/vf.hpp>



inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description beamoptions( "Beam options" );
    beamoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.01 ), "first h value to start convergence" )
    ( "beta", Feel::po::value<double>()->default_value( 1.0 ), "beta value in -Delta u + beta u = f" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "bctype", Feel::po::value<int>()->default_value( 1 ), "Dirichlet condition type(0=elimination,1=penalisation, 2=weak" )
    ( "scale", Feel::po::value<double>()->default_value( 10 ), "scale factor for mesh mover" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-matlab", "export matrix and vectors in matlab" )

    ;
    return beamoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "beam" ,
                           "beam" ,
                           "0.2",
                           "Linear elasticity model for a beam",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007-2012 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
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
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef fusion::vector<Lagrange<feOrder, Vectorial> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Beam( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        beta( this->vm()["beta"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        M_bctype( this->vm()["bctype"].template as<int>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers()
    {
        Log() << "[Beam] hsize = " << meshSize << "\n";
        Log() << "[Beam] beta = " << beta << "\n";
        Log() << "[Beam] bccoeff = " << bcCoeff << "\n";
        Log() << "[Beam] bctype = " <<  M_bctype << "\n";
        Log() << "[Beam] export = " << this->vm().count( "export" ) << "\n";

    }

    ~Beam()
    {
        std::map<std::string,std::pair<boost::timer,double> >::iterator it = timers.begin();
        std::map<std::string,std::pair<boost::timer,double> >::iterator en = timers.end();

        for ( ; it != en; ++it )
        {
            Log() << it->first << " : " << it->second.second << " s elapsed\n";
        }
    }

    /**
     * run the convergence test
     */
    void run();

private:

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

    std::map<std::string,std::pair<boost::timer,double> > timers;
}; // Beam

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
    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "beam-%1%" ) % nDim ).str() ,
                                                _shape="hypercube",
                                                _usenames=true,
                                                _xmin=0., _xmax=0.351,
                                                _ymin=0., _ymax=0.02,
                                                _h=meshSize ),
                                        _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                        _partitions=this->comm().size()  );
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
#if 0
    const double E = 21*1e5;
    const double sigma = 0.28;
#else
    const double E = 1.4e6;
    const double sigma = 0.4;
#endif
    const double mu = E/( 2*( 1+sigma ) );
    const double lambda = E*sigma/( ( 1+sigma )*( 1-2*sigma ) );
    const double density = 1e3;
    const double gravity = -2;//-density*0.05;
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

    form1( Xh, F, _init=true ) =   integrate( elements( mesh ), trans( gravity*oneY() )*id( v ) );
    F->close();

    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F1.m" );

    timers["assembly"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
    timers["assembly"].first.restart();
    auto deft = 0.5*( gradt( u )+trans( gradt( u ) ) );
    auto def = 0.5*( grad( v )+trans( grad( v ) ) );
    form2( Xh, Xh, D, _init=true ) =
        integrate( elements( mesh ),
                   lambda*divt( u )*div( v )  +
                   2*mu*trace( trans( deft )*def ) );

    if ( M_bctype == 1 ) // weak Dirichlet bc
    {
        auto Id = ( mat<nDim,nDim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) );
        form2( Xh, Xh, D ) +=
            integrate( markedfaces( mesh,1 ),
                       - trans( ( 2*mu*deft+lambda*trace( deft )*Id )*N() )*id( v )
                       - trans( ( 2*mu*def+lambda*trace( def )*Id )*N() )*idt( u )
                       + bcCoeff*trans( idt( u ) )*id( v )/hFace() );
    }

    if ( M_bctype == 0 )
        form2( Xh, Xh, D ) += on( markedfaces( mesh,( nDim==2 )?1:23 ), u, F, constant( 0 )*one() );

    if ( this->vm().count( "export-matlab" ) )
    {
        F->printMatlab( "F2.m" );
        D->printMatlab( "elas.m" );
    }

    timers["assembly"].second += timers["assembly"].first.elapsed();

    M_backend->solve( _matrix=D, _solution=u, _rhs=F );

    auto i1 = integrate( markedfaces( mesh,3 ), idv( u ) ).evaluate();
    std::cout << "deflection: " << i1/0.02 << "\n";
    v = vf::project( Xh, elements( Xh->mesh() ), P() );

    this->exportResults( 0, u, v );
#if 0
    MeshMover<mesh_type> meshmove;
    u.vec() *= this->vm()["scale"].template as<double>();
    meshmove.apply( Xh->mesh(), u );

    element_type w( Xh, "w" );

    for ( int i = 1; i < 5; ++i )
    {
        w = vf::project( Xh,
                         elements( Xh->mesh() ),
                         P() );
        this->exportResults( i, u, w );
        meshmove.apply( Xh->mesh(), u );
    }

    std::cout << "||v||_2 = " << v.l2Norm() << " (P() before move)\n";
    std::cout << "||w||_2 = " << w.l2Norm() << " (P() after move)\n";


    w.add( -1.0, v );

    std::cout << "||u||_2 = " << w.l2Norm() << " (displacement u)\n";
    std::cout << "||w-v||_2 = " << w.l2Norm() << " (displacement w-v\n";

    w.add( -1.0, u );
    std::cout << "||(w-v)-u||_2 = " << w.l2Norm() << " (should be 0, ie u=w-v)\n";
#endif
} // Beam::run

template<int nDim, int nOrder>
void
Beam<nDim,nOrder>::exportResults( double time, element_type const& u, element_type const &v  )
{
    timers["export"].first.restart();

    exporter->step( time )->setMesh( u.functionSpace()->mesh() );
    exporter->step( time )->add( "displ", u );
    exporter->step( time )->add( "P", v );
    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
} // Beam::export
} // Feel


int
main( int argc, char** argv )
{
    Environment env( argc, argv );
    const int nDim = 2;
    const int nOrder = 3;

    Feel::Beam<nDim,nOrder> beam( argc, argv, makeAbout(), makeOptions() );
    beam.run();
}




