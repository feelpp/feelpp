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

#include <feel/feel.hpp>
#include <feel/feelcore/units.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description beamoptions( "Beam options" );
    beamoptions.add_options()
        ( "E", Feel::po::value<double>()->default_value( 1.4e6 ), "Young modulus" )
        ( "nu", Feel::po::value<double>()->default_value( 0.4 ), "Poisson coefficient" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.01 ), "first h value to start convergence" )
        ( "beta", Feel::po::value<double>()->default_value( 1.0 ), "beta value in -Delta u + beta u = f" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
        ( "bctype", Feel::po::value<int>()->default_value( 1 ), "Dirichlet condition type(0=elimination,1=weak" )
        ( "scale", Feel::po::value<double>()->default_value( 10 ), "scale factor for mesh mover" )
        ;
    return beamoptions.add( Feel::feel_options() );
}


namespace Feel
{
/**
   \page LinearElasticity Linear Elasticity : a Beam example
   \author Christophe Prud'homme
   \date 2007-02-18

   This is the documentation for the beam example
*/
template<int nDim, int nOrder>
class Beam
    :
public Simget
{
    typedef Simget super;
public:

    // -- TYPEDEFS --
    static const uint16_type Dim = nDim;

    typedef double value_type;

    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef bases<Lagrange<nOrder, Vectorial> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Beam()
        :
        super(),
        meshSize( doption("hsize"  )),
        beta(     doption("beta"   )),
        bcCoeff(  doption("bccoeff")),
        M_bctype( ioption("bctype" )),
        exporter( Exporter<mesh_type>::New( this->about().appName() ) ),
        timers()
    {
        LOG(INFO) << "[Beam] hsize = " << meshSize << "\n";
        LOG(INFO) << "[Beam] beta = " << beta << "\n";
        LOG(INFO) << "[Beam] bccoeff = " << bcCoeff << "\n";
        LOG(INFO) << "[Beam] bctype = " <<  M_bctype << "\n";

    }

    ~Beam()
    {
        for( auto const& d: timers )
        {
            LOG(INFO) << d.first << " : " << d.second.second << " s elapsed\n";
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

    double meshSize;
    double beta;
    double bcCoeff;
    int M_bctype;

    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<mpi::timer,double> > timers;
}; // Beam

template<int nDim, int nOrder>
void
Beam<nDim,nOrder>::run()
{

    this->changeRepository( boost::format( "doc/manual/solid/%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % nOrder
                            % meshSize );
    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                        _desc=domain( _name=( boost::format( "beam-%1%" ) % nDim ).str() ,

                                                      _shape="hypercube",
                                                      _xmin=0., _xmax=0.351,
                                                      _ymin=0., _ymax=0.02,
                                                      _zmin=0., _zmax=0.02,
                                                      _h=meshSize ) );
    // add marker clamped to the mesh
    mesh->addMarkerName( "clamped",( nDim==2 )?1:19, (nDim==2)?1:2);
    mesh->addMarkerName( "tip",( nDim==2)?3:27, (nDim==2)?1:2);
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
    auto E = doption(_name="E")*pascal;
    const double nu = doption(_name="nu");

    auto mu = E/( 2*( 1+nu ) );
    auto lambda = E*nu/( ( 1+nu )*( 1-2*nu ) );
    auto density = 1e3;
    auto gravity = -2*newton/pow<Dim>(meter);//-density*0.05;
    LOG(INFO) << "lambda = " << lambda << "\n"
          << "mu     = " << mu << "\n"
          << "gravity= " << gravity << "\n";

    /*
     * Construction of the right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ y \f$ direction.
     */
    auto F = backend()->newVector( Xh );
    F->zero();
    timers["assembly"].first.restart();

    if ( Dim == 3 )
        form1( _test=Xh, _vector=F ) = integrate( elements( mesh ), trans( gravity.value()*oneZ() )*id( v ) );
    else
        form1( _test=Xh, _vector=F ) = integrate( elements( mesh ), trans( gravity.value()*oneY() )*id( v ) );

    timers["assembly"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    auto D = backend()->newMatrix( Xh, Xh );
    timers["assembly"].first.restart();
    auto deft = sym(gradt(u));
    auto def = sym(grad(u));
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D );
    a = integrate( elements( mesh ),
                   lambda.value()*divt( u )*div( v )  +
                   2.*mu.value()*trace( trans( deft )*def ) );

    if ( M_bctype == 1 ) // weak Dirichlet bc
    {
        auto Id = eye<nDim>();
        a += integrate( markedfaces( mesh, "clamped" ),
                        - trans( ( 2.*mu.value()*deft+lambda.value()*trace( deft )*Id )*N() )*id( v )
                        - trans( ( 2.*mu.value()*def+lambda.value()*trace( def )*Id )*N() )*idt( u )
                        + bcCoeff*std::max(2.*mu.value(),lambda.value())*trans( idt( u ) )*id( v )/hFace() );
    }

    if ( M_bctype == 0 )
        a += on( markedfaces( mesh, "clamped" ), u, F, zero<nDim,1>() );

    timers["assembly"].second += timers["assembly"].first.elapsed();

    backend(_rebuild=true)->solve( _matrix=D, _solution=u, _rhs=F );

    v = vf::project( Xh, elements( Xh->mesh() ), P() );
    this->exportResults( 0, u, v );

    auto i1 = mean( _range=markedfaces( mesh, "tip"  ), _expr=idv( u ) );
    LOG(INFO) << "deflection: " << i1 << "\n";

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
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="beam",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application beam;
    beam.add( new Feel::Beam<2,3>() );
    beam.add( new Feel::Beam<3,3>() );
    beam.run();
}
