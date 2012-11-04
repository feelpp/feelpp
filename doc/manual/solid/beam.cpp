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
#include <feel/feel.hpp>

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


namespace Feel
{
template<int nDim, int nOrder>
class Beam
    :
public Simget
{
    typedef Simget super;
public:

    // -- TYPEDEFS --
    static const uint16_type Dim = nDim;
    static const uint16_type feOrder = nOrder;
    static const uint16_type imOrder = nOrder;

    typedef double value_type;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef bases<Lagrange<feOrder, Vectorial> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Beam()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        beta( this->vm()["beta"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        M_bctype( this->vm()["bctype"].template as<int>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers()
    {
        LOG(INFO) << "[Beam] hsize = " << meshSize << "\n";
        LOG(INFO) << "[Beam] beta = " << beta << "\n";
        LOG(INFO) << "[Beam] bccoeff = " << bcCoeff << "\n";
        LOG(INFO) << "[Beam] bctype = " <<  M_bctype << "\n";
        LOG(INFO) << "[Beam] export = " << this->vm().count( "export" ) << "\n";

    }

    ~Beam()
    {
        std::map<std::string,std::pair<boost::timer,double> >::iterator it = timers.begin();
        std::map<std::string,std::pair<boost::timer,double> >::iterator en = timers.end();

        for ( ; it != en; ++it )
        {
            LOG(INFO) << it->first << " : " << it->second.second << " s elapsed\n";
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

    this->changeRepository( boost::format( "examples/solid/%1%/%2%/P%3%/h_%4%/" )
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
    LOG(INFO) << "lambda = " << lambda << "\n"
          << "mu     = " << mu << "\n"
          << "gravity= " << gravity << "\n";

    /*
     * Construction of the right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ y \f$ direction.
     */
    auto F = M_backend->newVector( Xh );
    F->zero();
    timers["assembly"].first.restart();

    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F0.m" );

    if ( Dim == 3 )
        form1( _test=Xh, _vector=F ) = integrate( elements( mesh ), trans( gravity*oneZ() )*id( v ) );
    else
        form1( _test=Xh, _vector=F ) = integrate( elements( mesh ), trans( gravity*oneY() )*id( v ) );

    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F1.m" );

    timers["assembly"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    auto D = M_backend->newMatrix( Xh, Xh );
    timers["assembly"].first.restart();
    auto deft = 0.5*( gradt( u )+trans( gradt( u ) ) );
    auto def = 0.5*( grad( v )+trans( grad( v ) ) );
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D );
    a = integrate( elements( mesh ),
                   lambda*divt( u )*div( v )  +
                   2*mu*trace( trans( deft )*def ) );

    if ( M_bctype == 1 ) // weak Dirichlet bc
    {
        auto Id = ( mat<nDim,nDim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) );
        a += integrate( markedfaces( mesh,1 ),
                        - trans( ( 2*mu*deft+lambda*trace( deft )*Id )*N() )*id( v )
                        - trans( ( 2*mu*def+lambda*trace( def )*Id )*N() )*idt( u )
                        + bcCoeff*trans( idt( u ) )*id( v )/hFace() );
    }

    if ( M_bctype == 0 )
        a += on( markedfaces( mesh,( nDim==2 )?1:23 ), u, F, constant( 0 )*one() );

    if ( this->vm().count( "export-matlab" ) )
    {
        F->printMatlab( "F2.m" );
        D->printMatlab( "elas.m" );
    }

    timers["assembly"].second += timers["assembly"].first.elapsed();

    M_backend->solve( _matrix=D, _solution=u, _rhs=F );

    v = vf::project( Xh, elements( Xh->mesh() ), P() );
    this->exportResults( 0, u, v );

    auto i1 = integrate( markedfaces( mesh,3 ), idv( u ) ).evaluate();
    LOG(INFO) << "deflection: " << i1/0.02 << "\n";



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

    LOG(INFO) << "||v||_2 = " << v.l2Norm() << " (P() before move)\n";
    LOG(INFO) << "||w||_2 = " << w.l2Norm() << " (P() after move)\n";


    w.add( -1.0, v );

    LOG(INFO) << "||u||_2 = " << w.l2Norm() << " (displacement u)\n";
    LOG(INFO) << "||w-v||_2 = " << w.l2Norm() << " (displacement w-v\n";

    w.add( -1.0, u );
    LOG(INFO) << "||(w-v)-u||_2 = " << w.l2Norm() << " (should be 0, ie u=w-v)\n";
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




