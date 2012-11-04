/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-03-15

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file dofpoints.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-03-15
 */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>

#if defined(__LP64__) || defined(_LP64)

#define hrtime   unsigned long
#define hrtime_t unsigned long

static inline hrtime gethrtime()
{
    hrtime t;
    unsigned int _a_, _d_;
    asm volatile( "rdtsc" : "=a" ( _a_ ), "=d" ( _d_ ) );
    t = ( ( unsigned long ) _a_ ) | ( ( ( unsigned long )_d_ ) << 32 );
    return t;
}

#else

#define hrtime   unsigned long long
#define hrtime_t unsigned long long

static inline hrtime gethrtime()
{
    hrtime t;
    asm volatile( "rdtsc" : "=A" ( t ) );
    return t;
}

#endif

template<typename MeshT, typename EltT>
FEELPP_DONT_INLINE double
f( MeshT& mesh, EltT const& u, bool use_tbb, int gs, std::string const& part )
{
    using namespace Feel;
    using namespace Feel::vf;

    //return integrate( _range=elements(mesh), _quad=_Q<10>(), _expr=trans(idv(u))*(sym(gradv(u))*idv(u)), _use_tbb=use_tbb ).evaluate()(0,0);
    //return integrate( _range=elements(mesh), _quad=_Q<10>(), _expr=(sym(hessv(u))*idv(u)), _use_tbb=use_tbb ).evaluate()(0,0);
    //return integrate( _range=elements(mesh), _quad=_Q<20>(), _expr=cos(idv(u))*sin(idv(u))+trace(P()*trans(P()))+exp(Px())*sin(Px()*Py()),
    return integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cos( Px() )*Py()*exp( Px() ),
                      _use_tbb=use_tbb, _grainsize=gs, _partitioner=part ).evaluate()( 0,0 );
}
template<typename SpaceT, typename EltT, typename MatrixT>
FEELPP_DONT_INLINE void
g( SpaceT& Xh, EltT const& u, MatrixT& M, bool use_tbb, int gs, std::string const& part )
{
    using namespace Feel;
    using namespace Feel::vf;

#if 0
    form2( _test=Xh, _trial=Xh, _matrix=M, _init=true ) =
        integrate( _range=elements( Xh->mesh() ), _expr=idt( u )*id( u ),
                   _use_tbb=use_tbb, _grainsize=gs, _partitioner=part  );
    M->close();
#else
    form1( _test=Xh, _vector=M, _init=true ) =
        integrate( _range=elements( Xh->mesh() ),
                   _expr=cos( Px() )*sin( Py() )*id( u ),
                   //_expr=id(u),
                   _use_tbb=use_tbb, _grainsize=gs, _partitioner=part  );
    M->close();
#endif

}

int main( int argc, char** argv )
{

    double hsize = 2;
    int nt = 2;
    int gs = 100;
    std::string part = "auto";
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options()
    ( "help", "produce help message" )
    ( "hsize", po::value<double>( &hsize )->default_value( 2 ), "h size" )
    ( "nt", po::value<int>( &nt )->default_value( 16 ), "nt" )
    ( "gs", po::value<int>( &gs )->default_value( 100 ), "grainsize" )
    ( "part", po::value<std::string>( &part )->default_value( "auto" ), "thread partitioner (auto,affinity,simple)" )
    ;
    desc.add( Feel::feel_options() );
    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );
#if defined( FEELPP_HAS_TBB )
    tbb::task_scheduler_init init( nt );
    std::cout << "is_active: " << init.is_active() << "\n";
    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<3,1> > mesh_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="simplex",
                                        _usenames=true,
                                        _shape="simplex",
                                        _dim=3,
                                        _order=1,
                                        _h=hsize ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    typedef FunctionSpace<mesh_type, Lagrange<10,Scalar> > space_type;
    auto Xh = space_type::New ( mesh  );
    auto u = vf::project( _space=Xh, _range=elements( mesh ), _expr=Px() );
    auto I = integrate( elements( mesh ), cst( 1. ), _use_tbb=false ).evaluate();

    tbb::tick_count t0 = tbb::tick_count::now();
    double r = f( mesh,u,false,gs,part );
    tbb::tick_count t1= tbb::tick_count::now();
    double singlethread_time2 = ( t1-t0 ).seconds();

    t0 = tbb::tick_count::now();
    r = f( mesh,u,true,gs,part );
    t1 = tbb::tick_count::now();
    double tbb_time2 = ( t1-t0 ).seconds();


    cout << "1T  tbbtime:  " << singlethread_time2 << " ticks" << endl;
    cout << "TBB  tbbtime: " << tbb_time2 << " ticks" << endl;
    cout << "speedup: " << singlethread_time2/tbb_time2 << " " << endl;

    u = vf::project( _space=Xh, _range=elements( mesh ), _expr=cst( 1.0 ) );

    auto backend = Backend<double>::build( vm );
    //auto M = backend->newMatrix( Xh, Xh );
    auto M = backend->newVector( Xh );
    t0 = tbb::tick_count::now();
    g( Xh, u, M, false, gs, part );
    t1 = tbb::tick_count::now();
    singlethread_time2 = ( t1-t0 ).seconds();
    M->printMatlab( "M1.m" );
    double i1 = inner_product( *M, u );
    std::cout << "inner1 = "  << i1 << "\n";

    t0 = tbb::tick_count::now();
    g( Xh, u, M, true, gs, part );
    t1 = tbb::tick_count::now();
    tbb_time2 = ( t1-t0 ).seconds();
    M->printMatlab( "M2.m" );
    double i2 = inner_product( *M, u );
    std::cout << "inner2 = "  << inner_product( *M, u ) << "\n";
    std::cout << "|i1=i2| = "  << std::setprecision( 16 ) << std::scientific << math::abs( i1-i2 ) << "\n";
    cout << "1T  tbbtime:  " << singlethread_time2 << " ticks" << endl;
    cout << "TBB  tbbtime: " << tbb_time2 << " ticks" << endl;
    cout << "speedup: " << singlethread_time2/tbb_time2 << " " << endl;
#endif

}
