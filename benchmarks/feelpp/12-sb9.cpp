/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <benchmark/benchmark.h>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelvf/vf.hpp>

#include <Eigen/Core>

#include <iomanip>
#include <sstream>

using namespace Feel;
using namespace Feel::vf;

namespace
{

using mesh_type = Mesh<Hypercube<3>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;
using shell_vector_type = Eigen::Matrix<double, 3, 1>;

enum class Sb9AssemblyMode
{
    Membrane,
    Shear,
    Full
};

const char*
modeName( Sb9AssemblyMode mode )
{
    switch ( mode )
    {
    case Sb9AssemblyMode::Membrane:
        return "membrane";
    case Sb9AssemblyMode::Shear:
        return "shear";
    default:
        return "full";
    }
}

mesh_ptrtype
createShellPatchMesh( int nxy )
{
    shell_vector_type const origin( 0.0, 0.0, 0.0 );
    shell_vector_type const tangent1( 2.0, 0.0, 0.0 );
    shell_vector_type const tangent2( 0.0, 1.0, 0.0 );
    shell_vector_type const thicknessVector( 0.0, 0.0, 0.15 );

    auto const p1 = origin - 0.5 * thicknessVector;
    auto const p2 = p1 + tangent1;
    auto const p3 = p2 + tangent2;
    auto const p4 = p1 + tangent2;

    std::ostringstream geoDesc;
    geoDesc << std::setprecision( 16 );
    geoDesc << "Mesh.RecombineAll = 1;\n"
            << "Point(1) = {" << p1( 0 ) << ", " << p1( 1 ) << ", " << p1( 2 ) << ", 1};\n"
            << "Point(2) = {" << p2( 0 ) << ", " << p2( 1 ) << ", " << p2( 2 ) << ", 1};\n"
            << "Point(3) = {" << p3( 0 ) << ", " << p3( 1 ) << ", " << p3( 2 ) << ", 1};\n"
            << "Point(4) = {" << p4( 0 ) << ", " << p4( 1 ) << ", " << p4( 2 ) << ", 1};\n"
            << "Line(1) = {1, 2};\n"
            << "Line(2) = {2, 3};\n"
            << "Line(3) = {3, 4};\n"
            << "Line(4) = {4, 1};\n"
            << "Line Loop(1) = {1, 2, 3, 4};\n"
            << "Plane Surface(1) = {1};\n"
            << "Transfinite Line {1, 3} = " << ( nxy + 1 ) << ";\n"
            << "Transfinite Line {2, 4} = " << ( nxy + 1 ) << ";\n"
            << "Transfinite Surface {1} = {1, 2, 3, 4};\n"
            << "Recombine Surface {1};\n"
            << "out[] = Extrude {" << thicknessVector( 0 ) << ", "
            << thicknessVector( 1 ) << ", " << thicknessVector( 2 ) << "} {\n"
            << "  Surface{1};\n"
            << "  Layers{1};\n"
            << "  Recombine;\n"
            << "};\n"
            << "Physical Volume(\"Shell\") = {out[1]};\n";

    Environment::changeRepository( _directory=boost::format( "benchmarks/feelpp/%1%/%2%/" )
                                   % Environment::about().appName()
                                   % ( "nxy_" + std::to_string( nxy ) ) );

    return createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="sb9_patch_" + std::to_string( nxy ) + ".geo",
                                      _desc=geoDesc.str(),
                                      _dim=3,
                                      _order=1,
                                      _h=1.0 ),
                           _force_rebuild=true );
}

template <Sb9AssemblyMode Mode>
void
BM_Sb9Assembly( benchmark::State& state )
{
    auto mesh = createShellPatchMesh( state.range( 0 ) );
    auto Uh = Pchv<1>( mesh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    constexpr double lambda = 2.3;
    constexpr double mu = 1.1;
    auto C = isotropic_stiffness<3>( lambda, mu );

    for ( auto _ : state )
    {
        auto a = form2( _trial=Uh, _test=Uh );

        if constexpr ( Mode == Sb9AssemblyMode::Membrane )
        {
            a = integrate( _range=elements( mesh ),
                           _expr=ddot( C,
                                       sb9MembraneBending( u, zeta() ),
                                       sb9MembraneBending( v, zeta() ) ) );
        }
        else if constexpr ( Mode == Sb9AssemblyMode::Shear )
        {
            a = integrate( _range=elements( mesh ),
                           _expr=ddot( C,
                                       sb9Shear( u, cst( 1.0 ) ),
                                       sb9Shear( v, cst( 1.0 ) ) ) );
        }
        else
        {
            a = integrate( _range=elements( mesh ),
                           _expr=ddot( C,
                                       sb9MembraneBending( u, zeta() ),
                                       sb9MembraneBending( v, zeta() ) ) +
                                 ddot( C,
                                       sb9Pinching( u ),
                                       sb9Pinching( v ) ) +
                                 ddot( C,
                                       sb9Shear( u, cst( 1.0 ) ),
                                       sb9Shear( v, cst( 1.0 ) ) ) );
        }

        a.close();
        benchmark::DoNotOptimize( a.matrixPtr() );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations() * mesh->numElements() );
    state.SetLabel( std::string( "sb9:" ) + modeName( Mode ) +
                    ":nxy=" + std::to_string( state.range( 0 ) ) +
                    ":nelems=" + std::to_string( mesh->numElements() ) );
}

} // namespace

BENCHMARK_TEMPLATE( BM_Sb9Assembly, Sb9AssemblyMode::Membrane )->Unit( benchmark::kMillisecond )->Arg( 4 )->Arg( 8 );
BENCHMARK_TEMPLATE( BM_Sb9Assembly, Sb9AssemblyMode::Shear )->Unit( benchmark::kMillisecond )->Arg( 4 )->Arg( 8 );
BENCHMARK_TEMPLATE( BM_Sb9Assembly, Sb9AssemblyMode::Full )->Unit( benchmark::kMillisecond )->Arg( 4 )->Arg( 8 );

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_sb9",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
