//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme
//! @date 25 Mar 2026
//!

#include <benchmark/benchmark.h>

#include <numbers>

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

enum class ConversionKind
{
    Voigt,
    Mandel,
    UnVoigt,
    UnMandel
};

template <ConversionKind Kind>
const char*
conversionName()
{
    if constexpr ( Kind == ConversionKind::Voigt )
        return "voigt";
    else if constexpr ( Kind == ConversionKind::Mandel )
        return "mandel";
    else if constexpr ( Kind == ConversionKind::UnVoigt )
        return "unvoigt";
    else
        return "unmandel";
}

template <int Dim>
auto
makeStructuredMesh( int n )
{
    using mesh_t = MeshStructured<Hypercube<Dim>>;

    auto discretisation = nl::json::array();
    for ( int axis = 0; axis < Dim; ++axis )
        discretisation.push_back( n );

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", discretisation } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();
    return mesh;
}

template <bool UseHelper, ConversionKind Kind, int Dim>
void
BM_VoigtConversion( benchmark::State& state )
{
    constexpr auto sqrt2 = std::numbers::sqrt2_v<double>;
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto mesh = makeStructuredMesh<Dim>( state.range( 0 ) );

    if constexpr ( Dim == 2 )
    {
        auto S = mat<2, 2>( cst( 1.0 ), Px(),
                            Px(), cst( 2.0 ) + Py() );
        auto v = vec( cst( 1.0 ), Px(), cst( 2.0 ) + Py() );
        auto m = vec( cst( 1.0 ), cst( sqrt2 )*Px(), cst( 2.0 ) + Py() );

        auto manualVoigt = vec( S( 0, 0 ), S( 0, 1 ), S( 1, 1 ) );
        auto manualMandel = vec( S( 0, 0 ), cst( sqrt2 )*S( 0, 1 ), S( 1, 1 ) );
        auto manualUnVoigt = mat<2, 2>( v( 0, 0 ), v( 1, 0 ),
                                        v( 1, 0 ), v( 2, 0 ) );
        auto manualUnMandel = mat<2, 2>( m( 0, 0 ), cst( invSqrt2 )*m( 1, 0 ),
                                         cst( invSqrt2 )*m( 1, 0 ), m( 2, 0 ) );

        for ( auto _ : state )
        {
            double value = 0.0;

            if constexpr ( UseHelper && Kind == ConversionKind::Voigt )
                value = integrate( _range=elements( mesh ), _expr=inner( voigt( S ), voigt( S ) ) ).evaluate()( 0, 0 );
            else if constexpr ( UseHelper && Kind == ConversionKind::Mandel )
                value = integrate( _range=elements( mesh ), _expr=inner( mandel( S ), mandel( S ) ) ).evaluate()( 0, 0 );
            else if constexpr ( UseHelper && Kind == ConversionKind::UnVoigt )
                value = integrate( _range=elements( mesh ), _expr=inner( unvoigt( v ), unvoigt( v ) ) ).evaluate()( 0, 0 );
            else if constexpr ( UseHelper && Kind == ConversionKind::UnMandel )
                value = integrate( _range=elements( mesh ), _expr=inner( unmandel( m ), unmandel( m ) ) ).evaluate()( 0, 0 );
            else if constexpr ( Kind == ConversionKind::Voigt )
                value = integrate( _range=elements( mesh ), _expr=inner( manualVoigt, manualVoigt ) ).evaluate()( 0, 0 );
            else if constexpr ( Kind == ConversionKind::Mandel )
                value = integrate( _range=elements( mesh ), _expr=inner( manualMandel, manualMandel ) ).evaluate()( 0, 0 );
            else if constexpr ( Kind == ConversionKind::UnVoigt )
                value = integrate( _range=elements( mesh ), _expr=inner( manualUnVoigt, manualUnVoigt ) ).evaluate()( 0, 0 );
            else
                value = integrate( _range=elements( mesh ), _expr=inner( manualUnMandel, manualUnMandel ) ).evaluate()( 0, 0 );

            benchmark::DoNotOptimize( value );
            benchmark::ClobberMemory();
        }
    }
    else
    {
        auto S = mat<3, 3>( cst( 1.0 ), Px(), Py(),
                            Px(), cst( 2.0 ) + Px(), cst( 3.0 ) + Pz(),
                            Py(), cst( 3.0 ) + Pz(), cst( 4.0 ) + Px() );

        auto v = vec( cst( 1.0 ), Px(), Py(),
                      cst( 2.0 ) + Px(), cst( 3.0 ) + Pz(), cst( 4.0 ) + Px() );

        auto m = vec( cst( 1.0 ), cst( sqrt2 )*Px(), cst( sqrt2 )*Py(),
                      cst( 2.0 ) + Px(), cst( sqrt2 )*( cst( 3.0 ) + Pz() ), cst( 4.0 ) + Px() );

        auto manualVoigt = vec( S( 0, 0 ), S( 0, 1 ), S( 0, 2 ), S( 1, 1 ), S( 1, 2 ), S( 2, 2 ) );
        auto manualMandel = vec( S( 0, 0 ), cst( sqrt2 )*S( 0, 1 ), cst( sqrt2 )*S( 0, 2 ),
                                 S( 1, 1 ), cst( sqrt2 )*S( 1, 2 ), S( 2, 2 ) );
        auto manualUnVoigt = mat<3, 3>( v( 0, 0 ), v( 1, 0 ), v( 2, 0 ),
                                        v( 1, 0 ), v( 3, 0 ), v( 4, 0 ),
                                        v( 2, 0 ), v( 4, 0 ), v( 5, 0 ) );
        auto manualUnMandel = mat<3, 3>( m( 0, 0 ), cst( invSqrt2 )*m( 1, 0 ), cst( invSqrt2 )*m( 2, 0 ),
                                         cst( invSqrt2 )*m( 1, 0 ), m( 3, 0 ), cst( invSqrt2 )*m( 4, 0 ),
                                         cst( invSqrt2 )*m( 2, 0 ), cst( invSqrt2 )*m( 4, 0 ), m( 5, 0 ) );

        for ( auto _ : state )
        {
            double value = 0.0;

            if constexpr ( UseHelper && Kind == ConversionKind::Voigt )
                value = integrate( _range=elements( mesh ), _expr=inner( voigt( S ), voigt( S ) ) ).evaluate()( 0, 0 );
            else if constexpr ( UseHelper && Kind == ConversionKind::Mandel )
                value = integrate( _range=elements( mesh ), _expr=inner( mandel( S ), mandel( S ) ) ).evaluate()( 0, 0 );
            else if constexpr ( UseHelper && Kind == ConversionKind::UnVoigt )
                value = integrate( _range=elements( mesh ), _expr=inner( unvoigt( v ), unvoigt( v ) ) ).evaluate()( 0, 0 );
            else if constexpr ( UseHelper && Kind == ConversionKind::UnMandel )
                value = integrate( _range=elements( mesh ), _expr=inner( unmandel( m ), unmandel( m ) ) ).evaluate()( 0, 0 );
            else if constexpr ( Kind == ConversionKind::Voigt )
                value = integrate( _range=elements( mesh ), _expr=inner( manualVoigt, manualVoigt ) ).evaluate()( 0, 0 );
            else if constexpr ( Kind == ConversionKind::Mandel )
                value = integrate( _range=elements( mesh ), _expr=inner( manualMandel, manualMandel ) ).evaluate()( 0, 0 );
            else if constexpr ( Kind == ConversionKind::UnVoigt )
                value = integrate( _range=elements( mesh ), _expr=inner( manualUnVoigt, manualUnVoigt ) ).evaluate()( 0, 0 );
            else
                value = integrate( _range=elements( mesh ), _expr=inner( manualUnMandel, manualUnMandel ) ).evaluate()( 0, 0 );

            benchmark::DoNotOptimize( value );
            benchmark::ClobberMemory();
        }
    }
    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( UseHelper ? "helper" : "manual" ) + ":" + conversionName<Kind>() + ":dim=" + std::to_string( Dim ) );
}

} // namespace

BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::Voigt, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::Voigt, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::Mandel, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::Mandel, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::UnVoigt, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::UnVoigt, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::UnMandel, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::UnMandel, 2 )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );

BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::Voigt, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::Voigt, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::Mandel, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::Mandel, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::UnVoigt, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::UnVoigt, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, false, ConversionKind::UnMandel, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );
BENCHMARK_TEMPLATE( BM_VoigtConversion, true, ConversionKind::UnMandel, 3 )->Unit( benchmark::kMillisecond )->Arg( 16 )->Arg( 24 );

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_voigt",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
