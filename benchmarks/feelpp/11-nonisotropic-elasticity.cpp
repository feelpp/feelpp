/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <benchmark/benchmark.h>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

#include "tensor-notation-bench-utils.hpp"

using namespace Feel;
using namespace Feel::vf;

namespace
{

enum class ConstitutiveFormulation
{
    Tensor,
    Mandel,
    Voigt
};

template <ConstitutiveFormulation Formulation>
const char*
formulationName()
{
    if constexpr ( Formulation == ConstitutiveFormulation::Tensor )
        return "tensor";
    else if constexpr ( Formulation == ConstitutiveFormulation::Mandel )
        return "mandel";
    else
        return "voigt";
}

template <typename TrialExprT, typename TestExprT>
auto
orthotropicTensorElasticityForm( TrialExprT const& u, TestExprT const& v )
{
    auto epsu = symm_grad( u );
    auto epsv = symm_grad( v );
    auto sigma = benchmark_detail::orthotropicTensorStress2D( epsu );

    return component<0, 0>( sigma ) * component<0, 0>( epsv ) +
           cst( 2.0 ) * component<0, 1>( sigma ) * component<0, 1>( epsv ) +
           component<1, 1>( sigma ) * component<1, 1>( epsv );
}

auto
orthotropicMandelMatrix()
{
    return benchmark_detail::orthotropicMandelMatrix2D();
}

auto
orthotropicVoigtMatrix()
{
    return benchmark_detail::orthotropicVoigtMatrix2D();
}

template <typename TrialExprT, typename TestExprT>
auto
orthotropicMandelElasticityForm( TrialExprT const& u, TestExprT const& v )
{
    return ddot( orthotropicMandelMatrix(), symm_grad( u ), symm_grad( v ) );
}

template <typename TrialExprT, typename TestExprT>
auto
orthotropicVoigtElasticityForm( TrialExprT const& u, TestExprT const& v )
{
    return ddot<SymmetricTensorNotation::Voigt>( orthotropicVoigtMatrix(),
                                                 symm_grad( u ),
                                                 symm_grad( v ) );
}

template <ConstitutiveFormulation Formulation>
void
BM_NonisotropicElasticityAssembly( benchmark::State& state )
{
    auto mesh = benchmark_detail::makeStructuredMesh<2>( state.range( 0 ) );
    auto Xh = Pchv<1>( mesh );
    auto u = trial( Xh, "u" );
    auto v = test( Xh, "v" );

    size_type const ndofs = Xh->nLocalDof();

    for ( auto _ : state )
    {
        auto a = form2( _trial=Xh, _test=Xh );

        if constexpr ( Formulation == ConstitutiveFormulation::Tensor )
        {
            a = integrate( _range=elements( mesh ),
                           _expr=orthotropicTensorElasticityForm( u, v ) );
        }
        else if constexpr ( Formulation == ConstitutiveFormulation::Mandel )
        {
            a = integrate( _range=elements( mesh ),
                           _expr=orthotropicMandelElasticityForm( u, v ) );
        }
        else
        {
            a = integrate( _range=elements( mesh ),
                           _expr=orthotropicVoigtElasticityForm( u, v ) );
        }

        a.close();
        benchmark::DoNotOptimize( a.matrixPtr() );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations() * mesh->numElements() );
    state.SetLabel( std::string( formulationName<Formulation>() ) +
                    ":orthotropic:dim=2:ndofs=" + std::to_string( ndofs ) );
}

} // namespace

BENCHMARK_TEMPLATE( BM_NonisotropicElasticityAssembly, ConstitutiveFormulation::Tensor )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 );
BENCHMARK_TEMPLATE( BM_NonisotropicElasticityAssembly, ConstitutiveFormulation::Mandel )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 );
BENCHMARK_TEMPLATE( BM_NonisotropicElasticityAssembly, ConstitutiveFormulation::Voigt )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 );

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_nonisotropic_elasticity",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
