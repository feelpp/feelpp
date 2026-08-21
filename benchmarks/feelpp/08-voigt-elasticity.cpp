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

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

enum class ElasticityFormulation
{
    Tensor,
    Mandel,
    Voigt
};

template <ElasticityFormulation Formulation>
const char*
formulationName()
{
    if constexpr ( Formulation == ElasticityFormulation::Tensor )
        return "tensor";
    else if constexpr ( Formulation == ElasticityFormulation::Mandel )
        return "mandel";
    else
        return "voigt";
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

template <int Dim, typename TrialExprT, typename TestExprT>
auto
tensorElasticityForm( TrialExprT const& u, TestExprT const& v, double lambda, double mu )
{
    auto epsu = symm_grad( u );
    auto epsv = symm_grad( v );
    return cst( lambda )*trace( epsu )*trace( epsv ) + cst( 2.0*mu )*inner( epsu, epsv );
}

template <int Dim, typename TrialExprT, typename TestExprT>
auto
mandelElasticityForm( TrialExprT const& u, TestExprT const& v, double lambda, double mu )
{
    auto C = isotropic_stiffness<Dim>( lambda, mu );
    return ddot( C, symm_grad( u ), symm_grad( v ) );
}

template <int Dim, typename TrialExprT, typename TestExprT>
auto
voigtElasticityForm( TrialExprT const& u, TestExprT const& v, double lambda, double mu )
{
    auto C = isotropic_stiffness<Dim, SymmetricTensorNotation::Voigt>( lambda, mu );
    return ddot<SymmetricTensorNotation::Voigt>( C, symm_grad( u ), symm_grad( v ) );
}

template <ElasticityFormulation Formulation, int Dim>
void
BM_ElasticityAssembly( benchmark::State& state )
{
    auto mesh = makeStructuredMesh<Dim>( state.range( 0 ) );
    auto Xh = Pchv<1>( mesh );
    auto u = trial( Xh, "u" );
    auto v = test( Xh, "v" );

    constexpr double E = 1.3e6;
    constexpr double nu = 0.31;
    constexpr double lambda = E*nu/( ( 1.0 + nu )*( 1.0 - 2.0*nu ) );
    constexpr double mu = E/( 2.0*( 1.0 + nu ) );

    size_type ndofs = Xh->nLocalDof();

    for ( auto _ : state )
    {
        auto a = form2( _trial=Xh, _test=Xh );

        if constexpr ( Formulation == ElasticityFormulation::Tensor )
        {
            a = integrate( _range=elements( mesh ),
                           _expr=tensorElasticityForm<Dim>( u, v, lambda, mu ) );
        }
        else if constexpr ( Formulation == ElasticityFormulation::Mandel )
        {
            a = integrate( _range=elements( mesh ),
                           _expr=mandelElasticityForm<Dim>( u, v, lambda, mu ) );
        }
        else
        {
            a = integrate( _range=elements( mesh ),
                           _expr=voigtElasticityForm<Dim>( u, v, lambda, mu ) );
        }

        a.close();
        benchmark::DoNotOptimize( a.matrixPtr() );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( formulationName<Formulation>() ) +
                    ":dim=" + std::to_string( Dim ) +
                    ":ndofs=" + std::to_string( ndofs ) );
}

} // namespace

BENCHMARK_TEMPLATE( BM_ElasticityAssembly, ElasticityFormulation::Tensor, 2 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 );
BENCHMARK_TEMPLATE( BM_ElasticityAssembly, ElasticityFormulation::Mandel, 2 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 );
BENCHMARK_TEMPLATE( BM_ElasticityAssembly, ElasticityFormulation::Voigt, 2 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 );

BENCHMARK_TEMPLATE( BM_ElasticityAssembly, ElasticityFormulation::Tensor, 3 )->Unit( benchmark::kMillisecond )->Arg( 5 )->Arg( 6 );
BENCHMARK_TEMPLATE( BM_ElasticityAssembly, ElasticityFormulation::Mandel, 3 )->Unit( benchmark::kMillisecond )->Arg( 5 )->Arg( 6 );
BENCHMARK_TEMPLATE( BM_ElasticityAssembly, ElasticityFormulation::Voigt, 3 )->Unit( benchmark::kMillisecond )->Arg( 5 )->Arg( 6 );

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_voigt_elasticity",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
