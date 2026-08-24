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
//! @author Feel++ Consortium
//! @date 03 Apr 2026
//!

#include <array>
#include <string>

#include <benchmark/benchmark.h>

#include <feel/feelvf/vf.hpp>

#include "tensor-notation-bench-utils.hpp"

using namespace Feel;
using namespace Feel::vf;

namespace
{

enum class ConstitutiveLaw
{
    Isotropic2D,
    Isotropic3D,
    Orthotropic2D
};

enum class ConstitutiveOperation
{
    StressAction,
    BilinearAction
};

enum class ConstitutiveVariant
{
    Manual,
    MandelTensor,
    MandelStorage,
    VoigtTensor,
    VoigtStorage
};

template <ConstitutiveLaw Law>
struct ConstitutiveLawTraits;

template <>
struct ConstitutiveLawTraits<ConstitutiveLaw::Isotropic2D>
{
    static constexpr int dim = 2;
    static constexpr double lambda = 2.7;
    static constexpr double mu = 1.9;
    static constexpr std::array<int, 2> benchmarkArgs = { 32, 64 };

    static std::string name() { return "isotropic-2d"; }

    template <typename TensorExprT>
    static auto manualStress( TensorExprT const& eps )
    {
        return benchmark_detail::isotropicTensorStress<dim>( lambda, mu, eps );
    }

    static auto mandelStiffness()
    {
        return isotropic_stiffness<dim>( lambda, mu );
    }

    static auto voigtStiffness()
    {
        return isotropic_stiffness<dim, SymmetricTensorNotation::Voigt>( lambda, mu );
    }
};

template <>
struct ConstitutiveLawTraits<ConstitutiveLaw::Isotropic3D>
{
    static constexpr int dim = 3;
    static constexpr double lambda = 1.8;
    static constexpr double mu = 2.4;
    static constexpr std::array<int, 2> benchmarkArgs = { 5, 6 };

    static std::string name() { return "isotropic-3d"; }

    template <typename TensorExprT>
    static auto manualStress( TensorExprT const& eps )
    {
        return benchmark_detail::isotropicTensorStress<dim>( lambda, mu, eps );
    }

    static auto mandelStiffness()
    {
        return isotropic_stiffness<dim>( lambda, mu );
    }

    static auto voigtStiffness()
    {
        return isotropic_stiffness<dim, SymmetricTensorNotation::Voigt>( lambda, mu );
    }
};

template <>
struct ConstitutiveLawTraits<ConstitutiveLaw::Orthotropic2D>
{
    static constexpr int dim = 2;
    static constexpr std::array<int, 2> benchmarkArgs = { 32, 64 };

    static std::string name() { return "orthotropic-2d"; }

    template <typename TensorExprT>
    static auto manualStress( TensorExprT const& eps )
    {
        return benchmark_detail::orthotropicTensorStress2D( eps );
    }

    static auto mandelStiffness()
    {
        return benchmark_detail::orthotropicMandelMatrix2D();
    }

    static auto voigtStiffness()
    {
        return benchmark_detail::orthotropicVoigtMatrix2D();
    }
};

template <ConstitutiveOperation Operation>
std::string
operationName()
{
    if constexpr ( Operation == ConstitutiveOperation::StressAction )
        return "stress";
    else
        return "bilinear";
}

template <ConstitutiveVariant Variant>
std::string
variantName()
{
    if constexpr ( Variant == ConstitutiveVariant::Manual )
        return "manual";
    else if constexpr ( Variant == ConstitutiveVariant::MandelTensor )
        return "mandel-tensor";
    else if constexpr ( Variant == ConstitutiveVariant::MandelStorage )
        return "mandel-storage";
    else if constexpr ( Variant == ConstitutiveVariant::VoigtTensor )
        return "voigt-tensor";
    else
        return "voigt-storage";
}

template <ConstitutiveLaw Law, ConstitutiveOperation Operation, ConstitutiveVariant Variant>
void
BM_TensorNotation( benchmark::State& state )
{
    using law_t = ConstitutiveLawTraits<Law>;

    auto mesh = benchmark_detail::makeStructuredMesh<law_t::dim>( state.range( 0 ) );
    auto eps = benchmark_detail::sampleSymmetricTensor<law_t::dim>();
    auto eta = benchmark_detail::sampleSymmetricTensor2<law_t::dim>();

    for ( auto _ : state )
    {
        auto value =
            [&]() -> double
            {
                if constexpr ( Operation == ConstitutiveOperation::StressAction )
                {
                    if constexpr ( Variant == ConstitutiveVariant::Manual )
                    {
                        auto sigma = law_t::manualStress( eps );
                        return integrate( _range=elements( mesh ), _expr=inner( sigma, sigma ) ).evaluate()( 0, 0 );
                    }
                    else if constexpr ( Variant == ConstitutiveVariant::MandelTensor )
                    {
                        auto sigma = contract( law_t::mandelStiffness(), eps );
                        return integrate( _range=elements( mesh ), _expr=inner( sigma, sigma ) ).evaluate()( 0, 0 );
                    }
                    else if constexpr ( Variant == ConstitutiveVariant::MandelStorage )
                    {
                        auto sigma = unmandel( contract( law_t::mandelStiffness(), mandel( eps ) ) );
                        return integrate( _range=elements( mesh ), _expr=inner( sigma, sigma ) ).evaluate()( 0, 0 );
                    }
                    else if constexpr ( Variant == ConstitutiveVariant::VoigtTensor )
                    {
                        auto sigma = contract<SymmetricTensorNotation::Voigt>( law_t::voigtStiffness(), eps );
                        return integrate( _range=elements( mesh ), _expr=inner( sigma, sigma ) ).evaluate()( 0, 0 );
                    }
                    else
                    {
                        auto sigma = unvoigt( contract<SymmetricTensorNotation::Voigt>( law_t::voigtStiffness(), voigt( eps ) ) );
                        return integrate( _range=elements( mesh ), _expr=inner( sigma, sigma ) ).evaluate()( 0, 0 );
                    }
                }
                else
                {
                    if constexpr ( Variant == ConstitutiveVariant::Manual )
                    {
                        return integrate( _range=elements( mesh ),
                                          _expr=inner( law_t::manualStress( eps ), eta ) ).evaluate()( 0, 0 );
                    }
                    else if constexpr ( Variant == ConstitutiveVariant::MandelTensor )
                    {
                        return integrate( _range=elements( mesh ),
                                          _expr=ddot( law_t::mandelStiffness(), eps, eta ) ).evaluate()( 0, 0 );
                    }
                    else if constexpr ( Variant == ConstitutiveVariant::MandelStorage )
                    {
                        return integrate( _range=elements( mesh ),
                                          _expr=ddot( law_t::mandelStiffness(), mandel( eps ), mandel( eta ) ) ).evaluate()( 0, 0 );
                    }
                    else if constexpr ( Variant == ConstitutiveVariant::VoigtTensor )
                    {
                        return integrate( _range=elements( mesh ),
                                          _expr=ddot<SymmetricTensorNotation::Voigt>( law_t::voigtStiffness(), eps, eta ) ).evaluate()( 0, 0 );
                    }
                    else
                    {
                        return integrate( _range=elements( mesh ),
                                          _expr=ddot<SymmetricTensorNotation::Voigt>( law_t::voigtStiffness(),
                                                                                       voigt( eps ),
                                                                                       voigt( eta ) ) ).evaluate()( 0, 0 );
                    }
                }
            }();

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations() * mesh->numElements() );
    state.SetLabel( variantName<Variant>() + ":" +
                    operationName<Operation>() + ":" +
                    law_t::name() );
}

template <ConstitutiveLaw Law, ConstitutiveOperation Operation, ConstitutiveVariant Variant>
void
registerTensorNotationBenchmark()
{
    using law_t = ConstitutiveLawTraits<Law>;

    auto* benchmarkHandle =
        benchmark::RegisterBenchmark( ( "tensor-notation/" + law_t::name() + "/" +
                                        operationName<Operation>() + "/" +
                                        variantName<Variant>() ).c_str(),
                                      &BM_TensorNotation<Law, Operation, Variant> );

    benchmarkHandle->Unit( benchmark::kMillisecond );
    for ( int arg : law_t::benchmarkArgs )
        benchmarkHandle->Arg( arg );
}

template <ConstitutiveLaw Law>
void
registerTensorNotationBenchmarksForLaw()
{
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::StressAction, ConstitutiveVariant::Manual>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::StressAction, ConstitutiveVariant::MandelTensor>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::StressAction, ConstitutiveVariant::MandelStorage>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::StressAction, ConstitutiveVariant::VoigtTensor>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::StressAction, ConstitutiveVariant::VoigtStorage>();

    registerTensorNotationBenchmark<Law, ConstitutiveOperation::BilinearAction, ConstitutiveVariant::Manual>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::BilinearAction, ConstitutiveVariant::MandelTensor>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::BilinearAction, ConstitutiveVariant::MandelStorage>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::BilinearAction, ConstitutiveVariant::VoigtTensor>();
    registerTensorNotationBenchmark<Law, ConstitutiveOperation::BilinearAction, ConstitutiveVariant::VoigtStorage>();
}

void
registerTensorNotationBenchmarks()
{
    registerTensorNotationBenchmarksForLaw<ConstitutiveLaw::Isotropic2D>();
    registerTensorNotationBenchmarksForLaw<ConstitutiveLaw::Isotropic3D>();
    registerTensorNotationBenchmarksForLaw<ConstitutiveLaw::Orthotropic2D>();
}

} // namespace

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_tensor_notation",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    registerTensorNotationBenchmarks();
    benchmark::RunSpecifiedBenchmarks();
}
