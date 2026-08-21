//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
#ifndef FEELPP_BENCHMARKS_TENSOR_NOTATION_BENCH_UTILS_HPP
#define FEELPP_BENCHMARKS_TENSOR_NOTATION_BENCH_UTILS_HPP

#include <numbers>

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel::benchmark_detail
{

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

template <int Dim, typename TensorExprT>
auto
isotropicTensorStress( double lambda, double mu, TensorExprT const& eps )
{
    using namespace Feel::vf;
    return cst( lambda )*trace( eps )*eye<Dim, Dim>() + cst( 2.0*mu )*eps;
}

template <int Dim>
auto
sampleSymmetricTensor()
{
    using namespace Feel::vf;

    if constexpr ( Dim == 2 )
    {
        return mat<2, 2>( 2.0*Px() + 3.0*Py(), -Px() + 4.0*Py(),
                          -Px() + 4.0*Py(), 5.0*Px() - 2.0*Py() );
    }
    else
    {
        return mat<3, 3>( 2.0*Px() + 3.0*Py() - Pz(), -Px() + 4.0*Py() + 2.0*Pz(), 3.0*Px() - Py() + Pz(),
                          -Px() + 4.0*Py() + 2.0*Pz(), 5.0*Px() - 2.0*Py() + 4.0*Pz(), 2.0*Px() + Py() - 3.0*Pz(),
                          3.0*Px() - Py() + Pz(), 2.0*Px() + Py() - 3.0*Pz(), -Px() + 2.0*Py() + 6.0*Pz() );
    }
}

template <int Dim>
auto
sampleSymmetricTensor2()
{
    using namespace Feel::vf;

    if constexpr ( Dim == 2 )
    {
        return mat<2, 2>( -Px() + 2.0*Py(), 3.0*Px() + Py(),
                          3.0*Px() + Py(), 4.0*Px() - Py() );
    }
    else
    {
        return mat<3, 3>( -Px() + 2.0*Py() + Pz(), 3.0*Px() + Py() - Pz(), Px() - 2.0*Py() + 4.0*Pz(),
                          3.0*Px() + Py() - Pz(), 4.0*Px() - Py() + 2.0*Pz(), -2.0*Px() + 3.0*Py() + Pz(),
                          Px() - 2.0*Py() + 4.0*Pz(), -2.0*Px() + 3.0*Py() + Pz(), 5.0*Px() + Py() - 3.0*Pz() );
    }
}

struct Orthotropic2DCoefficients
{
    static constexpr double c0000 = 8.0;
    static constexpr double c1111 = 6.0;
    static constexpr double c0011 = 1.5;
    static constexpr double c0001 = 0.6;
    static constexpr double c1101 = 0.4;
    static constexpr double c0101 = 2.2;
};

template <typename TensorExprT>
auto
orthotropicTensorStress2D( TensorExprT const& eps )
{
    using namespace Feel::vf;

    return mat<2, 2>(
        cst( Orthotropic2DCoefficients::c0000 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c0011 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c0001 ) * component<0, 1>( eps ),
        cst( Orthotropic2DCoefficients::c0001 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c1101 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ) * component<0, 1>( eps ),
        cst( Orthotropic2DCoefficients::c0001 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c1101 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ) * component<0, 1>( eps ),
        cst( Orthotropic2DCoefficients::c0011 ) * component<0, 0>( eps ) +
        cst( Orthotropic2DCoefficients::c1111 ) * component<1, 1>( eps ) +
        cst( 2.0 * Orthotropic2DCoefficients::c1101 ) * component<0, 1>( eps ) );
}

inline auto
orthotropicMandelMatrix2D()
{
    using namespace Feel::vf;
    constexpr double sqrt2 = std::numbers::sqrt2_v<double>;

    return mat<3, 3>(
        cst( Orthotropic2DCoefficients::c0000 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c0001 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c0001 ),
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( sqrt2 * Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c1111 ) );
}

inline auto
orthotropicVoigtMatrix2D()
{
    using namespace Feel::vf;

    return mat<3, 3>(
        cst( Orthotropic2DCoefficients::c0000 ),
        cst( 2.0 * Orthotropic2DCoefficients::c0001 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( Orthotropic2DCoefficients::c0001 ),
        cst( 2.0 * Orthotropic2DCoefficients::c0101 ),
        cst( Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c0011 ),
        cst( 2.0 * Orthotropic2DCoefficients::c1101 ),
        cst( Orthotropic2DCoefficients::c1111 ) );
}

} // namespace Feel::benchmark_detail

#endif
