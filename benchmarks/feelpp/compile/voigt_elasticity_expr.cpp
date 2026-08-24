/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <int Dim>
auto
makeMesh()
{
    if constexpr ( Dim == 2 )
        return unitSquare();
    else
        return unitCube();
}

template <int Dim>
int
instantiate()
{
    auto mesh = makeMesh<Dim>();
    auto Xh = Pchv<1>( mesh );
    auto u = trial( Xh, "u" );
    auto v = test( Xh, "v" );
    constexpr double E = 1.3e6;
    constexpr double nu = 0.31;
    constexpr double lambda = E*nu/( ( 1.0 + nu )*( 1.0 - 2.0*nu ) );
    constexpr double mu = E/( 2.0*( 1.0 + nu ) );

    auto a = form2( _trial=Xh, _test=Xh );
    a = integrate( _range=elements( mesh ),
                   _expr=ddot<SymmetricTensorNotation::Voigt>(
                       isotropic_stiffness<Dim, SymmetricTensorNotation::Voigt>( lambda, mu ),
                       symm_grad( u ),
                       symm_grad( v ) ) );
    return 0;
}

} // namespace

int
main()
{
    return instantiate<2>() + instantiate<3>();
}
