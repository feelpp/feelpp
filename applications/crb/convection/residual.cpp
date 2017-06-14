/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include "convection.hpp"

void
ConvectionCrb::updateR( const vector_ptrtype& X, vector_ptrtype& R)
{
    mesh_ptrtype mesh = Xh->mesh();
    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
    U = *X;

    auto u = U.template element<0>(); // velocity
    auto v = V.template element<0>();
    auto t = U.template element<3>(); // temperature
    auto s = V.template element<3>();

    R->zero();
    // -- NS Convection -- //
    form1( Xh, _vector=R ) =
        integrate ( elements(mesh),
                    trans( gradv(u)*idv(u) )*id( v ) );

    // -- Heat Convection -- //
    form1( Xh, _vector=R ) +=
        integrate ( elements(mesh),
                    id(s)*(gradv(t)*idv(u)) );

    R->close();
    // add the linear part
    R->add( -1, F );
    // add the bilinear part
    R->addVector( X, D );
}


typename ConvectionCrb::vector_ptrtype
ConvectionCrb::residual( const element_type& X )
{
    vector_ptrtype R ( M_backend->newVector( Xh ) );
    vector_ptrtype XX( M_backend->newVector( Xh ) );

    *XX = X;
    updateR( XX, R);
    return R;
}
