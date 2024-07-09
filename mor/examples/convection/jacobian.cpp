/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include "convection.hpp"

void ConvectionCrb ::updateJ( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    LOG(INFO) << "[updateJacobian] start\n";

    auto mesh = Xh->mesh();
    element_type U = Xh->element( "u" );
    element_type V = Xh->element( "v" );
    U = *X;
    auto u = U.template element<0>(); // velocity
    auto v = V.template element<0>();
    auto t = U.template element<3>(); // temperature
    auto s = V.template element<3>();

    if (!J)
        J = M_backend->newMatrix( _test=Xh, _trial=Xh );
    else
        J->zero();

    auto f2 = form2( _test=Xh,_trial=Xh, _matrix=J );

    // -- Fluid NL terms : 2 terms !!! -- //
     f2 += integrate( _range = elements(mesh),
                      _expr = trans( id(v) )*gradv(u)*idt(u)
                      +trans( id(v) )*gradt(u)*idv(u) );

    // -- Temperature NL terms : 2 terms !!! -- //
    f2 += integrate ( _range = elements(mesh),
                      _expr = id(s)*(gradv(t)*idt(u))
                      +id(s)*(gradt(t)*idv(u)) );

    if ( M_psiT )
    {
        auto norm_u=max(1e-12,sqrt(inner(idv(u),idv(u))));
        f2 += integrate( _range = elements(mesh),
                         _expr = norm_u/hMin()/(M_delta)*inner(idt(u),id(u)) );
        f2 += integrate( _range = elements(mesh),
                         _expr = norm_u/hMin()/(M_delta)*inner(idt(t),id(s)) );
    }


    J->close();
    J->addMatrix(1.,D);
}

typename ConvectionCrb::sparse_matrix_ptrtype
ConvectionCrb::jacobian( const element_type& X )
{
    sparse_matrix_ptrtype J;
    vector_ptrtype XX( M_backend->newVector( Xh ) );

    J = M_backend->newMatrix( _test=Xh, _trial=Xh );
    *XX = X;

    updateJ( XX,  J);
    return J;
}
