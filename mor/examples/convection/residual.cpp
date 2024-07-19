/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include "convection.hpp"

void
ConvectionCrb::updateR( const vector_ptrtype& X, vector_ptrtype& R)
{
    LOG(INFO) << "[updateResidual] start\n";
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
    form1( _test = Xh, _vector=R ) =
        integrate ( _range = elements(mesh),
                    _expr = trans( gradv(u)*idv(u) )*id( v ) );

    // -- Heat Convection -- //
    form1( _test = Xh, _vector=R ) +=
        integrate ( _range = elements(mesh),
                    _expr = id(s)*(gradv(t)*idv(u)) );

    R->close();
    // add the linear part
    R->add( -1, F );
    // add the bilinear part
    R->addVector( X, D );

    double new_rez = R->l2Norm();
    if ( M_psiT )
    {
        if ( M_rez==-1 )
            M_rez=new_rez;
        M_delta = M_delta*M_rez/new_rez;
        Feel::cout<<"psiT : new residual="<<new_rez<<", new delta="<<M_delta <<std::endl;
        M_rez = new_rez;
    }

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
