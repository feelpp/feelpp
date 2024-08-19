/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "convection.hpp"

ConvectionCrb::ConvectionCrb():
    super_type( name() ),
    M_psiT( false ),
    M_delta( 0 ),
    M_rez(-1),
    M_backend( backend() )
{}

void
ConvectionCrb::update( parameter_type const& mu )
{
    D->zero();
    for ( size_type q = 0; q < M_betaAqm.size(); ++q )
        for ( size_type m = 0; m < mMaxA(q); ++m )
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
    D->close();

    F->zero();
    for ( size_type q = 0; q < M_Fqm[0].size(); ++q )
        for ( size_type m = 0; m < mMaxF(0,q); ++m )
            F->add( M_betaFqm[0][q][m], M_Fqm[0][q][m] );
    F->close();
}

/**
 * \brief solve the model for parameter \p mu
 * \param mu the model parameter
 */
typename ConvectionCrb::element_type
ConvectionCrb::solve( parameter_type const& mu )
{
    this->solve( mu, pT );
    return *pT;
}

/**
 * Solve linear problem M u=f
 */
void
ConvectionCrb::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    M_backend->solve( _matrix=M,  _solution=u, _rhs=f );
}


/**
 * returns the scalar product of the std::shared_ptr vector x and
 * std::shared_ptr vector y
 */
double
ConvectionCrb::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
/**
 * returns the scalar product of the std::shared_ptr vector x and
 * std::shared_ptr vector y
 */
double
ConvectionCrb::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}


/**
 * specific interface for OpenTURNS
 *
 * \param X input vector of size N
 * \param N size of input vector X
 * \param Y input vector of size P
 * \param P size of input vector Y
 */
void
ConvectionCrb::run( const double * X, unsigned long N, double * Y, unsigned long P )
{}


typename ConvectionCrb ::sparse_matrix_ptrtype
ConvectionCrb::computeTrilinearForm( const element_type& X )
{
    auto mesh = Xh->mesh();
    element_type U = Xh->element( "U" );
    element_type V = Xh->element( "v" );
    U = X;

    auto u = U. element<0>(); // velocity
    auto v = V. element<0>();
    auto t = U. element<3>(); // temperature
    auto s = V. element<3>();

    // -- Fluid Convection -- //
    form2( _test=Xh, _trial=Xh, _matrix=M_A_tril ) =
        integrate ( _range = elements(mesh),
                    _expr = trans( id(v) )*gradt(v)*idv(u)  );

    // -- Heat Convection -- //
    form2( _test=Xh, _trial=Xh, _matrix=M_A_tril ) +=
        integrate ( _range = elements(mesh),
                    _expr = idv(t)*(gradt(s)*id(v)) );

    return M_A_tril;
}
