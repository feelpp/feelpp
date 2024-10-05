/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "../convection.hpp"

using namespace Feel::vf;

void ConvectionCrb::initModel()
{
    // -- MESH SETTING -- //
    auto mesh  =  loadMesh( _mesh=new mesh_type );

    // -- FUNCTION SPACES AND ELEMENTS SETTING -- //
    this->setFunctionSpaces( space_type::New(mesh) );
    pT = element_ptrtype( new element_type( Xh ) );

    Feel::cout << "Number of dof : "<< Xh->nDof() << std::endl;
    LOG(INFO) << "Number of dof : "<< Xh->nDof();
    element_type U( Xh, "u" );
    element_type V( Xh, "v" );

    element_u_type u = U. element<0>(); // velocity
    element_u_type v = V. element<0>();
    element_p_type p = U. element<1>(); // pressure
    element_p_type q = V. element<1>();
    element_l_type xi = U. element<2>(); // lagrange
    element_l_type eta = V. element<2>();
    element_t_type t = U. element<3>(); // temperature
    element_t_type s = V. element<3>();

    // -- PARAMETERS SETTING -- //
    Dmu->setDimension( 2 );
    auto mu_min = Dmu->element();
    mu_min << doption( "Grmin" ), doption( "Prmin" );
    auto mu_max = Dmu->element();
    mu_max << doption( "Grmax"), doption( "Prmax" );
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );

    // -- READ PARAMETERS -- //
    double gamma = doption("penalbc");
    double expansion = 1;

    // -- INITIALIZE MATRIX -- //
    // lhs
    M_Aqm.resize( Qa() );
    for (int ii=0; ii<Qa(); ii++)
    {
        M_Aqm[ii].resize(1);
        M_Aqm[ii][0] = M_backend->newMatrix( _test = Xh, _trial= Xh );
    }
    // rhs
    M_Fqm.resize( Nl() );
    for (int ii=0; ii<Nl(); ii++)
    {
        M_Fqm[ii].resize( Ql(ii) );
        for (int jj=0; jj<Ql(ii); jj++)
        {
            M_Fqm[ii][jj].resize(1);
            M_Fqm[ii][jj][0] = M_backend->newVector( Xh );
        }
    }
    // others
    D = M_backend->newMatrix(  _test = Xh, _trial= Xh);
    F = M_backend->newVector( Xh );
    M_A_tril = M_backend->newMatrix(  _test = Xh, _trial= Xh );

    // -- BUILD MATRIX -- //
    // fluid diffusion
    // M_Aqm[0] = C = grad(u)*grav(v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0], _init=true )
        = integrate( _range=elements( mesh ),
                     _expr=trace( gradt( u )*trans( grad( v ) ) )  );
    // pressure-velocity terms
    // M_Aqm[2] = -B = -p*div(v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0], _init=true )
        = integrate ( _range=elements( mesh ),
                      _expr = -idt( p )*div( v ) );
    // M_Aqm[2] += B^t = q*div(u)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate ( _range=elements( mesh ),
                       _expr = divt( u )*id( q ) );
    // heat diffusion
    // M_Aqm[1] = G = grad(t)*grad(s)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0], _init=true )
        = integrate( _range=elements( mesh ),
                     _expr=gradt( t )*trans( grad( s ) ) );

    // buyoancy forces c(theta,v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        +=integrate( _range = elements( mesh ),
                     _expr = -expansion*idt( t )*( trans( oneY() )*id( v ) ) );

    // multipliers for zero-mean pressure
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate ( _range = elements( mesh ),
                       _expr = id( q )*idt( xi )
                       +idt( p )*id( eta ) );

    // weak Dirichlet condition at the walls (u=0)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "wallInsulated" ),
                       _expr = -trans( gradt( u )*N() )*id( v )
                       -trans( grad( v )*N() )*idt( u )
                       +gamma*trans( idt( u ) )*id( v )/hFace() );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate ( _range = markedfaces( mesh, "wallInsulated" ),
                       _expr = -trans( -idt( p )*N() )*id( v )
                       -trans( -id( q )*N() )*idt( u ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "wallRight" ),
                       _expr = -trans( gradt( u )*N() )*id( v )
                       -trans( grad( v )*N() )*idt( u )
                       +gamma*trans( idt( u ) )*id( v )/hFace() );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate ( _range = markedfaces( mesh, "wallRight" ),
                       _expr = -trans( -idt( p )*N() )*id( v )
                       -trans( -id( q )*N() )*idt( u ) );
   form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "wallLeft" ),
                       _expr = -trans( gradt( u )*N() )*id( v )
                       -trans( grad( v )*N() )*idt( u )
                       +gamma*trans( idt( u ) )*id( v )/hFace() );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate ( _range = markedfaces( mesh, "wallLeft" ),
                       _expr = -trans( -idt( p )*N() )*id( v )
                       -trans( -id( q )*N() )*idt( u ) );


    // weak Dirichlet on temperature (T=0|left wall)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0] )
        += integrate ( _range = markedfaces( mesh, "wallLeft" ),
                       _expr = - gradt( t )*N()*id( s )
                       - grad( s )*N()*idt( t )
                       + gamma*idt( t )*id( s )/hFace() );

    M = M_backend->newMatrix(  _test = Xh, _trial= Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M, _init=true ) =
        integrate( _range = elements( mesh ),
                   _expr = trans( id( v ) )*idt( u ) + trace( grad( v )*trans( gradt( u ) ) )
                   + id( q )*idt( p ) + grad( q )*trans( gradt( p ) )
                   + id( s )*idt( t ) + grad( s )*trans( gradt( t ) ) );
    M->close();

    form1( _test = Xh, _vector=M_Fqm[0][0][0] ) =
        integrate (_range =  markedfaces( mesh, "wallRight"),
                   _expr =  id( s )  );
    M_Fqm[0][0][0]->close();

    //output : \int_{\Omega} T
    form1( _test=Xh, _vector=M_Fqm[1][0][0] ) =
        integrate( _range = markedfaces( mesh,"wallRight" ),
                   _expr = id(s) );//*(1.0/area) ) ;
    M_Fqm[1][0][0]->close();

    LOG(INFO) << "Natural Convection : Cavity Model is initilized";
}
