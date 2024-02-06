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

    A_OUT << "Number of dof : "<< Xh->nDof() << std::endl;
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
    mu_min << doption( "TinletMin" ), doption("UinletMin");
    auto mu_max = Dmu->element();
    mu_max << doption( "TinletMax" ), doption("UinletMax");
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );

    double penalbc = doption("penalbc");
    double nu = doption("nu");
    double rho = doption("rho");
    double k = doption("k");
    double expansion = 9.81*3.4112e-3;
    double Cp = 1000;

    M_psiT = boption("psiT");
    if ( M_psiT )
        M_delta = doption("delta0");

    // -- INITIALIZE MATRIX -- //
    // lhs
    M_Aqm.resize( Qa() );
    M_Aqm[0].resize(1);
    M_Aqm[0][0] = M_backend->newMatrix( _trial = Xh, _test =  Xh );

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
    D = M_backend->newMatrix( _trial = Xh, _test = Xh );
    F = M_backend->newVector( Xh );
    M_A_tril = M_backend->newMatrix( _trial = Xh, _test = Xh );

    // -- BUILD MATRIX -- //
    // fluid diffusion
    auto f2 = form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] );

    f2 += integrate( _range = elements( mesh ),
                    _expr = nu*trace( gradt( u )*trans( grad( v ) ) )  );
    // pressure-velocity terms
    f2 += integrate ( _range = elements( mesh ),
                      _expr = - 1./rho*idt( p )*div( v )
                      + 1./rho*divt( u )*id( q ) );
    // heat diffusion
    f2 += integrate( _range = elements( mesh ),
                     _expr = k*gradt( t )*trans( grad( s ) ) );

    // buyoancy forces c(theta,v)
    f2 +=integrate( _range = elements( mesh ),
                    _expr = -expansion*idt( t )*( trans( oneY() )*id( v ) ) );

    // multipliers for zero-mean pressure (not use in this case)
    f2 += integrate ( _range = elements( mesh ),
                      _expr = idt( xi )*id( eta ) );

    // -- BOUNDARY CONDITIONS -- //
    auto SigmaNt = ( -1./rho*idt( p )*N() + nu*gradt( u )*N() );
    auto SigmaN = ( -1./rho*id( q )*N() + nu*grad( v )*N() );

    // weak Dirichlet condition on the walls (u=0)
    f2 += integrate ( _range = markedfaces( mesh, "walls" ),
                      _expr = - trans( SigmaNt )*id(v)
                      - trans( SigmaN )*idt(u)
                      + nu*penalbc*trans( idt(u) )*id(v)/hFace() );
    f2 += integrate ( _range = markedfaces( mesh, "passengers" ),
                      _expr = - trans( SigmaNt )*id(v)
                      - trans( SigmaN )*idt(u)
                      + nu*penalbc*trans( idt(u) )*id(v)/hFace() );
    f2 += integrate ( _range = markedfaces( mesh, "floor" ),
                      _expr = - trans( SigmaNt )*id(v)
                      - trans( SigmaN )*idt(u)
                      + nu*penalbc*trans( idt(u) )*id(v)/hFace() );

    // weak Dirichlet on inlet (u=U0)
    f2 += integrate ( _range = markedfaces( mesh, "inlet" ),
                      _expr = - trans( SigmaNt )*id(v)
                      - trans( SigmaN )*idt(u)
                       + nu*penalbc*trans( idt(u) )*id(v)/hFace() );
    auto poiseuille = -(0.05-Px())*(Px()+0.05)*400*oneY();
    form1( _test = Xh, _vector=M_Fqm[0][1][0] )
        = integrate ( _range = markedfaces( mesh, "inlet"),
                      _expr = - trans( SigmaN )*poiseuille
                      + nu*penalbc/hFace()*trans(poiseuille)*id(v) );

    // weak Dirichlet on temperature inlet
    f2 += integrate ( _range = markedfaces( mesh, "inlet" ),
                      _expr = - k*gradt( t )*N()*id( s )
                      - k*grad( s )*N()*idt( t )
                      + k*penalbc*idt( t )*id( s )/hFace() );
    form1( _test = Xh, _vector=M_Fqm[0][0][0] )
        = integrate ( _range = markedfaces( mesh, "inlet"),
                      _expr = - k*grad(s)*N()
                      + k*penalbc*id(s)/hFace()  );


    // weak Dirichlet on temperature passengers
    double flux = doption( "passengers-flux" );
    form1( _test = Xh, _vector=M_Fqm[0][2][0] )
        = integrate ( _range = markedfaces( mesh, "passengers"),
                      _expr = flux*id(s)/rho/Cp );


    for (int i=0 ; i<Ql(0) ; i++ )
        M_Fqm[0][i][0]->close();


    M = M_backend->newMatrix( _trial = Xh, _test = Xh);
    form2( _test=Xh, _trial=Xh, _matrix=M ) =
        integrate( _range = elements( mesh ),
                   _expr = trans( id( v ) )*idt( u ) + trace( grad( v )*trans( gradt( u ) ) )
                   + id( q )*idt( p ) + grad( q )*trans( gradt( p ) )
                   + id( s )*idt( t ) + grad( s )*trans( gradt( t ) ) );
    M->close();

    // Output 1 : mean temperature in cabin
    double domain = integrate( _range = elements(mesh), _expr = cst(1.) ).evaluate()(0,0);
    form1( _test=Xh, _vector=M_Fqm[1][0][0] ) =
        integrate( _range = elements(mesh),
                   _expr = id(s)/domain );
    M_Fqm[1][0][0]->close();

    LOG(INFO) << "Natural Convection : Cabin Model is initilized";
}
